#!/usr/bin/env nextflow
/*
 * relBE MGE pipeline -- bait-anchored discovery of mobile genetic elements.
 *
 *   bait -> local BLAST -> Tier A (sequence) -> metadata -> Tier B (record)
 *        -> +/-flank windows -> IS homology + replicon typing
 *
 * The bait is a parameter, so this works for relBE or any other gene.
 * Portable: -profile local | firefly | gadi.
 */

def PY = { params.python_exec ?: 'python3' }

process BUILD_DB {
    publishDir "${projectDir}/blastdb", mode: 'copy'
    input:  path genome_dir
    output: path "${params.db_name}*", emit: db
    script:
    """
    BLAST_BIN='${params.blast_bin}' ${projectDir}/scripts/build_blastdb.sh \\
        ${genome_dir} ${params.db_name} '${params.db_name} (nextflow)'
    """
}

process BAIT_SEARCH {
    publishDir params.outdir, mode: 'copy'
    label 'python'
    input:  path bait
            path db_files
    output: path 'bait_hits.tsv', emit: hits
    script:
    """
    BLAST_BIN='${params.blast_bin}' ${PY()} ${projectDir}/scripts/bait_search.py \\
        --bait ${bait} --db ${params.db_name} --out bait_hits.tsv \\
        --subregions '${params.subregions}' \\
        --full-pident ${params.full_pident} --full-qcov ${params.full_qcov} \\
        --paralogue-max-pident ${params.paralogue_max_pident} \\
        --task ${params.blast_task} --evalue ${params.evalue}
    """
}

process FETCH_METADATA {
    publishDir params.outdir, mode: 'copy'
    label 'python'
    input:  path hits
    output: path 'record_metadata.tsv', emit: meta
    script:
    def off = params.offline_meta ? '--offline' : ''
    """
    ${PY()} ${projectDir}/scripts/fetch_metadata.py \\
        --accessions ${hits} --col sacc --out record_metadata.tsv \\
        --cache ${projectDir}/data/metadata_cache ${off}
    """
}

process SELECT_LOCI {
    publishDir params.outdir, mode: 'copy'
    label 'python'
    input:  path hits
            path meta
    output: path 'loci_selected.tsv', emit: sel
    script:
    """
    ${PY()} ${projectDir}/scripts/select_loci.py \\
        --hits ${hits} --meta ${meta} --out loci_selected.tsv \\
        --flank ${params.flank} --tier-a '${params.tier_a}' \\
        --replicon-type '${params.replicon_type}' --min-flank ${params.min_flank}
    """
}

process EXTRACT_LOCI {
    publishDir params.outdir, mode: 'copy', pattern: 'loci.tsv'
    publishDir params.loci_dir, mode: 'copy', pattern: '*.fa'
    label 'python'
    input:  path sel
    output: path '*.fa',      emit: loci
            path 'loci.tsv',  emit: table
    script:
    """
    ${PY()} ${projectDir}/scripts/extract_loci.py \\
        --selected ${sel} --fasta-dir ${params.genome_dir} \\
        --out-dir . --out-table loci.tsv --flank ${params.flank}
    """
}

process IS_SEARCH {
    publishDir "${params.outdir}/is", mode: 'copy'
    label 'python'
    input:  path locus
    output: path "*.is.tsv", emit: hits
    script:
    """
    BLAST_BIN='${params.blast_bin}' ${PY()} ${projectDir}/scripts/is_search.py \\
        --locus ${locus} --db ${params.isfinder_db} \\
        --out ${locus.baseName}.is.tsv \\
        --min-length ${params.is_min_length} --evalue ${params.is_evalue} \\
        --max-target ${params.is_max_target}
    """
}

process PLASMIDFINDER {
    publishDir "${params.outdir}/pf", mode: 'copy'
    label 'python'
    input:  tuple val(acc), path(replicon)
    output: path "*.pf.tsv", emit: hits
    script:
    """
    BLAST_BIN='${params.blast_bin}' ${PY()} ${projectDir}/scripts/plasmidfinder.py \\
        --replicon ${replicon} --db ${params.plasmidfinder_fa} \\
        --accession ${acc} --out ${acc}.pf.tsv \\
        --min-cov ${params.pf_min_cov} --min-id ${params.pf_min_id}
    """
}

process COLLATE {
    publishDir params.outdir, mode: 'copy'
    label 'python'
    input:  path is_tsv
            path pf_tsv
            path loci_tsv
    output: path 'mge_hits.tsv'
            path 'replicons.tsv'
            path 'summary.tsv'
    script:
    """
    ${PY()} ${projectDir}/scripts/collate.py \\
        --is-dir . --pf-dir . --loci ${loci_tsv} \\
        --out-mge mge_hits.tsv --out-rep replicons.tsv --out-summary summary.tsv
    """
}

workflow {
    bait    = Channel.fromPath(params.bait,       checkIfExists: true)
    gdir    = Channel.fromPath(params.genome_dir, checkIfExists: true, type: 'dir')

    db      = BUILD_DB(gdir)
    hits    = BAIT_SEARCH(bait, db.db)
    meta    = FETCH_METADATA(hits.hits)
    sel     = SELECT_LOCI(hits.hits, meta.meta)
    loci    = EXTRACT_LOCI(sel.sel)

    is_ch   = IS_SEARCH(loci.loci.flatten())

    // PlasmidFinder runs on the WHOLE PARENT REPLICON, not the window.
    rep_ch  = loci.table
                  .splitCsv(header: true, sep: '\t')
                  .map { r -> tuple(r.accession,
                                    file("${params.genome_dir}/${r.accession}.1.fna")) }
                  .filter { it[1].exists() }
                  .unique { it[0] }
    pf_ch   = PLASMIDFINDER(rep_ch)

    COLLATE(is_ch.hits.collect(), pf_ch.hits.collect(), loci.table)
}
