#!/usr/bin/env nextflow
/*
 * TAflank -- bait-anchored discovery of mobile genetic elements in
 * gene flanking regions.
 *
 *   bait -> local BLAST -> Tier A (sequence) -> metadata -> Tier B (record)
 *        -> +/-flank windows -> IS homology + replicon typing
 *
 * The bait is a parameter, so this works for relBE or any other gene.
 * Portable: -profile local | firefly | gadi.
 */

def PY = { params.python_exec ?: 'python3' }

// Nextflow parses `--tier_a ''` on the command line as the BOOLEAN true, so an
// "empty means keep everything" filter is unexpressible that way: the value
// arrives as the string "true" and matches nothing. Accept `all` (and a bare
// flag) as the explicit way to disable a filter.
def filterOf = { v ->
    (v == null || v == true || v.toString().toLowerCase() in ['all', 'any', '']) ? '' : v.toString()
}

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
    publishDir params.outdir, mode: 'copy', pattern: 'bait_hits.tsv'
    publishDir "${params.outdir}/blast", mode: 'copy', pattern: 'blast/*.tsv', saveAs: { it.split('/')[-1] }
    label 'python'
    cpus params.blast_cpus
    // db_files is empty for a pre-formatted database: nt_prok is far too large
    // to stage per task, so it is referenced by absolute path and left in place.
    input:  path bait
            tuple val(db), path(db_files)
    output: path 'bait_hits.tsv',  emit: hits
            path 'blast/*.tsv',     emit: staged
    script:
    def tax = params.taxonomy ? '--taxonomy' : ''
    """
    BLAST_BIN='${params.blast_bin}' ${PY()} ${projectDir}/scripts/bait_search.py \\
        --bait ${bait} --db ${db} --out bait_hits.tsv \\
        --subregions '${filterOf(params.subregions)}' \\
        --full-pident ${params.full_pident} --full-qcov ${params.full_qcov} \\
        --paralogue-max-pident ${params.paralogue_max_pident} \\
        --task ${params.blast_task} --evalue ${params.evalue} \\
        --max-target-seqs ${params.max_target_seqs} \\
        --num-threads ${task.cpus} ${tax} \\
        --raw-out blast/01_raw.tsv --kept-out blast/02_tier_a.tsv
    """
}

process IMPORT_BLAST {
    publishDir params.outdir, mode: 'copy', pattern: 'bait_hits.tsv'
    publishDir "${params.outdir}/blast", mode: 'copy', pattern: 'blast/*.tsv', saveAs: { it.split('/')[-1] }
    label 'python'
    // An NCBI web BLAST "Hit Table" download IS outfmt 6, so results from the
    // online search can be filtered and processed by exactly the same code as
    // a local one. It carries no slen/stitle, so metadata must come from
    // esummary -- params.meta_source is forced to 'ncbi' below.
    input:  path hit_table
            path bait
    output: path 'bait_hits.tsv', emit: hits
            path 'blast/*.tsv',   emit: staged
    script:
    """
    ${PY()} ${projectDir}/scripts/import_blast.py \\
        --hit-table ${hit_table} --bait ${bait} --out bait_hits.tsv \\
        --subregions '${filterOf(params.subregions)}' \\
        --full-pident ${params.full_pident} --full-qcov ${params.full_qcov} \\
        --paralogue-max-pident ${params.paralogue_max_pident} \\
        --raw-out blast/01_raw.tsv --kept-out blast/02_tier_a.tsv
    """
}

process FETCH_METADATA {
    publishDir params.outdir, mode: 'copy'
    label 'python'
    input:  path hits
    output: path 'record_metadata.tsv', emit: meta
    script:
    // 'blast'  -> derive from stitle/slen, NO network. Required on compute
    //             nodes without outbound access, and the only workable mode
    //             for a sweep where one esummary per hit is hopeless.
    // 'cache'  -> rebuild from the on-disk cache, no network.
    // 'ncbi'   -> batched esummary (default for small curated sets).
    def mode = params.meta_source == 'blast' ? '--from-blast'
             : params.meta_source == 'cache' ? '--offline' : ''
    """
    ${PY()} ${projectDir}/scripts/fetch_metadata.py \\
        --accessions ${hits} --col sacc --out record_metadata.tsv \\
        --cache ${params.meta_cache} ${mode}
    """
}

process SELECT_LOCI {
    publishDir params.outdir, mode: 'copy', pattern: '*.tsv'
    publishDir "${params.outdir}/blast", mode: 'copy', pattern: 'blast/*.tsv', saveAs: { it.split('/')[-1] }
    label 'python'
    input:  path hits
            path meta
    output: path 'loci_selected.tsv',  emit: sel
            path 'filter_summary.tsv',  emit: audit
            path 'blast/*.tsv',         emit: staged
    script:
    """
    ${PY()} ${projectDir}/scripts/select_loci.py \\
        --hits ${hits} --meta ${meta} --out loci_selected.tsv \\
        --flank ${params.flank} --tier-a '${filterOf(params.tier_a)}' \\
        --replicon-type '${filterOf(params.replicon_type)}' --min-flank ${params.min_flank} \\
        --summary-out filter_summary.tsv \\
        --blast-out blast/03_selected.tsv
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
    def src = params.genome_dir ? "--fasta-dir ${params.genome_dir}" : ''
    def bdb = params.use_blastdb_seq ? "--blastdb ${params.blastdb_path}" : ''
    // Imported web-BLAST hits have no local copy, so the window comes from NCBI.
    def efe = params.hit_table ? "--efetch --efetch-cache ${params.efetch_cache}" : ''
    """
    BLAST_BIN='${params.blast_bin}' ${PY()} ${projectDir}/scripts/extract_loci.py \\
        --selected ${sel} ${src} ${bdb} ${efe} \\
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
    input:  val acc
    output: path "*.pf.tsv", emit: hits
    script:
    // Whole parent replicon, never the window: Inc rep genes are rarely within
    // 20 kb of the bait. Source is a FASTA if we have one, else the BLAST db.
    def fa = file("${params.genome_dir}/${acc}.1.fna")
    def src = (params.genome_dir && fa.exists()) ? "--replicon ${fa}"
            : params.hit_table ? "--efetch --efetch-cache ${params.efetch_cache}"
            : "--blastdb ${params.blastdb_path}"
    """
    BLAST_BIN='${params.blast_bin}' ${PY()} ${projectDir}/scripts/plasmidfinder.py \\
        ${src} --db ${params.plasmidfinder_fa} \\
        --accession ${acc} --out ${acc}.pf.tsv \\
        --min-cov ${params.pf_min_cov} --min-id ${params.pf_min_id}
    """
}

process CARGO_SEARCH {
    publishDir "${params.outdir}/cargo", mode: 'copy'
    label 'python'
    input:  path locus
    output: path "*.cargo.tsv", emit: hits
    script:
    """
    BLAST_BIN='${params.blast_bin}' ${PY()} ${projectDir}/scripts/cargo_search.py \\
        --locus ${locus} --db ${params.cargo_fa} \\
        --out ${locus.baseName}.cargo.tsv \\
        --min-pident ${params.cargo_min_pident} --min-cov ${params.cargo_min_cov}
    """
}

process BAKTA {
    publishDir "${params.outdir}/bakta", mode: 'copy'
    label 'bakta'
    cpus 4
    input:  path locus
    output: tuple val(locus.baseName), path("${locus.baseName}.gff3"), emit: gff
            path "${locus.baseName}.faa",                            emit: faa
    // .faa is what Stage 3 clusters. Without it declared here,
    // publishDir copies only the GFF3 and the proteins die with the
    // work directory -- which on Gadi is \$PBS_JOBFS, gone at job end.
    when:   params.run_bakta
    script:
    """
    bakta --db ${params.bakta_db} --output . --prefix ${locus.baseName} \\
          --threads ${task.cpus} --keep-contig-headers --force \\
          ${params.bakta_args} ${locus}
    """
}

process GENE_BLOCKS {
    publishDir "${params.outdir}/blocks", mode: 'copy'
    label 'python'
    input:  tuple val(lid), path(locus), path(is_tsv), path(cargo_tsv), path(gff)
    output: path "*.blocks.tsv", emit: blocks
    script:
    // Annotation source, in order: Bakta GFF -> parent GenBank -> none.
    // The GenBank fallback exists so this merge is testable before Bakta is
    // available; Bakta cannot run on macOS at all (amrfinderplus has no macOS
    // build). Swapping it in changes the source, not the logic.
    def acc  = lid.split('__')[0]
    def gbk  = file("${params.gbk_dir}/${acc}.1.gbk")
    def gsrc = gff.name != 'NO_GFF' ? "--gff ${gff}"
             : (gbk.exists() ? "--gbk ${gbk}" : '')
    """
    ${PY()} ${projectDir}/scripts/gene_blocks.py \\
        --locus ${locus} --is ${is_tsv} --cargo ${cargo_tsv} ${gsrc} \\
        --out ${lid}.blocks.tsv
    """
}

process CLASSIFY_CONTEXT {
    publishDir params.outdir, mode: 'copy'
    label 'python'
    input:  path loci_tsv
            path is_tsv
            path cargo_tsv
    output: path 'context.tsv', emit: ctx
    script:
    """
    ${PY()} ${projectDir}/scripts/classify_context.py \\
        --loci ${loci_tsv} --is-dir . --cargo-dir . --out context.tsv \
        --ends '${params.vehicle_ends}'
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
        --out-mge mge_hits.tsv --out-rep replicons.tsv --out-summary summary.tsv \
        --aliases ${projectDir}/assets/replicon_aliases.tsv
    """
}

process RUN_REPORT {
    publishDir params.outdir, mode: 'copy'
    label 'python'
    // Reads only the PUBLISHED tables, so the report is reproducible from
    // results/ alone and never needs the pipeline re-run to regenerate.
    input:  val ready
    output: path 'run_report.txt'
            path 'run_report.tsv'
    script:
    """
    ${PY()} ${projectDir}/scripts/run_report.py \\
        --results ${params.outdir} \\
        --out run_report.txt --tsv run_report.tsv
    """
}

workflow {
    bait = Channel.fromPath(params.bait, checkIfExists: true)

    // Either build a db from local FASTAs, or point at a pre-formatted one
    // (update_blastdb.pl nt_prok ...). A pre-formatted db is referenced by
    // absolute path and never staged -- nt_prok is far too large to copy per
    // task, and it already contains every sequence we need.
    if (params.hit_table) {
        // Entry point for an NCBI web BLAST result: no database is searched,
        // so BUILD_DB and BAIT_SEARCH are both skipped.
        hits = IMPORT_BLAST(Channel.fromPath(params.hit_table, checkIfExists: true),
                            bait)
    } else if (params.blastdb) {
        db_ch = Channel.value(tuple(params.blastdb, []))
        hits  = BAIT_SEARCH(bait, db_ch)
    } else {
        gdir  = Channel.fromPath(params.genome_dir, checkIfExists: true, type: 'dir')
        db_ch = BUILD_DB(gdir).db.map { f -> tuple(params.db_name, f) }
        hits  = BAIT_SEARCH(bait, db_ch)
    }
    meta    = FETCH_METADATA(hits.hits)
    sel     = SELECT_LOCI(hits.hits, meta.meta)
    loci    = EXTRACT_LOCI(sel.sel)

    loci_fa = loci.loci.flatten()
    is_ch   = IS_SEARCH(loci_fa)

    // Tier C input. Deliberately NOT an annotation pass: ISfinder covers the
    // IS layer and a small curated set covers the cargo families the biology
    // turns on, so classification is cheap. Bakta is for the figure's complete
    // gene track and runs later, on far fewer loci.
    cargo_ch = CARGO_SEARCH(loci_fa)

    // Bakta and ISfinder are COMPLEMENTARY, not alternatives, and independent
    // of each other -- so they run in parallel and are merged afterwards, with
    // the ISfinder element name overriding annotation's generic family label.
    bakta_ch = params.run_bakta ? BAKTA(loci_fa).gff
                                : Channel.empty()

    keyed = loci_fa.map        { f -> tuple(f.baseName, f) }
        .join(is_ch.hits.map   { f -> tuple(f.baseName.replace('.is',''), f) })
        .join(cargo_ch.hits.map{ f -> tuple(f.baseName.replace('.cargo',''), f) })
        .join(bakta_ch, remainder: true)
        .map { lid, fa, is_t, cg, gff -> tuple(lid, fa, is_t, cg,
                                               gff ?: file('NO_GFF')) }
    blocks_ch = GENE_BLOCKS(keyed)

    // PlasmidFinder runs on the WHOLE PARENT REPLICON, not the window.
    rep_ch  = loci.table
                  .splitCsv(header: true, sep: '\t')
                  .map { r -> r.accession }
                  .unique()
    pf_ch   = PLASMIDFINDER(rep_ch)

    CLASSIFY_CONTEXT(loci.table, is_ch.hits.collect(), cargo_ch.hits.collect())

    collated = COLLATE(is_ch.hits.collect(), pf_ch.hits.collect(), loci.table)

    // Last: one readable account of what went in, what each filter removed,
    // and what came out at every layer. Gated on everything else finishing.
    RUN_REPORT(collated[0].mix(blocks_ch, CLASSIFY_CONTEXT.out.ctx).collect())
}
