#!/usr/bin/env bash
# Firefly capability probe -- READ ONLY, installs nothing, writes nothing.
# Run once on Firefly (VPN up) and paste the output back.
# It answers every question that gates the Stage 3 cluster move.
echo "=== identity / os ==="
hostname; whoami; uname -srm
( . /etc/os-release 2>/dev/null && echo "$PRETTY_NAME" ) || true

echo; echo "=== scheduler ==="
for s in sbatch qsub bsub; do printf '%-8s ' "$s"; command -v $s || echo MISSING; done
sinfo   -s 2>/dev/null | head -5
qstat   -Q 2>/dev/null | head -5

echo; echo "=== containers (THE decisive one) ==="
for c in apptainer singularity podman docker; do
  printf '%-12s ' "$c"
  if command -v $c >/dev/null 2>&1; then $c --version 2>&1 | head -1; else echo MISSING; fi
done

echo; echo "=== outbound internet from THIS node ==="
curl -sS -m 12 -o /dev/null -w 'ncbi      http=%{http_code}\n' \
  https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi 2>&1 || echo "ncbi      BLOCKED"
curl -sS -m 12 -o /dev/null -w 'quay.io   http=%{http_code}\n' \
  https://quay.io/v2/ 2>&1 || echo "quay.io   BLOCKED"
echo "proxy vars: ${http_proxy:-none} ${https_proxy:-none}"

echo; echo "=== java / nextflow ==="
command -v nextflow && nextflow -version 2>&1 | head -3 || echo "nextflow MISSING"
command -v java && java -version 2>&1 | head -1 || echo "java MISSING"
module avail 2>&1 | grep -iE "nextflow|java|openjdk|singularity|apptainer" | head -10

echo; echo "=== tools already present ==="
for t in blastn makeblastdb seqkit mmseqs mafft diamond hmmsearch prodigal mash; do
  printf '%-14s ' "$t"; command -v $t || echo MISSING
done

echo; echo "=== python ==="
python3 -V 2>&1
python3 -c "import Bio,pandas;print('biopython',Bio.__version__,'pandas',pandas.__version__)" 2>&1 | head -1

echo; echo "=== filesystems / quota ==="
echo "HOME=$HOME"; df -h "$HOME" 2>/dev/null | tail -1
for d in /data /scratch /project /g/data /srv; do [ -d "$d" ] && df -h "$d" 2>/dev/null | tail -1; done
ls -d /data/projects/* 2>/dev/null | head -5

echo; echo "=== cpu / mem on login node ==="
nproc 2>/dev/null; free -h 2>/dev/null | head -2
