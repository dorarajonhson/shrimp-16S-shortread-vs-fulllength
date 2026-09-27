#!/bin/bash
# Step 2b - Species-level classification of full-length 16S reads with Emu
#
# Emu estimates species-level relative abundances and read counts from the filtered
# Nanopore reads of step 1b, using the SILVA v138.1 database formatted for Emu
# (taxonomy names and identifiers from the DADA2 SILVA release).
#
# Usage:  bash workflow/02b_classify_nanopore_emu.sh <filtered_fastq_dir> <emu_silva_db_dir> <output_dir> [threads]

set -euo pipefail

in_dir=${1:?"give the folder with NanoFilt output"}
db_dir=${2:?"give the Emu SILVA database folder"}
out_dir=${3:?"give the output folder"}
threads=${4:-8}
mkdir -p "${out_dir}"

shopt -s nullglob
files=("${in_dir}"/*.filtered.fastq.gz)
if [ ${#files[@]} -eq 0 ]; then
  echo "No *.filtered.fastq.gz files found in ${in_dir}" >&2
  exit 1
fi

for fq in "${files[@]}"; do
  sample=$(basename "${fq}" .filtered.fastq.gz)
  emu abundance "${fq}" \
    --db "${db_dir}" \
    --type map-ont \
    --keep-counts \
    --output-dir "${out_dir}" \
    --output-basename "${sample}" \
    --threads "${threads}"
done

# one table with estimated read counts for all samples and all taxonomic ranks
emu combine-outputs "${out_dir}" species --counts

echo "Emu finished: ${#files[@]} samples in ${out_dir}"
