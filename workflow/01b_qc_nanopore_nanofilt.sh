#!/bin/bash
# Step 1b - Quality control of Oxford Nanopore full-length 16S reads with NanoFilt
#
# Reads were basecalled, adapter-trimmed and demultiplexed beforehand with Guppy 6.5.7
# in super-accuracy (SUP) mode (R10.4.1 flow cell, SQK-NBD114.96 kit).
# NanoFilt then keeps only near-complete 16S rRNA gene reads:
#   - length between 1,000 and 2,000 bp
#   - mean read quality (Q-score) of at least 10
#
# Usage:  bash workflow/01b_qc_nanopore_nanofilt.sh <demultiplexed_fastq_dir> <output_dir>
# One gzipped FASTQ per sample is expected: <sample>.fastq.gz

set -euo pipefail

in_dir=${1:?"give the folder with demultiplexed FASTQ files"}
out_dir=${2:?"give the output folder"}
mkdir -p "${out_dir}"

shopt -s nullglob
files=("${in_dir}"/*.fastq.gz)
if [ ${#files[@]} -eq 0 ]; then
  echo "No *.fastq.gz files found in ${in_dir}" >&2
  exit 1
fi

for fq in "${files[@]}"; do
  sample=$(basename "${fq}" .fastq.gz)
  gunzip -c "${fq}" \
    | NanoFilt --length 1000 --maxlength 2000 --quality 10 \
    | gzip > "${out_dir}/${sample}.filtered.fastq.gz"
done

echo "NanoFilt finished: ${#files[@]} samples in ${out_dir}"
