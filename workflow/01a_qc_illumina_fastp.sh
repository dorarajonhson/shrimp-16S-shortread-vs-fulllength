#!/bin/bash
# Step 1a - Quality control of Illumina V3-V4 paired-end reads with fastp
#
# Filters applied to every read pair:
#   - bases with Phred quality < 28 count as unqualified (--qualified_quality_phred 28)
#   - reads shorter than 200 bp after trimming are discarded (--length_required 200)
#   - adapters detected automatically for paired-end data
#   - poly-X tails of 6 bp or longer are trimmed
# An md5 checksum is written for every trimmed file.
#
# Usage:  bash workflow/01a_qc_illumina_fastp.sh <raw_fastq_dir> <output_dir> [threads]
# Input file names are expected as <sample>_1.<...>.gz and <sample>_2.<...>.gz
# Raw files are never modified or deleted.

set -euo pipefail

raw_dir=${1:?"give the folder with raw FASTQ files"}
out_dir=${2:?"give the output folder"}
threads=${3:-10}

mkdir -p "${out_dir}/reports"
checksums="${out_dir}/trimmed_fastq_checksums.md5"
: > "${checksums}"

shopt -s nullglob
r1_files=("${raw_dir}"/*_1.*.gz)
if [ ${#r1_files[@]} -eq 0 ]; then
  echo "No *_1.*.gz files found in ${raw_dir}" >&2
  exit 1
fi

# record which files went into the analysis
ls "${raw_dir}"/*.gz > "${out_dir}/fastq.list.txt"

for r1 in "${r1_files[@]}"; do
  r2=${r1/_1./_2.}
  if [ ! -f "${r2}" ]; then
    echo "Missing mate for ${r1}" >&2
    exit 1
  fi
  sample=$(basename "${r1}" | awk -F"_1." '{print $1}')

  fastp \
    --in1 "${r1}" \
    --in2 "${r2}" \
    --qualified_quality_phred 28 \
    --length_required 200 \
    --detect_adapter_for_pe \
    --trim_poly_x \
    --poly_x_min_len 6 \
    --thread "${threads}" \
    --html "${out_dir}/reports/${sample}.fastp.html" \
    --json "${out_dir}/reports/${sample}.fastp.json" \
    --out1 "${out_dir}/${sample}_R1.fastp-trim.fq.gz" \
    --out2 "${out_dir}/${sample}_R2.fastp-trim.fq.gz"

  md5sum "${out_dir}/${sample}_R1.fastp-trim.fq.gz" "${out_dir}/${sample}_R2.fastp-trim.fq.gz" >> "${checksums}"
done

echo "fastp finished: $(ls "${out_dir}"/*_R1.fastp-trim.fq.gz | wc -l) samples in ${out_dir}"
