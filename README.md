# Shrimp pond 16S rRNA microbiome: short-read vs full-length sequencing

Analysis workflow for:

> Rajonhson D.M. *et al.* (2024). Integrating short- and full-length 16S rRNA gene sequencing to elucidate microbiome profiles in Pacific white shrimp (*Litopenaeus vannamei*) ponds. *Microbiology Spectrum* 12(11): e00965-24. https://doi.org/10.1128/spectrum.00965-24

Raw reads: NCBI BioProject [PRJNA1087723](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1087723)

The same samples (shrimp intestine, pond water and sediment from two ponds, without and with probiotic supplementation) were sequenced two ways:

- **Short read**: Illumina MiSeq, 250 bp paired-end, V3–V4 region (primers 341F/806R)
- **Full length**: Oxford Nanopore, R10.4.1 flow cell, near-complete 16S rRNA gene (primers 27F/1492R)

The two platforms differ only in read processing (steps 1–2). From step 3 onward both datasets are converted to the same table format and go through the same analysis, so any difference in the results comes from the sequencing strategy and not from the analysis. The two sets of results are then compared.

## Workflow

```
            Illumina V3-V4                          Nanopore full-length
                  |                                          |
 1a  fastp: Q28, >=200 bp, adapters,       1b  NanoFilt: 1,000-2,000 bp, Q>=10
     poly-X >= 6 bp                                (after Guppy 6.5.7 SUP basecalling)
                  |                                          |
 2a  DADA2: primer trim, filter, error      2b  Emu: species-level abundance
     model, denoise, merge, chimera              (SILVA v138.1)
     removal, SILVA v138.1 taxonomy                          |
                  |                                          |
                  +-------------------+----------------------+
                                      |
                 3   phylum / genus / species count tables (identical format)
                                      |
                 4   downstream analysis in R (same code for both platforms)
                     rarefaction, taxonomic profiles, alpha diversity
                     (Shannon, observed features, Pielou; Wilcoxon, Kruskal-Wallis),
                     beta diversity (Bray-Curtis PCoA, PERMANOVA), DESeq2 differential
                     abundance, LEfSe, Spearman co-occurrence networks
                                      |
                 5   comparison of short-read and full-length results
```

| Step | Script | Tools |
|---|---|---|
| 1a | `workflow/01a_qc_illumina_fastp.sh` | fastp |
| 1b | `workflow/01b_qc_nanopore_nanofilt.sh` | NanoFilt |
| 2a | `workflow/02a_denoise_illumina_dada2.R` | DADA2, SILVA v138.1 |
| 2b | `workflow/02b_classify_nanopore_emu.sh` | Emu, SILVA v138.1 |
| 3 | `workflow/03_build_taxon_tables.R` | R |
| 4 | `analysis/short_read.Rmd`, `analysis/long_read.Rmd` | phyloseq, vegan, DESeq2, lefser, igraph |
| — | `R/` | helper functions used by the notebooks |

> **Note on step 2a.** The published short-read analysis ran DADA2 inside QIIME2 2022.2 (q2-cutadapt, q2-dada2, q2-classify-sklearn with SILVA 138). This repository gives a standalone version: fastp quality control followed by the DADA2 R package, with primers removed by fixed-length trimming and taxonomy assigned against SILVA v138.1. It follows the same denoising logic but is not a byte-for-byte rerun of the QIIME2 run.

## Running it

```bash
# 0. environment
conda env create -f environment.yml
conda activate shrimp-16s
Rscript -e 'remotes::install_github("twbattaglia/btools")'
Rscript -e 'remotes::install_github("NicolasH2/ggdendroplot")'

# 1. quality control
bash workflow/01a_qc_illumina_fastp.sh raw/illumina work/fastp 10
bash workflow/01b_qc_nanopore_nanofilt.sh raw/nanopore work/nanofilt

# 2. denoising / classification
Rscript workflow/02a_denoise_illumina_dada2.R work/fastp work/dada2 silva_nr99_v138.1_wSpecies_train_set.fa.gz 10
bash workflow/02b_classify_nanopore_emu.sh work/nanofilt emu_silva_db work/emu 8

# 3. count tables in a common format
Rscript workflow/03_build_taxon_tables.R dada2 work/dada2 data/short-read
Rscript workflow/03_build_taxon_tables.R emu work/emu/emu-combined-species-counts.tsv data/long-read

# 4. downstream analysis (figures and tables go to results/short-read and results/long-read)
Rscript -e 'rmarkdown::render("analysis/short_read.Rmd")'
Rscript -e 'rmarkdown::render("analysis/long_read.Rmd")'
```

**Sample metadata.** `data/metadata.csv` describes the 33 sequenced samples: shrimp intestine (7 from pond A, where 3 of 10 samples had too little DNA, and 10 from pond B), sediment (4 per pond) and water (4 per pond) and is used by both notebooks. Columns: `sample_id` (e.g. IA01 = intestine, pond A, sample 01), `shortread` (the Illumina sample name, e.g. IN-A1, used to rename short-read columns), `type` (intestine / sediment / water) and `pond` (A = no probiotic, B = probiotic supplementation).

File naming: Illumina FASTQ files are expected as `<shortread>_1.<...>.gz` / `<shortread>_2.<...>.gz` (e.g. `IN-A1_1.fq.gz`); Nanopore FASTQ files as `<sample_id>.fastq.gz` (e.g. `IA01.fastq.gz`).

**Databases.** DADA2-formatted SILVA v138.1: https://zenodo.org/records/4587955. Emu SILVA database: see the Emu documentation (https://github.com/treangenlab/emu).

## Comparing the two platforms

Each notebook writes its figures and tables to its own folder (`results/short-read/`, `results/long-read/`). Because both started from tables in the same format and went through the same code, the outputs can be compared directly: taxonomic composition at phylum, genus and species level, the taxa detected by one platform only, alpha and beta diversity, and differentially abundant taxa between ponds. The comparison reported in the article was made from these outputs.

## Licence and citation

Code released under the MIT licence. If you use it, please cite the article above (see `CITATION.cff`).

Contact: Dora M. Rajonhson — ORCID [0000-0003-3247-4510](https://orcid.org/0000-0003-3247-4510)
