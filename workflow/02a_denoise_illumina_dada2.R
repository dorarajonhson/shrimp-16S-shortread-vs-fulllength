# Step 2a - Denoising of Illumina V3-V4 reads with DADA2 (amplicon sequence variants)
#
# Input : fastp-trimmed read pairs from step 1a (<sample>_R1/_R2.fastp-trim.fq.gz)
# Output: ASV count table, taxonomy (SILVA v138.1) and a table tracking reads kept
#         at each step, all written to the output folder.
#
# Usage:  Rscript workflow/02a_denoise_illumina_dada2.R <trimmed_dir> <output_dir> <silva_train_set> [threads]
#   <silva_train_set>: silva_nr99_v138.1_wSpecies_train_set.fa.gz (DADA2-formatted SILVA,
#                      https://zenodo.org/records/4587955)

suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 02a_denoise_illumina_dada2.R <trimmed_dir> <output_dir> <silva_train_set> [threads]")
}
in_dir   <- args[1]
out_dir  <- args[2]
silva    <- args[3]
threads  <- if (length(args) >= 4) as.integer(args[4]) else TRUE

## ---- Parameters ---------------------------------------------------------------
## Primers 341F (CCTACGGGNGGCWGCAG, 17 nt) and 806R (GGACTACHVGGGTWTCTAAT, 20 nt)
## are removed by trimming their length from the 5' end of each read.
trim_left  <- c(17, 20)
## Truncation positions for 250 bp paired-end MiSeq reads. Reads shorter than
## truncLen are discarded, so check the quality profiles before changing these.
trunc_len  <- c(250, 231)
max_ee     <- c(2, 2)
trunc_q    <- 2

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "filtered"), showWarnings = FALSE)

## ---- Input files ----------------------------------------------------------------
fnFs <- sort(list.files(in_dir, pattern = "_R1\\.fastp-trim\\.fq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(in_dir, pattern = "_R2\\.fastp-trim\\.fq\\.gz$", full.names = TRUE))
if (length(fnFs) == 0 || length(fnFs) != length(fnRs)) {
  stop("Forward and reverse files missing or unpaired in ", in_dir)
}
sample.names <- sub("_R1\\.fastp-trim\\.fq\\.gz$", "", basename(fnFs))
stopifnot(identical(sample.names, sub("_R2\\.fastp-trim\\.fq\\.gz$", "", basename(fnRs))))

## Quality profiles of the first samples, for choosing truncation positions
pdf(file.path(out_dir, "quality_profiles.pdf"))
print(plotQualityProfile(fnFs[1:min(2, length(fnFs))]))
print(plotQualityProfile(fnRs[1:min(2, length(fnRs))]))
dev.off()

## ---- Filter and trim --------------------------------------------------------------
filtFs <- file.path(out_dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(out_dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     trimLeft = trim_left, truncLen = trunc_len,
                     maxN = 0, maxEE = max_ee, truncQ = trunc_q, rm.phix = TRUE,
                     compress = TRUE, verbose = TRUE, multithread = threads)
rownames(out) <- sample.names

## samples with no reads left after filtering are dropped
keep   <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[keep]
filtRs <- filtRs[keep]

## ---- Error model, denoising, merging -------------------------------------------------
set.seed(2024)
errF <- learnErrors(filtFs, multithread = threads)
errR <- learnErrors(filtRs, multithread = threads)
pdf(file.path(out_dir, "error_models.pdf"))
print(plotErrors(errF, nominalQ = TRUE))
print(plotErrors(errR, nominalQ = TRUE))
dev.off()

dadaFs <- dada(filtFs, err = errF, multithread = threads)
dadaRs <- dada(filtRs, err = errR, multithread = threads)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

## ---- ASV table and chimera removal --------------------------------------------------
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    minFoldParentOverAbundance = 5,
                                    multithread = threads, verbose = TRUE)
message(sprintf("Chimeras removed: %.1f%% of reads kept", 100 * sum(seqtab.nochim) / sum(seqtab)))

## ---- Taxonomy (SILVA v138.1, down to species) ---------------------------------------
taxa <- assignTaxonomy(seqtab.nochim, silva, multithread = threads, verbose = TRUE)

## ---- Read tracking --------------------------------------------------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(out[names(filtFs), , drop = FALSE],
               denoisedF = sapply(dadaFs, getN),
               denoisedR = sapply(dadaRs, getN),
               merged    = sapply(mergers, getN),
               nonchim   = rowSums(seqtab.nochim))
colnames(track)[1:2] <- c("input", "filtered")

## ---- Save ----------------------------------------------------------------------------------
saveRDS(seqtab.nochim, file.path(out_dir, "asv_table.rds"))
saveRDS(taxa,          file.path(out_dir, "asv_taxonomy.rds"))
write.csv(track, file.path(out_dir, "read_tracking.csv"))
message("DADA2 finished: ", ncol(seqtab.nochim), " ASVs in ", nrow(seqtab.nochim), " samples")
