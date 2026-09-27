# Step 3 - Absolute abundance tables at phylum, genus and species level
#
# Turns the output of step 2a (DADA2) or step 2b (Emu) into the three tables read by
# the analysis notebooks: level-2 (phylum), level-6 (genus) and level-7 (species).
# Rows are samples, columns are full lineages written as
#   d__Bacteria;p__Proteobacteria;c__...;g__Vibrio
# with "__" for a missing rank. Both platforms therefore enter the downstream analysis
# in exactly the same format.
#
# Usage:
#   Rscript workflow/03_build_taxon_tables.R dada2 <dada2_output_dir> data/short-read
#   Rscript workflow/03_build_taxon_tables.R emu   <emu-combined-species-counts.tsv> data/long-read

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3 || !args[1] %in% c("dada2", "emu")) {
  stop("Usage: Rscript 03_build_taxon_tables.R <dada2|emu> <input> <output_dir>")
}
mode    <- args[1]
input   <- args[2]
out_dir <- args[3]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ranks    <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
prefixes <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")

## ---- Read counts (samples x features) and taxonomy (features x 7 ranks) ------------
if (mode == "dada2") {
  counts <- readRDS(file.path(input, "asv_table.rds"))      # samples x ASVs
  tax    <- readRDS(file.path(input, "asv_taxonomy.rds"))   # ASVs x ranks
  tax    <- as.matrix(tax[colnames(counts), , drop = FALSE])
  if (ncol(tax) < 7) tax <- cbind(tax, matrix(NA, nrow(tax), 7 - ncol(tax)))
  tax    <- tax[, 1:7, drop = FALSE]
} else {
  emu <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE)
  lower <- tolower(colnames(emu))
  emu_ranks <- list(Domain  = c("superkingdom", "domain", "kingdom"),
                    Phylum  = "phylum", Class = "class", Order = "order",
                    Family  = "family", Genus = "genus", Species = "species")
  rank_cols <- sapply(emu_ranks, function(r) {
    hit <- which(lower %in% r)
    if (length(hit) == 0) NA_integer_ else hit[1]
  })
  if (is.na(rank_cols["Species"])) stop("No 'species' column found in ", input)
  info_cols   <- c(which(lower %in% c(unlist(emu_ranks), "tax_id", "clade", "subspecies", "species subgroup", "species group")))
  sample_cols <- setdiff(seq_along(emu), info_cols)
  tax <- sapply(names(emu_ranks), function(r) {
    if (is.na(rank_cols[r])) rep(NA_character_, nrow(emu)) else as.character(emu[[rank_cols[r]]])
  })
  counts <- t(as.matrix(emu[, sample_cols, drop = FALSE]))  # samples x taxa
  counts[is.na(counts)] <- 0
  counts <- round(counts)
}
colnames(tax) <- ranks

## ---- Build lineage strings and collapse counts ------------------------------------------
clean <- function(x) {
  x <- ifelse(is.na(x) | x == "", "", x)
  gsub(" ", "_", x)
}
lineage <- function(level) {
  parts <- sapply(seq_len(level), function(i) {
    v <- clean(tax[, i])
    ifelse(v == "", "__", paste0(prefixes[i], v))
  })
  if (is.null(dim(parts))) parts <- matrix(parts, nrow = 1)
  apply(parts, 1, paste, collapse = ";")
}

write_level <- function(level, file_level) {
  lin <- lineage(level)
  collapsed <- t(rowsum(t(counts), group = lin))            # samples x lineages
  collapsed <- collapsed[, colSums(collapsed) > 0, drop = FALSE]
  tab <- data.frame(index = rownames(collapsed), collapsed, check.names = FALSE)
  write.csv(tab, file.path(out_dir, sprintf("level-%d.csv", file_level)),
            row.names = FALSE, quote = TRUE)
  message(sprintf("level-%d: %d taxa x %d samples", file_level, ncol(collapsed), nrow(collapsed)))
}

write_level(2, 2)   # phylum
write_level(6, 6)   # genus
write_level(7, 7)   # species
