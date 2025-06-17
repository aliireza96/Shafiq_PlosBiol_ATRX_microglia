#!/usr/bin/env Rscript
# =========================================================================
# Title: Identification and Classification of Upregulated LTRs as Full-Length or Solo Using RetroTector
# Script: ltr_classification_solo_full_length.R
# Author: Alireza Ghahramani
# Contact: aghahram@uwo.ca
#
# Description: This pipeline filters significant upregulated LTRs, extracts ±15kb
# genomic FASTA, runs RetroTector, parses and merges output for downstream analysis.
################################################################################
# Load required packages
library(dplyr)
library(stringr)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
library(purrr)
library(tidyr)

################################################################################
# Step 1: Load and Filter LTRs for Upregulated Elements
################################################################################

# Read in CSV containing significant LTR info (adjust path as needed)
ltr_info <- read.csv("path/to/sig_ltrs.csv")

# Filter for upregulated LTRs (log2FoldChange > 0)
up_ltrs <- ltr_info %>%
  filter(log2FoldChange > 0)

# Helper: Convert CSV RepName to FASTA-compatible header
csv_to_fasta_header <- function(name) {
  parts <- str_split(name, "#")[[1]]
  paste(parts, collapse = "_")
}

# Generate a FASTA-formatted header for each upregulated LTR
up_ltrs <- up_ltrs %>%
  mutate(fasta_header = sapply(name, csv_to_fasta_header))

################################################################################
# Step 2: Extract ±15kb Genomic Sequence for Each LTR and Write FASTA Files
################################################################################

# Set extension length (e.g., 15kb)
EXT <- 15000

# Ensure coordinates stay within chromosome bounds
up_ltrs$ext_start <- pmax(up_ltrs$Start - EXT, 1)
up_ltrs$ext_end   <- up_ltrs$End + EXT

# Create a GRanges object for extraction
gr <- GRanges(
  seqnames = up_ltrs$Chr,
  ranges   = IRanges(start = up_ltrs$ext_start, end = up_ltrs$ext_end),
  strand   = if ("strand" %in% names(up_ltrs)) up_ltrs$Strand else "+"
)
names(gr) <- up_ltrs$fasta_header

# Extract sequences from the mm10 genome
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, gr)
names(seqs) <- names(gr)

# Write all sequences to a single multi-FASTA file (for batch analysis)
writeXStringSet(
  seqs,
  "sig_upregulated_LTRs_15kb_extended.fa"
)

# Write each sequence to its own individual FASTA file (for RetroTector input)
fasta_dir <- "LTRs_individual_fasta"
dir.create(fasta_dir, showWarnings = FALSE)
for (i in seq_along(seqs)) {
  header <- names(seqs)[i]
  seq <- seqs[i]
  file_name <- file.path(fasta_dir, paste0(header, ".fa"))
  writeXStringSet(seq, filepath = file_name)
}

################################################################################
# Step 3: [Run RetroTector outside of R as documented]
# Use SweepDNA/SweepScripts on the generated FASTA files.
################################################################################

# After running RetroTector, collect all 001SelectedChains.txt files from subfolders.

################################################################################
# Step 4: Parse RetroTector Output and Merge with LTR Metadata
################################################################################

# List all RetroTector output chain files (edit path as needed)
chain_files <- list.files(
  "path/tp/retrotector/ReTe1.0.1/Workplace",
  pattern = "001SelectedChains.txt$",
  recursive = TRUE,
  full.names = TRUE
)

# Helper: Parse info from folder name (handles underscores in repName robustly)
parse_folder_name <- function(folder_name) {
  parts <- unlist(strsplit(folder_name, "_"))
  n <- length(parts)
  end   <- parts[n]
  start <- parts[n-1]
  chr   <- parts[n-2]
  repClass  <- parts[n-3]
  repFamily <- parts[n-4]
  repName <- if (n > 5) paste(parts[1:(n-4)], collapse = "_") else parts[1]
  data.frame(
    folder = folder_name,
    repName = repName,
    repFamily = repFamily,
    repClass = repClass,
    chr = chr,
    start = as.numeric(start),
    end = as.numeric(end),
    stringsAsFactors = FALSE
  )
}

# Helper: Parse a single chain file for provirus classification
parse_chain_file <- function(file_path) {
  lines <- readLines(file_path)
  folder_name <- basename(dirname(file_path))
  info <- parse_folder_name(folder_name)
  
  chain_lines <- grep("^Chain", lines, value = TRUE)
  classification <- "Solo_LTR"
  domains <- NA
  
  if (length(chain_lines) > 0) {
    subgene_lines <- grep("^SubGene ", lines, value = TRUE)
    if (length(subgene_lines) > 0) {
      subgenes <- unique(str_extract(subgene_lines, "(?<=SubGene )\\w+"))
      domains <- paste(subgenes, collapse = ";")
      # Full-length: all main retroviral domains + LTR
      if (all(c("MA", "CA", "NC", "RT", "IN") %in% subgenes) && any(grepl("LTR", subgenes))) {
        classification <- "Full_length_provirus"
      } else {
        classification <- "Partial_provirus"
      }
    } else {
      classification <- "Partial_provirus"
      domains <- NA
    }
  }
  bind_cols(info, tibble(
    classification = classification,
    domains = domains
  ))
}

# Parse all chain files
chain_results <- map_dfr(chain_files, parse_chain_file)

# Save chain parsing results (optional, for record-keeping)
write.csv(chain_results, "retrotector_chain_classification.csv", row.names = FALSE)


