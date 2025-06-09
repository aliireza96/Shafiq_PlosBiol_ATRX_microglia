#!/usr/bin/env Rscript
# =========================================================================
# Script: ltr_classification_solo_full_length.R
# Author: Alireza Ghahramani
# Contact: aghahram@uwo.ca
#
# Purpose:
#   - Classify significant LTR elements into "Solo LTR" and "Full-Length"
#     based on naming conventions (e.g., "-int" in name = full-length)
#
#
# Inputs:
#   - A table of significant LTRs (e.g., output of DESeq2)
#     Required column: "TE_name" (character, e.g., "ERVB2_1-I_MM-int:ERVK:LTR")
#
# Outputs:
#   - Table with additional column: "LTR_type" ("Full-Length" or "Solo LTR")
#   - Counts of each category
#
# =========================================================================

# -------- USER PARAMETERS (EDIT THESE) --------
input_file <- "significant_LTRs.csv"      # Input: table with at least TE_name column
te_column  <- "TE_name"                   # Name of column with LTR names
output_file <- "LTR_solo_full_length_classification.csv"

# -------- 1. LOAD DATA --------
ltr_tbl <- read.csv(input_file, stringsAsFactors=FALSE)

# -------- 2. CLASSIFY LTRs --------
# Rule: If TE name contains "-int" => Full-Length, else => Solo LTR
ltr_tbl$LTR_type <- ifelse(grepl("-int", ltr_tbl[[te_column]], ignore.case = TRUE),
                           "Full-Length", "Solo LTR")

# -------- 3. EXPORT CLASSIFIED TABLE --------
write.csv(ltr_tbl, output_file, row.names=FALSE)

# -------- 4. SUMMARY COUNTS --------
table_counts <- table(ltr_tbl$LTR_type)
cat("LTR classification summary:\n")
print(table_counts)

# -------- 5. OPTIONAL: SAVE SUMMARY --------
write.csv(as.data.frame(table_counts), "LTR_type_counts_summary.csv", row.names=FALSE)

cat("\nClassification complete. Output: ", output_file, "\n")
