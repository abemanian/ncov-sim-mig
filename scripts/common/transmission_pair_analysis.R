#!/usr/bin/env Rscript
# Transmission Pair Sequence Identity Analysis
# Calculates the proportion of parent-child transmission pairs
# that have identical sequences after phastSim evolution

options(tidyverse.quiet = TRUE)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))

collect_args <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--seed",
    type = "integer",
    default = 66,
    help = "Simulation seed")
  parser$add_argument("--meta",
    type = "character",
    required = TRUE,
    help = "Path to metadata TSV file")
  parser$add_argument("--pairs",
    type = "character",
    required = TRUE,
    help = "Path to pairsnp output (idseq.tsv)")
  parser$add_argument("--output",
    type = "character",
    required = TRUE,
    help = "Path to output results file")
  parser$parse_args()
}

args <- collect_args()
sim_seed <- args$seed
meta_file <- args$meta
pairs_file <- args$pairs
output_file <- args$output

# Load metadata (has parent-child relationships via inf.by column)
meta <- read_tsv(meta_file, show_col_types = FALSE)

# Load pairsnp output (identical sequence pairs, distance = 0)
seq_pairs <- read_tsv(pairs_file,
                      col_names = c("X", "Y", "Dist"),
                      show_col_types = FALSE)

# Extract parent-child transmission pairs from metadata
transmission_pairs <- meta %>%
  filter(!grepl("NA", inf.by)) %>%  # Exclude root node (has NA parent)
  select(child = hosts.ID, parent = inf.by) %>%
  mutate(is_transmission = TRUE)

# Symmetrize sequence pairs (both X->Y and Y->X represent same pair)
seq_pairs_sym <- bind_rows(
  seq_pairs %>% select(X, Y),
  seq_pairs %>% select(X = Y, Y = X)
) %>%
  mutate(is_identical = TRUE) %>%
  distinct()

# Check which transmission pairs have identical sequences
results <- transmission_pairs %>%
  left_join(seq_pairs_sym, by = c("parent" = "X", "child" = "Y")) %>%
  mutate(identical_seq = ifelse(!is.na(is_identical),is_identical,FALSE)) %>%
  select(parent, child, identical_seq)

# Calculate summary statistics
summary_stats <- results %>%
  summarize(
    seed = sim_seed,
    total_transmission_pairs = n(),
    identical_sequences = sum(identical_seq),
    proportion = mean(identical_seq)
  )

# Print to console
cat(sprintf("Seed %d: %.1f%% (%d/%d) of transmission pairs have identical sequences\n",
            sim_seed,
            summary_stats$proportion * 100,
            summary_stats$identical_sequences,
            summary_stats$total_transmission_pairs))

# Append to output file
if (file.exists(output_file)) {
  write_tsv(summary_stats, file = output_file, append = TRUE)
} else {
  write_tsv(summary_stats, file = output_file)
}
