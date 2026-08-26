# carrying_capacity_RR.R
# Calculates RR matrices for carrying capacity simulation

options(tidyverse.quiet = TRUE)
library(tidyverse, quietly = TRUE)
library(argparse, quietly = TRUE)

source("scripts/common/calc_rr.R")
source("scripts/common/plot_pair_time.R")

collect_args <- function() {
  parser <- ArgumentParser()
  parser$add_argument('--seed', type = 'integer', default = 66, help = 'Simulation Seed')
  parser$parse_args()
}

args <- collect_args()
sim_seed <- args$seed

# Read and format data
carrying_meta <- read_tsv("./results/carrying_capacity/meta.tsv", show_col_types = FALSE)
carrying_pairs <- read_tsv("./results/carrying_capacity/idseq.tsv")
colnames(carrying_pairs) <- c("X", "Y", "Dist")
carrying_pairs <- filter(carrying_pairs, X %in% carrying_meta$hosts.ID) %>%
  filter(Y %in% carrying_meta$hosts.ID)

# Calculate overall RR
rr_total <- calc_rr(carrying_meta, carrying_pairs) %>%
  mutate(seed = sim_seed)

# Write RR to TSV (append if exists)
rr_fn <- "./results/carrying_capacity/rr.tsv"
if (file.exists(rr_fn)) {
  write_tsv(rr_total, file = rr_fn, append = TRUE)
} else {
  write_tsv(rr_total, file = rr_fn)
}

# Sliding window RR and normalized RR
nrr_bin_width <- 20
nrr_time_spacing <- 5
nrr_max_time <- max(carrying_meta$inf.time)
nrr_tpoints <- seq(nrr_bin_width, nrr_max_time, nrr_time_spacing)

rr_sliding <- map_dfr(nrr_tpoints, ~ calc_rr(carrying_meta, carrying_pairs,
                                             TIME_WINDOW = TRUE,
                                             TIME_START = .x - nrr_bin_width,
                                             TIME_END = .x + nrr_bin_width,
                                             time_point = .x)) %>%
  normalized_rr() %>%
  mutate(seed = sim_seed)

# Write sliding window RR/nRR to TSV (append if exists)
rr_sliding_fn <- "./results/carrying_capacity/rr_sliding.tsv"
if (file.exists(rr_sliding_fn)) {
  write_tsv(rr_sliding, file = rr_sliding_fn, append = TRUE)
} else {
  write_tsv(rr_sliding, file = rr_sliding_fn)
}