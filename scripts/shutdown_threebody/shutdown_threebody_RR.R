# shutdown_threebody_RR.R
# Calculates RR matrices for shutdown simulation, including phase-specific and sliding window time-based RR and nRR

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

shutdown_meta <- read_tsv("./results/shutdown_threebody/meta.tsv", show_col_types = FALSE)
shutdown_pairs <- read_tsv("./results/shutdown_threebody/idseq.tsv")
colnames(shutdown_pairs) <- c("X", "Y", "Dist")
shutdown_pairs <- filter(shutdown_pairs, X %in% shutdown_meta$hosts.ID) %>%
  filter(Y %in% shutdown_meta$hosts.ID)

# Define shutdown phases
t_shutdown_start <- 40
t_shutdown_end <- 120
t_max <- max(shutdown_meta$inf.time)

# Calculate RR for each phase with 10-day overlap
rr_pre <- calc_rr(shutdown_meta, shutdown_pairs, TIME_WINDOW = TRUE,
                  TIME_START = 0,
                  TIME_END = t_shutdown_start + 10)
rr_pre <- mutate(rr_pre, phase = "Pre-Lockdown")

rr_during <- calc_rr(shutdown_meta, shutdown_pairs, TIME_WINDOW = TRUE,
                     TIME_START = t_shutdown_start - 10,
                     TIME_END = t_shutdown_end + 10)
rr_during <- mutate(rr_during, phase = "During Lockdown")

rr_post <- calc_rr(shutdown_meta, shutdown_pairs, TIME_WINDOW = TRUE,
                   TIME_START = t_shutdown_end - 10,
                   TIME_END = t_max)
rr_post <- mutate(rr_post, phase = "Post-Lockdown")

# Combine and add seed
rr_phases <- bind_rows(rr_pre, rr_during, rr_post) %>%
  mutate(seed = sim_seed)


# Write phase RR to TSV (append if file exists)
rr_phase_fn <- "./results/shutdown_threebody/rr.tsv"
if (file.exists(rr_phase_fn)) {
  write_tsv(rr_phases, file = rr_phase_fn, append = TRUE)
} else {
  write_tsv(rr_phases, file = rr_phase_fn)
}

# Sliding window RR and normalized RR
nrr_bin_width <- 20
nrr_time_spacing <- 5
nrr_max_time <- 200
nrr_tpoints <- seq(nrr_bin_width, nrr_max_time, nrr_time_spacing)

rr_sliding <- map_dfr(nrr_tpoints, ~ calc_rr(shutdown_meta, shutdown_pairs,
                                             TIME_WINDOW = TRUE,
                                             TIME_START = .x - nrr_bin_width,
                                             TIME_END = .x + nrr_bin_width,
                                             time_point = .x)) %>%
  normalized_rr() %>%
  mutate(seed = sim_seed)


# Write sliding window RR/nRR to TSV (append if file exists)
rr_sliding_fn <- "./results/shutdown_threebody/rr_sliding.tsv"
if (file.exists(rr_sliding_fn)) {
  write_tsv(rr_sliding, file = rr_sliding_fn, append = TRUE)
} else {
  write_tsv(rr_sliding, file = rr_sliding_fn)
}
