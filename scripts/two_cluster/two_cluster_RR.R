# two_cluster_RR.R
# RR / nRR for the four-deme two-cluster simulation, at the deme level (4x4) and
# the cluster level (C1 / C2). Phase boundaries are fixed and known, so RR is
# written per phase as well as on a sliding window.

options(tidyverse.quiet = TRUE)
library(tidyverse, quietly = TRUE)
library(argparse, quietly = TRUE)

source("scripts/common/calc_rr.R")

collect_args <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--seed", type = "integer", default = 66, help = "Simulation Seed")
  parser$add_argument("--scenario", type = "character", default = "homogeneous",
                      help = "Scenario label carried into the output")
  parser$add_argument("--intervention", action = "store_true", default = FALSE,
                      help = "Replicate was run with the intervention (sim --vaccine)")
  parser$parse_args()
}

args <- collect_args()
sim_seed <- args$seed
scenario <- args$scenario
intervention <- args$intervention
# Output prefix, e.g. "homogeneous" or "homogeneous_intervention"
arm <- if (intervention) paste0(scenario, "_intervention") else scenario

DATA_DIR <- "./results/two_cluster"
PHASE_STARTS <- c(0, 100, 150, 200, 250)
PHASE_NAMES <- c("burn-in", "shutdown", "initial", "moderate", "full")

# Accumulators are scenario-tagged so sequential sweeps of different scenarios
# do not overwrite each other's results.
acc <- function(name) file.path(DATA_DIR, paste0(arm, "_", name))

append_tsv <- function(df, path) {
  if (file.exists(path)) {
    write_tsv(df, file = path, append = TRUE)
  } else {
    write_tsv(df, file = path)
  }
}

meta <- read_tsv(file.path(DATA_DIR, "meta.tsv"), show_col_types = FALSE)
pairs <- read_tsv(file.path(DATA_DIR, "idseq.tsv"), show_col_types = FALSE)
colnames(pairs) <- c("X", "Y", "Dist")
pairs <- pairs %>%
  filter(X %in% meta$hosts.ID) %>%
  filter(Y %in% meta$hosts.ID)

t_max <- max(meta$inf.time)
cluster_meta <- meta %>% mutate(inf.in = cluster)

# ---------------------------------------------------------------------------
# Overall RR
# ---------------------------------------------------------------------------
append_tsv(
  calc_rr(meta, pairs) %>% mutate(seed = sim_seed, scenario = scenario, intervention = intervention),
  acc("rr.tsv")
)
append_tsv(
  calc_rr(cluster_meta, pairs) %>% mutate(seed = sim_seed, scenario = scenario, intervention = intervention),
  acc("rr_cluster.tsv")
)

# ---------------------------------------------------------------------------
# Phase RR
# ---------------------------------------------------------------------------
phase_ends <- c(PHASE_STARTS[-1], t_max)
phase_rr <- function(m) {
  map_dfr(seq_along(PHASE_NAMES), function(i) {
    if (PHASE_STARTS[i] >= t_max) {
      return(NULL)
    }
    # time_point is set so normalized_rr() has a key to join the diagonal on;
    # pooling every pair in the phase makes nRR far more stable than taking a
    # median over narrow sliding windows, which collapses to zero whenever more
    # than half the windows happen to be empty.
    calc_rr(m, pairs,
            TIME_WINDOW = TRUE,
            TIME_START = PHASE_STARTS[i],
            TIME_END = min(phase_ends[i], t_max),
            time_point = i) %>%
      mutate(phase = PHASE_NAMES[i], phase_index = i)
  }) %>%
    normalized_rr() %>%
    mutate(seed = sim_seed, scenario = scenario, intervention = intervention)
}

append_tsv(phase_rr(meta), acc("rr_phase.tsv"))
append_tsv(phase_rr(cluster_meta), acc("rr_cluster_phase.tsv"))

# ---------------------------------------------------------------------------
# Sliding window RR / nRR
# ---------------------------------------------------------------------------
# Right-sided (trailing) window: tpoint t covers [t - NRR_WINDOW, t], inclusive,
# so each value uses only infections up to t and a change at time t0 cannot show
# before t0. Until 2026-11-12 the window was centred ([t - 20, t + 20]); the
# accumulators of the 2026-11-11 sweep were converted by relabelling t -> t + 20,
# which gives identical windows, rather than rerun.
NRR_WINDOW <- 40
NRR_TIME_SPACING <- 5
nrr_tpoints <- seq(NRR_WINDOW, t_max, NRR_TIME_SPACING)

sliding_rr <- function(m) {
  map_dfr(nrr_tpoints, ~ calc_rr(m, pairs,
                                 TIME_WINDOW = TRUE,
                                 TIME_START = .x - NRR_WINDOW,
                                 TIME_END = .x,
                                 time_point = .x)) %>%
    normalized_rr() %>%
    mutate(seed = sim_seed, scenario = scenario, intervention = intervention)
}

append_tsv(sliding_rr(meta), acc("rr_sliding.tsv"))
append_tsv(sliding_rr(cluster_meta), acc("rr_cluster_sliding.tsv"))

# Carry the realised schedule through, tagged by seed
schedule_fn <- file.path(DATA_DIR, "schedule.tsv")
if (file.exists(schedule_fn)) {
  append_tsv(
    read_tsv(schedule_fn, show_col_types = FALSE),
    acc("schedules.tsv")
  )
}
