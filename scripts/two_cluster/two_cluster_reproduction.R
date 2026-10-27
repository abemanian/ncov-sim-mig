# two_cluster_reproduction.R
# Effective reproduction number over time for one two-cluster replicate, read
# straight off the transmission chain in meta.tsv: R(t) is the mean number of
# onward infections caused by hosts infected in a window centred on t.
#
# Works on the full simulated meta, not the subsampled FASTA, so it does not
# depend on phastSim / pairsnp and is unaffected by SUBSAMPLE_N.

options(tidyverse.quiet = TRUE)
library(tidyverse, quietly = TRUE)
library(argparse, quietly = TRUE)

collect_args <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--seed", type = "integer", default = 66, help = "Simulation Seed")
  parser$add_argument("--scenario", type = "character", default = "homogeneous",
                      help = "Scenario label carried into the output")
  parser$parse_args()
}

args <- collect_args()
sim_seed <- args$seed
scenario <- args$scenario

DATA_DIR <- "./results/two_cluster"

# Total window length in time steps; must be odd so the window centres on t
RT_WINDOW <- 7

acc <- function(name) file.path(DATA_DIR, paste0(scenario, "_", name))

append_tsv <- function(df, path) {
  if (file.exists(path)) {
    write_tsv(df, file = path, append = TRUE)
  } else {
    write_tsv(df, file = path)
  }
}

meta <- read_tsv(file.path(DATA_DIR, "meta.tsv"), show_col_types = FALSE)

# Offspring count for every host, including the ones that infected nobody.
# Counting inf.by alone only sees hosts with at least one infectee, which drops
# the zeros and biases R(t) upward.
offspring <- meta %>%
  count(inf.by, name = "n") %>%
  rename(hosts.ID = inf.by)

hosts <- meta %>%
  select(hosts.ID, inf.time) %>%
  left_join(offspring, by = "hosts.ID") %>%
  mutate(n = coalesce(n, 0L))

rt_window <- function(t, df, window = RT_WINDOW) {
  if (window %% 2 == 0) {
    stop("Window must be odd!")
  }
  w_radius <- (window - 1) / 2
  in_window <- df$inf.time >= t - w_radius & df$inf.time <= t + w_radius
  # Empty windows give NaN rather than 0: no hosts means no estimate, not R = 0
  tibble(tpoint = t, Rt = mean(df$n[in_window]), n_hosts = sum(in_window))
}

# Hosts infected near the end of the run have their onward transmission cut off
# by MAX_TIME, so the tail of this series is biased low. It is written out in
# full; the plotting scripts crop it.
t_max <- max(meta$inf.time)
rt <- map_dfr(0:t_max, rt_window, df = hosts) %>%
  mutate(window = RT_WINDOW, seed = sim_seed, scenario = scenario)

append_tsv(rt, acc("rt.tsv"))
