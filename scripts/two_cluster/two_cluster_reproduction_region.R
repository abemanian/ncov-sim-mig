# two_cluster_reproduction_region.R
# Effective reproduction number over time, stratified by region, for one
# two-cluster replicate, read straight off the transmission chain in meta.tsv:
# R(t) is the mean number of onward infections caused by hosts infected in a
# window centred on t, zeros included, computed separately for each region.
# Uses the full meta, not the subsample, so it is unaffected by SUBSAMPLE_N.
#
# A host is assigned to the region it was infected in (meta.tsv `cluster`, from
# inf.in), and all of its offspring count towards that region even if they were
# infected after it moved. Movement is mostly within-cluster (only the A-X edge
# crosses), so this attribution rarely differs from the infector's location at
# transmission.
#
# Appends to <arm>_rt_region.tsv. (The pooled two_cluster_reproduction.R was
# removed on 2026-11-12; a pooled series is the n_hosts-weighted mean of the two
# regions' rows.)

options(tidyverse.quiet = TRUE)
library(tidyverse, quietly = TRUE)
library(argparse, quietly = TRUE)

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
arm <- if (intervention) paste0(scenario, "_intervention") else scenario

DATA_DIR <- "./results/two_cluster"

# Total window length in days; must be odd so the window centres on t
RT_WINDOW <- 7
stopifnot(RT_WINDOW %% 2 == 1)

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
  select(hosts.ID, inf.time, region = cluster) %>%
  left_join(offspring, by = "hosts.ID") %>%
  mutate(n = coalesce(n, 0L))

w_radius <- (RT_WINDOW - 1) / 2
t_max <- max(meta$inf.time)

# Hosts infected per (region, day) and their total offspring, then a
# centred rolling sum over the window. Empty windows give NaN, as in the pooled
# series: no hosts means no estimate, not R = 0.
rt <- hosts %>%
  count(region, inf.time, wt = n, name = "offspring") %>%
  full_join(count(hosts, region, inf.time, name = "n_hosts"), by = c("region", "inf.time")) %>%
  complete(region, inf.time = 0:t_max, fill = list(offspring = 0, n_hosts = 0)) %>%
  arrange(region, inf.time) %>%
  group_by(region) %>%
  mutate(
    win_off = map_dbl(inf.time, ~ sum(offspring[abs(inf.time - .x) <= w_radius])),
    win_n = map_dbl(inf.time, ~ sum(n_hosts[abs(inf.time - .x) <= w_radius]))
  ) %>%
  ungroup() %>%
  transmute(
    region, tpoint = inf.time, Rt = win_off / win_n, n_hosts = win_n,
    window = RT_WINDOW, seed = sim_seed, scenario = scenario,
    intervention = intervention
  )

append_tsv(rt, file.path(DATA_DIR, paste0(arm, "_rt_region.tsv")))
