# two_cluster_plot.R
# Summary plots across the seeds of a two-cluster batch.
#
# The comparison of interest is intra-cluster (A-B, A-C, B-C, X-Y, X-Z, Y-Z)
# against inter-cluster (A-X) recovery. Under the homogeneous schedule both
# should recover together; under the heterogeneous schedule A-X should lag.

options(tidyverse.quiet = TRUE)
library(tidyverse, quietly = TRUE)
library(viridis, quietly = TRUE)
library(argparse, quietly = TRUE)

collect_args <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--scenario", type = "character", default = "homogeneous",
                      help = "Scenario label used in output filenames")
  parser$parse_args()
}

args <- collect_args()
scenario <- args$scenario

DATA_DIR <- "results/two_cluster"

# Display name for plot titles; `scenario` itself stays as-is for filenames
scenario_label <- function(s) if (s == "none") "No Shutdown" else s
out_file <- function(name) file.path(DATA_DIR, paste0(scenario, "_", name))

DEME_ORDER <- c("A", "B", "C", "X", "Y", "Z")
PHASE_STARTS <- c(0, 100, 150, 200, 250)
PHASE_MARKS <- PHASE_STARTS[-1]
PHASE_NAMES <- c("burn-in", "shutdown", "initial", "moderate", "full")

# Label each unordered deme pair by whether it crosses the A-X bridge
pair_kind <- function(x, y) {
  c1 <- c("A", "B", "C")
  same <- (x %in% c1) == (y %in% c1)
  if_else(x == y, "self", if_else(same, "Intra-cluster", "Inter-cluster"))
}

add_marks <- function(p) {
  p + geom_vline(
    xintercept = PHASE_MARKS, linetype = "dashed",
    color = "grey30", linewidth = 0.4
  )
}

# ---------------------------------------------------------------------------
# Shared axis limits across scenarios
# ---------------------------------------------------------------------------
# Each scenario is plotted by a separate invocation of this script, so without
# a common scale the three sets of figures cannot be compared by eye: an arm
# whose nRR collapses would simply get a rescaled axis and look unchanged.
# Limits are pooled over whatever scenario files are present, so a full sweep
# yields identical axes for all three arms.
ALL_SCENARIOS <- c("none", "homogeneous", "heterogeneous")

pooled_column <- function(suffix, column) {
  vals <- map(ALL_SCENARIOS, function(s) {
    f <- file.path(DATA_DIR, paste0(s, "_", suffix))
    if (!file.exists(f)) {
      return(numeric(0))
    }
    d <- read_tsv(f, show_col_types = FALSE)
    if (!column %in% names(d)) {
      return(numeric(0))
    }
    d[[column]]
  })
  unlist(vals, use.names = FALSE)
}

# Positive values only: these axes are log-scaled, and nRR is exactly 0 whenever
# a window holds no identical pairs. The lower bound is a low quantile rather
# than the minimum so one near-zero value cannot stretch the axis flat.
log_limits <- function(x, lower_q = 0.01, pad = 1.6) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 2) {
    return(NULL)
  }
  c(unname(quantile(x, lower_q)) / pad, max(x) * pad)
}

NRR_SLIDING_LIM <- log_limits(pooled_column("rr_sliding.tsv", "nRR"))
RECOVERY_LIM <- NULL
RR_CLUSTER_LIM <- log_limits(pooled_column("rr_cluster_phase.tsv", "RR"))
RR_FILL_LIM <- log_limits(pooled_column("rr.tsv", "RR"), lower_q = 0.02, pad = 1.0)

fix_y <- function(p, lim) {
  if (is.null(lim)) p else p + coord_cartesian(ylim = lim)
}

# ---------------------------------------------------------------------------
# Sliding-window nRR: intra vs inter
# ---------------------------------------------------------------------------
rr_sliding <- read_tsv(out_file("rr_sliding.tsv"), show_col_types = FALSE) %>%
  mutate(kind = pair_kind(deme_X, deme_Y)) %>%
  filter(kind != "self")

p <- ggplot(rr_sliding, aes(x = tpoint, y = nRR, color = kind))
p <- add_marks(p) +
  geom_smooth(se = TRUE, method = "loess", formula = y ~ x, span = 0.3) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  scale_y_log10() +
  scale_color_manual(values = c("Intra-cluster" = "#38c", "Inter-cluster" = "#e94")) +
  labs(
    x = "Simulation Time", y = "nRR (log scale)", color = NULL,
    title = paste0("Intra- vs inter-cluster nRR (scenario: ", scenario_label(scenario), ")"),
    subtitle = "Phases: burn-in | shutdown | initial | moderate | full"
  ) +
  theme_bw()
p <- fix_y(p, NRR_SLIDING_LIM)
ggsave(out_file("nrr_intra_vs_inter.jpg"), plot = p, width = 10, height = 6, dpi = 192)

# Per-pair detail
p <- ggplot(
  rr_sliding %>% mutate(pair = paste0(deme_X, "-", deme_Y)),
  aes(x = tpoint, y = nRR, color = kind, group = pair)
)
p <- add_marks(p) +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x, span = 0.4, linewidth = 0.6) +
  scale_y_log10() +
  scale_color_manual(values = c("Intra-cluster" = "#38c", "Inter-cluster" = "#e94")) +
  labs(
    x = "Simulation Time", y = "nRR (log scale)", color = NULL,
    title = paste0("nRR by deme pair (scenario: ", scenario_label(scenario), ")")
  ) +
  theme_bw()
p <- fix_y(p, NRR_SLIDING_LIM)
ggsave(out_file("nrr_by_pair.jpg"), plot = p, width = 10, height = 6, dpi = 192)

# ---------------------------------------------------------------------------
# Recovery relative to each pair's own burn-in baseline
# ---------------------------------------------------------------------------
# This is the quantity the intra-vs-inter argument turns on: how far each pair
# type has returned toward its pre-shutdown level, phase by phase.
# Uses phase-pooled nRR rather than sliding windows: a narrow window that
# happens to contain no identical pairs contributes a hard zero, and the median
# of those zeros swamps the real signal for the rarer inter-cluster pairs.
rr_phase <- read_tsv(out_file("rr_phase.tsv"), show_col_types = FALSE) %>%
  mutate(kind = pair_kind(deme_X, deme_Y)) %>%
  filter(kind != "self")

baseline <- rr_phase %>%
  filter(phase == PHASE_NAMES[1]) %>%
  select(seed, deme_X, deme_Y, base = nRR)

recovery <- rr_phase %>%
  left_join(baseline, by = c("seed", "deme_X", "deme_Y")) %>%
  mutate(ratio = nRR / base) %>%
  filter(is.finite(ratio)) %>%
  group_by(seed, kind, phase) %>%
  summarize(ratio = median(ratio, na.rm = TRUE), .groups = "drop") %>%
  mutate(phase = factor(phase, levels = PHASE_NAMES))

recovery_lim <- log_limits(
  c(recovery$ratio, pooled_column("recovery_by_phase.tsv", "ratio"))
)

p <- ggplot(recovery, aes(x = phase, y = ratio, fill = kind)) +
  geom_boxplot(alpha = 0.75, outlier.size = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_log10() +
  scale_fill_manual(values = c("Intra-cluster" = "#38c", "Inter-cluster" = "#e94")) +
  labs(
    x = "Phase", y = "nRR relative to burn-in baseline (log scale)", fill = NULL,
    title = paste0("Recovery by phase (scenario: ", scenario_label(scenario), ")"),
    subtitle = "1.0 = back to pre-shutdown level"
  ) +
  theme_bw()
p <- fix_y(p, recovery_lim)
ggsave(out_file("recovery_by_phase.jpg"), plot = p, width = 10, height = 6, dpi = 192)

write_tsv(recovery, out_file("recovery_by_phase.tsv"))

# ---------------------------------------------------------------------------
# Phase RR at the cluster level
# ---------------------------------------------------------------------------
rr_cluster_phase <- read_tsv(
  out_file("rr_cluster_phase.tsv"),
  show_col_types = FALSE
) %>%
  mutate(phase = factor(phase, levels = PHASE_NAMES))

p <- rr_cluster_phase %>%
  mutate(kind = if_else(deme_X == deme_Y, "Within cluster", "Between clusters")) %>%
  ggplot(aes(x = phase, y = RR, fill = kind)) +
  geom_boxplot(alpha = 0.75, outlier.size = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_log10() +
  scale_fill_manual(values = c("Within cluster" = "#38c", "Between clusters" = "#e94")) +
  labs(
    x = "Phase", y = "RR (log scale)", fill = NULL,
    title = paste0("Cluster-level RR by phase (scenario: ", scenario_label(scenario), ")")
  ) +
  theme_bw()
p <- fix_y(p, RR_CLUSTER_LIM)
ggsave(out_file("rr_cluster_phase.jpg"), plot = p, width = 10, height = 6, dpi = 192)

# ---------------------------------------------------------------------------
# Deme-level RR heatmap
# ---------------------------------------------------------------------------
read_tsv(out_file("rr.tsv"), show_col_types = FALSE) %>%
  mutate(
    deme_X = factor(deme_X, levels = DEME_ORDER),
    deme_Y = factor(deme_Y, levels = DEME_ORDER)
  ) %>%
  group_by(deme_X, deme_Y) %>%
  summarize(RR = median(RR, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = deme_X, y = deme_Y, fill = RR)) +
  geom_tile() +
  scale_fill_gradient2(
    high = "#e94", low = "#38c", mid = "white",
    transform = "log10", midpoint = 1, limits = RR_FILL_LIM
  ) +
  labs(
    x = "Deme", y = "Deme",
    title = paste0("Median RR across seeds (scenario: ", scenario_label(scenario), ")")
  ) +
  theme_minimal()
ggsave(out_file("rr_heatmap.jpg"), width = 7, height = 6, dpi = 192)
