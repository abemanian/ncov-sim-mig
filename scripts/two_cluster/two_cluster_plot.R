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
SCENARIO_LABELS <- c(
  none = "No Shutdown",
  homogeneous = "Uniform Recovery",
  heterogeneous = "Staggered Recovery"
)
scenario_label <- function(s) unname(coalesce(SCENARIO_LABELS[s], s))
out_file <- function(name) file.path(DATA_DIR, paste0(scenario, "_", name))

DEME_ORDER <- c("A", "B", "C", "X", "Y", "Z")
PHASE_STARTS <- c(0, 100, 150, 200, 250)
PHASE_MARKS <- PHASE_STARTS[-1]
PHASE_NAMES <- c("burn-in", "shutdown", "initial", "moderate", "full")

# Time range drawn in the sliding-window figures. The left edge is the baseline
# window's centre -- everything earlier is the seeding period, which is dense in
# identical sequences for reasons unrelated to the intervention. The right edge
# stops short of MAX_TIME (350): the simulation ends abruptly, so hosts infected
# near the end have their onward transmission truncated, which depresses nRR in
# the final windows as an artefact of the stop rather than of connectivity.
# Filtering is applied per plot, never to the baseline, which stays independent
# of what the figures happen to display.
PLOT_TIME_MIN <- 90
PLOT_TIME_MAX <- 300

# Label each unordered deme pair by whether it crosses the A-X bridge
pair_kind <- function(x, y) {
  c1 <- c("A", "B", "C")
  same <- (x %in% c1) == (y %in% c1)
  if_else(x == y, "self", if_else(same, "Intra-cluster", "Inter-cluster"))
}

add_marks <- function(p) {
  marks <- PHASE_MARKS[
    PHASE_MARKS >= PLOT_TIME_MIN & PHASE_MARKS <= PLOT_TIME_MAX
  ]
  p + geom_vline(
    xintercept = marks, linetype = "dashed",
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

# Recovery ratios are exactly 0 whenever a phase holds no identical pairs for a
# pair of demes -- a real result, not missing data, and the one a log axis cannot
# draw. A pseudo-log axis is log-like above sigma and linear through 0, so those
# zeros land on the axis instead of being silently dropped. sigma sits just below
# the smallest non-zero ratio observed (~0.011) so it does not compress real data.
PSEUDO_SIGMA <- 0.005
RECOVERY_BREAKS <- c(0, 0.01, 0.03, 0.1, 0.3, 1, 3)

recovery_y_scale <- function() {
  scale_y_continuous(
    transform = scales::transform_pseudo_log(sigma = PSEUDO_SIGMA, base = 10),
    breaks = RECOVERY_BREAKS,
    labels = format(RECOVERY_BREAKS, trim = TRUE, drop0trailing = TRUE)
  )
}

# Pooled over the plotted window only; including the censored tail would stretch
# the shared axis to accommodate values the figures never show.
pooled_sliding_nrr <- function() {
  vals <- map(ALL_SCENARIOS, function(s) {
    f <- file.path(DATA_DIR, paste0(s, "_rr_sliding.tsv"))
    if (!file.exists(f)) {
      return(numeric(0))
    }
    read_tsv(f, show_col_types = FALSE) %>%
      filter(between(tpoint, PLOT_TIME_MIN, PLOT_TIME_MAX)) %>%
      pull(nRR)
  })
  unlist(vals, use.names = FALSE)
}

NRR_SLIDING_LIM <- log_limits(pooled_sliding_nrr())
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

rr_sliding_win <- rr_sliding %>%
  filter(between(tpoint, PLOT_TIME_MIN, PLOT_TIME_MAX))

p <- ggplot(rr_sliding_win, aes(x = tpoint, y = nRR, color = kind))
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
  rr_sliding_win %>% mutate(pair = paste0(deme_X, "-", deme_Y)),
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

# Baseline is the sliding window centred on BASELINE_TPOINT, not the burn-in
# phase. The ±20 half-width is the point: inter-cluster identical pairs are rare,
# and any narrower window sends a large fraction of them to zero, which is both
# unstable and biased (the pairs that vanish are the least-mixed ones). A window
# bracketing the intervention start also avoids the early burn-in, when the
# outbreak is still confined to the seed deme.
BASELINE_TPOINT <- 90

baseline <- rr_sliding %>%
  filter(tpoint == BASELINE_TPOINT) %>%
  select(seed, deme_X, deme_Y, base = nRR)

# A pair with no identical sequences in the baseline window yields base = 0 and a
# non-finite ratio. Those rows are dropped, which silently shrinks the sample for
# exactly the rarest pairs, so report how many were lost rather than letting it
# pass unnoticed.
recovery <- rr_phase %>%
  left_join(baseline, by = c("seed", "deme_X", "deme_Y")) %>%
  mutate(ratio = nRR / base)

dropped <- recovery %>%
  filter(!is.finite(ratio)) %>%
  count(kind, name = "n_dropped")
if (nrow(dropped) > 0) {
  message(
    "Dropped rows with an empty baseline window: ",
    paste0(dropped$kind, "=", dropped$n_dropped, collapse = ", "),
    " (of ", nrow(recovery), " total)"
  )
}

# Burn-in is excluded: the seeding period is intrinsically dense in identical
# sequences (one introduction, no standing diversity), so it is not comparable to
# the later phases and would sit high for reasons that have nothing to do with
# the intervention.
RECOVERY_PHASES <- PHASE_NAMES[PHASE_NAMES != "burn-in"]

recovery <- recovery %>%
  filter(is.finite(ratio)) %>%
  filter(phase %in% RECOVERY_PHASES) %>%
  group_by(seed, kind, phase) %>%
  summarize(ratio = median(ratio, na.rm = TRUE), .groups = "drop") %>%
  mutate(phase = factor(phase, levels = RECOVERY_PHASES))

p <- ggplot(recovery, aes(x = phase, y = ratio, fill = kind)) +
  geom_boxplot(alpha = 0.75, outlier.size = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  recovery_y_scale() +
  scale_fill_manual(values = c("Intra-cluster" = "#38c", "Inter-cluster" = "#e94")) +
  labs(
    x = "Phase", y = "nRR relative to baseline (pseudo-log scale)", fill = NULL,
    title = paste0("Recovery by phase (scenario: ", scenario_label(scenario), ")"),
    subtitle = paste0(
      "1.0 = back to baseline (sliding window at t = ", BASELINE_TPOINT, ")"
    )
  ) +
  theme_bw()
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

# ---------------------------------------------------------------------------
# Effective reproduction number
# ---------------------------------------------------------------------------
# One trace per seed, with a loess fit pooled over all of them. Drawn from t = 0,
# with t < RT_SEEDING_END shaded as the seeding period: R(t) spikes to ~10 there
# and runs off the top of the fixed axis, which is expected rather than a signal.
# Named "Seeding" rather than "burn-in" because the burn-in *phase* in
# PHASE_NAMES is the longer t < 100. Note the spike also widens the loess
# confidence band across the whole series, since loess uses one residual scale.
# The right edge is cropped at PLOT_TIME_MAX, because hosts infected near
# MAX_TIME lose onward infections to the stop and R(t) falls as an artefact. As
# with the nRR figures, that crop is applied to the data so the censored tail
# cannot pull the loess fit down inside the visible range.
RT_SPAN <- 0.2
RT_SEEDING_END <- 50
# Fixed rather than pooled so all three arms (and the combined figure) share one
# axis; coord_cartesian() clips the view without dropping data from the fit.
# Mirrored, with RT_SEEDING_END, in two_cluster_combined_plot.R.
RT_LIM <- c(0, 3)

rt_file <- out_file("rt.tsv")
if (file.exists(rt_file)) {
  rt <- read_tsv(rt_file, show_col_types = FALSE) %>%
    filter(tpoint <= PLOT_TIME_MAX, is.finite(Rt))

  marks <- PHASE_MARKS[PHASE_MARKS <= PLOT_TIME_MAX]
  p <- ggplot(rt, aes(x = tpoint, y = Rt)) +
    annotate("rect", xmin = 0, xmax = RT_SEEDING_END, ymin = -Inf, ymax = Inf,
             fill = "#e94", alpha = 0.15) +
    annotate("text", x = RT_SEEDING_END / 2, y = RT_LIM[1] + 0.12,
             label = "Seeding", color = "#b65f1c", size = 3.5) +
    geom_vline(xintercept = marks, linetype = "dashed", color = "grey30", linewidth = 0.4) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_line(aes(group = seed), color = "grey50", alpha = 0.35, linewidth = 0.3) +
    geom_smooth(se = TRUE, method = "loess", formula = y ~ x, span = RT_SPAN,
                color = "#38c", fill = "#38c") +
    coord_cartesian(ylim = RT_LIM) +
    labs(
      x = "Simulation Time (days)", y = expression(R[t]),
      title = paste0("Effective reproduction number (scenario: ", scenario_label(scenario), ")"),
      subtitle = paste0(
        "Mean offspring of hosts infected in a ", rt$window[1],
        "-day window; one line per seed"
      )
    ) +
    theme_bw()
  ggsave(out_file("rt.jpg"), plot = p, width = 10, height = 6, dpi = 192)
} else {
  message("No ", rt_file, "; skipping R(t) plot (run two_cluster_reproduction.R)")
}
