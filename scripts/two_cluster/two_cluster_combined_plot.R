# two_cluster_combined_plot.R
# Lays out all six arms (three scenarios x with/without the intervention) on one
# shared axis per figure: rows are scenarios, columns are the intervention.
#
# Reads the per-seed accumulators the batch writes (rr_phase, rr_sliding,
# rt_region), so run two_cluster_batch.sh for every arm first. Deme pairs are
# coloured by four groups (Within Region 1, Within Region 2, A-X bridge, Other
# across region) -- see "RR and nRR by deme pair" below.

options(tidyverse.quiet = TRUE)
library(tidyverse, quietly = TRUE)

DATA_DIR <- "results/two_cluster"
# The main combined figures go to figs/pdf/, alongside the movement design figures
FIG_PDF_DIR <- "figs/pdf"
dir.create(FIG_PDF_DIR, recursive = TRUE, showWarnings = FALSE)
# Burn-in is excluded: the seeding period is intrinsically dense in identical
# sequences, so it sits high for reasons unrelated to the intervention.
PHASE_NAMES <- c("shutdown", "initial", "moderate", "full")

# Ratios of exactly 0 are a real result -- a phase with no identical pairs for
# that deme pair -- and a log axis cannot draw them. A pseudo-log axis is
# log-like above sigma and linear through 0. sigma sits just below the smallest
# non-zero ratio observed (~0.011) so real data is not compressed.
PSEUDO_SIGMA <- 0.005
RECOVERY_BREAKS <- c(0, 0.01, 0.03, 0.1, 0.3, 1, 3)

# Panel order top to bottom: null model first, then the two policies
SCENARIO_LABELS <- c(
  none = "No Shutdown",
  homogeneous = "Uniform Recovery",
  heterogeneous = "Staggered Recovery"
)

# Columns left to right: without, then with, the intervention (sim --vaccine)
INTERVENTION_LABELS <- c("No Intervention", "Intervention")

# Reads <scenario>[_intervention]_<suffix> for all six arms; missing arms are
# skipped with a warning so a partial sweep still plots
read_arms <- function(suffix) {
  arms <- expand_grid(scenario = names(SCENARIO_LABELS), intervention = c(FALSE, TRUE))
  pmap_dfr(arms, function(scenario, intervention) {
    arm <- if (intervention) paste0(scenario, "_intervention") else scenario
    f <- file.path(DATA_DIR, paste0(arm, "_", suffix))
    if (!file.exists(f)) {
      warning("Missing ", f, "; skipping ", arm)
      return(NULL)
    }
    read_tsv(f, show_col_types = FALSE) %>%
      mutate(scenario = !!scenario, intervention = !!intervention)
  })
}

add_panels <- function(d) {
  d %>% mutate(
    panel = factor(unname(SCENARIO_LABELS[scenario]), levels = unname(SCENARIO_LABELS)),
    arm_col = factor(INTERVENTION_LABELS[intervention + 1], levels = INTERVENTION_LABELS)
  )
}

# ---------------------------------------------------------------------------
# Effective reproduction number by region
# ---------------------------------------------------------------------------
# Reads the per-seed accumulators from two_cluster_reproduction_region.R; one
# loess fit per region, pooled over seeds. The intervention acts on transmission
# in region 1 only, from t = 175, so Region 1 dropping while Region 2 holds is
# the check that it works. (A pooled, all-region R(t) figure was dropped on
# 2026-11-12: it only averaged these two curves.)
# PLOT_TIME_MAX and RT_SPAN mirror two_cluster_plot.R: the tail past 300 is
# depressed by the MAX_TIME stop, and is dropped before fitting, not just hidden.
PLOT_TIME_MAX <- 300
PHASE_MARKS <- c(100, 150, 200, 250)
RT_SPAN <- 0.2
# t < RT_SEEDING_END is shaded as the seeding period; its spike runs off the top
# of the fixed axis, which clips the view without dropping data from the fit.
RT_SEEDING_END <- 50
RT_LIM <- c(0, 3)
REGION_LABELS <- c(C1 = "Region 1 (A-B)", C2 = "Region 2 (X-Y)")

rt_region <- read_arms("rt_region.tsv")
if (nrow(rt_region) == 0) {
  message("No R(t) by region tables found under ", DATA_DIR, "; skipping R(t) by region")
} else {
  rt_region <- rt_region %>%
    filter(tpoint <= PLOT_TIME_MAX, is.finite(Rt)) %>%
    mutate(region = factor(unname(REGION_LABELS[region]), levels = unname(REGION_LABELS))) %>%
    add_panels()

  p <- ggplot(rt_region, aes(x = tpoint, y = Rt, color = region, fill = region)) +
    annotate("rect", xmin = 0, xmax = RT_SEEDING_END, ymin = -Inf, ymax = Inf,
             fill = "grey60", alpha = 0.15) +
    annotate("text", x = RT_SEEDING_END / 2, y = RT_LIM[1] + 0.15,
             label = "Seeding", color = "grey35", size = 2.5) +
    geom_vline(xintercept = PHASE_MARKS, linetype = "dashed",
               color = "grey40", linewidth = 0.3) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
    # Intervention start (TIME_VACC in two_cluster_sim.R)
    geom_vline(data = tibble(arm_col = factor("Intervention", levels = INTERVENTION_LABELS)),
               aes(xintercept = 175), color = "#b65f1c", linewidth = 0.4,
               inherit.aes = FALSE) +
    geom_smooth(se = TRUE, method = "loess", formula = y ~ x, span = RT_SPAN,
                linewidth = 0.6, alpha = 0.25) +
    scale_color_manual(values = c("#38c", "#e94"), aesthetics = c("color", "fill")) +
    facet_grid(rows = vars(panel), cols = vars(arm_col)) +
    coord_cartesian(ylim = RT_LIM) +
    labs(x = "Simulation Time (days)", y = expression(R[t]), color = NULL, fill = NULL) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold")
    )

  out <- file.path(FIG_PDF_DIR, "rt_region_combined.pdf")
  ggsave(out, plot = p, width = 9, height = 8, units = "in", device = cairo_pdf)
  message("Wrote ", out)
}

# ---------------------------------------------------------------------------
# RR and nRR by deme pair
# ---------------------------------------------------------------------------
# Every figure below is computed per deme pair -- no pooling across pairs --
# and coloured by the group the pair belongs to: Within Region 1, Within
# Region 2, the A-X bridge (the only inter-cluster edge), or Other across
# region. Same-deme pairs appear only in the RR figures (nRR is 1 there by
# definition). Each unordered pair is drawn once.
#
# The three non-bridge across-region pairs (A-Y, B-X, B-Y) need a chain that
# crosses A-X and then moves on, so they are the sparsest: expect exact zeros
# in the shutdown and initial phases (drawn at 0 on the pseudo-log axes) and
# the odd pair dropped for an empty baseline window.
C1_DEMES <- c("A", "B")
DEME_ORDER <- c("A", "B", "X", "Y")
POOLS <- c("Within Region 1", "Within Region 2", "A-X bridge", "Other across region")
SELF_POOLS <- c("Same deme, Region 1", "Same deme, Region 2")
POOL_COLORS <- c(
  "Within Region 1" = "#38c", "Within Region 2" = "#3a9a6a",
  "A-X bridge" = "#8b0000", "Other across region" = "#e94",
  "Same deme, Region 1" = "grey35", "Same deme, Region 2" = "grey65"
)
# Sliding windows are right-sided, [t - 40, t]. These three must match
# two_cluster_plot.R; see the comments there.
PLOT_TIME_MIN <- 100
SLIDING_TIME_MAX <- 320
BASELINE_TPOINT <- 100 # window [60, 100], the 40 days before the shutdown

pair_pool <- function(x, y) {
  in1 <- function(d) d %in% C1_DEMES
  case_when(
    x == y & in1(x) ~ "Same deme, Region 1",
    x == y ~ "Same deme, Region 2",
    (x == "A" & y == "X") | (x == "X" & y == "A") ~ "A-X bridge",
    in1(x) & in1(y) ~ "Within Region 1",
    !in1(x) & !in1(y) ~ "Within Region 2",
    TRUE ~ "Other across region"
  )
}

# One orientation per unordered pair, tagged with its group
by_pair <- function(d) {
  d %>%
    filter(match(deme_X, DEME_ORDER) <= match(deme_Y, DEME_ORDER)) %>%
    mutate(pool = pair_pool(deme_X, deme_Y), pair = paste0(deme_X, "-", deme_Y))
}

sliding <- read_arms("rr_sliding.tsv")
rr_phase <- read_arms("rr_phase.tsv")
if (nrow(sliding) == 0 || nrow(rr_phase) == 0) {
  stop("No sliding-window or phase RR tables found under ", DATA_DIR)
}
sliding <- by_pair(sliding)
rr_phase <- by_pair(filter(rr_phase, phase %in% PHASE_NAMES))

pair_theme <- function() {
  theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold")
    )
}

# --- Traces ----------------------------------------------------------------
# Spaghetti plot: one line per deme pair per seed, straight from the
# sliding-window table (RR and nRR exactly as calc_rr() / normalized_rr()
# compute them), so up to 15 pairs x 20 seeds per panel. Windows with no
# identical pairs for a pair are real zeros; the pseudo-log axis (as in the
# recovery figures) draws them at 0 instead of dropping them. Cropped to
# PLOT_TIME_MIN..SLIDING_TIME_MAX: the right edge is censored by the MAX_TIME stop.
# Each point is the window ending at tpoint, so a phase change at t0 shows from t0.
# On top, one opaque line per group: the median at each tpoint over every pair
# x seed in the group, zeros included.
TRACE_ALPHA <- 0.15
TRACE_BREAKS <- c(0, 0.01, 0.1, 1, 10)

plot_pair_traces <- function(value, y_label, pools, out_name) {
  d <- sliding %>%
    filter(between(tpoint, PLOT_TIME_MIN, SLIDING_TIME_MAX), pool %in% pools) %>%
    mutate(y = .data[[value]]) %>%
    filter(is.finite(y)) %>%
    mutate(pool = factor(pool, levels = pools)) %>%
    add_panels()
  med <- d %>%
    group_by(panel, arm_col, pool, tpoint) %>%
    summarize(y = median(y), .groups = "drop")

  marks <- PHASE_MARKS[PHASE_MARKS >= PLOT_TIME_MIN & PHASE_MARKS <= SLIDING_TIME_MAX]
  p <- ggplot(d, aes(x = tpoint, y = y, color = pool, group = interaction(pair, seed))) +
    geom_vline(xintercept = marks, linetype = "dashed", color = "grey40", linewidth = 0.3) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
    geom_line(linewidth = 0.25, alpha = TRACE_ALPHA) +
    geom_line(data = med, aes(group = pool), linewidth = 0.8) +
    # Intervention start (TIME_VACC in two_cluster_sim.R)
    geom_vline(data = tibble(arm_col = factor("Intervention", levels = INTERVENTION_LABELS)),
               aes(xintercept = 175), color = "#b65f1c", linewidth = 0.4,
               inherit.aes = FALSE) +
    scale_y_continuous(
      transform = scales::transform_pseudo_log(sigma = PSEUDO_SIGMA, base = 10),
      breaks = TRACE_BREAKS,
      labels = format(TRACE_BREAKS, trim = TRUE, drop0trailing = TRUE, scientific = FALSE)
    ) +
    scale_color_manual(values = POOL_COLORS[pools]) +
    facet_grid(rows = vars(panel), cols = vars(arm_col)) +
    labs(x = "Simulation Time (days)", y = y_label, color = NULL) +
    pair_theme() +
    guides(color = guide_legend(nrow = if (length(pools) > 4) 2 else 1,
                                override.aes = list(alpha = 1, linewidth = 0.8)))
  out <- file.path(FIG_PDF_DIR, out_name)
  ggsave(out, plot = p, width = 9, height = 8, units = "in", device = cairo_pdf)
  message("Wrote ", out)
}

plot_pair_traces("RR", "RR (pseudo-log scale)", c(POOLS, SELF_POOLS), "rr_by_pair_combined.pdf")
plot_pair_traces("nRR", "nRR (pseudo-log scale)", POOLS, "nrr_by_pair_combined.pdf")

# --- Recovery by phase -------------------------------------------------------
# Per seed and deme pair: the phase value over that pair's value in the
# sliding window at BASELINE_TPOINT (same seed). Each box pools every pair in
# the group across all 20 seeds (1 pair x 20 seeds = up to 20 points for each
# within-region group and A-X; 3 x 20 for Other across region). Pairs with an empty
# baseline window have no defined ratio and are dropped; the number kept is
# printed above each box. Burn-in is excluded.
recovery_for <- function(value) {
  base <- sliding %>%
    filter(tpoint == BASELINE_TPOINT) %>%
    select(scenario, intervention, seed, pair, base = all_of(value))
  rr_phase %>%
    filter(pool %in% POOLS) %>%
    left_join(base, by = c("scenario", "intervention", "seed", "pair")) %>%
    mutate(ratio = .data[[value]] / base) %>%
    filter(is.finite(ratio)) %>%
    mutate(
      phase = factor(phase, levels = PHASE_NAMES),
      pool = factor(pool, levels = POOLS)
    ) %>%
    select(scenario, intervention, seed, pair, pool, phase, ratio) %>%
    add_panels()
}

plot_recovery <- function(recovery, y_label, out_name) {
  n_kept <- recovery %>% count(panel, arm_col, phase, pool)
  dodge <- position_dodge(width = 0.85)
  p <- ggplot(recovery, aes(x = phase, y = ratio, fill = pool)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    geom_boxplot(alpha = 0.5, outlier.size = 0.4, outlier.alpha = 0.4,
                 linewidth = 0.3, position = dodge) +
    geom_text(data = n_kept, aes(y = Inf, label = n, color = pool), position = dodge,
              vjust = 1.5, size = 2, show.legend = FALSE) +
    facet_grid(rows = vars(panel), cols = vars(arm_col)) +
    scale_y_continuous(
      transform = scales::transform_pseudo_log(sigma = PSEUDO_SIGMA, base = 10),
      breaks = RECOVERY_BREAKS,
      labels = format(RECOVERY_BREAKS, trim = TRUE, drop0trailing = TRUE)
    ) +
    scale_fill_manual(values = POOL_COLORS[POOLS], aesthetics = c("fill", "color")) +
    labs(x = "Phase", y = y_label, fill = NULL) +
    pair_theme()
  out <- file.path(FIG_PDF_DIR, out_name)
  ggsave(out, plot = p, width = 9, height = 8, units = "in", device = cairo_pdf)
  message("Wrote ", out)
}

plot_recovery(recovery_for("nRR"), "nRR relative to baseline (pseudo-log scale)",
              "recovery_by_phase_combined.pdf")
plot_recovery(recovery_for("RR"), "RR relative to baseline (pseudo-log scale)",
              "recovery_by_phase_rr_combined.pdf")
