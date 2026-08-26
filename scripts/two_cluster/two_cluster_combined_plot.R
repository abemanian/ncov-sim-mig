# two_cluster_combined_plot.R
# Stacks the per-scenario recovery_by_phase panels into a single figure, so the
# three arms can be read against one shared axis instead of across three files.
#
# Reads the recovery tables written by two_cluster_plot.R, so run that for every
# scenario first (two_cluster_batch.sh does this at the end of each sweep).

options(tidyverse.quiet = TRUE)
library(tidyverse, quietly = TRUE)

DATA_DIR <- "results/two_cluster"
# Burn-in is excluded: the seeding period is intrinsically dense in identical
# sequences (one introduction, no standing diversity), so it sits high for
# reasons unrelated to the intervention. two_cluster_plot.R already drops it;
# filtering again here keeps the figure correct against older tables.
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

recovery <- imap_dfr(SCENARIO_LABELS, function(label, scen) {
  f <- file.path(DATA_DIR, paste0(scen, "_recovery_by_phase.tsv"))
  if (!file.exists(f)) {
    warning("Missing ", f, "; skipping ", scen)
    return(NULL)
  }
  read_tsv(f, show_col_types = FALSE) %>% mutate(scenario = scen)
})

if (nrow(recovery) == 0) {
  stop("No recovery tables found under ", DATA_DIR)
}

recovery <- recovery %>%
  filter(phase %in% PHASE_NAMES) %>%
  mutate(
    phase = factor(phase, levels = PHASE_NAMES),
    panel = factor(
      unname(SCENARIO_LABELS[scenario]),
      levels = unname(SCENARIO_LABELS)
    )
  )

p <- ggplot(recovery, aes(x = phase, y = ratio, fill = kind)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_boxplot(alpha = 0.75, outlier.size = 0.5, linewidth = 0.35) +
  facet_wrap(vars(panel), ncol = 1) +
  scale_y_continuous(
    transform = scales::transform_pseudo_log(sigma = PSEUDO_SIGMA, base = 10),
    breaks = RECOVERY_BREAKS,
    labels = format(RECOVERY_BREAKS, trim = TRUE, drop0trailing = TRUE)
  ) +
  scale_fill_manual(
    values = c("Intra-cluster" = "#38c", "Inter-cluster" = "#e94")
  ) +
  labs(
    x = "Phase", y = "nRR relative to baseline (pseudo-log scale)", fill = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold")
  )

out <- file.path(DATA_DIR, "recovery_by_phase_combined.pdf")
ggsave(out, plot = p, width = 6, height = 8, units = "in", device = cairo_pdf)
message("Wrote ", out)
