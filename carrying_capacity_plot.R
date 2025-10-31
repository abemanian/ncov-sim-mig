library(tidyverse)

collect_args <- function() {
  parser <- argparse::ArgumentParser()
  parser$add_argument("--scenario",
  type = "character",
  default = "basic",
  help = "Scenario name: basic, seasonal_cyclic, seasonal_threshold, isolated, immune")
  parser$parse_args()
}

args <- collect_args()
scenario <- args$scenario

# Load RR data
rr_total <- read_tsv("carrying_three/rr.tsv")

# Boxplot of RR
rr_boxplot <- ggplot(rr_total, aes(x = deme_X, y = RR, fill = deme_Y)) +
  geom_hline(yintercept = 1, linewidth = 1, linetype = "dashed") +
  geom_boxplot(linewidth = 0.3) +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Deme", y = "Identical Sequence RR", fill = "Deme") +
  theme(legend.position = "bottom")

ggsave(
  filename = paste0("carrying_three/", scenario, "_rr_boxplot.jpg"),
  plot = rr_boxplot,
  dpi = 192,
  height = 7,
  width = 7,
  units = "in"
)

# Load sliding window RR data
rr_sliding <- read_tsv("carrying_three/rr_sliding.tsv") %>%
  mutate(seed = as.factor(seed))

# RR over time plot
rr_scatter_smooth <- ggplot(
  rr_sliding,
  aes(
    x = tpoint,
    y = RR + 1E-6,
    color = deme_Y,
    group = interaction(seed, deme_Y)
  )
) +
  geom_line(alpha = 0.2, linewidth = 0.2) +
  geom_smooth(aes(group = deme_Y),
    se = TRUE,
    method = "loess",
    formula = y ~ x,
    linewidth = 1
  ) +
  theme_bw() +
  labs(
    x = "Time",
    y = "RR",
    color = "Deme"
  ) +
  scale_y_log10() +
  coord_cartesian(ylim = c(1E-6, 10)) +
  facet_wrap(vars(deme_X), ncol = 1) +
  scale_color_manual(values = c(
    "A" = "#F8766D",
    "O" = "#00BFC4",
    "B" = "#7CAE00"
  ))

ggsave(
  filename = paste0("carrying_three/", scenario, "_rr_scatter_smooth.jpg"),
  plot = rr_scatter_smooth,
  dpi = 192,
  height = 9,
  width = 7,
  units = "in"
)

# Only plot off-diagonal (between demes)
rr_sliding_offdiag <- rr_sliding %>% filter(deme_X != deme_Y)

# Scatter and smooth fit for nRR over time
nrr_scatter_smooth <- ggplot(
  rr_sliding_offdiag,
  aes(
    x = tpoint, y = nRR, color = deme_Y,
    group = interaction(seed, deme_Y)
  )
) +
  geom_line(alpha = 0.2, linewidth = 0.2) +
  geom_smooth(aes(group = deme_Y),
    se = TRUE,
    method = "loess",
    formula = y ~ x,
    linewidth = 1
  ) +
  theme_bw() +
  labs(
    x = "Time",
    y = "nRR",
    color = "Deme"
  ) +
  facet_wrap(vars(deme_X), ncol = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_manual(values = c(
    "A" = "#F8766D",
    "O" = "#00BFC4",
    "B" = "#7CAE00"
  ))

ggsave(
  filename = paste0("carrying_three/", scenario, "_nrr_scatter_smooth.jpg"),
  plot = nrr_scatter_smooth,
  dpi = 192,
  height = 9,
  width = 7,
  units = "in"
)

# Load transmission pair analysis results
transmission_pairs <- read_tsv("carrying_three/transmission_pairs.tsv")

# Summary statistics
cat("\n=== Transmission Pair Sequence Identity Summary ===\n")
cat(sprintf("Mean proportion: %.3f\n", mean(transmission_pairs$proportion)))
cat(sprintf("Median proportion: %.3f\n", median(transmission_pairs$proportion)))
cat(sprintf("SD: %.3f\n", sd(transmission_pairs$proportion)))
cat(sprintf("Min: %.3f, Max: %.3f\n",
            min(transmission_pairs$proportion),
            max(transmission_pairs$proportion)))
cat(sprintf("Mean transmission pairs per simulation: %.1f\n",
            mean(transmission_pairs$total_transmission_pairs)))
cat(sprintf("Mean identical sequences per simulation: %.1f\n",
            mean(transmission_pairs$identical_sequences)))

# Histogram of proportion with identical sequences
transmission_hist <- ggplot(transmission_pairs, aes(x = proportion)) +
  geom_histogram(bins = 10, 
    fill = "#00BFC4", 
    color = "black",
    alpha = 0.7
  ) +
  geom_vline(aes(xintercept = mean(proportion)),
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Proportion of parent-child pairs with identical sequences",
    y = "Count (across simulations)",
    title = "Transmission Pair Sequence Identity"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = paste0("carrying_three/", scenario, "_transmission_pairs_hist.jpg"),
  plot = transmission_hist,
  dpi = 192,
  height = 6,
  width = 8,
  units = "in"
)

# Time series of proportion across seeds
transmission_timeseries <- ggplot(transmission_pairs, aes(x = seed, y = proportion)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, color = "red") +
  labs(
    x = "Simulation Seed",
    y = "Proportion with identical sequences",
    title = "Transmission Pair Sequence Identity Across Seeds"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 1))

ggsave(
  filename = paste0("carrying_three/", scenario, "_transmission_pairs_timeseries.jpg"),
  plot = transmission_timeseries,
  dpi = 192,
  height = 6,
  width = 10,
  units = "in"
)
