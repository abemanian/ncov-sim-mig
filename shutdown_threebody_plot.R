library(tidyverse)

# Load phase RR data and set phase order
rr_phase <- read_tsv("shutdown_three/rr.tsv") %>%
  mutate(phase = factor(phase, 
                       levels = c("Pre-Lockdown", "During Lockdown", "Post-Lockdown"),
                       ordered = TRUE))

# Boxplot faceted by phase
rr_boxplot <- ggplot(rr_phase, aes(x = deme_X, y = RR, fill = deme_Y)) +
  geom_hline(yintercept = 1, linewidth = 1, linetype = "dashed") +
  geom_boxplot(linewidth = 0.3) +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Deme", y = "Identical Sequence RR", fill = "Deme") +
  theme(legend.position = "bottom") +
  facet_wrap(vars(phase), ncol = 1)

ggsave(filename = "shutdown_three/rr_boxplot.jpg",
  plot = rr_boxplot,
  dpi = 192,
  height = 12,
  width = 7,
  units = "in")

# Load sliding window nRR data
rr_sliding <- read_tsv("shutdown_three/rr_sliding.tsv")
rr_scatter_smooth <- ggplot(rr_sliding,
       aes(x = tpoint, y = RR + 1E-3, color = deme_Y)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(se = TRUE, method = "loess", formula = y ~ x, linewidth = 1) +
  theme_bw() +
  labs(
    x = "Date",
    y = "RR",
    color = "Deme"
  ) +
  geom_vline(xintercept = 40, linetype = "dashed") +
  geom_vline(xintercept = 120, linetype = "dashed") +
  scale_y_log10() +
  coord_cartesian(ylim = c(0.001, 10)) +
  facet_wrap(vars(deme_X), ncol = 1) +
  scale_color_manual(values = c("A" = "#F8766D", "O" = "#00BFC4", "B" = "#7CAE00"))
ggsave(filename = "shutdown_three/rr_scatter_smooth.jpg",
  plot = rr_scatter_smooth,
  dpi = 192,
  height = 9,
  width = 7,
  units = "in")

# Only plot off-diagonal (between demes)
rr_sliding_offdiag <- rr_sliding %>% filter(deme_X != deme_Y)

# Scatter and smooth fit for nRR over time, colored by deme_Y, faceted by deme_X
nrr_scatter_smooth <- ggplot(rr_sliding_offdiag,
       aes(x = tpoint, y = nRR, color = deme_Y)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(se = TRUE, method = "loess", formula = y ~ x, linewidth = 1) +
  theme_bw() +
  labs(
    x = "Date",
    y = "nRR",
    color = "Deme"
  ) +
  geom_vline(xintercept = 40, linetype = "dashed") +
  geom_vline(xintercept = 120, linetype = "dashed") +
  facet_wrap(vars(deme_X), ncol = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_manual(values = c("A" = "#F8766D", "O" = "#00BFC4", "B" = "#7CAE00"))

ggsave(filename = "shutdown_three/nrr_scatter_smooth.jpg",
  plot = nrr_scatter_smooth,
  dpi = 192,
  height = 9,
  width = 7,
  units = "in")
