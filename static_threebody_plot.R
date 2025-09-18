library(tidyverse)

df_rr <- read_tsv("static_three/rr.tsv")

rr_boxplot <- ggplot(df_rr,aes(x = deme_X, y = RR, fill = deme_Y)) + 
  geom_hline(yintercept = 1, linewidth = 1, linetype = "dashed") +
  geom_boxplot(linewidth = 0.3) +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Deme", y = "Identical Sequence RR", fill = "Deme") +
  theme(legend.position = "bottom")

ggsave(filename = "static_three/rr_boxplot.jpg",
  plot = rr_boxplot,
  dpi = 192,
  height = 7,
  width = 5,
  units = "in")
