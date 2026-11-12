# two_cluster_network_plot.R
# Diagram of the movement structure: two deme pairs joined by the single A-X
# bridge, with the baseline destination probability on each directed edge.
# Also writes figs/pdf/movement_design.pdf: the same diagram stacked over the movement
# schedule of each scenario.
#
# NOTE: BASE_WEIGHTS, PHASE_STARTS, PHASE_ALL and PHASE_BRIDGE_HETERO below must match the matrix of the same name in
# two_cluster_sim.R. It is duplicated rather than sourced because sourcing that
# script would run a simulation.

options(tidyverse.quiet = TRUE)
library(tidyverse, quietly = TRUE)
library(patchwork, quietly = TRUE)

# Written to figs/pdf/, not results/: these figures are a fixed property of the
# model rather than the output of any one sweep.
FIG_PDF_DIR <- "figs/pdf"
dir.create(FIG_PDF_DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(FIG_PDF_DIR, "movement_network.pdf")
OUT_DESIGN <- file.path(FIG_PDF_DIR, "movement_design.pdf")

DEMES <- c("A", "B", "X", "Y")
PROB_MOVE <- 0.10

BASE_WEIGHTS <- matrix(
  c(
    0.00, 0.70, 0.30, 0.00, # A -> B .70  X .30
    1.00, 0.00, 0.00, 0.00, # B -> A 1.0
    0.30, 0.00, 0.00, 0.70, # X -> A .30  Y .70
    0.00, 0.00, 1.00, 0.00  # Y -> X 1.0
  ),
  nrow = 4, byrow = TRUE,
  dimnames = list(DEMES, DEMES)
)

# One row: B-A on the left, X-Y on the right, joined across the middle by A-X
nodes <- tribble(
  ~name, ~x,  ~y,
  "B",   0.0,  0.00,
  "A",   1.0,  0.00,
  "X",   2.0,  0.00,
  "Y",   3.0,  0.00
) %>%
  mutate(cluster = if_else(name %in% c("A", "B"), "C1", "C2"))

edges <- as.data.frame(as.table(BASE_WEIGHTS)) %>%
  set_names(c("from", "to", "prob")) %>%
  filter(prob > 0) %>%
  mutate(across(c(from, to), as.character)) %>%
  left_join(nodes %>% select(from = name, x0 = x, y0 = y), by = "from") %>%
  left_join(nodes %>% select(to = name, x1 = x, y1 = y), by = "to") %>%
  mutate(bridge = (from == "A" & to == "X") | (from == "X" & to == "A"))

NODE_R <- 0.26    # trim edges to the node circle so arrowheads sit outside it
OFFSET <- 0.055   # each direction runs parallel to the axis, this far to one side
LABEL_GAP <- 0.11 # label distance beyond its arrow

# Straight, parallel forward/reverse arrows between each pair, as in a
# reversible-reaction scheme. Offsetting to the left of the direction of travel
# puts opposite directions on opposite sides, so they and their labels never
# collide.
edges <- edges %>%
  mutate(
    dx = x1 - x0, dy = y1 - y0,
    len = sqrt(dx^2 + dy^2),
    ux = dx / len, uy = dy / len,
    px = -uy, py = ux,
    sx = x0 + ux * NODE_R + px * OFFSET,
    sy = y0 + uy * NODE_R + py * OFFSET,
    ex = x1 - ux * NODE_R + px * OFFSET,
    ey = y1 - uy * NODE_R + py * OFFSET
  )

edge_labels <- edges %>%
  mutate(
    lx = (x0 + x1) / 2 + px * (OFFSET + LABEL_GAP),
    ly = (y0 + y1) / 2 + py * (OFFSET + LABEL_GAP)
  )

FONT_PT <- 16     # every label, in points

EDGE_COLORS <- c("Within Region" = "#38c", "A-X bridge" = "#e94")
edges <- edges %>%
  mutate(kind = if_else(bridge, "A-X bridge", "Within Region"))
edge_labels <- edge_labels %>%
  mutate(kind = if_else(bridge, "A-X bridge", "Within Region"))

p <- ggplot() +
  geom_segment(
    data = edges,
    aes(x = sx, y = sy, xend = ex, yend = ey, color = kind),
    linewidth = 0.8,
    arrow = arrow(length = unit(0.12, "in"), type = "closed", angle = 20)
  ) +
  geom_text(
    data = edge_labels,
    aes(x = lx, y = ly, label = format(prob, nsmall = 2), color = kind),
    size = FONT_PT / .pt, show.legend = FALSE
  ) +
  geom_point(
    data = nodes, aes(x = x, y = y, fill = cluster),
    shape = 21, size = 24, stroke = 0.8, color = "grey25"
  ) +
  geom_text(data = nodes, aes(x = x, y = y, label = name),
            size = FONT_PT / .pt, fontface = "bold") +
  scale_color_manual(values = EDGE_COLORS, name = NULL) +
  # Green/purple for cluster membership, deliberately clear of the blue/orange
  # used for edge type, so node fill and edge color cannot be read as the same
  # encoding. Kept pale so the bold deme labels stay legible on top.
  scale_fill_manual(
    values = c(C1 = "#d3e8cd", C2 = "#e0d5ef"), guide = "none"
  ) +
  coord_fixed(clip = "off") +
  theme_void(base_size = FONT_PT) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = FONT_PT),
    legend.key.width = unit(0.5, "in"),
    legend.margin = margin(t = 14),
    plot.margin = margin(6, 10, 4, 10)
  )

ggsave(OUT, plot = p, width = 9, height = 2.6, units = "in", device = cairo_pdf)
message("Wrote ", OUT)

# ---------------------------------------------------------------------------
# Movement schedule by scenario, stacked under the diagram
# ---------------------------------------------------------------------------
# Mirrors the phase constants in two_cluster_sim.R. The intervention arms share
# these schedules (the intervention changes transmission, not movement), so
# there are three panels, not six.
MAX_TIME <- 350
PHASE_STARTS <- c(0, 100, 150, 200, 250)
PHASE_NAMES <- c("burn-in", "shutdown", "initial", "moderate", "full")
PHASE_ALL <- c(1.00, 0.10, 0.33, 0.66, 1.00)
PHASE_BRIDGE_HETERO <- c(1.00, 0.10, 0.10, 0.33, 0.50)

# Same display names as two_cluster_combined_plot.R
SCENARIO_LABELS <- c(
  none = "No Shutdown", homogeneous = "Uniform Recovery",
  heterogeneous = "Staggered Recovery"
)

schedule <- tibble(
  phase = PHASE_NAMES, start = PHASE_STARTS, end = c(PHASE_STARTS[-1], MAX_TIME),
  none_w = 1, none_b = 1,
  homogeneous_w = PHASE_ALL, homogeneous_b = PHASE_ALL,
  heterogeneous_w = PHASE_ALL, heterogeneous_b = PHASE_BRIDGE_HETERO
) %>%
  pivot_longer(-c(phase, start, end), names_to = c("scenario", "edge"),
               names_sep = "_", values_to = "mult") %>%
  mutate(
    scenario = factor(SCENARIO_LABELS[scenario], levels = SCENARIO_LABELS),
    kind = if_else(edge == "b", "A-X bridge", "Within Region")
  )

# Step path per scenario x edge: horizontal run per phase, vertical jumps between
sched_path <- schedule %>%
  arrange(scenario, kind, start) %>%
  group_by(scenario, kind) %>%
  reframe(t = c(rbind(start, end)), mult = rep(mult, each = 2))

phase_labels <- tibble(
  t = (PHASE_STARTS + c(PHASE_STARTS[-1], MAX_TIME)) / 2, label = PHASE_NAMES,
  scenario = factor(SCENARIO_LABELS[1], levels = SCENARIO_LABELS)
)

p_sched <- ggplot(sched_path, aes(x = t, y = mult, color = kind)) +
  geom_vline(xintercept = PHASE_STARTS[-1], linetype = "dashed", color = "grey60",
             linewidth = 0.4) +
  # Within drawn wide and solid, bridge narrower and dashed on top, so the two
  # stay distinguishable where they coincide (every phase under Uniform Recovery)
  geom_path(data = filter(sched_path, kind == "Within Region"), linewidth = 2.2) +
  geom_path(data = filter(sched_path, kind == "A-X bridge"), linewidth = 1.2,
            linetype = "22") +
  geom_text(data = phase_labels, aes(x = t, y = 1.28, label = label),
            inherit.aes = FALSE, size = FONT_PT / .pt, color = "grey30") +
  facet_wrap(~scenario, ncol = 1) +
  scale_color_manual(values = EDGE_COLORS, guide = "none") +
  scale_x_continuous(breaks = seq(0, MAX_TIME, 50), expand = expansion(mult = 0.01)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.4)) +
  coord_cartesian(clip = "off") +
  labs(x = "Simulation Time (days)", y = "Fraction of baseline movement") +
  theme_bw(base_size = FONT_PT) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = FONT_PT),
    axis.text = element_text(size = FONT_PT),
    axis.title = element_text(size = FONT_PT)
  )

design <- p / p_sched + plot_layout(heights = c(0.8, 3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = FONT_PT, face = "bold"))
ggsave(OUT_DESIGN, plot = design, width = 9, height = 10, units = "in", device = cairo_pdf)
message("Wrote ", OUT_DESIGN)
