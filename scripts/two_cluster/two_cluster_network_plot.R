# two_cluster_network_plot.R
# Diagram of the movement structure: two triangles joined by the single A-X
# bridge, with the baseline destination probability on each directed edge.
#
# NOTE: BASE_WEIGHTS below must match the matrix of the same name in
# two_cluster_sim.R. It is duplicated rather than sourced because sourcing that
# script would run a simulation.

options(tidyverse.quiet = TRUE)
library(tidyverse, quietly = TRUE)

# Written to the repository root, not results/: this diagram is a fixed property
# of the model rather than the output of any one sweep, so it does not belong
# with the per-run artefacts that get cleared between batches.
OUT <- "movement_network.pdf"

DEMES <- c("A", "B", "C", "X", "Y", "Z")
PROB_MOVE <- 0.10

BASE_WEIGHTS <- matrix(
  c(
    0.00, 0.35, 0.35, 0.30, 0.00, 0.00, # A -> B .35  C .35  X .30
    0.50, 0.00, 0.50, 0.00, 0.00, 0.00, # B -> A .50  C .50
    0.50, 0.50, 0.00, 0.00, 0.00, 0.00, # C -> A .50  B .50
    0.30, 0.00, 0.00, 0.00, 0.35, 0.35, # X -> A .30  Y .35  Z .35
    0.00, 0.00, 0.00, 0.50, 0.00, 0.50, # Y -> X .50  Z .50
    0.00, 0.00, 0.00, 0.50, 0.50, 0.00  # Z -> X .50  Y .50
  ),
  nrow = 6, byrow = TRUE,
  dimnames = list(DEMES, DEMES)
)

# Triangle ABC on the left, XYZ on the right, joined across the middle by A-X
nodes <- tribble(
  ~name, ~x,  ~y,
  "B",   0.0,  0.85,
  "C",   0.0, -0.85,
  "A",   1.0,  0.00,
  "X",   2.0,  0.00,
  "Y",   3.0,  0.85,
  "Z",   3.0, -0.85
) %>%
  mutate(cluster = if_else(name %in% c("A", "B", "C"), "C1", "C2"))

edges <- as.data.frame(as.table(BASE_WEIGHTS)) %>%
  set_names(c("from", "to", "prob")) %>%
  filter(prob > 0) %>%
  mutate(across(c(from, to), as.character)) %>%
  left_join(nodes %>% select(from = name, x0 = x, y0 = y), by = "from") %>%
  left_join(nodes %>% select(to = name, x1 = x, y1 = y), by = "to") %>%
  mutate(bridge = (from == "A" & to == "X") | (from == "X" & to == "A"))

NODE_R <- 0.155   # trim edges to the node circle so arrowheads sit outside it
CURVE <- 0.16     # perpendicular offset of the arc's control point

# Quadratic Bezier per directed edge. Drawing the arc explicitly (rather than
# with geom_curve) means the midpoint is known, so labels can be placed on the
# arc instead of guessed at relative to the straight chord.
arc <- function(x0, y0, x1, y1, offset, n = 60) {
  dx <- x1 - x0
  dy <- y1 - y0
  len <- sqrt(dx^2 + dy^2)
  ux <- dx / len
  uy <- dy / len
  # perpendicular; reversing the edge reverses it too, so opposite directions
  # bow to opposite sides and their labels never collide
  px <- -uy
  py <- ux

  sx <- x0 + ux * NODE_R
  sy <- y0 + uy * NODE_R
  ex <- x1 - ux * NODE_R
  ey <- y1 - uy * NODE_R

  cx <- (sx + ex) / 2 + px * offset * 2
  cy <- (sy + ey) / 2 + py * offset * 2

  t <- seq(0, 1, length.out = n)
  tibble(
    x = (1 - t)^2 * sx + 2 * (1 - t) * t * cx + t^2 * ex,
    y = (1 - t)^2 * sy + 2 * (1 - t) * t * cy + t^2 * ey
  )
}

edge_paths <- edges %>%
  mutate(id = row_number()) %>%
  rowwise() %>%
  mutate(pts = list(arc(x0, y0, x1, y1, CURVE))) %>%
  ungroup() %>%
  select(id, bridge, pts) %>%
  unnest(pts)

# Label sits on the arc's midpoint, nudged a little further out
edge_labels <- edges %>%
  mutate(id = row_number()) %>%
  rowwise() %>%
  mutate(
    mid = list(arc(x0, y0, x1, y1, CURVE, n = 3) %>% slice(2))
  ) %>%
  ungroup() %>%
  unnest(mid) %>%
  mutate(
    dx = x1 - x0, dy = y1 - y0,
    len = sqrt(dx^2 + dy^2),
    lx = x + (-dy / len) * 0.10,
    ly = y + (dx / len) * 0.10
  )

EDGE_COLORS <- c("Within cluster" = "#38c", "A-X bridge" = "#e94")
edge_paths <- edge_paths %>%
  mutate(kind = if_else(bridge, "A-X bridge", "Within cluster"))
edge_labels <- edge_labels %>%
  mutate(kind = if_else(bridge, "A-X bridge", "Within cluster"))

p <- ggplot() +
  geom_path(
    data = edge_paths,
    aes(x = x, y = y, group = id, color = kind),
    linewidth = 0.4,
    arrow = arrow(length = unit(0.05, "in"), type = "closed", angle = 18)
  ) +
  geom_text(
    data = edge_labels,
    aes(x = lx, y = ly, label = format(prob, nsmall = 2), color = kind),
    size = 2.1, show.legend = FALSE
  ) +
  geom_point(
    data = nodes, aes(x = x, y = y, fill = cluster),
    shape = 21, size = 8, stroke = 0.5, color = "grey25"
  ) +
  geom_text(data = nodes, aes(x = x, y = y, label = name),
            size = 3.1, fontface = "bold") +
  scale_color_manual(values = EDGE_COLORS, name = NULL) +
  # Green/purple for cluster membership, deliberately clear of the blue/orange
  # used for edge type, so node fill and edge color cannot be read as the same
  # encoding. Kept pale so the bold deme labels stay legible on top.
  scale_fill_manual(
    values = c(C1 = "#d3e8cd", C2 = "#e0d5ef"), guide = "none"
  ) +
  coord_fixed(clip = "off") +
  theme_void(base_size = 8) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -4),
    plot.margin = margin(6, 10, 4, 10)
  )

ggsave(OUT, plot = p, width = 6, height = 3, units = "in", device = cairo_pdf)
message("Wrote ", OUT)
