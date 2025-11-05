options(tidyverse.quiet = TRUE)
qlibrary <- function(pkg, ...) {
  suppressPackageStartupMessages(library(pkg, 
  quietly = TRUE, 
  character.only = TRUE, ...))
}

qlibrary("nosoi")
qlibrary("tidyverse")
qlibrary("viridis")
qlibrary("igraph")
qlibrary("ggraph")
qlibrary("ape")
qlibrary("tidytree")
qlibrary("treeio")
qlibrary("ggtree")

collect_args <- function() {
  parser <- argparse::ArgumentParser()
  parser$add_argument("--seed",
  type = "integer",
  default = 66,
  help = "Simulation Seed")
  parser$add_argument("--draw",
  action = "store_true",
  default = FALSE,
  help = "Draw plots")
  parser$add_argument("--scenario",
  type = "character",
  default = "basic",
  help = "Scenario name: basic, seasonal_cyclic, seasonal_threshold, isolated, immune")
  parser$parse_args()
}

args <- collect_args()
sim_seed <- args$seed
draw_plots <- args$draw
scenario <- args$scenario

# Constants
INCUB_MEAN <- 3
INCUB_SD <- 0.5
PROB_EXIT <- 1.0/7.0
CONTACT_MEAN <- 1.5
PROB_MOVE <- 0.02

# Set scenario flags based on scenario name
SEASONAL_CYCLIC <- (scenario == "seasonal_cyclic")
SEASONAL_THRESHOLD <- (scenario == "seasonal_threshold")
ISOLATED_SCENARIO <- (scenario == "isolated")
IMMUNE_SCENARIO <- (scenario == "immune")

# Seasonal Cyclic Parameters
# If true, global transmission is scaled cyclically over time
CYCLE_SCALE <- 0.7
CYCLE_PERIOD <- 200

# Seasonal threshold parameters
# If true, global transmission is scaled down during a specific time window
THRESHOLD_SCALE <- 0.5
THRESHOLD_TIME <- 200
THRESHOLD_WIDTH <- 200

# Isolation Scenario Parameters
# If true, deme A will have reduced movement after a certain time point
ISOLATED_TIME <- 200
ISOLATED_MOVEMENT <- 0

# Carrying capacity parameters
MAX_CAPACITY <- 500
BASE_INF_ALPHA <- 4
BASE_INF_BETA <- 8

# Immunization Scenario Parameters
# For this scenario, deme A will have a fraction of
# the population immunized at a certain time point
IMMUNIZATION_TIME <- 150
IMMUNIZATION_FRACTION <- 0.8

# Movement matrix
MOVE_MATRIX <- matrix(
  c(
    0, 0.8, 0.2,
    0.5, 0, 0.5,
    0.2, 0.8, 0
  ),
  nrow = 3, byrow = TRUE,
  dimnames = list(
    c("A", "O", "B"),
    c("A", "O", "B")
  )
)

# Main functions
p_Exit_fct <- function(t) {
  PROB_EXIT
}

p_Move_fct <- function(t, prestime, current.in) {
  if (ISOLATED_SCENARIO && prestime >= ISOLATED_TIME) {
    if (current.in == "A") {
      p <- ISOLATED_MOVEMENT
    }
    if (current.in == "B") {
      p <- PROB_MOVE
    }
    if (current.in == "O") {
      p <- PROB_MOVE
    }
  } else {
    p <- PROB_MOVE
  }
  return(p) # Explicit return due to the if logic
}

n_Contact_fct <- function(t) rpois(1, lambda = CONTACT_MEAN)

# Transmission probability with carrying capacity effect
# To use host.count, need to set hostCount.pTrans = TRUE and diff.pTrans = TRUE
calc_seasonal_factor <- function(t) {
  s <- 1
  if (SEASONAL_CYCLIC) {
    s <- s * (CYCLE_SCALE * cos(2 * pi * t / CYCLE_PERIOD) + 1)
  }
  if (SEASONAL_THRESHOLD) {
    s <- s * ifelse(t >= THRESHOLD_TIME & t <= THRESHOLD_TIME + THRESHOLD_WIDTH,
      THRESHOLD_SCALE,
      1
    )
  }
  return(s)
}
calc_capacity_factor <- function(n) {
  1 / (1 + exp((n - MAX_CAPACITY / 2) / (MAX_CAPACITY / 10)))
}

p_Trans_fct <- function(t, prestime, current.in, host.count, t_incub) {
  # Logistic decline in transmission as population grows
  if (current.in == "O") {
    place_scale <- 1
  }
  if (current.in == "A") {
    place_scale <- 1
    if (IMMUNE_SCENARIO && prestime >= IMMUNIZATION_TIME) {
      place_scale <- place_scale * (1 - IMMUNIZATION_FRACTION)
    }
  }
  if (current.in == "B") {
    place_scale <- 1
  }
  seasonal_factor <- calc_seasonal_factor(prestime) # Seasonal variation
  capacity_factor <- calc_capacity_factor(host.count)
  p <- rbeta(1, BASE_INF_ALPHA, BASE_INF_BETA) * 
  capacity_factor * 
  seasonal_factor * 
  place_scale
  ifelse(t < t_incub, 0, p)
}



d_Incub <- function(n, mean = INCUB_MEAN, sd = INCUB_SD) {
  shape <- (mean / sd)^2
  rate <- shape / mean
  rgamma(n, shape = shape, rate = rate)
}

param_pTrans <- list(t_incub = d_Incub)

# Outbreak simulation parameters
set.seed(sim_seed)
MAX_OUTBREAK <- 50000
MAX_TIME <- 1000
MAX_TRIES <- 50
MIN_INFECTED <- 50

simulate_until_size <- function(tries, min_inf, ...) {
  for (k in 1:tries) {
    sim <- nosoi::nosoiSim(
      type = "single",
      popStructure = "discrete",
      structure.matrix = MOVE_MATRIX,
      length.sim = MAX_TIME,
      max.infected = MAX_OUTBREAK,
      init.individuals = 1,
      init.structure = "O",
      nContact = n_Contact_fct,
      param.nContact = list(),
      timeDep.nContact = FALSE,
      pExit = p_Exit_fct,
      param.pExit = NULL,
      timeDep.pExit = FALSE,
      pTrans = p_Trans_fct,
      param.pTrans = param_pTrans,
      timeDep.pTrans = TRUE,
      diff.pTrans = TRUE,
      hostCount.pTrans = TRUE,
      pMove = p_Move_fct,
      param.pMove = list(),
      timeDep.pMove = TRUE,
      diff.pMove = TRUE,
      prefix.host = "H",
      print.progress = TRUE
    )
    n_inf <- nrow(nosoi::getTableHosts(sim))
    if (n_inf >= min_inf) {
      return(sim)
    }
  }
  stop("Failed to reach target size after ", tries, " attempts.")
}

threebody_sim <- simulate_until_size(tries = MAX_TRIES, min_inf = MIN_INFECTED)

sim_summary <- summary(threebody_sim)$dynamics %>%
  group_by(t) %>%
  summarize(n = sum(Count))

data_sim <- getTableHosts(threebody_sim)

# Ensure output directory exists
if (!dir.exists("carrying_three")) dir.create("carrying_three")

# Full transmission tree
tree_sim <- getTransmissionTree(threebody_sim)
write.tree(tree_sim@phylo, "carrying_three/tree.nwk")
write_tsv(data_sim, "carrying_three/meta.tsv")


# Plots
if (draw_plots) {
  # Incidence over time
  ggplot(sim_summary, aes(x = t, y = n)) +
    geom_point(color = "firebrick") +
    geom_line(alpha = 0.6) +
    labs(x = "Simulation Time", y = "Active Infected") +
    theme_bw()
  ggsave(
    filename = paste0("carrying_three/", scenario, "_incidence.jpg"),
    dpi = 192,
    height = 5,
    width = 15,
    units = "in"
  )

  t <- seq(0, 400, by = 1)
  seasonal_values <- calc_seasonal_factor(t)
  df_seasonal <- data.frame(time = t, seasonal_factor = seasonal_values)
  ggplot(df_seasonal, aes(x = time, y = seasonal_factor)) +
    geom_line(color = "blue", size = 1) +
    labs(x = "Time", y = "Seasonal Scaling") +
    theme_bw()
  ggsave(
    filename = paste0("carrying_three/", scenario, "_seasonal_factor.jpg"),
    dpi = 192,
    height = 5,
    width = 15,
    units = "in"
  )

  n <- seq(0, MAX_CAPACITY, by = 1)
  capacity_values <- calc_capacity_factor(n)
  df_capacity <- data.frame(population = n, capacity_factor = capacity_values)
  ggplot(df_capacity, aes(x = population, y = capacity_factor)) +
    geom_line(color = "red", size = 1) +
    labs(x = "Population Size", y = "Capacity Scaling") +
    theme_bw()
  ggsave(
    filename = paste0("carrying_three/", scenario, "_capacity_factor.jpg"),
    dpi = 192,
    height = 5,
    width = 15,
    units = "in"
  )

  # Transmission tree plot
  graph_sim <- graph_from_data_frame(
    data_sim %>%
      filter(!grepl("NA", inf.by)) %>%
      select("inf.by", "hosts.ID"),
    directed = TRUE,
    vertices = data_sim
  )

  # Manual layout by location and time with relative Y positioning
  # Calculate Y position relative to parent node
  manual_layout <- data_sim %>%
    mutate(
      x = inf.time,
      y_delta = case_when(
        inf.in == "A" ~ 1,
        inf.in == "O" ~ 0,
        inf.in == "B" ~ -1
      )
    )

  # Initialize y positions
  y_positions <- numeric(nrow(manual_layout))
  names(y_positions) <- manual_layout$hosts.ID

  # Set root node (first infection) to y = 0
  y_positions[1] <- 0

  # Calculate y for each subsequent node based on parent's y + delta
  for (i in 2:nrow(manual_layout)) {
    parent_id <- manual_layout$inf.by[i]
    if (!is.na(parent_id) && parent_id %in% names(y_positions)) {
      y_positions[i] <- y_positions[parent_id] + manual_layout$y_delta[i]
    } else {
      y_positions[i] <- manual_layout$y_delta[i]
    }
  }

  manual_layout <- manual_layout %>%
    mutate(y = y_positions) %>%
    select(x, y)

  ggraph(graph_sim, layout = manual_layout) +
    geom_edge_link(alpha = 0.2, color = "grey80") +
    geom_node_point(aes(color = inf.in), size = 3,alpha = 0.7) +
    labs(x = "Infection Time", 
      y = "Cumulative Deme Position",
      color = "Current Deme") +
    theme_minimal() +
    theme(
      legend.position = "bottom"
    )
  ggsave(
    filename = paste0("carrying_three/", scenario, "_manual_layout.jpg"),
    dpi = 192,
    height = 8,
    width = 15,
    units = "in"
  )

  # Heatmap of connections between demes over time
  max_time <- max(data_sim$inf.time)
  time_breaks <- seq(0, ceiling(max_time / 20) * 20, by = 20)

  transmission_counts <- data_sim %>%
    filter(!grepl("NA", inf.by)) %>%
    left_join(data_sim %>% select(hosts.ID, inf.in),
              by = c("inf.by" = "hosts.ID")) %>%
    rename(source.in = inf.in.y, target.in = inf.in.x) %>%
    mutate(time_bin = cut(inf.time, breaks = time_breaks, include.lowest = TRUE)) %>%
    count(time_bin, source.in, target.in)

  ggplot(transmission_counts, aes(x = time_bin, y = interaction(source.in, target.in), fill = n)) +
    geom_tile() +
    scale_fill_viridis(name = "Transmissions", na.value = "gray80") +
    labs(x = "Time Period", y = "Source → Target Deme") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(
    filename = paste0("carrying_three/", scenario, "_transmission_heatmap.jpg"),
    dpi = 192,
    height = 8,
    width = 15,
    units = "in"
  )

  #Making phylo tree for visualization
  # Phylodynamic sampling: ~20 sequences every 10 time steps
  max_time <- max(data_sim$inf.time)
  time_bins <- seq(0, ceiling(max_time / 10) * 10, by = 10)
  samples_per_bin <- 50

  sampled_hosts <- data_sim %>%
    mutate(time_bin = cut(inf.time, breaks = time_bins, include.lowest = TRUE)) %>%
    filter(!is.na(time_bin)) %>%
    group_by(time_bin) %>%
    {
      bind_rows(lapply(split(., .$time_bin), function(df) {
        sample_n(df, min(samples_per_bin, nrow(df)))
      }))
    }

  # Create phylodynamic tree with sampled tips
  # Prune the full tree to keep only sampled tips
  sampled_tips <- sampled_hosts$hosts.ID
  tree_phylodyn <- ape::keep.tip(tree_sim@phylo, sampled_tips)

  # Plot phylodynamic tree colored by deme
  # Prepare metadata with label column for ggtree joining
  tree_metadata <- sampled_hosts %>%
    rename(label = hosts.ID)

  ggtree(tree_phylodyn, layout = "ape", alpha = 0.1) %<+% tree_metadata +
    geom_tippoint(aes(color = inf.in), size = 3, alpha = 0.8) +
    scale_color_discrete(name = "Infection Deme") +
    theme_tree()
  ggsave(
    filename = paste0("carrying_three/", scenario, "_phylodyn_tree.jpg"),
    dpi = 192,
    height = 10,
    width = 20,
    units = "in"
  )
}
