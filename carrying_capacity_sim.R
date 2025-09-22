options(tidyverse.quiet = TRUE)
qlibrary <- function(pkg, ...) {
  suppressPackageStartupMessages(library(pkg, quietly = TRUE, character.only = TRUE, ...))
}

qlibrary("nosoi")
qlibrary("tidyverse")
qlibrary('viridis')
qlibrary('igraph')
qlibrary('ggraph')
qlibrary('ape')
qlibrary('tidytree')
qlibrary('treeio')
qlibrary('ggtree')
qlibrary('argparse')

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--seed', type = 'integer', default = 66, help = 'Simulation Seed')
  parser$parse_args()
}

args <- collect_args()
sim_seed <- args$seed

DRAW_PLOTS <- FALSE

# Constants
INCUB_MEAN <- 3
INCUB_SD <- 0.5
PROB_EXIT <- 0.1
CONTACT_MEAN <- 1
PROB_MOVE <- 0.1

SEASONAL_CYCLIC <- FALSE
CYCLE_SCALE <- 0.9
CYCLE_PERIOD <- 200

SEASONAL_THRESHOLD <- FALSE
THRESHOLD_SCALE <- 0.5
THRESHOLD_TIME <- 150
THRESHOLD_WIDTH <- 200


# Carrying capacity parameters
MAX_CAPACITY <- 300
BASE_INF_ALPHA <- 4
BASE_INF_BETA <- 8

#Immunization Scenario Parameters
#For this scenario, deme A will have a fraction of
# the population immunized at a certain time point
IMMUNE_SCENARIO <- TRUE
IMMUNIZATION_TIME <- 150
IMMUNIZATION_FRACTION <- 0.8

# Movement matrix
MOVE_MATRIX <- matrix(c(
  0,0.9,0.1,
  0.5,0,0.5,
  0.1,0.9,0),
  nrow = 3, byrow = TRUE,
  dimnames=list(c("A","O","B"),
                c("A","O","B"))
)

# Main functions
p_Exit_fct <- function(t){PROB_EXIT}
p_Move_fct <- function(t){PROB_MOVE}
n_Contact_fct <- function(t) rpois(1, lambda = CONTACT_MEAN)

# Transmission probability with carrying capacity effect
# To use host.count, need to set hostCount.pTrans = TRUE and diff.pTrans = TRUE
calc_seasonal_factor <- function(t) {
  s <- 1
  if(SEASONAL_CYCLIC){
    s <- s * (CYCLE_SCALE * cos(2 * pi * t / CYCLE_PERIOD) + 1)
  }
  if(SEASONAL_THRESHOLD){
    s <- s * ifelse(t >= THRESHOLD_TIME & t <= THRESHOLD_TIME + THRESHOLD_WIDTH,
                    THRESHOLD_SCALE,
                    1)
  }
  return(s)
}
calc_capacity_factor <- function(n) {
  1 / (1 + exp((n - MAX_CAPACITY/2)/(MAX_CAPACITY/10)))
}

p_Trans_fct <- function(t, prestime, current.in, host.count, t_incub) {
  # Logistic decline in transmission as population grows
  if(current.in == "O") {place_scale = 1}
  if(current.in == "A") {
    place_scale = 1
    if (IMMUNE_SCENARIO && prestime >= IMMUNIZATION_TIME) {
      place_scale <- place_scale * (1 - IMMUNIZATION_FRACTION)
    }
  }
  if(current.in == "B") {place_scale = 1}
  seasonal_factor <- calc_seasonal_factor(prestime) # Seasonal variation
  capacity_factor <- calc_capacity_factor(host.count)
  p <- rbeta(1, BASE_INF_ALPHA, BASE_INF_BETA) * capacity_factor * seasonal_factor * place_scale
  ifelse(t < t_incub, 0, p)
}



d_Incub <- function(n, mean=INCUB_MEAN, sd=INCUB_SD){
  shape <- (mean/sd)^2
  rate  <- shape/mean
  rgamma(n, shape = shape, rate = rate)
}

param_pTrans <- list(t_incub = d_Incub)

# Outbreak simulation parameters
set.seed(sim_seed)
MAX_OUTBREAK <- 10000
MAX_TIME <- 1000
MAX_TRIES <- 50
MIN_INFECTED <- 300

simulate_until_size <- function(tries, min_inf, ...) {
  for (k in 1:tries) {
    sim <- nosoiSim(type="single",
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
                    timeDep.pMove = FALSE,
                    diff.pMove = FALSE,
                    prefix.host = "H",
                    print.progress = TRUE)
    n_inf <- nrow(getTableHosts(sim))
    if (n_inf >= min_inf) return(sim)
  }
  stop("Failed to reach target size after ", tries, " attempts.")
}

threebody_sim <- simulate_until_size(tries=MAX_TRIES, min_inf=MIN_INFECTED)

sim_summary <- summary(threebody_sim)$dynamics %>%
  group_by(t) %>%
  summarize(n = sum(Count))

data_sim <- getTableHosts(threebody_sim)

# Plots
if (DRAW_PLOTS) {
  # Incidence over time
  ggplot(sim_summary, aes(x = t, y = n)) +
    geom_point(color = "firebrick") +
    geom_line(alpha = 0.6) +
    labs(x = "Simulation Time", y = "New Cases") +
    theme_bw()
  ggsave(filename = "carrying_three/cases_over_time.jpg",
    dpi = 192,
    height = 5,
    width = 15,
    units = "in")
  
  t <- seq(0, 400, by=1)
  seasonal_values <- calc_seasonal_factor(t)
  df_seasonal <- data.frame(time = t, seasonal_factor = seasonal_values)
  ggplot(df_seasonal, aes(x = time, y = seasonal_factor)) +
    geom_line(color = "blue", size = 1) +
    labs(x = "Time", y = "Seasonal Scaling") +
    theme_bw()
  ggsave(filename = "carrying_three/seasonal_factor.jpg",
    dpi = 192,
    height = 5,
    width = 15,
    units = "in")
  
  n <- seq(0, MAX_CAPACITY, by=1)
  capacity_values <- calc_capacity_factor(n)
  df_capacity <- data.frame(population = n, capacity_factor = capacity_values)
  ggplot(df_capacity, aes(x = population, y = capacity_factor)) +
    geom_line(color = "red", size = 1) +
    labs(x = "Population Size", y = "Capacity Scaling") +
    theme_bw()
  ggsave(filename = "carrying_three/capacity_factor.jpg",
    dpi = 192,
    height = 5,
    width = 15,
    units = "in")

  # Transmission tree plot
  graph_sim <- graph_from_data_frame(data_sim %>%
                                       filter(!grepl("NA", inf.by)) %>%
                                       select("inf.by", "hosts.ID"),
                                     directed = TRUE,
                                     vertices = data_sim)
  
  ggraph(graph_sim, layout = "kk") +
    geom_edge_link(alpha = 0.2) +
    geom_node_point(size = 0.2, aes(color = inf.time)) +
    scale_color_viridis(name = "Infection Time", option = "plasma") +
    theme_void()
  ggsave(filename = "carrying_three/transmission_tree.jpg",
    dpi = 192,
    height = 10,
    width = 10,
    units = "in")
}

# Ensure output directory exists
if (!dir.exists("carrying_three")) dir.create("carrying_three")

tree_sim <- getTransmissionTree(threebody_sim)
write.tree(tree_sim@phylo, "carrying_three/tree.nwk")
write_tsv(data_sim, "carrying_three/meta.tsv")