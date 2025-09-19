
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
qlibrary('treeio') #Via bioconductor
qlibrary('ggtree') # Via bioconductor
qlibrary('argparse')


collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--seed', type = 'integer', default = 66, help = 'Simulation Seed')
  parser$parse_args()
}

args <- collect_args()
sim_seed <- args$seed

INF_ALPHA_PRE <- 4 
INF_BETA_PRE <-6
INF_ALPHA_DURING <- 2
INF_BETA_DURING <- 8
INF_ALPHA_POST <- 3
INF_BETA_POST <- 7

INCUB_MEAN <- 3
INCUB_SD <- 0.5
PROB_EXIT <- 0.1

CONTACT_PRE <- 1.5
CONTACT_DURING <- 0.5
CONTACT_POST <- 0.9

PROB_MOVE_NORMAL <- 0.1
PROB_MOVE_SHUTDOWN <- 0.001
T_START_SHUTDOWN <- 40
T_STOP_SHUTDOWN <- 120

DRAW_PLOTS <- FALSE

#For Nosoi, diagonal should be zero as there is a separate movement probability
MOVE_MATRIX <- matrix(c(
  0,0.9,0.1,
  0.5,0,0.5,
  0.1,0.9,0),
  nrow = 3, byrow = TRUE,
  dimnames=list(c("A","O","B"),
                c("A","O","B"))
)


#Main functions necessary for each agent are the exit probability, number of contacts, and transmission probability
p_Exit_fct <- function(t){PROB_EXIT} #Probability host leaves sim (death or recovery)
p_Move_fct <- function(t,prestime, current.in){ #Probability host moves to a different deme
  if(prestime >= T_START_SHUTDOWN && prestime < T_STOP_SHUTDOWN){
    return(PROB_MOVE_SHUTDOWN)
  }else{
    return(PROB_MOVE_NORMAL)
  } 
}
n_Contact_fct <- function(t,prestime,current.in){
  contact_mean <- case_when(
    prestime <= T_START_SHUTDOWN ~ CONTACT_PRE,
    prestime <= T_STOP_SHUTDOWN ~ CONTACT_DURING,
    TRUE ~ CONTACT_POST
  )
  rpois(1, lambda = contact_mean) #Number of contacts for host
}
#Time dependent, host-independent force of infection
#Host-dependent incubation
p_Trans_fct <- function(t, prestime,t_incub) {
  p <- if (prestime < T_START_SHUTDOWN) {
    rbeta(1, INF_ALPHA_PRE, INF_BETA_PRE)
  } else if (prestime < T_STOP_SHUTDOWN) {
    rbeta(1, INF_ALPHA_DURING,  INF_BETA_DURING)
  } else {
    rbeta(1, INF_ALPHA_POST, INF_BETA_POST)
  }
  ifelse(t < t_incub, 0,p)  # no transmission until incubation has elapsed
}

d_Incub <- function(n, mean=INCUB_MEAN, sd=INCUB_SD){
  shape <- (mean/sd)^2
  rate  <- shape/mean
  rgamma(n, shape = shape, rate = rate)
}

param_pTrans <- list(
  t_incub = d_Incub
)

# d_Inf <- function(x){
#   rbeta(x,shape1 = INF_ALPHA,shape2=INF_BETA)
# } #Infection probability distribution
#param_pTrans = list(p_inf = d_Inf, t_incub=d_Incub) #Wrapper to pass back to p_Trans_fct



# Infection probability distribution plot data
td_Inf <- tibble(x = seq(0,1,0.002)) %>%
  mutate(pre_y    = dbeta(x, INF_ALPHA_PRE, INF_BETA_PRE),
         during_y = dbeta(x, INF_ALPHA_DURING, INF_BETA_DURING),
         post_y   = dbeta(x, INF_ALPHA_POST, INF_BETA_POST))

# Outbreak simulation wrapper
set.seed(sim_seed)
MAX_OUTBREAK <- 20000
MAX_TIME <- 1000
MAX_TRIES <- 50
MIN_INFECTED <- 100

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
                    timeDep.nContact = TRUE,
                    pExit = p_Exit_fct,
                    param.pExit = NULL,
                    timeDep.pExit = FALSE,
                    pTrans = p_Trans_fct,
                    param.pTrans = param_pTrans,
                    timeDep.pTrans = TRUE,
                    pMove = p_Move_fct,
                    param.pMove = list(),
                    timeDep.pMove = TRUE,
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
  ggplot(td_Inf, aes(x = x)) + 
    geom_area(aes(y = pre_y, fill = "Pre-Lockdown"), alpha = 0.4, color = "black", linewidth = 1) +
    geom_area(aes(y = during_y, fill = "During"), alpha = 0.4, color = "black", linewidth = 1) +
    geom_area(aes(y = post_y, fill = "Post-Lockdown"), alpha = 0.4, color = "black", linewidth = 1) +
    scale_fill_manual(
      name   = "Phase",
      values = c("Pre-Lockdown" = "firebrick1",
                 "During"       = "royalblue",
                 "Post-Lockdown"= "forestgreen"),
      breaks = c("Pre-Lockdown", "During", "Post-Lockdown")
    ) +
    labs(x = "Infection Probability", y = "Density") +
    theme_bw()

  ggplot(sim_summary, aes(x = t, y = n)) +
    geom_point(color = "firebrick") +
    geom_line(alpha = 0.6) +
    geom_vline(xintercept = T_START_SHUTDOWN, color = "royalblue", lty = "dashed", linewidth = 1) +
    geom_vline(xintercept = T_STOP_SHUTDOWN, color = "royalblue", lty = "dashed", linewidth = 1) +
    scale_y_log10() +
    labs(x = "Simulation Time", y = "New Cases") +
    theme_bw()

  graph_sim <- graph_from_data_frame(data_sim %>%
                                       filter(!grepl("NA", inf.by)) %>%
                                       select("inf.by", "hosts.ID"),
                                     directed = TRUE,
                                     vertices = data_sim)
  ggraph(graph_sim, layout = "kk") +
    geom_edge_link(alpha = 0.2) +
    geom_node_point(size = 3, aes(color = inf.time)) +
    scale_color_viridis(name = "Infection Time", option = "plasma") +
    theme_void()
  ggraph(graph_sim, layout = "kk") +
    geom_edge_link(alpha = 0.2) +
    geom_node_point(size = 3, aes(color = inf.in), alpha = 0.4) +
    scale_color_discrete(name = "Infection Deme") +
    theme_void()

  tree_sim <- getTransmissionTree(threebody_sim)
  tr <- tree_sim@phylo

  # ---- sampling parameters ----
  bin_size     <- 10     # time units per bin
  max_per_bin  <- 40     # cap per bin
  per_deme     <- FALSE  # TRUE = cap per (bin x deme) instead of per bin

  tip_meta <- data_sim %>%
    transmute(
      label = as.character(hosts.ID),
      inf_in = as.factor(inf.in),
      t_inf = as.numeric(inf.time),
      bin = floor(t_inf / bin_size) * bin_size
    ) %>%
    filter(!is.na(t_inf), label %in% tr$tip.label)

  if (!per_deme) {
    tree_sample <- tip_meta %>%
      group_by(bin) %>%
      group_modify(~ {
        .x %>% slice_sample(n = min(max_per_bin, nrow(.x)), replace = FALSE)
      }) %>%
      ungroup()
  } else {
    tree_sample <- tip_meta %>%
      group_by(bin, inf_in) %>%
      group_modify(~ {
        .x %>% slice_sample(n = min(max_per_bin, nrow(.x)), replace = FALSE)
      }) %>%
      ungroup()
  }

  tree_labels <- tree_sample$label
  tr <- tr %>% keep.tip(tree_labels)
  sample_meta <- tree_sample %>% select(label, inf_in, t_inf, bin)

  ggtree(tr, size = 0.2, alpha = 0.4) %<+% sample_meta +
    geom_tippoint(aes(color = inf_in), size = 3, alpha = 0.6) +
    theme_tree2() +
    labs(color = "Infection Deme")
}


# Ensure output directory exists
if (!dir.exists("shutdown_three")) dir.create("shutdown_three")

tree_sim <- getTransmissionTree(threebody_sim)
write.tree(tree_sim@phylo, "shutdown_three/tree.nwk")
write_tsv(data_sim, "shutdown_three/meta.tsv")
