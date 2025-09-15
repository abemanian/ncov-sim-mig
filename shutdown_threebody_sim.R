library(nosoi)
library(tidyverse)
library(viridis)
library(igraph)
library(ggraph)
library(ape)
library(tidytree)
library(treeio) #Via bioconductor
library(ggtree) # Via bioconductor

INF_ALPHA_PRE <- 3
INF_BETA_PRE <- 8
INF_ALPHA_DURING <- 2
INF_BETA_DURING <- 8
INF_ALPHA_POST <- 2
INF_BETA_POST <- 8

INCUB_MEAN <- 7
INCUB_SD <- 1
PROB_EXIT <- 0.1
CONTACT_MEAN <- 1

PROB_MOVE_NORMAL <- 0.1
PROB_MOVE_SHUTDOWN <- 0.001
T_START_SHUTDOWN <- 40
T_STOP_SHUTDOWN <- 80

MAKE_GRAPH <- FALSE

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
n_Contact_fct <- function(t) rpois(1, lambda = CONTACT_MEAN) #Number of contacts for host

#Time dependent, host-independent force of infection
#Host-dependent incubation
p_Trans_fct <- function(t, prestime) {
  p <- if (prestime < T_START_SHUTDOWN) {
    rbeta(1, INF_ALPHA_PRE,  INF_BETA_PRE)
  } else if (prestime < T_STOP_SHUTDOWN) {
    rbeta(1, INF_ALPHA_DURING,  INF_BETA_DURING)
  } else {
    rbeta(1, INF_ALPHA_POST, INF_BETA_POST)
  }
  ifelse(t < 3, 0, p)  # no transmission until incubation has elapsed
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


td_Inf <- tibble(x = seq(0,1,0.002)) |>
  mutate(pre_y    = dbeta(x, INF_ALPHA_PRE, INF_BETA_PRE),
         during_y = dbeta(x, INF_ALPHA_DURING, INF_BETA_DURING),
         post_y   = dbeta(x, INF_ALPHA_POST, INF_BETA_POST))

ggplot(td_Inf,aes(x=x)) + 
  geom_area(aes(y=pre_y),alpha=0.4,linewidth=1,fill="firebrick1",color="black") + 
  geom_area(aes(y=during_y),alpha=0.4,linewidth=1,fill="royalblue",color="black") + 
  geom_area(aes(y=post_y),alpha=0.4,linewidth=1,fill="forestgreen",color="black") + 
  labs(x = "Infection Probability",
       y = "Density") +
  theme_bw() 


set.seed(66)
MAX_OUTBREAK <- 20000
MAX_TIME <- 1000
threebody_sim <- nosoiSim(type="single", #Single vs two host species
                      popStructure = "discrete", #Tells nosoi there will be a structure
                      structure.matrix=MOVE_MATRIX,

                      length.sim = MAX_TIME, #Maximum number of timesteps
                      max.infected = MAX_OUTBREAK, #Maximum total number the simulation will infect
                      init.individuals = 1,
                      init.structure = "O",
                      
                      nContact = n_Contact_fct, #Main function for agent contacts
                      param.nContact = NULL, #Any required parameters for main function 
                      timeDep.nContact = FALSE, #Does the main function depend on the simulation time
                      
                      pExit = p_Exit_fct, #Same structure for exit probability
                      param.pExit = NULL,
                      timeDep.pExit = FALSE,
                      
                      pTrans = p_Trans_fct, #Finally for the transmission probability
                      param.pTrans = param_pTrans, #Pass the distributions we created for it
                      timeDep.pTrans = TRUE, #Even though pTrans does depend on the individual time, it doesn't depend on the absolute simulation time
                      
                      pMove = p_Move_fct,
                      param.pMove = list(),
                      timeDep.pMove = TRUE,
                      diff.pMove = FALSE,
                      
                      prefix.host = "H",
                      print.progress = TRUE
)

sim_summary <- summary(threebody_sim)$dynamics %>% 
  group_by(t) %>%
  summarize(n=sum(Count))

ggplot(sim_summary, aes(x=t,y=n)) +
  geom_point(color="firebrick") +
  geom_line(alpha=0.6) +
  geom_vline(xintercept = 40,color="royalblue",lty="dashed",linewidth = 1) +
  geom_vline(xintercept = 80,color="royalblue",lty="dashed",linewidth = 1) +
  scale_y_log10() +
  labs(x="Simulation Time",y="New Cases") +
  theme_bw()

data_sim<-getTableHosts(threebody_sim)
if(MAKE_GRAPH == TRUE){
  graph_sim <- graph_from_data_frame(data_sim %>% 
                                       filter(!grepl("NA",inf.by)) %>% 
                                       select("inf.by","hosts.ID"), #Remove NA nodes for the root infections
                                directed=T,
                                vertices=data_sim)
  ggraph(graph_sim, layout = "kk") +
    geom_edge_link(alpha = 0.2) +
    geom_node_point(size=3,aes(color = inf.time)) +   # vertex attr
  #  geom_node_text(aes(label=name),color="#FFF") +
    scale_color_viridis(name = "Infection Time",option="plasma") +
    theme_void()
  ggraph(graph_sim, layout = "kk") +
    geom_edge_link(alpha = 0.2) +
    geom_node_point(size=3,aes(color = inf.in),alpha=0.4) +   # vertex attr
  #  geom_node_text(aes(label=name),color="#FFF") +
    scale_color_discrete(name = "Infection Deme") +
    theme_void()
}
tree_sim <- getTransmissionTree(threebody_sim)
tr <- tree_sim@phylo
#Tree visualization

# ---- sampling parameters ----
bin_size     <- 10     # time units per bin
max_per_bin  <- 40     # cap per bin
per_deme     <- FALSE  # TRUE = cap per (bin x deme) instead of per bin

tip_meta <- data_sim %>%
  transmute(
    label = as.character(hosts.ID),
    inf_in = as.factor(inf.in),
    t_inf = as.numeric(inf.time),
    bin = floor(t_inf/bin_size) * bin_size
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
sample_meta <- tree_sample %>% select(label,inf_in,t_inf,bin)

ggtree(tr,size = 0.2, alpha = 0.4) %<+% sample_meta +
  geom_tippoint(aes(color=inf_in),size=3,alpha=0.6) +
  theme_tree2() +
  labs(color="Infection Deme")
  write.tree(tree_sim@phylo,"shutdown_three/tree.nwk")
write_tsv(data_sim,"shutdown_three/meta.tsv")
