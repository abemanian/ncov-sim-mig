library(nosoi)
library(tidyverse)
library(viridis)
library(igraph)
library(ggraph)
library(ape)
library(tidytree)
library(treeio) #Via bioconductor
library(ggtree) # Via bioconductor

INF_ALPHA = 2
INF_BETA = 4

INCUB_MEAN = 5
INCUB_SD = 3

PROB_EXIT = 0.07

CONTACT_MEAN = 1

#Main functions necessary for each agent are the exit probability, number of contacts, and transmission probability
p_Exit_fct <- function(t){PROB_EXIT} #Probability host leaves sim (death or recovery)
n_Contact_fct <- function(t) rpois(1, lambda = CONTACT_MEAN) #Number of contacts for host
p_Trans_fct <- function(t,p_inf,t_incub){ #Transmission probability, dependent on incubation period for timing 
  ifelse(t <  t_incub,0,p_inf)
}


d_Incub <- function(n, mean=INCUB_MEAN, sd=INCUB_SD){
  shape <- (mean/sd)^2
  rate  <- shape/mean
  rgamma(n, shape = shape, rate = rate)
}

d_Inf <- function(x){rbeta(x,shape1 = INF_ALPHA,shape2=INF_BETA)} #Infection probability distribution

td_Inf <- tibble(
  x = seq(0,1,0.002)
) %>% mutate(y = dbeta(x,shape1 = INF_ALPHA,shape2 = INF_BETA))

ggplot(td_Inf,aes(x=x,y=y)) + 
  geom_area(alpha=0.4,linewidth=1,fill="firebrick1",color="black") + 
  labs(title = paste0("Infection Probability Distribution: Alpha = ",INF_ALPHA," Beta = ",INF_BETA),
       x = "Infection Probability",
       y = "Density") +
  theme_bw() 

param_pTrans = list(p_inf = d_Inf, t_incub=d_Incub) #Wrapper to pass back to p_Trans_fct

set.seed(25)
MAX_OUTBREAK <- 50
MAX_TIME <- 100
basic_sim <- nosoiSim(type="single", #Single vs two host species
                      popStructure = "none", #Homogenous pop
                      length.sim = MAX_TIME, #Maximum number of timesteps
                      max.infected = MAX_OUTBREAK, #Maximum total number the simulation will infect
                      init.individuals = 1,
                      nContact = n_Contact_fct, #Main function for agent contacts
                      param.nContact = NA, #Any required parameters for main function 
                      timeDep.nContact = FALSE, #Does the main function depend on the simulation time
                      pExit = p_Exit_fct, #Same structure for exit probability
                      param.pExit = NA,
                      timeDep.pExit = FALSE,
                      pTrans = p_Trans_fct, #Finally for the transmission probability
                      param.pTrans = param_pTrans, #Pass the distributions we created for it
                      timeDep.pTrans = FALSE, #Even though pTrans does depend on the individual time, it doesn't depend on the absolute simulation time
                      prefix.host = "H",
                      print.progress = TRUE
)

sim_summary <- summary(basic_sim)

ggplot(sim_summary$dynamics, aes(x=t,y=Count)) +
  geom_point() +
  geom_line() +
  theme_bw()

data_sim<-getTableHosts(basic_sim)
graph_sim <- graph_from_data_frame(data_sim %>% 
                                     filter(!grepl("NA",inf.by)) %>% 
                                     select("inf.by","hosts.ID"), #Remove NA nodes for the root infections
                              directed=T,
                              vertices=data_sim)
ggraph(graph_sim, layout = "kk") +
  geom_edge_link(alpha = 0.4,size=1) +
  geom_node_point(size=20,aes(color = inf.time)) +   # vertex attr
  geom_node_label(aes(label=name),fill="white") +
  scale_color_viridis(name = "Infection Time",option="plasma") +
  theme_void()

ggsave("basic_out/basic_sim.jpg",width=12,height=12,dpi=192)

tree_sim <- getTransmissionTree(basic_sim)
ggtree(tree_sim@phylo) + geom_tiplab() + theme_tree2()
write.tree(tree_sim@phylo,"basic_out/basic.nwk")
ggsave("basic_out/basic_tree.jpg",width=8,height=10,dpi=192)

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
  parser$parse_args()
}

args <- collect_args()
sim_seed <- args$seed
draw_plots <- args$draw
