library(nosoi)
library(tidyverse)
library(viridis)
library(igraph)
library(ggraph)
library(ape)
library(tidytree)
library(treeio) #Via bioconductor
library(ggtree) # Via bioconductor

INF_ALPHA = 5
INF_BETA = 5

INCUB_MEAN = 3
INCUB_SD = 0.5
PROB_EXIT = 0.1
CONTACT_MEAN = 1

PROB_MOVE = 0.1
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
p_Move_fct <- function(t){PROB_MOVE} #Probability host moves to a different deme
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

set.seed(66)
MAX_OUTBREAK <- 10000
MAX_TIME <- 1000
threebody_sim <- nosoiSim(type="single", #Single vs two host species
                      popStructure = "discrete", #Tells nosoi there will be a structure
                      structure.matrix=MOVE_MATRIX,

                      length.sim = MAX_TIME, #Maximum number of timesteps
                      max.infected = MAX_OUTBREAK, #Maximum total number the simulation will infect
                      init.individuals = 1,
                      init.structure = "O",
                      
                      nContact = n_Contact_fct, #Main function for agent contacts
                      param.nContact = NA, #Any required parameters for main function 
                      timeDep.nContact = FALSE, #Does the main function depend on the simulation time
                      
                      pExit = p_Exit_fct, #Same structure for exit probability
                      param.pExit = NA,
                      timeDep.pExit = FALSE,
                      
                      pTrans = p_Trans_fct, #Finally for the transmission probability
                      param.pTrans = param_pTrans, #Pass the distributions we created for it
                      timeDep.pTrans = FALSE, #Even though pTrans does depend on the individual time, it doesn't depend on the absolute simulation time
                      
                      pMove = p_Move_fct,
                      param.pMove = NA,
                      timeDep.pMove = FALSE,
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
  scale_y_log10() +
  labs(x="Simulation Time",y="New Cases") +
  theme_bw()
data_sim<-getTableHosts(threebody_sim)
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

tree_sim <- getTransmissionTree(threebody_sim)
tip_meta <- data_sim %>%
  transmute(
    label = as.character(hosts.ID),
    inf_in = as.factor(inf.in)
  )
tr <- tree_sim@phylo
ggtree(tr) %<+% tip_meta +
  geom_tippoint(aes(color=inf_in),size=2) +
  theme_tree2() +
  labs(color="Infection Deme")
  
write.tree(tree_sim@phylo,"static_three/static_threebody.nwk")
write_tsv(data_sim,"static_three/static_threebody_meta.tsv")
