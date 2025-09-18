#Takes simulation pairs and calculates the RR for the specific seed

options(tidyverse.quiet = TRUE)
library(tidyverse,quietly = T)
library(argparse,quietly = T)

source("calc_rr.R")
source("plot_pair_time.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--seed', type = 'integer', default = 66, help = 'Simulation Seed')
  parser$parse_args()
}

args <- collect_args()
sim_seed <- args$seed

static_meta <- read_tsv("./static_three/static_threebody_meta.tsv") %>% 
  group_by(inf.in) 
static_pairs <- read_tsv("./static_three/idseq.tsv")

colnames(static_pairs) <- c("X","Y","Dist")
static_pairs <- filter(static_pairs,X %in% static_meta$hosts.ID) %>%
  filter(Y %in% static_meta$hosts.ID)

static_rr <- calc_rr(static_meta,static_pairs) %>%
  mutate(seed = sim_seed)

rr_fn <- "./static_three/rr.tsv"
if(file.exists("./static_three/rr.tsv")){
  write_tsv(static_rr,file = rr_fn,append = TRUE)
} else{
  write_tsv(static_rr,file = rr_fn)
}
