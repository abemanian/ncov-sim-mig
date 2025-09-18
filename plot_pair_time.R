require(tidyverse)

plot_pair_time <- function(meta,pairs){
  deme_list <- meta %>% select(hosts.ID,inf.in)
  time_list <- meta %>% select(hosts.ID,inf.time)
  
  pairs <- pairs %>% 
    left_join(deme_list,join_by(X == hosts.ID)) %>%
    rename(deme_X = inf.in) %>%
    left_join(deme_list,join_by(Y == hosts.ID)) %>%
    rename(deme_Y = inf.in) %>%
    left_join(time_list,join_by(X == hosts.ID)) %>%
    rename(date_X = inf.time) %>%
    left_join(time_list,join_by(Y == hosts.ID)) %>%
    rename(date_Y = inf.time) %>%
    mutate(deltaT = abs(date_Y - date_X))
  
  time_hist <- ggplot(pairs,aes(x=deltaT)) +
    geom_histogram(binwidth = 1,color="black",linewidth=0.5,fill="royalblue",alpha = 0.5) +
    theme_bw() +
    labs(x = "Difference in Case Date",y = "Pair Counts")
  
  return(time_hist)
}