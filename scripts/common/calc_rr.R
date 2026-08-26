require(tidyverse)

calc_rr <- function(meta,pairs,
                    TIME_WINDOW = FALSE,
                    TIME_START = NULL,
                    TIME_END = NULL,
                    time_point = NULL,
                    sim_seed = NULL){
  
  reflected_pairs <- pairs %>%
    rename(X2 = X) %>%
    rename(X = Y) %>%
    rename(Y = X2)
  pairs_sym <- rbind(pairs,reflected_pairs)
  deme_list <- meta %>% select(hosts.ID,inf.in)
  time_list <- meta %>% select(hosts.ID,inf.time)
  
  pairs_sym <- pairs_sym %>% 
    left_join(deme_list,join_by(X == hosts.ID)) %>%
    rename(deme_X = inf.in) %>%
    left_join(deme_list,join_by(Y == hosts.ID)) %>%
    rename(deme_Y = inf.in) %>%
    left_join(time_list,join_by(X == hosts.ID)) %>%
    rename(date_X = inf.time) %>%
    left_join(time_list,join_by(Y == hosts.ID)) %>%
    rename(date_Y = inf.time) %>%
    mutate(deltaT = abs(date_Y - date_X))
  
  if(TIME_WINDOW){
    pairs_sym <- pairs_sym %>%
      filter(
        between(date_X, TIME_START, TIME_END), 
        between(date_Y, TIME_START, TIME_END)
      )
  }
  
  zero_grid <- expand_grid(deme_X = deme_list$inf.in %>% unique(),
                           deme_Y = deme_list$inf.in %>% unique())
  
  rr <- pairs_sym %>%
    group_by(deme_X,deme_Y) %>%
    summarize(n_XY = n() %>% as.numeric(), .groups = "drop") %>%
    right_join(zero_grid, by = c("deme_X", "deme_Y")) %>%
    mutate(n_XY = ifelse(is.na(n_XY),0,n_XY)) %>%
    group_by(deme_X) %>%
    mutate(n_X = sum(n_XY)) %>%
    ungroup() %>%
    group_by(deme_Y) %>%
    mutate(n_Y = sum(n_XY)) %>%
    ungroup() %>%
    mutate(N = sum(n_XY)) %>%
    mutate(RR = n_XY * N/(n_X * n_Y))
  
  #Add a time point for additional analyses like plotting and nRR
  if(TIME_WINDOW){
    rr <- rr %>% mutate(tpoint = time_point)
  }
  if(!is.null(sim_seed)){
    rr <- rr %>% mutate(seed = sim_seed)
  }
  return(rr)
}

normalized_rr <- function(df_rr){
  # First: filter to the diagonal (RR(x,x)) entries
  rr_diag <- df_rr %>%
    filter(deme_X == deme_Y) %>%
    select(tpoint, deme = deme_X, RR_diag = RR)
  # Join to get RR(x,x) and RR(y,y) for each row
  df_rr <- df_rr %>%
    left_join(rr_diag, by = c("tpoint", "deme_X" = "deme")) %>%
    rename(RR_xx = RR_diag) %>%
    left_join(rr_diag, by = c("tpoint", "deme_Y" = "deme")) %>%
    rename(RR_yy = RR_diag) %>%
    mutate(nRR = RR / rowMeans(across(c(RR_xx, RR_yy)), na.rm = TRUE))
  return(df_rr)
}
