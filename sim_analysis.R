library(tidyverse)

static_meta <- read_tsv("./static_three/static_threebody_meta.tsv") %>% 
  group_by(inf.in) 

static_pairs <- read_tsv("./static_three/idseq.tsv")
colnames(static_pairs) <- c("X","Y","Dist")
static_pairs <- filter(static_pairs,X %in% static_meta$hosts.ID) %>%
  filter(Y %in% static_meta$hosts.ID)
reflected_pairs <- static_pairs %>%
  rename(X2 = X) %>%
  rename(X = Y) %>%
  rename(Y = X2)

static_pairs_sym <- rbind(static_pairs,reflected_pairs)

deme_list <- static_meta %>% select(hosts.ID,inf.in)
time_list <- static_meta %>% select(hosts.ID,inf.time)

static_pairs_sym <- static_pairs_sym %>% 
  left_join(deme_list,join_by(X == hosts.ID)) %>%
  rename(deme_X = inf.in) %>%
  left_join(deme_list,join_by(Y == hosts.ID)) %>%
  rename(deme_Y = inf.in) %>%
  left_join(time_list,join_by(X == hosts.ID)) %>%
  rename(date_X = inf.time) %>%
  left_join(time_list,join_by(Y == hosts.ID)) %>%
  rename(date_Y = inf.time) %>%
  mutate(deltaT = abs(date_Y - date_X))

ggplot(static_pairs_sym,aes(x=deltaT)) +
  geom_histogram(binwidth = 1,color="black",linewidth=0.5,fill="royalblue",alpha = 0.5) +
  theme_bw() +
  labs(x = "Difference in Case Date",y = "Pair Counts")

rr_static <- static_pairs_sym %>%
  group_by(deme_X,deme_Y) %>%
  summarize(n_XY = n() %>% as.numeric()) %>%
  group_by(deme_X) %>%
  mutate(n_X = sum(n_XY)) %>%
  ungroup() %>%
  group_by(deme_Y) %>%
  mutate(n_Y = sum(n_XY)) %>%
  ungroup() %>%
  mutate(N = sum(n_XY)) %>%
  mutate(RR = n_XY * N/(n_X * n_Y))

shutdown_meta <- read_tsv("./shutdown_three/meta.tsv") %>% 
  group_by(inf.in) 

shutdown_pairs <- read_tsv("./shutdown_three/idseq.tsv")
colnames(shutdown_pairs) <- c("X","Y","Dist")
shutdown_pairs <- filter(shutdown_pairs,X %in% shutdown_meta$hosts.ID) %>%
  filter(Y %in% shutdown_meta$hosts.ID)
reflected_pairs <- shutdown_pairs %>%
  rename(X2 = X) %>%
  rename(X = Y) %>%
  rename(Y = X2)

shutdown_pairs_sym <- rbind(shutdown_pairs,reflected_pairs)

deme_list <- shutdown_meta %>% select(hosts.ID,inf.in)
time_list <- shutdown_meta %>% select(hosts.ID,inf.time)


shutdown_pairs_full <- shutdown_pairs_sym %>% 
  left_join(deme_list,join_by(X == hosts.ID)) %>%
  rename(deme_X = inf.in) %>%
  left_join(deme_list,join_by(Y == hosts.ID)) %>%
  rename(deme_Y = inf.in) %>%
  left_join(time_list,join_by(X == hosts.ID)) %>%
  rename(date_X = inf.time) %>%
  left_join(time_list,join_by(Y == hosts.ID)) %>%
  rename(date_Y = inf.time) %>%
  mutate(deltaT = abs(date_Y - date_X))

T_SHUTDOWN_START <- 40
T_SHUTDOWN_END <- 80

rr_shutdown_pre <- shutdown_pairs_full %>% 
  filter(date_X < T_SHUTDOWN_START & date_Y < T_SHUTDOWN_START) %>%
  group_by(deme_X,deme_Y) %>%
  summarize(n_XY = n() %>% as.numeric()) %>%
  group_by(deme_X) %>%
  mutate(n_X = sum(n_XY)) %>%
  ungroup() %>%
  group_by(deme_Y) %>%
  mutate(n_Y = sum(n_XY)) %>%
  ungroup() %>%
  mutate(N = sum(n_XY)) %>%
  mutate(RR = n_XY * N/(n_X * n_Y))

rr_shutdown_during <- shutdown_pairs_full %>% 
  filter(date_X >= T_SHUTDOWN_START & date_Y >= T_SHUTDOWN_START) %>%
  filter(date_X < T_SHUTDOWN_END & date_Y < T_SHUTDOWN_END )%>%
  group_by(deme_X,deme_Y) %>%
  summarize(n_XY = n() %>% as.numeric()) %>%
  group_by(deme_X) %>%
  mutate(n_X = sum(n_XY)) %>%
  ungroup() %>%
  group_by(deme_Y) %>%
  mutate(n_Y = sum(n_XY)) %>%
  ungroup() %>%
  mutate(N = sum(n_XY)) %>%
  mutate(RR = n_XY * N/(n_X * n_Y))

rr_shutdown_post <- shutdown_pairs_full %>% 
  filter(date_X >= T_SHUTDOWN_END & date_Y >= T_SHUTDOWN_END )%>%
  group_by(deme_X,deme_Y) %>%
  summarize(n_XY = n() %>% as.numeric()) %>%
  group_by(deme_X) %>%
  mutate(n_X = sum(n_XY)) %>%
  ungroup() %>%
  group_by(deme_Y) %>%
  mutate(n_Y = sum(n_XY)) %>%
  ungroup() %>%
  mutate(N = sum(n_XY)) %>%
  mutate(RR = n_XY * N/(n_X * n_Y))

nRR_BIN_WIDTH <- 10
nRR_MAX_TIME <- 150
nRR_tpoints <- seq(nRR_BIN_WIDTH,nRR_MAX_TIME,nRR_BIN_WIDTH)

# helper: compute RR within the window (t - width, t]
compute_rr_window <- function(df, tpoint, width) {
  df %>%
    filter(
      date_X >=  tpoint - width, date_X <= tpoint + width,
      date_Y >=  tpoint - width, date_Y <= tpoint + width
    ) %>%
    group_by(deme_X, deme_Y) %>%
    summarise(n_XY = n() %>% as.numeric(), .groups = "drop_last") %>%
    group_by(deme_X) %>%
    mutate(n_X = sum(n_XY)) %>%
    ungroup() %>%
    group_by(deme_Y) %>%
    mutate(n_Y = sum(n_XY)) %>%
    ungroup() %>%
    mutate(
      N  = sum(n_XY),
      RR = (n_XY * N) / (n_X * n_Y)
    ) %>%
    # attach window metadata
    mutate(
      tpoint    = tpoint,
      bin_start = tpoint - width,
      bin_end   = tpoint + width
    ) %>%
    # guard against zero marginals
    filter(n_X > 0, n_Y > 0) %>%
    select(tpoint, bin_start, bin_end, deme_X, deme_Y, n_XY, n_X, n_Y, N, RR)
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

# run over all time points and stack
rr_sliding <- map_dfr(nRR_tpoints, ~ compute_rr_window(shutdown_pairs_full, .x, nRR_BIN_WIDTH)) %>%
  normalized_rr()

ggplot(rr_sliding%>%filter(deme_X != deme_Y),
       aes(x = tpoint,
           y = nRR,
           color= deme_Y)) +
  geom_point() +
  geom_line() + 
  theme_bw() +
  facet_wrap(vars(deme_X),ncol=1)
