library(tidyverse)

source("calc_rr.R")
source("plot_pair_time.R")


static_meta <- read_tsv("./static_three/static_threebody_meta.tsv") %>% 
  group_by(inf.in) 
static_pairs <- read_tsv("./static_three/idseq.tsv")

colnames(static_pairs) <- c("X","Y","Dist")
static_pairs <- filter(static_pairs,X %in% static_meta$hosts.ID) %>%
  filter(Y %in% static_meta$hosts.ID)

plot_pair_time(static_meta,static_pairs)
static_rr <- calc_rr(static_meta,static_pairs)

ggplot(static_rr) +
  geom_tile(aes(x=deme_X,y=deme_Y,fill=RR)) +
  theme_minimal() +
  scale_fill_gradient2(high="#e94",low="#38c",
                       transform = "log10",midpoint=1) +
  labs(x = "Deme", y = "Deme")


shutdown_meta <- read_tsv("./shutdown_three/meta.tsv") %>% 
  group_by(inf.in) 
shutdown_pairs <- read_tsv("./shutdown_three/idseq.tsv")
colnames(shutdown_pairs) <- c("X","Y","Dist")
shutdown_pairs <- filter(shutdown_pairs,X %in% shutdown_meta$hosts.ID) %>%
  filter(Y %in% shutdown_meta$hosts.ID)

plot_pair_time(shutdown_meta,shutdown_pairs)
T_SHUTDOWN_START <- 40
T_SHUTDOWN_END <- 120
T_MAX <- max(shutdown_meta$inf.time)
rr_pre <-calc_rr(shutdown_meta,shutdown_pairs,TIME_WINDOW = TRUE,
                 TIME_START = 0,
                 TIME_END =  T_SHUTDOWN_START+10)
rr_during <-calc_rr(shutdown_meta,shutdown_pairs,TIME_WINDOW = TRUE,
                 TIME_START = T_SHUTDOWN_START-10,
                 TIME_END =  T_SHUTDOWN_END+10)
rr_post <-calc_rr(shutdown_meta,shutdown_pairs,TIME_WINDOW = TRUE,
                    TIME_START = T_SHUTDOWN_END-10,
                    TIME_END =  T_MAX)

ggplot(rr_pre) +
  geom_tile(aes(x=deme_X,y=deme_Y,fill=RR)) +
  theme_minimal() +
  scale_fill_gradient2(high="#e94",low="#38c",
                       transform = "log10",midpoint=1,
                       limits=c(0.05,5.0)) +
  labs(x = "Deme", y = "Deme",title="Pre-Lockdown")

ggplot(rr_during) +
  geom_tile(aes(x=deme_X,y=deme_Y,fill=RR)) +
  theme_minimal() +
  scale_fill_gradient2(high="#e94",low="#38c",
                       transform = "log10",midpoint=1,
                       limits=c(0.05,5.0)) +
  labs(x = "Deme", y = "Deme",title="During Lockdown")

ggplot(rr_post) +
  geom_tile(aes(x=deme_X,y=deme_Y,fill=RR)) +
  theme_minimal() +
  scale_fill_gradient2(high="#e94",low="#38c",
                       transform = "log10",midpoint=1,
                       limits=c(0.05,5.0))+
  labs(x = "Deme", y = "Deme",title="Post Lockdown")

nRR_BIN_WIDTH <- 20
nRR_TIME_SPACING <- 5
nRR_MAX_TIME <- 200
nRR_tpoints <- seq(nRR_BIN_WIDTH,nRR_MAX_TIME,nRR_TIME_SPACING)



# run over all time points and stack
rr_sliding <- map_dfr(nRR_tpoints, ~ calc_rr(shutdown_meta,shutdown_pairs,
                                             TIME_WINDOW = TRUE,
                                             TIME_START = .x - nRR_BIN_WIDTH,
                                             TIME_END =  .x + nRR_BIN_WIDTH,
                                             time_point = .x)) %>%
  normalized_rr()

ggplot(rr_sliding%>%filter(deme_X != deme_Y),
       aes(x = tpoint,
           y = nRR,
           color= deme_Y)) +
  geom_point() +
  geom_line() + 
  theme_bw() +
  labs(
    x = "Date",
    y = "nRR",
    color = "Deme"
  ) +
  geom_vline(xintercept = T_SHUTDOWN_START,linetype = "dashed") +
  geom_vline(xintercept = T_SHUTDOWN_END,linetype = "dashed") +
  facet_wrap(vars(deme_X),ncol=1)

