
# Libraries ---------------------------------------------------------------

library(bendr)
library(dplyr)

# Data Generation ---------------------------------------------------------

set.seed(3431)
R <- 5
d <- seq(from=0,to=R,by=0.01)
f_1 <- function(x) 2*((1/2)*dbeta(x/R,1,8) + (1/2)*dbeta(x/R,6,1))
f_2 <- function(x) 2*((1/5)*dbeta(x/R,3,2) + (2/3)*dbeta(x/R,3,1) + (2/15)*dbeta(x/R,1,1))
f_3 <- function(x) 2*((1/2)*dbeta(x/R,8,2) + (1/2)*dbeta(x/R,30,50))
num_schools <- 50

schools_1 <- rnhpp(nsim = num_schools,
                   lambda = function(y) f_1(y),
                   interval = c(0,R),seed = 3431,
                   max =max(f_1(d)))
schools_2 <- rnhpp(nsim = num_schools,
                   lambda = function(y) f_2(y),
                   interval = c(0,R),
                   seed = 3431,
                   max = max(f_2(d)))
schools_3 <- rnhpp(nsim = num_schools,
                   lambda = function(y) f_3(y) ,
                   interval = c(0,R),
                   seed = 3431,
                   max = max(f_3(d)))

school_data <- as_tibble(schools_1) %>%
    rename(school_id = sim_id,
           distances=event_times) %>%
    mutate(Intensity =1)
school_data <- rbind(school_data,
                     as_tibble(schools_2) %>%
                         mutate(school_id = sim_id + num_schools,
                                Intensity = 2) %>%
    rename(distances=event_times) %>%
        select(-sim_id)) %>%
    rbind(.,as_tibble(schools_3) %>%
              mutate(school_id = sim_id + 2*num_schools,
                     Intensity = 3) %>%
              rename(distances=event_times) %>%
              select(-sim_id))

school_data <- rbenvo::benvo(subject_data = school_data %>%
                  select(-distances) %>%
                  mutate(school_id=as.integer(school_id)) %>%
                      distinct(school_id),
              sub_bef_data = list(FFR=school_data %>%
                                  mutate(school_id=as.integer(school_id)) %>%
                                  rename(Distance=distances)))


# Write Data --------------------------------------------------------------

usethis::use_data(school_data, overwrite = TRUE)

