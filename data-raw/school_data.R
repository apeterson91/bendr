
# Libraries ---------------------------------------------------------------

require(rndpp)
require(dplyr)
require(magrittr)

# Data Generation ---------------------------------------------------------

set.seed(3431)
R <- 10
d <- seq(from=0,to=R,by=0.01)
f_1 <- function(x) ((1/2)*dbeta(x/R,1,8) + (1/2)*dbeta(x/R,6,1))
f_2 <- function(x) ((1/5)*dbeta(x/R,3,2) + (2/3)*dbeta(x/R,3,1) + (2/15)*dbeta(x/R,1,1))
f_3 <- function(x) ((1/2)*dbeta(x/R,8,2) + (1/2)*dbeta(x/R,30,50))
num_schools <- 50
avg_FFR <- 11
schools_1 <- rnhpp(nsim = num_schools,lambda = function(y) avg_FFR*f_1(y),
                   interval = c(0,R),seed = 3431,
                   max =max(f_1(d)))
schools_2 <- rnhpp(nsim = num_schools,lambda = function(y) avg_FFR*f_2(y),interval = c(0,R),seed = 3431,max = max(f_2(d)))
schools_3 <- rnhpp(nsim = num_schools,lambda = function(y) avg_FFR*f_3(y) ,interval = c(0,R),seed = 3431,max = max(f_3(d)))
school_data <- as_tibble(schools_1) %>% rename(school_id = sim_id) %>% rename(distances=event_times)
school_data <- rbind(school_data,as_tibble(schools_2) %>%
                         mutate(school_id = sim_id + num_schools) %>%
    rename(distances=event_times) %>%
        select(-sim_id)) %>%
    mutate(density = ifelse(school_id<=num_schools,1,2))


# Write Data --------------------------------------------------------------

usethis::use_data(school_data, overwrite = TRUE)

