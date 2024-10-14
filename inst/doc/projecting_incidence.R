## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  collapse = TRUE,
  comment = "#>"
)

## ----packages, include=TRUE---------------------------------------------------
library("epichains")
library("dplyr")
library("ggplot2")
library("lubridate")

## ----data---------------------------------------------------------------------
data("covid19_sa", package = "epichains")

## ----view_data----------------------------------------------------------------
seed_cases <- covid19_sa[1:5, ]
head(seed_cases)

## ----linelist_gen, message=FALSE----------------------------------------------
days_since_index <- as.integer(seed_cases$date - min(seed_cases$date))
days_since_index

## ----t0_setup-----------------------------------------------------------------
t0 <- rep(days_since_index, seed_cases$cases)
t0

## ----generation_time_setup, message=FALSE-------------------------------------
mu <- 4.7
sgma <- 2.9

log_mean <- log((mu^2) / (sqrt(sgma^2 + mu^2))) # log mean
log_sd <- sqrt(log(1 + (sgma / mu)^2)) # log sd

#' serial interval function
generation_time <- function(n) {
  gt <- rlnorm(n, meanlog = log_mean, sdlog = log_sd)
  return(gt)
}

## ----nbinom_args_setup, message=FALSE-----------------------------------------
mu <- 2.5
size <- 0.58

## ----time_args_setup, message=FALSE-------------------------------------------
#' Date to end simulation
projection_window <- 21
tf <- max(days_since_index) + projection_window
tf

## ----sim_reps_setup-----------------------------------------------------------
#' Number of simulations
sim_rep <- 100

## ----stat_threshold_setup-----------------------------------------------------
#' Maximum chain size allowed
stat_threshold <- 1000

## ----run_simulations, message=FALSE-------------------------------------------
set.seed(1234)
sim_chain_sizes <- lapply(
  seq_len(sim_rep),
  function(sim) {
    simulate_chains(
      n_chains = length(t0),
      offspring_dist = rnbinom,
      mu = mu,
      size = size,
      statistic = "size",
      stat_threshold = stat_threshold,
      generation_time = generation_time,
      t0 = t0,
      tf = tf
    ) %>%
      mutate(sim = sim)
  }
)

sim_output <- bind_rows(sim_chain_sizes)

## ----view_sim_output----------------------------------------------------------
head(sim_output)

## ----post_process_output------------------------------------------------------
# Daily number of cases for each simulation
incidence_ts <- sim_output %>%
  mutate(day = ceiling(time)) %>%
  count(sim, day, name = "cases") %>%
  as_tibble()

head(incidence_ts)

## ----add_dates----------------------------------------------------------------
# Get start date from the observed data
index_date <- min(seed_cases$date)
index_date

# Add a dates column to each simulation result
incidence_ts_by_date <- incidence_ts %>%
  group_by(sim) %>%
  mutate(date = index_date + days(seq(0, n() - 1))) %>%
  ungroup()

head(incidence_ts_by_date)

## ----aggregate_simulations----------------------------------------------------
# Median daily number of cases aggregated across all simulations
median_daily_cases <- incidence_ts_by_date %>%
  group_by(date) %>%
  summarise(median_cases = median(cases)) %>%
  ungroup() %>%
  arrange(date)

head(median_daily_cases)

## ----viz, fig.cap ="COVID-19 incidence in South Africa projected over a two week window in 2020. The light gray lines represent the individual simulations, the red line represents the median daily cases across all simulations, the black connected dots represent the observed data, and the dashed vertical line marks the beginning of the projection.", fig.width=6.0, fig.height=6----
# since all simulations may end at a different date, we will find the minimum
# final date for all simulations for the purposes of visualisation.
final_date <- incidence_ts_by_date %>%
  group_by(sim) %>%
  summarise(final_date = max(date), .groups = "drop") %>%
  summarise(min_final_date = min(final_date)) %>%
  pull(min_final_date)

incidence_ts_by_date <- incidence_ts_by_date %>%
  filter(date <= final_date)

median_daily_cases <- median_daily_cases %>%
  filter(date <= final_date)

ggplot(data = incidence_ts_by_date) +
  geom_line(
    aes(
      x = date,
      y = cases,
      group = sim
    ),
    color = "grey",
    linewidth = 0.2,
    alpha = 0.25
  ) +
  geom_line(
    data = median_daily_cases,
    aes(
      x = date,
      y = median_cases
    ),
    color = "tomato3",
    linewidth = 1.8
  ) +
  geom_point(
    data = covid19_sa,
    aes(
      x = date,
      y = cases
    ),
    color = "black",
    size = 1.75,
    shape = 21
  ) +
  geom_line(
    data = covid19_sa,
    aes(
      x = date,
      y = cases
    ),
    color = "black",
    linewidth = 1
  ) +
  scale_x_continuous(
    breaks = seq(
      min(incidence_ts_by_date$date),
      max(incidence_ts_by_date$date),
      5
    ),
    labels = seq(
      min(incidence_ts_by_date$date),
      max(incidence_ts_by_date$date),
      5
    )
  ) +
  scale_y_continuous(
    breaks = seq(
      0,
      max(incidence_ts_by_date$cases),
      30
    ),
    labels = seq(
      0,
      max(incidence_ts_by_date$cases),
      30
    )
  ) +
  geom_vline(
    mapping = aes(xintercept = max(seed_cases$date)),
    linetype = "dashed"
  ) +
  labs(x = "Date", y = "Daily cases")

