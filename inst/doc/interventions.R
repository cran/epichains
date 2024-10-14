## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  collapse = TRUE,
  comment = "#>"
)

## ----load_libraries-----------------------------------------------------------
## main package
library("epichains")
## for plotting
library("ggplot2")
## for truncating the offspring distribution later
library("truncdist")

## ----simulate_chains----------------------------------------------------------
sims <- simulate_chain_stats(
  n_chains = 200, offspring_dist = rnbinom, stat_threshold = 99, mu = 1.2,
  size = 0.5, statistic = "size"
)

## ----uncontrolled_chains_plot-------------------------------------------------
sims[is.infinite(sims)] <- 100 # Replace infections > 99 with 100 for plotting.
ggplot(data.frame(x = sims), aes(x = x)) +
  geom_histogram(breaks = seq(0, 100, by = 5), closed = "left") +
  scale_x_continuous(
    breaks = c(0, 25, 50, 75, 100),
    labels = c(0, 25, 50, 75, ">99")
  ) +
  theme_bw()

## ----simulate_chains_pop_control----------------------------------------------
sims <- simulate_chain_stats(
  n_chains = 200, offspring_dist = rnbinom, stat_threshold = 99, mu = 0.9,
  size = 0.5, statistic = "size"
)
sims[is.infinite(sims)] <- 100 # Replace infections > 99 with 100 for plotting.
ggplot(data.frame(x = sims), aes(x = x)) +
  geom_histogram(breaks = seq(0, 100, by = 5), closed = "left") +
  scale_x_continuous(
    breaks = c(0, 25, 50, 75, 100),
    labels = c(0, 25, 50, 75, ">99")
  ) +
  theme_bw()

## ----nbinom_ind_control-------------------------------------------------------
rnbinom_ind <- function(n, ..., control = 0) {
  ## initialise number of offspring to 0
  offspring <- rep(0L, n)
  ## for each individual, decide whether they transmit further
  transmits <- rbinom(n = n, prob = 1 - control, size = 1)
  ## check if anyone transmits further
  if (any(transmits == 1L)) {
    ## for those that transmit, sample from negative binomial with given
    ## parameters
    offspring[which(transmits == 1L)] <- rnbinom(n = n, ...)
  }
  return(offspring)
}

## ----simulate_chains_ind_control----------------------------------------------
sims <- simulate_chain_stats(
  n_chains = 200, offspring_dist = rnbinom_ind, stat_threshold = 99, mu = 1.2,
  size = 0.5, control = 0.25, statistic = "size"
)
sims[is.infinite(sims)] <- 100 # Replace infections > 99 with 100 for plotting.
ggplot(data.frame(x = sims), aes(x = x)) +
  geom_histogram(breaks = seq(0, 100, by = 5), closed = "left") +
  scale_x_continuous(
    breaks = c(0, 25, 50, 75, 100),
    labels = c(0, 25, 50, 75, ">99")
  ) +
  theme_bw()

## ----negbin_truncated---------------------------------------------------------
rnbinom_truncated <- function(n, ..., max = Inf) {
  return(rtrunc(n = n, spec = "nbinom", b = max, ...))
}

## ----simulate_chains_truncated------------------------------------------------
sims <- simulate_chain_stats(
  n_chains = 200, offspring_dist = rnbinom_truncated, stat_threshold = 99,
  mu = 1.2, size = 0.5, max = 10, statistic = "size"
)
sims[is.infinite(sims)] <- 100 # Replace infections > 99 with 100 for plotting.
ggplot(data.frame(x = sims), aes(x = x)) +
  geom_histogram(breaks = seq(0, 100, by = 5), closed = "left") +
  scale_x_continuous(
    breaks = c(0, 25, 50, 75, 100),
    labels = c(0, 25, 50, 75, ">99")
  ) +
  theme_bw()

## ----truncate_gen_int---------------------------------------------------------
control <- 1 - pgamma(6, shape = 25, rate = 5)
signif(control, 2)

