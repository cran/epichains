## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  collapse = TRUE,
  comment = "#>",
  fig.height = 6,
  fig.width = 6
)

## ----out.width = "75%", echo = FALSE, fig.cap = "An example of a transmission chain starting with a single case C1. Cases are represented by blue circles and arrows indicate who infected whom. The chain grows through generations Gen 1, Gen 2, and Gen 3, producing cases C2, C3, C4, C5, and C6. The chain ends at generation Gen 3 with cases C4, C5, and C6. The size of C1's chain is 6, including C1 (that is, the sum of all blue circles) and the length is 3, which includes Gen 1 (maximum number of generations reached by C1's chain).", fig.alt = "The figure shows an example of a transmission chain starting with a single case C1. The complete chain is depicted as blue circles linked by orange directed arrows showing the cases produced in each generation. Each generation is marked by a brown rectangle and labelled Gen 1, Gen 2, and Gen 3. The chain grows through generations Gen 1, Gen 2, and Gen 3. Case C1 starts in Gen 1 and produces cases C2, and C3 in Gen 2. In Gen 3, C2 produces C4 and C5, and C3 produces C6. The chain ends in Gen 3. The chain size of C1 is 6, which includes C1 (that is the sum of all the blue circles, representing cases) and the length is 3, which includes Gen 1 (maximum number of generations reached by C1's chain)."----
knitr::include_graphics(file.path("img", "transmission_chain_example.png"))

## ----load_packages------------------------------------------------------------
library("epichains")
library("epicontacts")

## ----estimate_likelihoods-----------------------------------------------------
set.seed(121)
# example of observed chain sizes
# randomly generate 20 chains of size between 1 to 10
chain_sizes <- sample(1:10, 20, replace = TRUE)
chain_sizes
# estimate loglikelihood of the observed chain sizes for given lambda
likelihood_eg <- likelihood(
  chains = chain_sizes,
  statistic = "size",
  offspring_dist = rpois,
  nsim_obs = 100,
  lambda = 0.5
)
# Print the estimate
likelihood_eg

## ----estimate_likelihoods2----------------------------------------------------
set.seed(121)
# example of observed chain sizes
# randomly generate 20 chains of size between 1 to 10
chain_sizes <- sample(1:10, 20, replace = TRUE)
chain_sizes

# estimate loglikelihood of the observed chain sizes
likelihood_ind_eg <- likelihood(
  chains = chain_sizes,
  statistic = "size",
  offspring_dist = rpois,
  nsim_obs = 100,
  lambda = 0.5,
  individual = TRUE
)
# Print the estimate
likelihood_ind_eg

## ----likelihoods_imperfect_obs------------------------------------------------
set.seed(121)
# example of observed chain sizes; randomly generate 20 chains of size 1 to 10
chain_sizes <- sample(1:10, 20, replace = TRUE)
# get their likelihood
liks <- likelihood(
  chains = chain_sizes,
  statistic = "size",
  offspring_dist = rpois,
  obs_prob = 0.3,
  lambda = 0.5,
  nsim_obs = 10
)
liks

## ----likelihoods_by_simulation------------------------------------------------
set.seed(121)
# example of observed chain sizes; randomly generate 20 chains of size 1 to 10
chain_sizes <- sample(1:10, 20, replace = TRUE)
# get their likelihood
liks <- likelihood(
  chains = chain_sizes,
  offspring_dist = rbinom,
  statistic = "size",
  size = 1,
  prob = 0.9,
  nsim_offspring = 250
)
liks

## ----sim_chains_pois----------------------------------------------------------
set.seed(32)
# Define generation time
generation_time_fn <- function(n) {
  gt <- rep(3, n)
  return(gt)
}

sim_chains <- simulate_chains(
  n_chains = 10,
  statistic = "size",
  offspring_dist = rpois,
  stat_threshold = 10,
  generation_time = generation_time_fn,
  lambda = 0.9
)

head(sim_chains)

## ----sim_chains_with_finite_pop-----------------------------------------------
set.seed(32)
# Define generation time
generation_time_fn <- function(n) {
  gt <- rep(3, n)
  return(gt)
}

sim_chains_with_pop <- simulate_chains(
  pop = 1000,
  n_chains = 10,
  percent_immune = 0.1,
  offspring_dist = rpois,
  statistic = "size",
  lambda = 1,
  generation_time = generation_time_fn
)

head(sim_chains_with_pop)

## ----sim_summary_pois---------------------------------------------------------
set.seed(123)

simulate_chain_stats_eg <- simulate_chain_stats(
  n_chains = 10,
  statistic = "size",
  offspring_dist = rpois,
  stat_threshold = 10,
  lambda = 0.9
)

# Print the results
simulate_chain_stats_eg

## ----summary_eg, include=TRUE,echo=TRUE---------------------------------------
set.seed(32)
# Define generation time
generation_time_fn <- function(n) {
  gt <- rep(3, n)
  return(gt)
}

sim_chains <- simulate_chains(
  n_chains = 10,
  statistic = "size",
  offspring_dist = rpois,
  stat_threshold = 10,
  generation_time = generation_time_fn,
  lambda = 0.9
)

summary(sim_chains)

# Example with simulate_chain_stats()
set.seed(32)

simulate_chain_stats_eg <- simulate_chain_stats(
  n_chains = 10,
  statistic = "size",
  offspring_dist = rpois,
  stat_threshold = 10,
  lambda = 0.9
)

# Get summaries
summary(simulate_chain_stats_eg)

## ----compare_sim_funcs, include=TRUE,echo=TRUE--------------------------------
setequal(
  # summary of output from simulate_chains()
  summary(sim_chains),
  # results from simulate_chain_stats()
  simulate_chain_stats_eg
)

## ----aggregation_eg, include=TRUE,echo=TRUE-----------------------------------
# Example with simulate_chains()
set.seed(32)

# Define generation time
generation_time_fn <- function(n) {
  gt <- rep(3, n)
  return(gt)
}

sim_chains <- simulate_chains(
  n_chains = 10,
  statistic = "size",
  offspring_dist = rpois,
  stat_threshold = 10,
  generation_time = generation_time_fn,
  lambda = 0.9
)

aggregate(sim_chains, by = "time")

## ----longest-chain------------------------------------------------------------
# Get the biggest chain
longest_chain <- sim_chains[sim_chains$chain == which.max(
  unname(table(sim_chains$chain))
), ]

# Convert to data.frame to view the whole data
as.data.frame(longest_chain)

## ----plot_chain---------------------------------------------------------------
# Create an `<epicontacts>` object
epc <- make_epicontacts(
  linelist = longest_chain,
  contacts = longest_chain,
  id = "infectee",
  from = "infector",
  to = "infectee",
  directed = TRUE
)

# Plot the chain
plot(epc)

## ----plot_eg------------------------------------------------------------------
# Run simulation with simulate_chains()
set.seed(32)
# Define generation time
generation_time_fn <- function(n) {
  gt <- rep(3, n)
  return(gt)
}

sim_chains <- simulate_chains(
  n_chains = 10,
  statistic = "size",
  offspring_dist = rpois,
  stat_threshold = 1000,
  generation_time = generation_time_fn,
  lambda = 2
)

# Aggregate cases over time
sim_chains_aggreg <- aggregate(sim_chains, by = "time")

# Plot cases over time
plot(
  sim_chains_aggreg,
  type = "b",
  col = "red",
  lwd = 2,
  xlab = "Time",
  ylab = "Cases"
)

