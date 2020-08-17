# Démaragge du cluster =============================================================================



## Packages
library(purrr)
library(dplyr)
library(tidyr)
library(tibble)
library(furrr)



## Functions
source("fonctions.R")

if(! dir.exists("sim")) dir.create("sim")

plan(multisession, seed = 1611)

#tic()

# Determiner les parametres =========================================================================

# on determine ici le scenario pour toute la simulation, trop gros de tout faire d'un coup
type_scenario <- c("none", "interaction", "both")

parameters <- list(
  scenario = c("none", "interaction", "both"),
  ttt = c(-2L, -5L),
  sig = c(10L),
  tau = c(10L),
  n_sim = 1:100,
  comp = c("AC", "BC", "AB"),
  essai = 1:10
)


parameters <- parameters  %>% cross_df() %>% 
  filter(
    (comp != "AB") |(comp == "AB" & essai <= 0),## on a moins d'essais pour la comp directe
    scenario %in% type_scenario,
    # sig == 0.01,
    # tau == 0.01,
    # ttt == -0.5,
    n_sim < 200) %>% 
  mutate(
    scenario = factor(scenario),
    comp = factor(comp))


# Generation des donnees ============================================================================


sim_df <- parameters %>% 
  mutate(
    poisson_df = future_pmap(list(scenario = scenario, ttt = ttt / 10, sig = 1 / sig, tau = 1 / tau, comp = comp), survsim2)  
  )


# Decoupage en Chunk  ================================================================================

n_chunk <- 8

chunk <- cut(sim_df$n_sim, n_chunk, labels = paste0("chunk",1:n_chunk))

sim_df <- split(sim_df, chunk)

walk2(sim_df, levels(chunk), ~ saveRDS(.x,  file = paste0("sim/sim_df-",.y)))



