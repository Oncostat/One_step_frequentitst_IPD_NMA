# Démaragge du cluster =============================================================================



## Packages
library(purrr)
library(dplyr)
library(tidyr)
library(tibble)
library(furrr)

# library(tictoc)

library(survival)
library(splines)
library(lme4)
library(netmeta)



## Functions
source("fonctions.R")
source("sim_parallel_3c.R")

plan(multisession)

# Preparation des data =================================================================================

## Definir le chunk a utiliser -------------------------------------------------------------------------

chunk <- "chunk1"

## Recuperation des donnees ----------------------------------------------------------------------------

sim_df <- readRDS(file = paste0("sim/sim_df-", chunk))

res <- sim_df %>% group_by(n_sim) %>% nest() %>% transmute(res = future_map(data, do_sim))

saveRDS(res, file = paste0("sim/", chunk, "-results.RDS"))
