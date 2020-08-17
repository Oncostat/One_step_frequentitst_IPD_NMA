
library(tidyverse)

nbre <- 1:4
type <- "config2"

chunks <- paste0("simulations/", type, "/chunk", nbre, "-results.RDS")
age_chunks <- paste0("simulations/", type, "/chunk", nbre, "-age.RDS")

objets <- map_dfr(chunks, readRDS) %>% unnest_wider(res)
ages <- map_dfr(age_chunks, readRDS)


objets <- inner_join(objets, ages, by = "n_sim")

nettoyage <- function(x, y) {
  ungroup(x) %>% mutate(
  ttt = ttt / 10, 
  sig = 1/ sig, 
  tau = 1/ tau,
  n_sim = y
)
}

fusion_age <- function(x, y, z) {
  
  out <- inner_join(x, y, by = c("scenario", "ttt", "sig", "tau", "comp", "essai")) %>% mutate(n_sim = z)

  return(out)
}



objets <- objets %>% mutate(
  ad_models = pmap(list(ad_models, age, n_sim), fusion_age),
  age = NULL,
  ad_nets = map2(ad_nets, n_sim, nettoyage),
  ipd_2steps_models = map2(ipd_2steps_models, n_sim, nettoyage),
  ipd_nets = map2(ipd_nets, n_sim, nettoyage),
  ipd_ma_models = map2(ipd_ma_models, n_sim, nettoyage),
  ipd_1step_models = map2(ipd_1step_models, n_sim, nettoyage) 
)

noms <- c("ad_models", "ad_nets", "ipd_2steps_models", "ipd_nets", "ipd_ma_models", "ipd_1step_models")

walk2(objets[, -1], noms, 
      ~ write_rds(
        bind_rows(.x),
        path = paste0("simulations/",type, "/", .y, ".RDS"))
      )

