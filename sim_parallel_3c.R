
do_sim <- function(sim_df) {

  
## Preparation de l'output
  
results <- list(
  ad_models = NULL,
  ad_nets = NULL,
  ipd_2steps_models = NULL,
  ipd_nets = NULL,
  ipd_ma_models = NULL,
  ipd_1step_models = NULL
)  
  
## Transformation pour le Poisson -----------------------------------------------------------------------

sim_df$poisson_df <- map(sim_df$poisson_df, survsim4)




# Exploitation en essais separe =====================================================================

## Cox/Poisson non ajuste (AD) ------------------------------------------------------------------------------

ad_models <- sim_df %>% select(- poisson_df)

ad_models$poisson <-  map(sim_df$poisson_df, ~ glm2(status ~  ns(periode, 3) + treatment + offset(offset) , family = "poisson", data = .))

results$ad_models <- ad_models




## Net meta non ajuste  -------------------------------------------------------------------------------
ad_nets <- ad_models %>% mutate(
  TE = map_dbl(poisson,  ~ .[[1]]["treatment"]),
  seTE = map_dbl(poisson, ~ sqrt(.[[2]]["treatment", "treatment"])),
  treat1 = substr(comp, 1, 1),
  treat2 = substr(comp, 2, 2),
  poisson = NULL,
  studlab = paste0(comp, essai),
  essai = NULL,
  comp = NULL
) %>% group_by(scenario, ttt, sig, tau) %>% nest()


ad_nets <- ad_nets %>% mutate(net = map(data, ~ netmeta2(.$TE, .$seTE, .$treat1, .$treat2, .$studlab, data = ., sm = "HR")))


results$ad_nets <- ad_nets



## Cox/Poisson ajuste (IPD)   -------------------------------------------------------------------------------

ipd_2steps_models <- sim_df %>% select(- poisson_df)

ipd_2steps_models$res <- map(sim_df$poisson_df, poisson_lin_cat_glm)

ipd_2steps_models <- ipd_2steps_models %>% mutate(
  poisson_lin = map(res,"lin"),
  poisson_cat = map(res,"cat"),
  res = NULL
)

results$ipd_2steps_models <- ipd_2steps_models





## Net meta ajuste  -------------------------------------------------------------------------------

idx <- map_dbl(ipd_2steps_models$poisson_cat, ~ length(.[[1]])) ## cas ou une categorie n'existe pas

ipd_2steps_models <- ipd_2steps_models %>% filter(idx == 11) %>% mutate(
  TE = map_dbl(poisson_lin,  ~ .[[1]]["treatment"]),
  seTE = map_dbl(poisson_lin, ~ sqrt(.[[2]]["treatment", "treatment"])),
  TE_q1 = map_dbl(poisson_cat,  ~ .[[1]]["treatment"]),
  TE_q2 = map_dbl(poisson_cat,  ~ .[[1]]["treatment"] + .[[1]]["age_clq2"]),
  TE_q3 = map_dbl(poisson_cat,  ~ .[[1]]["treatment"] + .[[1]]["age_clq3"]),
  TE_q4 = map_dbl(poisson_cat,  ~ .[[1]]["treatment"] + .[[1]]["age_clq4"]),
  seTE_q1 = map_dbl(poisson_cat, ~ sqrt(.[[2]]["treatment", "treatment"])),
  seTE_q2 = map_dbl(poisson_cat, ~ sqrt(.[[2]]["treatment", "treatment"] + .[[2]]["age_clq2", "age_clq2"] + .[[2]]["treatment", "age_clq2"]*2)),
  seTE_q3 = map_dbl(poisson_cat, ~ sqrt(.[[2]]["treatment", "treatment"] + .[[2]]["age_clq3", "age_clq3"] + .[[2]]["treatment", "age_clq3"]*2)),
  seTE_q4 = map_dbl(poisson_cat, ~ sqrt(.[[2]]["treatment", "treatment"] + .[[2]]["age_clq4", "age_clq4"] + .[[2]]["treatment", "age_clq4"]*2)),
  treat1 = substr(comp, 1, 1),
  treat2 = substr(comp, 2, 2),
  poisson = NULL,
  studlab = paste0(comp, essai),
  essai = NULL,
  comp = NULL
) %>% group_by(scenario, ttt, sig, tau) %>% nest()


ipd_nets <- ipd_2steps_models %>% mutate(
  nets = map(data, ~netmeta_ipd( 
    TE = .$TE, seTE = .$seTE, 
    TE_q1 = .$TE_q1, seTE_q1 = .$seTE_q1,
    TE_q2 = .$TE_q2, seTE_q2 = .$seTE_q2,
    TE_q3 = .$TE_q3, seTE_q3 = .$seTE_q3,
    TE_q4 = .$TE_q4, seTE_q4 = .$seTE_q4,
    treat1 = .$treat1, treat2 = .$treat2, studlab = .$studlab)),
  data = NULL
)

ipd_nets <- ipd_nets %>% mutate(
  net_lin = map(nets, "lin"),
  net_q1 = map(nets, "q1"),
  net_q2 = map(nets, "q2"),
  net_q3 = map(nets, "q3"),
  net_q4 = map(nets, "q4"),
  nets = NULL
)

results$ipd_nets <- ipd_nets




# Exploitation en MA ================================================================================

## Fusion pour former les MA ------------------------------------------------------------------------

ipd_ma  <- sim_df  %>% group_by(scenario, ttt, sig, tau, comp) %>% summarise(ma = list(bind_rows(poisson_df, .id = "essai")))



## calcul des modeles

ipd_ma_models <- ipd_ma %>% select(- ma)

ipd_ma_models$res <- map(ipd_ma$ma, poisson_lin_cat_glmer_ma)

ipd_ma_models <- ipd_ma_models %>% mutate(
  poisson_lin = map(res,"lin"),
  poisson_cat = map(res,"cat"),
  res = NULL
)


results$ipd_ma_models <- ipd_ma_models





# Exploitation en NMA ===============================================================================

## Fusion pour former les NMA -----------------------------------------------------------------------
ipd_nma  <- ipd_ma %>% fuse_nma()


## Creation du modele aleatoire (parametrisé selon Freeman)


ipd_1step_models <- ipd_nma %>% select(- nma)

ipd_1step_models$res <- map(ipd_nma$nma, poisson_lin_cat_glmer_nma)


ipd_1step_models <- ipd_1step_models %>% mutate(
  poisson_lin = map(res,"lin"),
  poisson_cat = map(res,"cat"),
  res = NULL
)


results$ipd_1step_models <- ipd_1step_models

## retour des résultats

return(results)

}