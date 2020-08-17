# Fonction de génération de survie -----------------------------------------------------------------------

  #Définie dans 02-Simulate survival
survsim <- function(n, intercept, scale, random_intercept = NULL,
                    covars = NULL, params = NULL, inter = NULL, 
                    betas = NULL, random = NULL, 
                    cens_intercept = NULL, cens_scale = 1, suivi_max = NULL) {
  
  
  # Generation des co-variables --------------------------------
  
  # fonction sample maison pour plus de cohérence
  sample2 <- function(ncat, prob = NULL, n) as.integer(sample.int(n = ncat, prob = prob, replace = TRUE, size = n) - 1)
  covars <- gsub( "sample", "sample2", covars, fixed = TRUE)
  
  
  # génération des variables aléatoires
  if (! is.null(params)) {
    x <- covars %>% set_names(paste0("X", 1:length(covars))) %>% invoke_map_df(params, n = n)
  }
  
  x_out <- x
  
  
  # gestion des interactions
  if(! is.null(inter)) {
    # int_final <- inter %>% set_names(paste0("INTER", 1:length(inter))) %>% map_df(~ apply(x[.], 1, prod))
    int_final <- inter %>% set_names(paste0("INTER", 1:length(inter))) %>% 
      map_df(~ apply(x[.], 1, function(y) (y[1] - params[[1]]$mean)*y[2]))
    x <- cbind(x, int_final)
  }
  
  # Generation des effets aleatoires -----------------------------
  if (! is.null(random)) {
    re <- map_dbl(random, ~ betas[.[1]] + rnorm(1, mean = 0, sd = .[2]))
    betas[map_dbl(random, 1)] <- re 
  }
  
  # Génération des predicteurs linéaires ------------------------
  x$lp <- apply(x, 1, function(y) sum(y * betas))
  
  # Generation des délais de survies ---------------------------
  ## Parametrisation lineaire
  
  if(! is.null(random_intercept)) intercept <- intercept + rnorm(1, 0, random_intercept)
  
  W <- log(-log(runif(n, 0, 1)))
  
  time <- exp(intercept + scale * (W - x$lp))
  
  ## Parametrisation classique
  # shape <- 1 / scale
  # scale <- 1 / exp(intercept)
  
  # W <- runif(n, 0, 1)
  # 
  # time <- ((- log(W)) /(scale * exp(x$lp)) )^ (1 / shape)
  
  # Génération des censures ------------------------------------  
  if(!is.null(suivi_max) | !is.null(cens_intercept)) { 
    
    if(is.null(suivi_max)) suivi_max<- max(time) + 1
    
    if(! is.null(cens_intercept)) {
      C <- rweibull(n = n, shape = 1/cens_scale, scale = exp(cens_intercept))
      C <- pmin(C, suivi_max)
    } else {
      C <- suivi_max
    }
    
    status <- (C >= time) + 0
    time <- pmin(time, C)
    
  } else {
    status <- 1
  }
  
  # optimisation de l'espace mémoire
  d <- as_tibble(cbind(time, status, x_out))
  d$time <- as.integer(ceiling(time)) ## le ceiling s'assure que les survies ne peuvent être à 0
  d$status <- as.logical(d$status)
  
  # renvoie du tibble final -----------------------------------
  
  return(d)
}


# Fonction pour la génération des data -------------------------------------------------------------------

survsim2 <- function(ttt, sig, tau, comp, n = 200, scenario) {
  
  # Determination des paramètres de simulation ===================================
  
  # Initialisation pour le cas "None"
  age <- 60
  interaction <- 0
  ttt <- ttt
  effet_age <- 0.011  ##  (on vise effet de l'age HR = 1.25 de HR en passant du 10eme au 90 eme percentile)
  
  # Differenciation des parametres si scenario particulier
  
    ## Age
    if(scenario %in% c("age", "both")) {
      age <- ifelse(comp == "AC",  60 - 8, 60 + 8)
      }
    
    ## Interaction
    if(scenario %in% c("interaction", "both") & comp == "AC")  {
      interaction <-  -(ttt * 0.25) / 8 # variation de 25% de l'effet par écart type
      
      effet_age <- effet_age - 1/2 * (interaction) ## puisque 1 malade sur deux est traité
      
      ## il faut prendre en compte que l'effet est modifié par l'interaction pour avoir un effet moyen
      if (scenario == "both") {
              ttt <- ttt + ttt * 0.25 ## comme l'age moyen va être différent on augmente l'effet traitement pour avoir celui d'une population de 52 ans donc 25% plus grand
      }
    }
    
    ## Essai AB  
    if(comp == "AB") {
        if(scenario %in% c("interaction", "both"))  {
          interaction <- (-ttt * 0.25) / ( 8)
          effet_age <- effet_age - 1/2 * (interaction) ## puisque 1 malade sur deux est traité
        }
        age <- 60       
        ttt <- 0

      }
  
  
  # Génération des données ======================================================
    
  out <- survsim(n = n, intercept = 8, scale = 0.73, random_intercept = sig,
                 covars = c("rnorm", "sample"), 
                 params = list(list(mean = age, sd = 8), list(ncat = 2, prob = c(0.5,0.5))),
                 inter = list(c(1, 2)),
                 betas = c(effet_age, ttt, interaction),  
                 random = list(c(2, tau)),
                 cens_intercept = 8.5, cens_scale = 0.35,
                 suivi_max = 365 * 8)
  names(out)[3:4] <- c("age", "treatment")
  return(out)
}


survsim3 <- function(ttt, sig, tau, comp, scenario) {
  
  out <- survsim2(ttt = ttt, sig = sig, tau = tau, comp = comp, scenario = scenario)
  
  ## on reclasse et on prepare le poisson
  out <- cut_center(out)
  out <- survSplit2(Surv(time, status) ~ ., data = out, cut = seq(0, 2548, by = 182*2), center_periode = TRUE)
  
  return(out)
}


survsim4 <- function(df) {
  
  df <- cut_center(df)
  df <- survSplit2(Surv(time, status) ~ ., data = df, cut = seq(0, 2548, by = 182*2), center_periode = TRUE)
  
  return(df)
}

cut_center <- function(x) {
  x$age_cl <- cut(x$age, breaks = c(1, 55, 60, 65, 150), labels = paste0("q", 1:4))
  x$age <- (x$age - 60)## on centre par la vraie valeur de la population pour que les résultats soient comparables
  x$treatment <- x$treatment -0.5
  return(x)
}


fuse_nma <- function(sim_data) {
  out <- 
    sim_data %>%  
    transmute(nma = map2(ma, comp, ~ bind_cols(.x, "comp" = rep(.y, nrow(.x))))) %>% 
    summarize(nma = list(bind_rows(nma))) %>% 
    mutate(nma = map(nma, create_treatment))
  
  return(out)
}


create_treatment <- function(nma) {
  
  ## Codage type Freeman et Carpenter
  nma$ttt1 <- 1 * (nma$comp == "AC") * (nma$treatment) + 1 * (nma$comp == "AB") * (nma$treatment)
  nma$ttt2 <- 1 * (nma$comp == "BC") * (nma$treatment) - 1 * (nma$comp == "AB") * (nma$treatment)
  
  nma$essai <- paste0(nma$essai, nma$comp)
  return(nma)
}


# Fonction de model maison pour ne renvoyer qu'une partie des résultats

coxph2 <- function(...) {
  mod <- coxph(...,  x = FALSE, y = FALSE, eps = 10e-6)
  out <- list(coef = coef(mod), var = diag(vcov(mod)))
  return(out)
}

glm2 <- function(...) {
  mod <- glm(... , model = FALSE, x = FALSE, y = FALSE)
  out <- list(coef = coef(mod), var = vcov(mod))
  return(out)  
}


coxme2 <- function(...) {
  mod <- coxme(...,  x = FALSE, y = FALSE)
  out <- list(fe = fixef(mod), fe_se = vcov(mod), re = VarCorr(mod))
  return(out)  
}

glmer2 <- function(...) {
  mod <- glmer(..., nAGQ = 0,
                 control = glmerControl(
                   calc.derivs = FALSE,
                   check.nobs.vs.rankZ = "ignore",
                   check.nobs.vs.nlev = "ignore",
                   check.nlev.gtreq.5 = "ignore",
                   check.nlev.gtr.1 = "ignore",
                   check.nobs.vs.nRE="ignore",
                   check.rankX = "ignore",
                   check.scaleX = "ignore",
                   check.formula.LHS = "ignore",
                   check.response.not.const = "ignore"
                 )
               )
  out <- list(fe = fixef(mod), fe_se = vcov(mod), re = VarCorr(mod))
  return(out)    
}

netmeta2 <- function(...) {
  mod <- netmeta(...)
  out <- list(TE.fixed = mod$TE.fixed, seTE.fixed = mod$seTE.fixed, TE.random = mod$TE.random, seTE.random = mod$seTE.random)
}


### Fonctions pour simplifier l'écriture dans la simulation

netmeta_ipd <- function(TE, TE_q1, TE_q2, TE_q3, TE_q4, seTE, seTE_q1, seTE_q2, seTE_q3, seTE_q4, treat1, treat2, studlab) {
  lin <- netmeta2(TE = TE, seTE = seTE, treat1 = treat1, treat2 = treat2, studlab = studlab,  sm = "HR")
  q1  <- netmeta2(TE = TE_q1, seTE = seTE_q1, treat1 = treat1, treat2 = treat2, studlab = studlab,  sm = "HR")
  q2  <- netmeta2(TE = TE_q2, seTE = seTE_q2, treat1 = treat1, treat2 = treat2, studlab = studlab,  sm = "HR")
  q3  <- netmeta2(TE = TE_q3, seTE = seTE_q3, treat1 = treat1, treat2 = treat2, studlab = studlab,  sm = "HR")
  q4  <- netmeta2(TE = TE_q4, seTE = seTE_q4, treat1 = treat1, treat2 = treat2, studlab = studlab,  sm = "HR")
  
  out <- list(lin = lin, q1 = q1, q2 = q2, q3 = q3, q4 = q4)
  return(out)
}

poisson_lin_cat_glm <- function(donnees) {
  lin <- glm2(status ~ ns(periode, 3) + offset(offset) + treatment * age, family = "poisson", data = donnees)
  cat <- glm2(status ~ ns(periode, 3) + offset(offset) + treatment * age_cl, family = "poisson", data = donnees)
  return(list(lin = lin, cat = cat))
}

poisson_lin_cat_glmer_ma <- function(donnees) {
  lin <- glmer2(status ~ ns(periode, 3) + offset(offset)  + treatment * age  + (1 | essai) + (0 + treatment|essai), family = "poisson" , data = donnees)
  cat <- glmer2(status ~ ns(periode, 3) + offset(offset)  + treatment * age_cl + (1 | essai) + (0 + treatment|essai), family = "poisson", data = donnees)
  return(list(lin = lin, cat = cat))
}

poisson_lin_cat_glmer_nma <- function(donnees) {
  lin <- glmer2(status ~ ns(periode, 3) + offset(offset) + age  + ttt1 * age + ttt2 * age + (1 | essai) + (0 + ttt1|essai) + (0 + ttt2|essai), data = donnees, family = "poisson")
  cat <- glmer2(status ~ ns(periode, 3) + offset(offset) + age_cl  + ttt1 * age_cl + ttt2 * age_cl + (1 | essai) + (0 + ttt1|essai) + (0 + ttt2|essai), data = donnees, family = "poisson")
  return(list(lin = lin, cat = cat))
}


poisson_lin_cat_glmer_nma2 <- function(donnees) {
  lin <- glmer2(status ~ ns(periode, 3) + offset(offset) + age  + ttt1 * age + ttt2 * age + (1 | comp/essai) + (0 + ttt1 + ttt1:age|essai) + (0 + ttt2 + ttt2:age|essai), data = donnees, family = "poisson")
  cat <- glmer2(status ~ ns(periode, 3) + offset(offset) + age_cl  + ttt1 * age_cl + ttt2 * age_cl + (1 | essai) + (0 + ttt1|essai) + (0 + ttt2|essai), data = donnees, family = "poisson")
  return(list(lin = lin, cat = cat))
}


# Fonction de split maison ---------------------------------------------------------------------------------

survSplit2 <- function (formula, data, subset, na.action = na.pass, cut, start = "tstart", 
                        id, zero = 0, episode, end = "tstop", event = "event", center_periode = FALSE) 
{
  Call <- match.call()
  # if (missing(formula) || is.data.frame(formula)) {
  #   if (missing(data)) {
  #     if (!missing(formula)) {
  #       names(Call)[[2]] <- "data"
  #       data <- formula
  #     }
  #     else stop("a data frame is required")
  #   }
  #   if (missing(end) || missing(event)) 
  #     stop("either a formula or the end and event arguments are required")
  #   if (!(is.character(event) && length(event) == 1 && event %in% 
  #         names(data))) 
  #     stop("'event' must be a variable name in the data set")
  #   if (!(is.character(end) && length(end) == 1 && end %in% 
  #         names(data))) 
  #     stop("'end' must be a variable name in the data set")
  #   if (!(is.character(start) && length(start) == 1)) 
  #     stop("'start' must be a variable name")
  #   if (start %in% names(data)) 
  #     temp <- paste(start, end, event, sep = ",")
  #   else temp <- paste(end, event, sep = ",")
  #   formula <- as.formula(paste("Surv(", temp, ")~ ."))
  # }
  # else if (missing(formula)) 
  #   stop("either a formula or the end and event arguments are required")
  indx <- match(c("data", "weights", "subset"), names(Call), 
                nomatch = 0)
  temp <- Call[c(1L, indx)]
  temp$formula <- formula
  temp$na.action <- na.action
  temp[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(temp)
  Y <- model.response(mf)
  states <- attr(Y, "states")
  # if (!is.Surv(Y)) 
  #   stop("the model must have a Surv object as the response")
  # if (!(attr(Y, "type") %in% c("right", "mright", "counting", 
  #                              "mcounting"))) 
  #   stop(paste("not valid for", attr(Y, "type"), "censored survival data"))
  nY <- ncol(Y)
  if (nY == 2) 
    Y <- cbind(zero, Y)
  temp <- (Y[, 1] >= Y[, 2])
  # if (any(temp & !is.na(temp))) 
  #   stop("start time must be < stop time")
  # if (!is.numeric(cut) || any(!is.finite(cut))) 
  #   stop("cut must be a vector of finite numbers")
  cut <- sort(cut)
  ntimes <- length(cut)
  n <- nrow(data)
  # if (!missing(id)) {
  #   if (!is.character(id)) 
  #     stop("id must be a variable name")
  #   if (id %in% names(mf)) 
  #     stop("the suggested id name is already present")
  #   id <- make.names(id)
  #   if (id %in% names(mf)) 
  #     stop("the suggested id name is already present")
  #   mf[[id]] <- 1:nrow(mf)
  # }
  storage.mode(Y) <- "double"
  index <- .Call(survival:::Csurvsplit, Y[, 1], Y[, 2], as.double(cut))
  newdata <- mf[index$row, -1, drop = FALSE]
  row.names(newdata) <- NULL
  attr(newdata, "terms") <- NULL
  status <- Y[index$row, 3]
  status[index$censor] <- 0
  if (!is.null(states)) 
    status <- factor(status, labels = c("censor", states))
  if (class(formula[[2]]) == "call" && formula[[2]][[1]] == 
      as.name("Surv")) {
    temp <- match.call(Surv, formula[[2]])
    if (nY == 2) {
      if (missing(end) && !is.null(temp[["time"]]) && 
          is.name(temp[["time"]])) 
        end <- as.character(temp[["time"]])
      if (missing(event) && !is.null(temp$time2) && is.name(temp$time2)) 
        event <- as.character(temp$time2)
      if (missing(event) && !is.null(temp$event) && is.name(temp$event)) 
        event <- as.character(temp$event)
    }
    else {
      if (missing(end) && !is.null(temp[["time"]]) && 
          is.name(temp["time"])) 
        start <- as.character(temp[["time"]])
      if (missing(end) && !is.null(temp$time2) && is.name(temp$time2)) 
        end <- as.character(temp$time2)
      if (missing(event) && !is.null(temp$event) && is.name(temp$event)) 
        event <- as.character(temp$event)
      if (missing(start) && !is.null(temp$time) && is.name(temp$time)) 
        start <- as.character(temp$time)
    }
    newdata[[start]] <- as.integer(index$start)
    newdata[[end]] <- as.integer(index$end)
    newdata[[event]] <- as.logical(status)
  }
  else {
    # if (class(formula[[2]]) != "name") 
    #   stop("left hand side not recognized")
    temp <- as.character(formula[[2]])
    newdata[temp] <- Surv(index$start, index$end, status)
  }
  # if (!missing(episode)) {
  #   # if (!is.character(episode)) 
  #   #   stop("episode must be a character string")
  #   newdata[[make.names(episode)]] <- index$interval + 1
  # }
  # newdata$periode <- cut[newdata[[episode]] - 1]
  if (! center_periode)  newdata$periode <- cut[index$interval]
  
  newdata$periode <- as.integer(index$interval - ceiling(ntimes/2))

  newdata$offset <- log(newdata[[as.character(formula[[2]])[2]]] - newdata$tstart + 0.5 )
  newdata
}



# Manipulation de matrice de variance/covariance
covi <- function(mat, variables, pond) {
  
  for(i in seq_along(pond)) {
    if (pond[i] != 0) {
      mat[variables[i], ] <- mat[variables[i], ] * pond[i]
      mat[, variables[i]] <- mat[, variables[i]] * pond[i]
    }
  }
  
  sum(mat[variables, variables])
}
