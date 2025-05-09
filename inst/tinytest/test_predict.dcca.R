data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
# use vegan::rda in step 2
divide <- FALSE # divide by site.totals if TRUE

Y <- dune_trait_env$comm[, -1]  # must delete "Sites"
# delete "Species", "Species_abbr" from traits and
# use all remaining variables due to formulaTraits = ~. (the default)
traits <- dune_trait_env$traits
envir <- dune_trait_env$envir

# test a normal model
mod1a <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
               formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
               response = Y, 
               dataEnv = envir,
               dataTraits = traits,
               divideBySiteTotals = divide, 
               verbose = FALSE)

# checks that a factor or character in newdata does need 
# to have all levels of the training data
CWMfit1a <- fitted(mod1a, type = "CWM")
dat <- mod1a$data$dataEnv[17, ]
dat$Mag <- levels(dat$Mag)[dat$Mag] # factor to chr!
CWMfit1a17 <- predict(mod1a, type = "CWM", newdata = dat)

expect_equal(CWMfit1a17[1, ], CWMfit1a[17, ])

dat1 <- mod1a$data$dataEnv
dat1$Mag <- factor(dat1$Mag, levels = rev(levels(dat1$Mag)))
expect_equal(CWMfit1a, predict(mod1a, type = "traitsFromEnv", newdata = dat1))
# check predict scores (lc and lc_traits) with full scores
sc <- scores(mod1a, choices = 1:2, display = "lc", scaling = "sites")
sc17  <- predict(mod1a, type = "lc", rank = 2, scaling = "sites", newdata = dat)
expect_equal(sc17[1, ], sc[17, ])

sc <- scores(mod1a, choices = 1:2, display = "lc_traits", scaling = "sites")
dat2 <- mod1a$data$dataTraits[3, , drop = FALSE]
sc3 <- predict(mod1a, type = "lc_traits", rank = 2, scaling = "sites", 
               newdata = dat2)
expect_equal(sc3[1, ], sc[3, ])

# show that prediction of the mean gives mean environment 
#(for factors, weighted by weights of species, weight = sum over factor level)
mm_t <- t(mod1a$c_traits_normed0[, "Avg"])[, -5]
meanpredict_an <- predict(mod1a, type = "envFromTraits", 
                          newdata = data.frame(t(mm_t), Lifespan = "annual"))
meanpredict_per <- predict(mod1a, type = "envFromTraits",
                           newdata = data.frame(t(mm_t), Lifespan = "perennial"))
a1 <- sum(mod1a$weights$columns[mod1a$data$dataTraits[, "Lifespan"] == "annual"])
a2 <- sum(mod1a$weights$columns[mod1a$data$dataTraits[, "Lifespan"] == "perennial"])
mm_e <- (a1 * meanpredict_an[1, ] + a2 * meanpredict_per[1, ]) / (a1 + a2)
expect_equivalent(mm_e, 
                  attr(scores(mod1a, display = "bp", normed = FALSE), 
                       which = "mean"))

expect_warning(meanpredict_an2 <- predict(mod1a, type = "envFromTraits", 
                                          newdata = as.data.frame(t(mm_t))), 
               "newdata does not")

# missing variables are set to their mean (an3 and an4) and 
# in an3: variable mm_t is not in the model and is thus ignored.
expect_warning(meanpredict_an3 <- 
                 predict(mod1a, type = "envFromTraits", 
                         newdata = data.frame(mm_t, Lifespan = "annual" )), 
               "newdata does not")

expect_warning(meanpredict_an4 <- 
                 predict(mod1a, type = "envFromTraits",
                         newdata = data.frame(Lifespan = "annual")), 
               "newdata does not")

expect_equivalent(meanpredict_an2[1, ], meanpredict_an[1, ]) 
expect_equivalent(meanpredict_an3[1, ], meanpredict_an[1, ]) 
expect_equivalent(meanpredict_an4[1, ], meanpredict_an[1, ])

# five types of full rank prediction:
predDivF1a_reg <- coef(mod1a, type = "env2traits")
predDivF1a_reg_traits <- coef(mod1a, type = "traits2env")
predDivF1a_env <- predict(mod1a, type = "env")
predDivF1a_envFitted <- fitted(mod1a, type = "SNC")
expect_equal(predDivF1a_env, predDivF1a_envFitted)

predDivF1a_traits <- fitted(mod1a, type = "CWM")
predDivF1a_traits <- predict(mod1a, type = "traits")
predDivF1a_response <- predict(mod1a, type = "response",
                               newdata = list(mod1a$data$dataTraits,
                                              mod1a$data$dataEnv))

expect_equal_to_reference(predDivF1a_reg, "predDivF1a_reg") 
expect_equal_to_reference(predDivF1a_reg_traits, "predDivF1a_reg_traits") 
expect_equal_to_reference(predDivF1a_env, "predDivF1a_env") 
expect_equal_to_reference(predDivF1a_traits, "predDivF1a_traits") 
expect_equal_to_reference(predDivF1a_response, "predDivF1a_response") 

# five type of rank 1 predDivFiction
predDivF1a_reg1 <- coef(mod1a, type = "env2traits", rank = 1)
predDivF1a_reg_traits1 <- coef(mod1a, type = "traits2env", rank = 1)
predDivF1a_env1 <- predict(mod1a, type = "env", rank = 1)
predDivF1a_traits1 <- predict(mod1a, type = "traits", rank = 1)
predDivF1a_response1 <- predict(mod1a, type = "response", rank = 1,
                                newdata = list(mod1a$data$dataTraits,
                                               mod1a$data$dataEnv))

expect_equal_to_reference(predDivF1a_reg1, "predDivF1a_reg1") 
expect_equal_to_reference(predDivF1a_reg_traits1, "predDivF1a_reg_traits1") 
expect_equal_to_reference(predDivF1a_env1, "predDivF1a_env1") 
expect_equal_to_reference(predDivF1a_traits1, "predDivF1a_traits1") 
expect_equal_to_reference(predDivF1a_response1, "predDivF1a_response1") 

expect_equal(predict(mod1a, type = "env", rank = 1), 
             fitted(mod1a, type = "SNC", rank = 1))
expect_equal(predict(mod1a, type = "traits", rank = 1), 
             fitted(mod1a, type = "CWM", rank = 1))
expect_equal(
  predict(mod1a, type = "response", rank = 1, weights = mod1a$weights) *
    sum(mod1a$data$Y), 
  fitted(mod1a, type = "response", rank = 1))

# test factors only
envir$fUse <- factor(envir$Use)
traits$fSLA <- cut(traits$SLA, 3)
mod_fact<- dc_CA(formulaEnv = ~ fUse + Mag,
                 formulaTraits = ~ fSLA + Lifespan ,
                 response = Y, 
                 dataEnv = envir, 
                 dataTraits = traits,
                 divideBySiteTotals = divide, 
                 verbose = FALSE)

expect_equal(colnames(predict(mod_fact, type = "env")), 
             c("fUse1", "fUse2", "fUse3", "MagSF", "MagBF", "MagHF", "MagNM"))

expect_equal(colnames(predict(mod_fact, type = "traits")),
             c( "fSLA(5.77,15]", "fSLA(15,24.2]", "fSLA(24.2,33.4]", 
                "Lifespanannual", "Lifespanperennial"))

traits$lifequant <- 1 * traits$Lifespan %in% "perennial"
mod_factb<- dc_CA(formulaEnv = ~ fUse + Mag,
                  formulaTraits = ~ fSLA + lifequant ,
                  response = Y, 
                  dataEnv = envir, 
                  dataTraits = traits,
                  divideBySiteTotals = divide,
                  verbose = FALSE)

expect_silent(mod_fact)
expect_equal(mod_fact$eigenvalues, mod_factb$eigenvalues)
expect_equivalent(mod_fact$c_traits_normed0, mod_factb$c_traits_normed0)

# test quant variable only, 1 quant trait

modq11<- dc_CA(formulaEnv = ~ Manure,
               formulaTraits = ~ SLA,
               response = Y, 
               dataEnv = envir, 
               dataTraits = traits,
               divideBySiteTotals = divide,
               verbose = FALSE)

modq11_predDivFenv <- predict(modq11, type = "env")
expect_equal_to_reference(modq11_predDivFenv, "modq11_predDivFenv")

modq11_predDivFtraits <- predict(modq11, type = "traits")
expect_equal_to_reference(modq11_predDivFtraits, "modq11_predDivFtraits")

# check how douconca manages collinear models A
envir$Sites <- factor(dune_trait_env$envir$Sites)
mod_dccaFF <- dc_CA(formulaEnv = ~ Sites,
                    formulaTraits = ~ Species,
                    response = dune_trait_env$comm[, -1],  # must delete "Sites"
                    dataEnv = envir,
                    dataTraits = traits,
                    divideBySiteTotals = divide,
                    verbose = FALSE)

Yfit <- fitted(mod_dccaFF, type = "response")
expect_equivalent(Yfit, as.matrix(mod_dccaFF$data$Y)) # full rank fit!

betaT <- rep(NA, 28)
for (k in 1:28) {
  betaT[k] <- coef(lm(Yfit[, k] ~ mod_dccaFF$data$Y[, k]), 
                   weights = mod_dccaFF$weights$rows)[2]
}
expect_equal(betaT, rep(1, length(betaT)))

beta <- rep(NA, 20)
for (k in 1:20) {
  beta[k] <- coef(lm(Yfit[k, ] ~ mod_dccaFF$data$Y[k, ]), 
                  weights = mod_dccaFF$weights$columns)[2]
}
expect_equal(beta, rep(1,length(beta)))

# CWM perfect fit

envir$Sites <- factor(dune_trait_env$envir$Sites)
mod_dccaCWM <- dc_CA(formulaEnv = ~ Sites,
                     formulaTraits = ~ F + R + N + L,
                     response = dune_trait_env$comm[, -1],  # must delete "Sites"
                     dataEnv = envir,
                     dataTraits = traits,
                     divideBySiteTotals = divide,
                     verbose = FALSE)

expect_equal(mod_dccaCWM$inertia["traits_explain", 1], 
             mod_dccaCWM$inertia["constraintsTE", 1])
expect_equal(mod_dccaCWM$inertia["total", 1], 
             mod_dccaCWM$inertia["env_explain", 1])

CWMSNC <- fCWM_SNC(formulaEnv = ~ Sites,
                   formulaTraits = ~ F + R + N + L,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv = envir,
                   dataTraits = traits,
                   divideBySiteTotals = divide)

CWMfit <- fitted(mod_dccaCWM, type = "CWM")
expect_equivalent(CWMfit, CWMSNC$CWM)# full rank fit! PASSED without the factor sqrt(20)

# SNC perfect fit 

mod_dccaSNC <- dc_CA(formulaEnv = ~ Moist+ Mag,
                     formulaTraits = ~ Species,
                     response = dune_trait_env$comm[, -1],  # must delete "Sites"
                     dataEnv = envir,
                     dataTraits = traits,
                     divideBySiteTotals = divide,
                     verbose = FALSE)

expect_equal(mod_dccaSNC$inertia["env_explain", 1], 
             mod_dccaSNC$inertia["constraintsTE", 1])
expect_equal(mod_dccaSNC$inertia["total", 1], 
             mod_dccaSNC$inertia["traits_explain", 1])

CWMSNCb <- fCWM_SNC(formulaEnv = ~ Moist + Mag,
                    formulaTraits = ~ Species,
                    response = dune_trait_env$comm[, -1],  # must delete "Sites"
                    dataEnv = envir,
                    dataTraits = traits,
                    divideBySiteTotals = divide)

SNCfit <- fitted(mod_dccaSNC, type = "SNC")
expect_equivalent(SNCfit[, -2], CWMSNCb$SNC)# Passed Without the factor sqrt(28)!

# check how douconca manages collinear models

#  A11 collinear variable
envir$Sites <- factor(dune_trait_env$envir$Sites)
envir$A11 <- envir$A1

mod_DivF_dccaA1 <- dc_CA(formulaEnv = ~ A1 + Manure + Moist,
                         formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                         response = Y, 
                         dataEnv = envir,
                         dataTraits = traits,
                         divideBySiteTotals = divide,
                         verbose = FALSE)
mod_DivF_dccaA11 <- 
    dc_CA(formulaEnv = ~ A1 + A11 + Manure + Moist,
          formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
          response = Y, 
          dataEnv = envir,
          dataTraits = traits,
          divideBySiteTotals = divide,
          verbose = FALSE)
expect_equal(predict(mod_DivF_dccaA1, type = "response"),
             predict(mod_DivF_dccaA11, type = "response"))

# SLA11 collinear
traits$SLA11 <- traits$SLA
mod_DivF_dccaSLA <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
                          formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                          response = Y, 
                          dataEnv = envir, 
                          dataTraits = traits,
                          divideBySiteTotals = divide,
                          verbose = FALSE)
expect_warning(mod_DivF_dccaSLA11 <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
                              formulaTraits = ~ SLA + SLA11+ Height + LDMC + Seedmass + Lifespan,
                              response = Y, 
                              dataEnv = envir,
                              dataTraits = traits,
                              divideBySiteTotals = divide,
                              verbose = FALSE),"Collinearity detected in CWM-model")

expect_equal(predict(mod_DivF_dccaSLA, type = "response"),
             predict(mod_DivF_dccaSLA11, type = "response"))
