data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
# use vegan::rda in step 2
divide <- TRUE # divide by site.totals if TRUE

Y <- dune_trait_env$comm[, -1]  # must delete "Sites"
# delete "Species", "Species_abbr" from traits and
# use all remaining variables due to formulaTraits = ~. (the default)
traits <- dune_trait_env$traits
envir <- dune_trait_env$envir

# test a normal model
mod_DivT1a <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
                    formulaTraits = ~ SLA + Height + LDMC + Seedmass +Lifespan,
                    response = Y, 
                    dataEnv = envir, 
                    dataTraits = traits,
                    divideBySiteTotals = divide,
                    verbose = FALSE)

# Because singular vectors have arbitrary sign, 
# the absolute value of scores is tested.
# The test is an implicit test on the singular values and eigen values as well.
scores_mod_DivT1a_abs <- sapply(X = scores(mod_DivT1a), FUN = abs)
expect_equal_to_reference(scores_mod_DivT1a_abs, "scores_mod_DivT1a_abs")

mod_DivT1a1 <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
                     formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                     response = Y, 
                     dataEnv = envir, 
                     dataTraits = traits, use_vegan_cca = TRUE,
                     divideBySiteTotals = divide, 
                     verbose = FALSE)

scores_mod_DivT1a1_abs <- sapply(X = scores(mod_DivT1a1), FUN = abs)
expect_equal(scores_mod_DivT1a_abs,scores_mod_DivT1a1_abs)


# test all variables(~.) in formulas
mod_DivT1b <- dc_CA(formulaEnv = ~ .,
                    formulaTraits = ~ .,
                    response = Y, 
                    dataEnv = envir[, -c(1, 7:10)], 
                    dataTraits = traits[, -c(1, 2, 8:11)],
                    divideBySiteTotals = divide, 
                    verbose = FALSE)

expect_equal(mod_DivT1a$eigenvalues, mod_DivT1b$eigenvalues)
expect_equal(mod_DivT1a$c_traits_normed0, mod_DivT1b$c_traits_normed0)

# model with covariates

mod_DivTc <- dc_CA(formulaEnv = ~ A1 + Moist + Manure + Use + Condition(Mag),
                   formulaTraits = ~ SLA + Height + LDMC + Condition(Seedmass) + Lifespan,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv = dune_trait_env$envir,
                   dataTraits = dune_trait_env$traits,
                   divideBySiteTotals = divide,
                   verbose = FALSE)

Y <- mod_DivTc$data$Y

mod_DivTc_vegan <- vegan::cca(formula = Y ~ A1 + Moist + Manure + Use + Condition(Mag),
                              data =  mod_DivTc$data$dataEnv)
expect_equal(mod_DivTc$inertia["env_explain", 1],
             sum(vegan::eigenvals(mod_DivTc_vegan, model = "constrained")))

scores_mod_DivTc_abs <- sapply(X = scores(mod_DivTc), FUN = abs)
expect_equal_to_reference(scores_mod_DivTc_abs, "scores_mod_DivTc_abs")

mod_DivTcV <- dc_CA(formulaEnv = ~ A1 + Moist + Manure + Use + Condition(Mag),
                    formulaTraits = ~ SLA + Height + LDMC + Condition(Seedmass) + Lifespan,
                    response = dune_trait_env$comm[, -1],  # must delete "Sites"
                    dataEnv = dune_trait_env$envir,
                    dataTraits = dune_trait_env$traits,
                    divideBySiteTotals = divide, use_vegan_cca = TRUE,
                    verbose = FALSE)
scores_mod_DivTcV_abs <- sapply(X = scores(mod_DivTcV), FUN = abs)
expect_equal(scores_mod_DivTc_abs,scores_mod_DivTcV_abs)

# test factors only
envir$fUse <- factor(envir$Use)
traits$fSLA <- cut(traits$SLA, 3)

mod_DivT_fact<- dc_CA(formulaEnv = ~ fUse + Mag,
                      formulaTraits = ~ fSLA + Lifespan,
                      response = Y, 
                      dataEnv = envir, 
                      dataTraits = traits,
                      divideBySiteTotals = divide, 
                      verbose = FALSE)

scores_mod_DivT_fact_abs <- sapply(X = scores(mod_DivT_fact), FUN = abs)
expect_equal_to_reference(scores_mod_DivT_fact_abs, "scores_mod_DivT_fact_abs")

traits$lifequant <- 1 * traits$Lifespan %in% "perennial"

mod_DivT_factb <- dc_CA(formulaEnv = ~ fUse + Mag,
                        formulaTraits = ~ fSLA + lifequant,
                        response = Y, 
                        dataEnv = envir, 
                        dataTraits = traits,
                        divideBySiteTotals = divide,
                        verbose = FALSE)

expect_silent(mod_DivT_fact)
expect_equal(mod_DivT_fact$eigenvalues, mod_DivT_factb$eigenvalues)
expect_equivalent(mod_DivT_fact$c_traits_normed0, mod_DivT_factb$c_traits_normed0)

# test quant variable only, 1 quant trait

mod_DivTq11<- dc_CA(formulaEnv = ~ Manure,
                    formulaTraits = ~ SLA,
                    response = Y, 
                    dataEnv = envir, 
                    dataTraits = traits,
                    divideBySiteTotals = divide,
                    verbose = FALSE)

scores_mod_DivTq11_abs <- sapply(X = scores(mod_DivTq11), FUN = abs)
expect_equal_to_reference(scores_mod_DivTq11_abs, "scores_mod_DivTq11_abs")

# check how douconca manages collinear models 

#  A11 collinear variable
envir$Sites <- factor(dune_trait_env$envir$Sites)
envir$A11 <- envir$A1

mod_DivT_dccaA1 <- dc_CA(formulaEnv = ~ A1 + Manure+Moist,
                         formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                         response = Y, dataEnv = envir, dataTraits = traits,
                         divideBySiteTotals = divide,
                         verbose = FALSE)
#expect_stdout(expect_warning(
  mod_DivT_dccaA11 <- dc_CA(formulaEnv = ~ A1 + A11 + Manure + Moist,
                            formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                            response = Y, 
                            dataEnv = envir, 
                            dataTraits = traits,
                            divideBySiteTotals = divide,
                            verbose = FALSE)#, "singular"))

expect_equal(mod_DivT_dccaA1$eigenvalues, mod_DivT_dccaA11$eigenvalues)

expect_equivalent(mod_DivT_dccaA1$inertia, 
                  mod_DivT_dccaA11$inertia)

expect_equal(abs(scores(mod_DivT_dccaA1, display = "tval_traits")), 
             abs(scores(mod_DivT_dccaA11, display = "tval_traits")))

expect_warning(scores(mod_DivT_dccaA11, display = "reg"), 
               "Collinearity detected")
expect_warning(scores(mod_DivT_dccaA11, display = "tval"), 
               "Collinearity detected")
# note that the VIFs differ due to NA in VIF for manure! so also unvalid tvalues

expect_warning(scores_A11 <- scores(mod_DivT_dccaA11, display = "reg"), 
               "Collinearity detected in CWM-model.")

expect_equal(abs(scores(mod_DivT_dccaA1, display = "reg"))[1:3, -3, drop = FALSE], 
             abs(scores_A11)[-2, -3, drop = FALSE])

# SLA11 collinear
traits$SLA11 <- traits$SLA
mod_DivT_dccaSLA <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
                          formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                          response = Y, 
                          dataEnv = envir,
                          dataTraits = traits,
                          divideBySiteTotals = divide,
                          verbose = FALSE)

expect_warning(
  mod_DivT_dccaSLA11 <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
                              formulaTraits = ~ SLA + SLA11+ Height + LDMC + Seedmass + Lifespan,
                              response = Y, 
                              dataEnv = envir,
                              dataTraits = traits,
                              divideBySiteTotals = divide,
                              verbose = FALSE), "Collinearity")

expect_equal(mod_DivT_dccaSLA$eigenvalues, mod_DivT_dccaSLA11$eigenvalues)
expect_equal(mod_DivT_dccaSLA$inertia, mod_DivT_dccaSLA11$inertia)

expect_equal(abs(scores(mod_DivT_dccaSLA, display = "tval")), 
             abs(scores(mod_DivT_dccaSLA11, display = "tval")))

expect_warning(scores_SLA11 <- scores(mod_DivT_dccaSLA11, 
                                      display = "reg_traits"), 
               "Collinearity")

expect_equal(abs(scores(mod_DivT_dccaSLA, display = "reg_traits"))[, -3], 
             abs(scores_SLA11)[-2, -3])

expect_warning(tval_SLA11 <- scores(mod_DivT_dccaSLA11, display = "tval_traits"),
               "Collinearity detected")

expect_true(all(is.na(tval_SLA11)))

# full fit
expect_stdout(mod_DivT_dcca_near_singular_sites <- 
                dc_CA(formulaEnv = ~ Sites,
                      formulaTraits = ~ SLA + Seedmass,
                      response = Y,
                      dataEnv = envir,
                      dataTraits = traits,
                      divideBySiteTotals = divide,
                      verbose = FALSE))

expect_true(all(scores(mod_DivT_dcca_near_singular_sites, display = "tval") == 0))

expect_stdout(mod_DivT_dcca_near_singular_species <- 
                dc_CA(formulaEnv = ~ A1 + Manure + Condition(Mag),
                      formulaTraits = ~ Species,
                      response = Y,
                      dataEnv = envir,
                      dataTraits = traits,
                      divideBySiteTotals = divide,
                      verbose = FALSE))

expect_message(dc_CA(formulaEnv = ~ A1 + Manure + Condition(Mag),
                     formulaTraits = ~ Species,
                     response = Y,
                     dataEnv = envir,
                     dataTraits = traits, use_vegan_cca = TRUE,
                     divideBySiteTotals = divide,
                     verbose = FALSE), 
               "The model is overfitted")

expect_true(all(is.nan(scores(mod_DivT_dcca_near_singular_species, 
                              display = "tval_traits"))))

#expect_stdout(expect_warning(
  mod_DivT_dcca_singular2 <- 
    dc_CA(formulaEnv = ~ Sites + A1,
          formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
          response = Y,
          dataEnv = envir,
          dataTraits = traits,
          divideBySiteTotals = divide, 
          verbose = FALSE)#, "singular"))

expect_warning(scores(mod_DivT_dcca_singular2),
               "Collinearity detected in CWM-model")		

expect_warning(
  mod_divT_dcca_singular3 <- 
    dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
          formulaTraits = ~ Species_abbr + SLA,
          response = Y,
          dataEnv = envir,
          dataTraits = traits,
          divideBySiteTotals = divide,
          verbose = FALSE),"Collinearity detected in CWM-model")

expect_warning(scores(mod_divT_dcca_singular3),
               "Collinearity detected in SNC-model")

expect_inherits(mod_DivT_dccaA11, "dccav")

expect_warning(scores_dccaA11 <- scores(mod_DivT_dccaA11),
               "Collinearity detected")

scores_dccavA11_abs <- sapply(X = scores_dccaA11, FUN = abs)
expect_equal_to_reference(scores_dccavA11_abs, "scores_dccav_abs")

expect_warning(expect_stdout(dcca_print <- print(mod_DivT_dccaA11)),
               "Collinearity detected")

expect_equal(names(dcca_print), 
             c("CCAonTraits", "formulaTraits", "formulaEnv", "data", "call", 
               "weights", "Nobs", "CWMs_orthonormal_traits", "RDAonEnv", 
               "eigenvalues", "c_traits_normed0", "inertia", "site_axes", 
               "species_axes", "c_env_normed", "c_traits_normed"))

# From dcCA dune

envData <- dune_trait_env$envir
traitData <- dune_trait_env$traits
names(envData)[3:4] <- c("Moisture", "Management")

mod_funct_traits <- dc_CA(
  formulaEnv = ~ Moisture + Management,
  formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
  response = Y, 
  envData, 
  traitData, 
  verbose = FALSE)

mod_ecol_traits <- dc_CA(
  formulaEnv = ~ Moisture + Management,
  formulaTraits = ~ F + R + N + L,
  response = Y, 
  envData, 
  traitData, 
  verbose = FALSE)

Y <- Y / rowSums(Y)
mod_cca <- vegan::cca(Y ~ Moisture + Management, data = envData)
mod_CCA <- dc_CA(formulaEnv = ~ Moisture + Management,
                 formulaTraits = ~ Species,
                 response = Y, 
                 envData, 
                 traitData, 
                 verbose = FALSE)

expect_equivalent(mod_CCA$eigenvalues, 
                  vegan::eigenvals(mod_cca, model = "constrained"), 
                  tolerance = 1.e-6)

set.seed(1457)
an_ecol <- anova(mod_ecol_traits, by = "axis") 
an_func <- anova(mod_funct_traits, by = "axis")

an_ecol_max <- an_ecol$maxP
an_func_max <- an_func$maxP

expect_equivalent(unlist(an_ecol_max[3, -4, drop =  TRUE]), 
                  c(1, 0.007899116, 0.964, 0.964), 
                  tolerance = 1.e-7)

expect_equivalent(unlist(an_func_max[2, , drop = TRUE]), 
                  c(1, 0.09752322, 0.10085067, 3.4687262, 0.412), 
                  tolerance = 1.e-7)

