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

# test a normal mod_DivF
mod_DivF1a <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
                    formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                    response = Y, 
                    dataEnv = envir, 
                    dataTraits = traits,
                    divideBySiteTotals = divide, 
                    verbose = FALSE)

# Because singular vectors have arbitrary sign, 
# the absolute value of scores is tested.
# The test is an implicit test on the singular values and eigen values as well.
scores_mod_DivF1a_abs <- sapply(X = scores(mod_DivF1a), FUN = abs)
expect_equal_to_reference(scores_mod_DivF1a_abs, "scores_mod_DivF1a_abs")

# test all variables(~.) in formulas
mod_DivF1b <- dc_CA(formulaEnv = ~ .,
                    formulaTraits = ~ . ,
                    response = Y, 
                    dataEnv = envir[, -c(1, 7:10)], 
                    dataTraits = traits[, -c(1, 2, 8:11)],
                    divideBySiteTotals = divide,
                    verbose = FALSE)
expect_equal(mod_DivF1a$eigenvalues, mod_DivF1b$eigenvalues)
expect_equal(mod_DivF1a$c_traits_normed0, mod_DivF1b$c_traits_normed0)

# model with covariates

mod_DivFc <- dc_CA(formulaEnv = ~ A1 + Moist + Manure + Use + Condition(Mag),
                   formulaTraits = ~ SLA + Height + LDMC + Condition(Seedmass) + Lifespan,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv = dune_trait_env$envir,
                   dataTraits = dune_trait_env$traits,
                   divideBySiteTotals = divide,
                   verbose = FALSE)

Y <- mod_DivFc$data$Y
mod_DivFc_vegan <- vegan::cca(formula = Y ~ A1 + Moist + Manure + Use + Condition(Mag),
                              data = mod_DivFc$data$dataEnv)
expect_equal(mod_DivFc$inertia["env_explain", 1],
             sum(vegan::eigenvals(mod_DivFc_vegan, model = "constrained")))


scores_mod_DivFc_abs <- sapply(X = scores(mod_DivFc), FUN = abs)
expect_equal_to_reference(scores_mod_DivFc_abs, "scores_mod_DivFc_abs")

# test factors only
envir$fUse <- factor(envir$Use)
traits$fSLA <- cut(traits$SLA, 3)
mod_DivF_fact <- dc_CA(formulaEnv = ~ fUse + Mag,
                       formulaTraits = ~ fSLA + Lifespan ,
                       response = Y, 
                       dataEnv = envir, 
                       dataTraits = traits,
                       divideBySiteTotals = divide, 
                       verbose = FALSE)

scores_mod_DivF_fact_abs <- sapply(X = scores(mod_DivF_fact), FUN = abs)
expect_equal_to_reference(scores_mod_DivF_fact_abs, "scores_mod_DivF_fact_abs")

traits$lifequant <- 1 * traits$Lifespan %in% "perennial"
mod_DivF_factb <- dc_CA(formulaEnv = ~fUse + Mag,
                        formulaTraits = ~fSLA + lifequant ,
                        response = Y, 
                        dataEnv = envir, 
                        dataTraits = traits,
                        divideBySiteTotals = divide, 
                        verbose = FALSE)

expect_silent(mod_DivF_fact)
expect_equal(mod_DivF_fact$eigenvalues, mod_DivF_factb$eigenvalues)
expect_equivalent(mod_DivF_fact$c_traits_normed0, mod_DivF_factb$c_traits_normed0)

# test quant variable only, 1 quant trait

mod_DivFq11<- dc_CA(formulaEnv = ~ Manure,
                    formulaTraits = ~ SLA,
                    response = Y, 
                    dataEnv =envir, 
                    dataTraits = traits,
                    divideBySiteTotals = divide,
                    verbose = FALSE)

scores_mod_DivFq11_abs <- sapply(X = scores(mod_DivFq11), FUN = abs)
expect_equal_to_reference(scores_mod_DivFq11_abs, "scores_mod_DivFq11_abs")

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

expect_stdout(expect_warning(
  mod_DivF_dccaA11 <- 
    dc_CA(formulaEnv = ~ A1 + A11+Manure+Moist,
          formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
          response = Y, 
          dataEnv = envir, 
          dataTraits = traits,
          divideBySiteTotals = divide,
          verbose = FALSE), "singular"))

expect_equal(mod_DivF_dccaA1$eigenvalues, mod_DivF_dccaA11$eigenvalues)
expect_equivalent(mod_DivF_dccaA1$inertia[-3, , drop = FALSE], 
                  mod_DivF_dccaA11$inertia)
expect_equal(abs(scores(mod_DivF_dccaA1, display = "tval_traits")), 
             abs(scores(mod_DivF_dccaA11, display = "tval_traits")))

expect_warning(scores(mod_DivF_dccaA11, display = "reg"), 
               "Collinearity detected")
expect_warning(scores(mod_DivF_dccaA11, display = "tval"), 
               "Collinearity detected")
# note that the VIFs differ due to NA in VIF for manure! so also unvalid tvalues

expect_warning(scores_A11 <- scores(mod_DivF_dccaA11, display = "reg"),
               "Collinearity detected in CWM-model.")

expect_equal(abs(scores(mod_DivF_dccaA1, display = "reg"))[1:3, -3, drop = FALSE], 
             abs(scores_A11)[-2, -3, drop = FALSE])

# SLA11 collinear
traits$SLA11 <- traits$SLA

mod_DivF_dccaSLA <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
                          formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                          response = Y, 
                          dataEnv = envir, 
                          dataTraits = traits,
                          divideBySiteTotals = divide,
                          verbose = FALSE)

expect_message(
  mod_DivF_dccaSLA11 <- 
    dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
          formulaTraits = ~ SLA + SLA11+ Height + LDMC + Seedmass + Lifespan,
          response = Y, 
          dataEnv = envir, 
          dataTraits = traits,
          divideBySiteTotals = divide,
          verbose = FALSE))

expect_equal(mod_DivF_dccaSLA$eigenvalues, mod_DivF_dccaSLA11$eigenvalues)
expect_equal(mod_DivF_dccaSLA$inertia, mod_DivF_dccaSLA11$inertia)
expect_equal(abs(scores(mod_DivF_dccaSLA, display = "tval")), 
             abs(scores(mod_DivF_dccaSLA11, display = "tval")))
expect_warning(scores_SLA11 <- scores(mod_DivF_dccaSLA11, display = "reg_traits"), 
               "Collinearity")
expect_equal(abs(scores(mod_DivF_dccaSLA, display = "reg_traits"))[, -3], 
             abs(scores_SLA11)[-2, -3])

expect_warning(tval_SLA11 <- scores(mod_DivF_dccaSLA11, display = "tval_traits"),
               "Collinearity detected")
expect_true(all(is.na(tval_SLA11)))

# full fit
expect_stdout( mod_DivF_dcca_near_singular_sites <- 
                 dc_CA(formulaEnv = ~ Sites,
                       formulaTraits = ~ SLA + Seedmass,
                       response = Y,
                       dataEnv = envir,
                       dataTraits = traits,
                       divideBySiteTotals = divide,
                       verbose = FALSE),)
expect_true(all(scores(mod_DivF_dcca_near_singular_sites, display = "tval") == 0))

expect_stdout(mod_DivF_dcca_near_singular_species <- 
                dc_CA(formulaEnv = ~ A1 + Manure + Condition(Mag),
                      formulaTraits = ~Species,
                      response = Y,
                      dataEnv = envir,
                      dataTraits = traits,
                      divideBySiteTotals = divide,
                      verbose = FALSE), )

expect_message(dc_CA(formulaEnv = ~ A1 + Manure + Condition(Mag),
                     formulaTraits = ~ Species,
                     response = Y,
                     dataEnv = envir,
                     dataTraits = traits,
                     divideBySiteTotals = divide,
                     verbose = FALSE), 
               "The model is overfitted")

expect_true(all(is.nan(scores(mod_DivF_dcca_near_singular_species, 
                              display = "tval_traits"))))

expect_stdout(expect_warning(
  mod_DivF_dcca_singular2 <- 
    dc_CA(formulaEnv = ~ Sites + A1,
          formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
          response = Y,  # must delete "Sites"
          dataEnv = envir,
          dataTraits = traits,
          divideBySiteTotals = divide, 
          verbose = FALSE), "singular"))

expect_warning(scores(mod_DivF_dcca_singular2),
               "Collinearity detected in CWM-model")		

expect_stdout(expect_message(
  mod_DivF_dcca_singular3 <- 
    dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
          formulaTraits = ~ Species_abbr + SLA,
          response = Y,
          dataEnv = envir,
          dataTraits = traits,
          divideBySiteTotals = divide,
          verbose = FALSE)))

expect_warning(scores(mod_DivF_dcca_singular3), 
               "Collinearity detected in SNC-model")
expect_inherits(mod_DivF_dccaA11, "dcca")

expect_warning(scores_dccaA11 <- scores(mod_DivF_dccaA11),
               "Collinearity detected")

scores_dccaA11DivF_abs <- sapply(X = scores_dccaA11, FUN = abs)
expect_equal_to_reference(scores_dccaA11DivF_abs, "scores_dccaA11DivF_abs")

expect_warning(expect_stdout(dcca_print <- print(mod_DivF_dccaA11)),
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
  divideBySiteTotals = FALSE, 
  verbose = FALSE)

mod_ecol_traits <- dc_CA(
  formulaEnv = ~ Moisture + Management,
  formulaTraits = ~ F + R + N + L,
  response = Y, 
  envData, 
  traitData,
  divideBySiteTotals = FALSE,
  verbose = FALSE)

mod_cca <- vegan::cca(Y ~ Moisture + Management, data = envData)
mod_CCA <- dc_CA(formulaEnv = ~ Moisture + Management,
                 formulaTraits = ~ Species,
                 response = Y, 
                 envData, 
                 traitData,
                 divideBySiteTotals = FALSE, 
                 verbose = FALSE)

expect_equivalent(mod_CCA$eigenvalues, 
                  as.numeric(vegan::eigenvals(mod_cca, model = "constrained")), 
                  tolerance = 1e-6)

set.seed(1457)
an_ecol <- anova(mod_ecol_traits, by = "axis")
an_func <- anova(mod_funct_traits, by = "axis")
an_ecol_max <- an_ecol$maxP
an_func_max <- an_func$maxP

expect_equivalent(unlist(an_ecol_max[3, , drop =  TRUE]), 
                  c(1, 0.008964815, 0.948, 0.92, 0.948), 
                  tolerance = 1.e-7)

expect_equivalent(unlist(an_func_max[2, , drop =  TRUE]), 
                  c(1, 0.08028837, 0.08730463, 2.84752775, 0.542), 
                  tolerance = 1.e-7)
