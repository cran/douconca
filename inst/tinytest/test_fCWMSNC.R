data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
# use vegan::rda in step 2
divide <- FALSE # divide by site.totals if TRUE

Y <-dune_trait_env$comm[, -1]  # must delete "Sites"
# delete "Species", "Species_abbr" from traits and
# use all remaining variables due to formulaTraits = ~. (the default)
traits <- dune_trait_env$traits
envir <- dune_trait_env$envir

dcca_mod_DivF <- 
  dc_CA(formulaEnv = ~ A1 + Moist + Manure + Use + Condition(Mag),
        formulaTraits = ~ SLA + Height + LDMC + Condition(Seedmass) + Lifespan,
        response = Y,  # must delete "Sites"
        dataEnv = envir,
        dataTraits = traits,
        divideBySiteTotals = divide,
        verbose = FALSE)

CWMSNCa <- fCWM_SNC(formulaEnv = dcca_mod_DivF$formulaEnv,
                    formulaTraits = dcca_mod_DivF$formulaTraits,
                    response = Y,  # must delete "Sites"
                    dataEnv = envir,
                    dataTraits = traits,
                    divideBySiteTotals = divide,
                    verbose = FALSE)

expect_inherits(CWMSNCa, "list")
expect_equal(names(CWMSNCa), 
             c("CWM", "SNC", "formulaEnv", "formulaTraits", "inertia", "weights", 
               "call", "data"))
dcca_mod_DivF0 <- dcca_mod_DivF
dcca_mod_DivF0$data$Y <-NULL
expect_equal(dcca_mod_DivF0[c("weights" ,"data")], 
             CWMSNCa[c( "weights","data")])

CWMSNCb <- fCWM_SNC(formulaEnv = dcca_mod_DivF$formulaEnv,
                    formulaTraits = dcca_mod_DivF$formulaTraits,
                    response = Y,  # must delete "Sites"
                    dataEnv = envir,
                    dataTraits = traits,
                    divideBySiteTotals = divide,
                    minimal_output = FALSE,
                    verbose = FALSE)

expect_inherits(CWMSNCb, "list")
expect_equal(names(CWMSNCb), 
             c("CWM", "SNC","formulaEnv", "formulaTraits", "inertia",
               "weights", "call", "data", "CWMs_orthonormal_traits", 
			   "SNCs_orthonormal_env", "trans2ortho", "T_ortho", "E_ortho"))													

dcca_mod_DivF2 <- dc_CA(formulaEnv = dcca_mod_DivF$formulaEnv,
                        response = CWMSNCa,  
                        verbose = FALSE)

expect_equal(dcca_mod_DivF$eigenvalues, dcca_mod_DivF2$eigenvalues)
# Note: the axes of dcca_mod_DivF2 may have switched sign compared to dcca_mod_DivF
# formulaTrait is taken from CWMSNCa
expect_equal(sapply(X = scores(dcca_mod_DivF, choices = 
                                 seq_along(dcca_mod_DivF$eigenvalues)), 
                    FUN = abs),
             sapply(X = scores(dcca_mod_DivF2, 
                               choices = seq_along(dcca_mod_DivF$eigenvalues)), 
                    FUN = abs))

# dc-CA without species data, without SNC
CWMSNCc <- CWMSNCa
CWMSNCc$SNC <- NULL
dcca_mod_DivF3 <- dc_CA(formulaEnv = dcca_mod_DivF$formulaEnv,
                        response = CWMSNCc,  
                        verbose = FALSE)

expect_equal(dcca_mod_DivF$eigenvalues, dcca_mod_DivF3$eigenvalues)

expect_warning(XT3 <- scores(dcca_mod_DivF3), "SNC analysis is not available")
expect_equal(sapply(X = scores(dcca_mod_DivF)[names(XT3)], FUN = abs),
             sapply(X = XT3, FUN = abs))

# example 1 of no weights specified
CWMSNCd <- CWMSNCb
CWMSNCd$weights$rows <- NULL

expect_warning(dcca_mod_DivF4 <- dc_CA(formulaEnv = dcca_mod_DivF$formulaEnv,
                                       response = CWMSNCd, 
                                       divideBySiteTotals = divide,
                                       verbose = FALSE),
               "no site weights supplied")

expect_equivalent(as.numeric(dcca_mod_DivF4$eigenvalues), 
                  c(0.0785279, 0.0400522, 0.02339385, 0.00417814), 1.0e-7)
expect_warning(scores(dcca_mod_DivF4),
               "The eigenvalues of the CWM and SNC analyses differ")

# example 2 of no weights specified
CWMSNCd <- CWMSNCb
CWMSNCd$weights <- NULL
expect_warning(dcca_mod_DivF4 <- dc_CA(formulaEnv = dcca_mod_DivF$formulaEnv,
                                       response = CWMSNCd, 
                                       verbose = FALSE),
               "no weights supplied")

expect_equivalent(as.numeric(dcca_mod_DivF4$eigenvalues), 
                  c(0.053673468, 0.024680073, 0.013251034, 0.001575069), 1.e-7) 
expect_warning(scores(dcca_mod_DivF4), 
               "The eigenvalues of the CWM and SNC analyses differ.")

# example of weights specified via dataTraits and dataEnv
CWMSNCe <- CWMSNCd
CWMSNCe$dataTraits <- CWMSNCa$data$dataTraits
CWMSNCe$dataTraits$weight <- CWMSNCa$weights$columns
CWMSNCe$data$dataTraits <- NULL
CWMSNCe$weights <- NULL 
envir$weight <- CWMSNCa$weights$rows

expect_warning(dcca_mod_DivF5 <- dc_CA(formulaEnv = dcca_mod_DivF$formulaEnv,
                                       dataEnv = envir,
                                       response = CWMSNCe,  
                                       verbose = FALSE),
               "species weights taken from")

expect_equal(dcca_mod_DivF5$eigenvalues, dcca_mod_DivF$eigenvalues)

dcca_mod_DivF6<- dc_CA(formulaEnv = dcca_mod_DivF5$formulaEnv,
                       formulaTraits = dcca_mod_DivF5$formulaTraits,
                       response = dcca_mod_DivF$data$Y,  
                       dataEnv = dcca_mod_DivF5$data$dataEnv,
                       dataTraits = dcca_mod_DivF5$data$dataTraits,
                       divideBySiteTotals = divide,
                       verbose = FALSE)

expect_equivalent(sapply(X = scores(dcca_mod_DivF6), FUN = abs), 
                  sapply(X = scores(dcca_mod_DivF5), FUN = abs))

# A minimal specification with a non-trivial trait model, giving 3 warnings
envir$weight <- NULL
CWMSNCf <- list(CWM = as.data.frame(CWMSNCa$CWM),
                weights = list(columns = 100 * dcca_mod_DivF$weights$columns),
                dataTraits = traits,
                formulaTraits = ~ SLA + Height + LDMC + Condition(Seedmass) + Lifespan)

# Without trait covariates and only traits in dataTraits, 
# even formulaTraits can be deleted from the list.
# For the same trait model, different environmental predictors can be used
expect_warning(dcca_mod_DivF7 <- dc_CA(response = CWMSNCf,
                                       dataEnv = envir,
                                       formulaEnv = ~ Moist,
                                       verbose = FALSE),
               "no site weights supplied")

# test predict with dc_CA from CWMSNC

# Ten 'new' sites with a subset of the variables in mod 
# X_lot will be ignored as it is not part of the model
newEnv <- dune_trait_env$envir[1:10, c("A1", "Mag", "Manure", "X_lot")]
newEnv[2, "A1"] <- 3.0
rownames(newEnv) <- paste0("NewSite", 1:10)
pred.traits <- predict(dcca_mod_DivF, type = "traitsFromEnv", newdata = newEnv)

# Eight 'new' species with a subset of traits that are included in the model 
# Variable "L" will be ignored as it is not in the model 
newTraits <- dune_trait_env$traits[1:8, c("Species","SLA", "LDMC", "L")]
newTraits[3, "SLA"]<- 18
rownames(newTraits) <- paste("Species", LETTERS[1:8])
pred.env <- predict(dcca_mod_DivF, type = "envFromTraits", newdata = newTraits)

pred.resp <- predict(dcca_mod_DivF, type = "response", 
                     newdata = list(newTraits, newEnv),
                     weights = list(species = rep(1:2, 4), sites = rep(1, 10)))

pred.traits2 <- predict(dcca_mod_DivF2, type = "traitsFromEnv", newdata = newEnv)

# Eight 'new' species with a subset of traits that are included in the model 
# Variable "L" will be ignored as it is not in the model 
newTraits <- dune_trait_env$traits[1:8, c("Species", "SLA", "LDMC", "L")]
newTraits[3, "SLA"] <- 18
rownames(newTraits) <- paste("Species",LETTERS[1:8])
pred.env2 <- predict(dcca_mod_DivF2, type = "envFromTraits", newdata = newTraits)

pred.resp2 <- predict(dcca_mod_DivF2, type = "response", 
                      newdata = list(newTraits, newEnv),
                      weights = list(species = rep(1:2, 4), sites = rep(1, 10))) 

expect_equal(pred.env, pred.env2)
expect_equal(pred.traits, pred.traits2)
expect_equal(pred.resp, pred.resp2)

# fitted and predict for type = response
r <- 2
fitted.resp <- fitted(dcca_mod_DivF, type = "response", rank = r)
fitted.resp2 <- fitted(dcca_mod_DivF2, type = "response", rank = r)
expect_equal(fitted.resp, fitted.resp2 * sum(dcca_mod_DivF$data$Y))

pred0.resp <- predict(dcca_mod_DivF, type = "response")
pred22.resp2 <- predict(dcca_mod_DivF2, type = "response")
expect_equal(pred0.resp, pred22.resp2)

all_reg <- coef(dcca_mod_DivF, type = "a")
all_reg2 <- coef(dcca_mod_DivF2, type = "a")
expect_equal(all_reg, all_reg2)

