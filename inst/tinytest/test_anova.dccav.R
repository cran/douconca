# test anova.dcca, anova_species, anova_sites
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
modDivT1a <- dc_CA(formulaEnv = ~A1 + Moist + Mag + Use + Manure,
                   formulaTraits = ~SLA + Height + LDMC + Seedmass + Lifespan ,
                   response = Y, 
                   dataEnv = envir, 
                   dataTraits = traits,
                   divideBySiteTotals = divide, 
                   verbose = FALSE)

set.seed(123)
modDivT1a_an <- anova(modDivT1a)
expect_equal_to_reference(modDivT1a_an, "modDivT1a_an")

set.seed(123)
expect_equivalent(anova_species(modDivT1a)$table, modDivT1a_an$species)
# equivalent only for very significant tests; the seed for anova_sites cannot 
# easily be set identical to that for the sites test in anova.dcca.
expect_equivalent(anova_sites(modDivT1a)$table, modDivT1a_an$sites)

# test quant variable only, 1 quant trait
modDivTq11<- dc_CA(formulaEnv = ~ Manure,
                   formulaTraits = ~ SLA,
                   response = Y, 
                   dataEnv =envir, 
                   dataTraits = traits,
                   divideBySiteTotals = divide, 
                   verbose = FALSE)

set.seed(123)
modDivTq11_an <- anova(modDivTq11)
expect_equal_to_reference(modDivTq11_an, "modDivTq11_an")

set.seed(123)
# test of the by axis of 1 single predictor. 
expect_equivalent(anova(modDivTq11, by ="axis")$species, 
                  modDivTq11_an$species)

# check how douconca manages collinear models
#  A11 collinear variable
envir$Sites <- factor(dune_trait_env$envir$Sites)
envir$A11 <- envir$A1

expect_silent(modDivT_dccaA1 <- 
                dc_CA(formulaEnv = ~ A1,
                      formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                      response = Y, dataEnv = envir, dataTraits = traits,
                      divideBySiteTotals = divide,
                      verbose = FALSE))
expect_stdout( modDivT_dccaA11 <- 
                 dc_CA(formulaEnv = ~ A1 + A11,
                       formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                       response = Y, dataEnv = envir, dataTraits = traits,
                       divideBySiteTotals = divide,
                       verbose = FALSE))

set.seed(237)
anA1 <- anova(modDivT_dccaA1)

set.seed(237)
anA11<- anova(modDivT_dccaA11)

expect_equivalent(anA1, anA11)

set.seed(237)
expect_equivalent(anova_species(modDivT_dccaA11)$table[, 1:5], 
                  anA11$species[, 1:5])

# full fit
modDivT_dcca_near_singular_species <- 
  dc_CA(formulaEnv = ~ A1 + Manure + Condition(Mag),
        formulaTraits = ~ Species,
        response = Y,
        dataEnv = envir,
        dataTraits = traits,
        divideBySiteTotals = FALSE,
        verbose = FALSE)

set.seed(37)
anova_dccaDivT <- anova(modDivT_dccaA11)
expect_equal_to_reference(anova_dccaDivT, "anova_dccaDivT")

set.seed(159)
an_env <- anova_sites(modDivT_dcca_near_singular_species)

set.seed(159)
an_env2 <- anova(wrda(modDivT_dcca_near_singular_species$formulaEnv,
                      response=modDivT_dcca_near_singular_species$CWMs_orthonormal_traits,
                      data = modDivT_dcca_near_singular_species$data$dataEnv,
                      weights = modDivT_dcca_near_singular_species$weights$rows))
expect_equivalent(an_env$table, an_env2$table)
