data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
divide <- FALSE # divide by site.totals if TRUE

dune_trait_env$envir$Sites <- factor(dune_trait_env$envir$Sites)
dune_trait_env$envir$A11 <- dune_trait_env$envir$A1

# collinear variable
mod_dcca <- dc_CA(formulaEnv = ~ A1 + A11,
                  formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                  response = dune_trait_env$comm[, -1],  # must delete "Sites"
                  dataEnv = dune_trait_env$envir,
                  dataTraits = dune_trait_env$traits,
                  divideBySiteTotals = divide,
                  verbose = FALSE)

# full fit
mod_dcca1 <- dc_CA(formulaEnv = ~ Sites,
                   formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv = dune_trait_env$envir,
                   dataTraits = dune_trait_env$traits,
                   divideBySiteTotals = divide,
                   verbose = FALSE)

mod_dcca2 <- dc_CA(formulaEnv = ~ Sites + A1,
                   formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv = dune_trait_env$envir,
                   dataTraits = dune_trait_env$traits,
                   divideBySiteTotals = divide, 
                   verbose = FALSE)

expect_warning(scores(mod_dcca2),"collinearity detected in CWM-model")																			
dune_trait_env$traits$Species <- factor(dune_trait_env$traits$Species)
dune_trait_env$traits$Species_abbr <- factor(dune_trait_env$traits$Species_abbr)

expect_message(mod_dcca3 <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
                                  formulaTraits = ~ Species_abbr + SLA,
                                  response = dune_trait_env$comm[, -1],
                                  dataEnv = dune_trait_env$envir,
                                  dataTraits = dune_trait_env$traits,
                                  divideBySiteTotals = divide,
                                  verbose = FALSE),
               "Some constraints or conditions were aliased because they were redundant")

expect_warning(scores(mod_dcca3),"collinearity detected in SNC-model")																			
expect_inherits(mod_dcca, "dcca")
expect_equal_to_reference(mod_dcca, "mod_dcca")

expect_warning(scores_dcca <- scores(mod_dcca),
               "collinearity detected")
expect_equal_to_reference(scores_dcca, "scores_dcca")

set.seed(37)
expect_message(anova_dcca <- anova(mod_dcca),
               "Some constraints or conditions were aliased because they were redundant")
expect_equal_to_reference(anova_dcca, "anova_dcca")

expect_message(anova_dcca_species <- anova_species(mod_dcca),
               "Some constraints or conditions were aliased because they were redundant")
expect_equivalent(anova_dcca_species$table[, 1:4], anova_dcca$species[, 1:4])

anova_dcca_sites <- anova_sites(mod_dcca)
expect_equivalent(anova_dcca_sites$table[, 1:4], anova_dcca$sites[, 1:4])

expect_equivalent(anova_dcca_species$eigenvalues, mod_dcca$eigenvalues)
expect_equivalent(anova_dcca_sites$eigenvalues, mod_dcca$eigenvalues)

expect_warning(expect_stdout(dcca_print <- print(mod_dcca)),
               "collinearity detected")
expect_equal(names(dcca_print), 
             c("CCAonTraits", "formulaTraits", "formulaEnv", "data", "call", 
               "weights", "Nobs", "CWMs_orthonormal_traits", "RDAonEnv", 
               "eigenvalues", "c_traits_normed0", "inertia", "site_axes", 
               "species_axes", "c_env_normed", "c_traits_normed"))

