data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
Y <- dune_trait_env$comm[, -1]  # must delete "Sites"  

mod_cca00 <- cca0(formula = Y ~ A1 + Moist + Mag + Use + Condition(Manure),
                  data = dune_trait_env$envir)

expect_error(cca0(formula = response ~ A1 + Moist + Mag + Use + Condition(Manure),
                  data = dune_trait_env$envir), "object 'response' not found")

mod_cca0 <- cca0(formula = ~ A1 + Moist + Mag + Use + Condition(Manure),
                 response = Y, 
                 data = dune_trait_env$envir)

# compare with vegan::cca
mod_vegan <- vegan::cca(formula = Y ~ A1 + Moist + Mag + Use + Condition(Manure), 
                        data = dune_trait_env$envir)

expect_equivalent(mod_cca0$CCA$eig, mod_vegan$CCA$eig)
expect_equal(mod_cca0$tot.chi, mod_vegan$tot.chi)
expect_equal(mod_cca0$CCA$tot.chi, mod_vegan$CCA$tot.chi)

expect_equivalent(abs(scores(mod_cca0,display = "species",scaling="sym")),
                  abs(scores(mod_vegan,display = "species",scaling="sym")))

expect_equivalent(abs(scores(mod_cca00, display = "species", scaling = "sym")),
                  abs(scores(mod_vegan, display = "species", scaling = "sym")))
expect_stdout(cca0_print <- print(mod_cca0))
expect_equal(names(cca0_print), 
             c("call", "method", "tot.chi", "formula", "site_axes", 
               "species_axes", "Nobs", "eigenvalues", "weights", "sumY","data", "Ybar", 
               "pCCA", "CCA", "CA","inertia" ,"c_env_normed"))


set.seed(129)

anova_cca0 <- anova(mod_cca0)
anova_byaxis_cca0 <- anova(mod_cca0, by = "axis")

expect_equivalent_to_reference(anova_cca0, "anova_cca0")
expect_equivalent_to_reference(anova_byaxis_cca0, "anova_byaxis_cca0")

expect_equivalent(anova_cca0$eigenvalues, mod_cca0$eigenvalues)

mod_cca1 <- cca0(formula = ~ A1 + Moist + Mag + Use + Condition(Manure),
                 response = response,  data = dune_trait_env$envir, 
                 cca_object = mod_cca0, object4QR = mod_cca0)

expect_equivalent(mod_cca1[-1],mod_cca0[-1])

