data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
response <- dune_trait_env$comm[, -1]  

w <- rep(1, 20)
w[1:10] <- 4 
w[17:20] <- 0.5

mod_wrda0 <- wrda(formula = ~ A1 + Moist + Manure + Use + Condition(Mag),
                  response = response, 
                  data = dune_trait_env$envir, 
                  weights = w)

# if response is specified, the left-hand-side of formula is not used.
mod_wrda00 <- wrda(formula = Y ~ A1 + Moist + Manure + Use + Condition(Mag),
                   response = response, 
                   data = dune_trait_env$envir, 
                   weights = w)

expect_error(wrda(formula = Y ~ A1 + Moist + Manure + Use + Condition(Mag),
                  data = dune_trait_env$envir, 
                  weights = w), "object 'Y' not found")

mod_wrda <- wrda(formula = response ~ A1 + Moist + Manure + Use + Condition(Mag),
                 data = dune_trait_env$envir, 
                 weights = w)

expect_equal(abs(mod_wrda00$site_axes$site_scores$site_scores_unconstrained),
             abs(mod_wrda0$site_axes$site_scores$site_scores_unconstrained))
expect_equal(abs(mod_wrda$site_axes$site_scores$site_scores_unconstrained),
             abs(mod_wrda0$site_axes$site_scores$site_scores_unconstrained)) 

expect_inherits(mod_wrda, "wrda")

scores_wrda <- scores(mod_wrda, display = "all")

scores_sub <- scores(mod_wrda, which_cor = c("A1", "Manure"), display = "cor")
expect_equivalent(scores_wrda$correlation[c("A1", "Manure"), ], scores_sub)

scores_wrda_abs <- sapply(X = scores_wrda, FUN = abs)
expect_equal_to_reference(scores_wrda_abs, "scores_wrda_abs")


expect_equivalent(mod_wrda$eig, anova(mod_wrda)$eig)

expect_equal(anova(mod_wrda, by = "axis")$table$Variance,
             c(13.4759764625272, 3.51996709660156, 2.83281997217771, 
               1.90283508422467, 24.4005894936307))

set.seed(129)

anova_wrda <- anova(mod_wrda)
anova_byaxis_wrda <- anova(mod_wrda, by = "axis")

expect_equivalent_to_reference(anova_wrda, "anova_wrda")
expect_equivalent_to_reference(anova_byaxis_wrda, "anova_byaxis_wrda")

# The default is equal weights, which allows checking against vegan
mod_wrda_ew <- wrda(formula = response ~ A1 + Moist + Mag + Use + Condition(Manure),
                    response = response, 
                    data = dune_trait_env$envir)

# compare with vegan::rda
mod_vegan <- vegan::rda(formula = response ~ A1 + Moist + Mag + Use + Condition(Manure), 
                        data = dune_trait_env$envir)

expect_equivalent(mod_wrda_ew$CCA$eig, mod_vegan$CCA$eig)
expect_equal(mod_wrda_ew$tot.chi, mod_vegan$tot.chi)
expect_equal(mod_wrda_ew$CCA$tot.chi, mod_vegan$CCA$tot.chi)
expect_equivalent(abs(mod_wrda_ew$CCA$u), abs(mod_vegan$CCA$u))

const <- 6.385427

expect_equivalent(abs(mod_wrda_ew$site_axes$site_scores$lc_env_scores),
                  abs(const *scores(mod_vegan, 
                                    choices = seq_len(ncol(mod_wrda_ew$site_axes$site_scores$lc_env_scores)),
                                    display= "lc", scaling = "sites")), tol = 1.0e-6)

expect_stdout(wrda_print <- print(mod_wrda))

expect_equal(names(wrda_print), 
             c("call", "method", "tot.chi", "formula", "site_axes", 
               "species_axes", "Nobs", "eigenvalues", "weights", "data", "Ybar", 
               "pCCA", "CCA", "CA","inertia" ,"c_env_normed"))

mod_wrda1 <- wrda(formula = mod_wrda$formula,
                  response = response,  data =  mod_wrda$data, 
                  cca_object = mod_wrda, object4QR = mod_wrda)

expect_equivalent(mod_wrda1[-1], mod_wrda[-1])

