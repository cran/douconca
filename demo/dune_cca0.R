data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
abun <- dune_trait_env$comm[, -1]  # must delete "Sites"

mod <- cca0(formula = abun ~ A1 + Moist + Mag + Use + Manure,
            data = dune_trait_env$envir)

mod # Proportions equal to those Canoco 5.15

scores(mod, which_cor = c("A1", "X_lot"), display = "cor")

set.seed(123)
anova(mod)
anova(mod, by = "axis")

mod2 <- vegan::cca(abun ~ A1 + Moist + Mag + Use + Manure,
                   data = dune_trait_env$envir)
anova(mod2, by = "axis")

dat <- dune_trait_env$envir
dat$Mag <- "SF"
predict(mod, type = "lc", newdata = dat)
