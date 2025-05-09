data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
# use vegan::rda in step 2
divide <- FALSE # divide by site.totals if TRUE
f2c <- TRUE #factor2categories 

Y <- dune_trait_env$comm[, -1]  # must delete "Sites"
# delete "Species", "Species_abbr" from traits and
# use all remaining variables due to formulaTraits = ~. (the default)
traits <- dune_trait_env$traits
envir <- dune_trait_env$envir
# test a normal mod_DivF
mod <- dc_CA(formulaEnv = Y ~ A1 + Moist + Mag + Use + Manure,
             formulaTraits = ~ F + R + N + L,
             dataEnv = envir, 
             dataTraits = traits,
             divideBySiteTotals = divide, 
             verbose = FALSE)

consider =  c("A1", "Moist", "Use", "Mag")

expect_error(FS(mod,consider =  consider, permutations = 999), 
             "none of the variables")

set.seed(123)
slc <- FS(mod, select = "e", consider = consider, permutations = 999,
         factor2categories = f2c)

set.seed(123)
slc2 <- FS(mod, select = "e",consider = rev(consider), permutations = 999,
          initial_model = "Moist",factor2categories = f2c)

set.seed(123)
slc3 <- FS(mod, select = "e", consider = consider, permutations = 999, 
          initial_model = "Moist + MagNM", factor2categories = f2c)

mod1 <- dc_CA(formulaEnv = Y ~ Moist,
              formulaTraits = mod$formulaTraits,
              dataEnv = envir,
              dataTraits = traits,
              divideBySiteTotals = divide, 
              verbose = FALSE)

set.seed(123)
an1 <- anova_sites(mod1, permutations = 999)
envir1 <- cbind(envir, model.matrix(~ 0 + Mag, data = envir))

mod2 <- dc_CA(formulaEnv = Y ~ MagNM + Condition(Moist),
              formulaTraits = mod$formulaTraits,
              dataEnv = envir1, 
              dataTraits = traits,
              divideBySiteTotals = divide, 
              verbose = FALSE)

set.seed(123)
an2 <- anova_sites(mod2, permutations = 999)
mod3 <- dc_CA(formulaEnv =  Y ~ A1 + Condition(Moist + MagNM),
              formulaTraits = mod$formulaTraits,
              dataEnv = envir1, 
              dataTraits = traits,
              divideBySiteTotals = divide, 
              verbose = FALSE)

set.seed(123)
an3 <- anova_sites(mod3, permutations = 999)
expect_equivalent(c(an1$table$`Pr(>F)`[1], 
                    an2$table$`Pr(>F)`[1],
                    an3$table$`Pr(>F)`[1]),
                  c(slc$finalWithOneExtra[1, "pval"],
                    slc2$finalWithOneExtra[1, "pval"],
                    slc3$final[1, "pval"]))

expect_equal(slc$finalWithOneExtra[, 1:2],
             rbind(slc$finalWithOneExtra[1, , drop = FALSE],
                   slc2$finalWithOneExtra[1, , drop = FALSE],
                   slc3$final_Non_Significant[1, , drop = FALSE])[, 1:2])
