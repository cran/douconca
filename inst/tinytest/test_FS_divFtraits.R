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
mod <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
             formulaTraits = ~ L,
             response = Y, 
             dataEnv = envir, 
             dataTraits = traits,
             divideBySiteTotals = divide, 
             verbose = FALSE)

consider =  c("F", "R", "L")

set.seed(123)
slc <- FS(mod, consider = consider, permutations = 999)

set.seed(123)
slc2 <- FS(mod, consider = c("R", "L"), permutations = 999, initial_model = "F")

set.seed(123)
slc3 <- FS(mod, consider = "L", permutations = 999, initial_model = "F + R")

mod1 <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
              formulaTraits = ~ F,
              response = Y, 
              dataEnv = envir, 
              dataTraits = traits,
              divideBySiteTotals = divide, 
              verbose = FALSE)

set.seed(123)
an1 <- anova_species(mod1)
mod2 <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
              formulaTraits = ~ R + Condition(F),
              response = Y, 
              dataEnv = envir, 
              dataTraits = traits,
              divideBySiteTotals = divide, 
              verbose = FALSE)

set.seed(123)
an2 <- anova_species(mod2)
mod3 <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
              formulaTraits = ~ L + Condition(F + R),
              response = Y, 
              dataEnv = envir, 
              dataTraits = traits,
              divideBySiteTotals = divide, 
              verbose = FALSE)

set.seed(123)
an3 <- anova_species(mod3)

expect_equivalent(c(an1$table$`Pr(>F)`[1], 
                    an2$table$`Pr(>F)`[1],
                    an3$table$`Pr(>F)`[1]),
                  c(slc$finalWithOneExtra[1, "pval"],
                    slc2$finalWithOneExtra[1, "pval"],
                    slc3$final[1, "pval"]))

expect_equal(slc$finalWithOneExtra[, 1:2],
             rbind(slc$finalWithOneExtra[1, , drop = FALSE],
                   slc2$finalWithOneExtra[1, , drop = FALSE],
                   slc3$final[1, , drop = FALSE])[, 1:2])
