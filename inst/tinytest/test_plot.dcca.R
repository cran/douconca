data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites

Y <- dune_trait_env$comm[, -1] # must delete "Sites"
traits <- dune_trait_env$traits
envir <- dune_trait_env$envir

mod_dcca <- dc_CA(formulaEnv = ~ A1 + Moist + Use + Manure + Mag,
                  formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                  response = Y,
                  dataEnv = envir,
                  dataTraits = traits,
                  verbose = FALSE)

plotDat <- getPlotdata(mod_dcca)
expect_equivalent_to_reference(plotDat, "plotdat_dcca")

p <- plot(mod_dcca, verbose = FALSE)

expect_inherits(p, "list")
expect_equal(names(p), c("plot", "nameList", "separateplots"))

expect_equal(names(p[[2]]), c("newnames", "weightnames", "centroidnames"))

expect_equal(names(p[[3]]), c("CWM_SNC", "traits", "env", "species"))

expect_inherits(p[[3]][[1]], "ggplot")
expect_inherits(p[[3]][[2]], "ggplot")
expect_inherits(p[[3]][[3]], "ggplot")
expect_inherits(p[[3]][[4]], "ggplot")


p2 <- plot(mod_dcca, gradient_description = c("corre", "tval"),
           envfactor = NA, verbose = FALSE)

expect_equal(p2[[3]]$traits$labels$title, "correlation")
expect_equal(p2[[3]]$env$labels$title, "t-value")

newnames <-list(traits = c("SLA", "Height", "LDMC", "Seedmass", 
                           "annual", "perennial"),
                env = c("A1 horizon", "Moisture", "Type of use", "Manure", 
                        "MAG SF", "MAG BF", "MAG HF", "MAG NM"))

p3 <- plot(mod_dcca, gradient_description = c("corre", "weights"), 
           newnames = newnames, verbose = FALSE)

expect_equal(p3[[3]]$traits$labels$title, "correlation")
expect_equal(p3[[3]]$env$labels$title, "weight")

expect_equal(p3$nameList$newnames$traits, newnames$traits)
expect_equal(p3$nameList$newnames$env, newnames$env)

# test factors only
envir$fUse <- factor(envir$Use)
traits$fSLA <- cut(traits$SLA, 3)

mod_fact2 <- dc_CA(formulaEnv = ~ fUse + Mag,
                   formulaTraits = ~ fSLA + Lifespan + Seedmass,
                   response = Y, 
                   dataEnv = envir, 
                   dataTraits = traits,
                   divideBySiteTotals = FALSE, 
                   verbose = FALSE)

pp2 <- plot(mod_fact2)
pp2b <- plot(mod_fact2, envfactor = "Mag", traitfactor = "fSLA")
