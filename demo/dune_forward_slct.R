data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
abun <- dune_trait_env$comm[, -1] # must delete "Sites"

mod <- dc_CA(formulaEnv = abun ~ Moist + Mag,
             formulaTraits = ~ F + R + N + L,
             dataEnv = dune_trait_env$envir,
             dataTraits = dune_trait_env$traits,
             verbose = FALSE)

# selection of traits with environmental model of mod (~ Moist+Mag)
out1 <- FS(mod, consider = c("F", "R", "N", "L"), 
           select = "traits", verbose = FALSE) 

names(out1)
out1$finalWithOneExtra
out1$model_final

# selection of environmental variables with trait model of mod (~ F + R + N + L)
out2 <- FS(mod, consider =  c("A1", "Moist", "Mag", "Use", "Manure"), 
           select= "env", verbose = FALSE) 

names(out2)
out2$finalWithOneExtra
out2$model_final

# selection of environmental variables without a trait model 
# i.e. with a single constraint
mod3 <- cca0(mod$data$Y ~ Moist, data = mod$data$dataEnv)
out3 <- FS(mod3, consider = c("A1", "Moist", "Mag", "Use", "Manure"), 
           threshold_P = 0.05)

out3$finalWithOneExtra
out3$model_final

# selection of traits without an environmental model 
#                         i.e. with a single constraint
tY <- t(mod$data$Y)

mod4 <- cca0(tY ~ L, data = mod$data$dataTraits)

names(mod$data$dataTraits)
out4 <- FS(mod4, 
           consider =  c("SLA", "Height", "LDMC", "Seedmass", "Lifespan", 
                         "F", "R", "N", "L"))

out4$finalWithOneExtra
out4$model_final