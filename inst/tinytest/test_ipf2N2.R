data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
Y <- dune_trait_env$comm[, -1] # must delete "Sites"
Y_N2 <- ipf2N2(Y, updateN2 = FALSE, N2N_N2_species = FALSE)
Y_iter0 <-ipf2N2(Y,max_iter = 0)
expect_equal(attr(Y_N2, "iter"), 61)

# show that column margins of the transform matrix are
# equal to the Hill N2 values
expect_equal(range(colSums(Y_N2) / apply(X = Y, MARGIN = 2, FUN = fN2)),
             c(1, 1), tol = 1.0e-6)
expect_equal(range(rowSums(Y_N2) / apply(X = Y, MARGIN = 1, FUN = fN2)),
             c(1.051043, 1.051043), tol = 1.0e-6)

Y_N2i <- ipf2N2(Y, updateN2 = TRUE, N2N_N2_species = FALSE)

expect_equal(range(colSums(Y_N2i) / apply(X = Y_N2i, MARGIN = 2, FUN = fN2)),
             rep(1,2), tol = 1.0e-6)
expect_equal(range(rowSums(Y_N2i) / apply(X = Y_N2i, MARGIN = 1, FUN = fN2)),
             c(1, 1), tol = 1.0e-6)

# the default version:
Y_N2N_N2i <- ipf2N2(Y)

expect_equal(attr(Y_N2N_N2i, "iter"), 29)

N2 <- apply(X = Y_N2N_N2i, MARGIN = 2, FUN = fN2)
N <- nrow(Y)
expect_equal(range(colSums(Y_N2N_N2i) / (N2 * (1 - N2/N))),
             c(1,1), tol = 1.0e-6)

N2_sites <- apply(X = Y_N2N_N2i, MARGIN = 1, FUN = fN2)
R <- rowSums(Y_N2N_N2i)
expect_equal(range(N2_sites / sum(N2_sites) / (R / sum(R))), 
             c(0.9731012, 1.0054688), tol = 1.0e-6)

expect_equal(range(N2_sites /R), 
             c(1.454648, 1.503033), tol = 1.0e-6)
mod0 <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
              formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
              response = Y,  
              dataEnv = dune_trait_env$envir,
              dataTraits = dune_trait_env$traits, 
              divide = FALSE,
              verbose = FALSE)

expect_message(mod1 <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
                             formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                             response = Y_N2N_N2i,  
                             dataEnv = dune_trait_env$envir,
                             dataTraits = dune_trait_env$traits, 
                             verbose = FALSE), "Argument divideBySiteTotals set to FALSE")

expect_equivalent(mod1$eigenvalues > mod0$eigenvalues, 
                  c(rep(TRUE, 4), FALSE))

Y_PA <- 1 * (Y > 0)
Y_PA_N2 <- ipf2N2(Y_PA, N2N_N2_species = FALSE)

expect_equivalent(Y_PA, Y_PA_N2) 

Y_PA_N2i <- ipf2N2(Y_PA, N2N_N2_species = TRUE)

N_occ <- colSums(Y_PA) # number of occurrences of species
N <- nrow(Y_PA)

expect_equal(cor(colSums(Y_PA_N2i), N_occ * (N - N_occ)),
             0.982612338165359, tol = 1.0e-6)

N2spp <- douconca:::fN2N_N2(Y_PA_N2i, 2, N2N_N2 = TRUE)
# the columns sum and N2(1-N2/N)=  N_occ * (N - N_occ)

expect_equivalent(N2spp, colSums(Y_PA_N2i)) 

mod2 <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
              formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
              response = Y_PA,  
              dataEnv = dune_trait_env$envir,
              dataTraits = dune_trait_env$traits,
              divideBySiteTotals = FALSE,
              verbose = FALSE)

mod3 <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
              formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
              response = Y_PA_N2i,  
              dataEnv = dune_trait_env$envir,
              dataTraits = dune_trait_env$traits,
              verbose = FALSE)

expect_equivalent(mod3$eigenvalues > mod2$eigenvalues, rep(TRUE, 5))
