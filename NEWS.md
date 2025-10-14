# douconca 1.2.4

* Update of function ipf2N2 for informative pre-processing of abundance data.
The row and column marginals are set equal to Hill N2 or, the column
marginals to N2(1-N2/N), the effective number of informative species.
For max_iter=0 only the species marginal is adapted to N2 or N2(1-N2/N) 
without further adjustment to the abundance table.
This is useful if the function did not converge or gives a warning
indicating very unequal site totals.

# douconca 1.2.3

* Forward selection of traits and of environmental variables added (function FS()).
* Function ipf2N2 for informative pre-processing of abundance data. The 
row and column marginals are set equal to Hill N2 or, the column
marginals to 2N2(N-N2)/N, the effective number of informative species.
informative species
* More efficiency for large data sets by addition of a new cca function (cca0).
* An anova method for cca0 to enable residual predictor permutation.
* Improved stability for 'exceptional' data sets.
* The response can now be supplied as left-hand side of the environmental formula,
instead of by the response argument.

# douconca 1.2.2

* New coef.dcca() and fitted.dcca() functions with predict.dcca() adapted.
The function coef() can give fourth-corner correlations and regression 
coefficients. 
* Patch release with extended test files and associated small corrections,
for example, SDS (standard deviation of predictors)
was in v1.2.1 a constant factor too large with the default of the argument
divideBySiteTotals (the regression weights and t-values were correct).

# douconca 1.2.1

* Patch release addressing check errors on several CRAN build machines.

# douconca 1.2.0

* Intial CRAN release

# douconca 1.1.6

* An issue with collinear predictors in v1.1.5 has been resolved.

# douconca 1.1.5

* The package can now do general dc-CA, instead of the vegan-based version with 
equal site weights only. For users of the previous version, the function 
dc_CA_vegan has been replaced by the more general function dc_CA. 
The default gives the same analysis. By specifying
the argument `divideBySiteTotals = FALSE`, obtain the original dc-CA analysis
with unequal site weights.
* The `plot_dcCA` function is now a method: `plot.`
* General dc-CA required weighted redundancy analysis. For this, a new function
`wrda` has been added, with methods for print, scores and anova.
* A `predict` function has been added.
* A dc-CA can be computed from community-weighted means (CWMs) with
trait and environment data with species and site weights. See the new function 
`fCWM_SNC`. This is of interest, for example, to make a dc-CA analysis 
reproducible when the abundance data cannot be made public, and
it may also allow to perform dcCA with intra-species trait variation. 
The user needs to be able to compute meaningful CWMs in this case and supply 
trait data that reflect the (species-weighted) inter-trait covariance.
* Several functions are updated. In particular, there are corrections to
the anova function.

# douconca 1.1.2

* The `scores.dccav` function is corrected concerning intra-set correlations for
traits and environmental variables.
* The plotting functions are updated to avoid ggplot2 warnings on color and 
size.
* The fitted straight lines in the plots use the implicit weights 
(they did already, but the help said they did not).

