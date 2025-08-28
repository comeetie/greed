# greed 0.6.2

* only small change in doc to comply with CRAN notes 

# greed 0.6.1

* fixe url problems in doc

# greed 0.6

* Finalization of the Models and Algorithms API with user friendly constructors and final model acronyms see greed::available_algorithms() and greed::available_models()
* Generic functions to access fit results slots ICL(), K(), clustering(), coefs(), prior() and analyze them cut(), plot()
* New CombinedModels() model to combine multiple observational models with the same discrete latent variable and conditional independence assumptions
* New Lca() model for Latent class analysis
* New interface for Mixture of regression models (MoR) with the use of a formula object to define the regression model to use.
* 7 Vignettes for the main package / model use cases on the website
* New datasets (mushroom, youngpeople survey, ...)
* Better tests for all models swap and merge moves
* Optimization of MoM model swap moves to avoid unnecessary matrix transpositions
* Bug correction on Lbm merge matrix update from the cluster proportions part 
* Bug correction on as.sparse matrix to avoid problems with empty rows / columns


# greed 0.5.1

* New Gaussian mixture model diaggmm for diagonal mixture models
* New generic function coef to extract MAP estimate of the models parameters
* New plot function gmmpairs to explore gmm fit results
* Better input checking for mvmreg and gmm
* Better input checking for greed_cond
* Correction of compilation problems on solaris
* Correction of pointer problem coming from shed_row/shed_col
* Added a `NEWS.md` file to track changes to the package.