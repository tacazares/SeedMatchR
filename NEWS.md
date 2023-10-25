# SeedMatchR 1.1.1

* Vignettes are no longer rendered. When evaluating, there was an issue with internet connectivity that caused the build checks to fail on CRAN. These are now provided as un-evaluated vignettes. Visit the documentation page at https://tacazares.github.io/SeedMatchR/. 

# SeedMatchR 1.1.0

* Addressed some notes and errors from CRAN check https://www.r-project.org/nosvn/R.check/r-release-macos-arm64/SeedMatchR-00check.html
* Removed the function `dtat_histogram.R`. This was not needed. We no longer
perform this bootstrapping method for testing. We use the built-in R functions 
`ks.test` and `wilcoxen.test`. 
* Removed `twosamples` package dependency. Added `stats` package dependency. 

# SeedMatchR 1.0.1

* Example data is saved to a temp directory
* Example data is returned from a function instead of writing to env
* `message()` is used instead of `print()`

# SeedMatchR 1.0.0

* First release prepared for CRAN submission

# SeedMatchR 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
