I am taking back the maintenance of the package, after Ho Tsung-wu kindly accepted to adopt the package, and kindly accepted that I take the maintenance back.

## Test environments
* local R installation, R 4.0.3 on Ubuntu 20.04
* win-builder (devel, release, old-release)
* R-hub standard platforms

## R CMD check results

0 errors | 0 warnings | 1 note

* There is only one note about *(possibly) invalid URLs: https://cran.r-project.org/package=tsDyn/vignettes/tsDyn.pdf*. This is URL for the vignette, which is cited in the package, and will be valid once the package is uploaded.

## Downstream dependencies

I checked downstream dependencies with `revdepcheck::revdep_check` and all passed the test.