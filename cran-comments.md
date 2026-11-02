## Submission

This is a resubmission of a package archived on 2026-08-21. It fixes the two issues from the CRAN checks:

* NOTE "Found calls to structure() using deprecated special names" in data/barry.R and data/zeroyld.R (fixed in 11.0.6).
* Differences in tests/compare_tsDyn_vars.Rout.save under the additional MKL checks (printed digits and an all.equal() tolerance, both in the test only).

## Test environments

* local R installation, R 4.6.1 on Ubuntu 24.04
* R-devel (2026-11-02 r90827) on Ubuntu 24.04, with the reference BLAS, OpenBLAS 0.3.26 and Intel MKL 2026.2 (test output identical under all three)
* win-builder (devel, release, old-release)

## R CMD check results

0 errors | 0 warnings | 1 note

* The NOTE is "Package was archived on CRAN" (new submission).

## Reverse dependencies

svars (Suggests) checked with the new version, no change in status.
