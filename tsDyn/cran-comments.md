

## Test environments

* local R installation, R 4.2.2 on Ubuntu 20.04
* win-builder (devel, release, old-release)

## R CMD check results

- winbuilder R-oldrelease has a NOTE, _Possibly mis-spelled words in DESCRIPTION_, but the words are all listed in inst/WORDLIST.
- there is a possible NOTE on ATLAS due to some differences in numerical precision. Unfortunately, as ATLAS is not replicable on r-hub to the best of my knowledge, it is not possible for me to verify/address this issue, see: https://stat.ethz.ch/pipermail/r-devel/2019-December/078852.html