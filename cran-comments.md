

## Test environments

* local R installation, R 4.2.2 on Ubuntu 20.04
* win-builder (devel, release, old-release)

## R CMD check results

- winbuilder R-oldrelease has a NOTE, _Possibly mis-spelled words in DESCRIPTION_, but the words are all listed in inst/WORDLIST.
- winbuilder R-release and devel have a NOTE, _Found the following (possibly) invalid DOIs:_, but this is a false positive. 
- there is a NOTE on ATLAS due to some differences in numerical precision. Unfortunaltey, as ATLAS is not replicable, it is not possible to address this issue, see: https://stat.ethz.ch/pipermail/r-devel/2019-December/078852.html