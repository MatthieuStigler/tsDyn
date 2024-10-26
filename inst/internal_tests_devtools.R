
## just once:
usethis::use_package("generics", "Imports")
usethis::use_package("broom", "Suggests")
usethis::use_package("dplyr", "Suggests")
usethis::use_package("stringr", "Suggests")

usethis::use_package("purrr", "Suggests")
usethis::use_package("tibble", "Suggests")
usethis::use_package("tidyr", "Suggests")

usethis::use_build_ignore(".github")
# usethis::use_news_md()
# usethis::use_cran_comments(open = FALSE)


## URL checks:
urlchecker::url_check(".")
# curlGetHeaders("https://faculty.chicagobooth.edu/ruey-s-tsay/research/analysis-of-financial-time-series-2nd-edition")


##
devtools::check(manual = TRUE,
                remote = TRUE,
                incoming = TRUE)


devtools::check_win_devel()
devtools::check_win_release()
devtools::check_win_oldrelease()

## Rhb  v2
# git remote set-url origin https://github.com/MatthieuStigler/tsDyn.git
rhub::rhub_setup()
rhub::rhub_doctor()

rhub::rhub_platforms()

rhub::rhub_check(platforms = "linux")
rhub::rhub_check(platforms = c("windows", "macos"))


## Online Checks
# _R_CHECK_FORCE_SUGGESTS_
rhub:::default_cran_check_platforms(devtools:::as.package(".")$path)
rhub::platforms()
my_platf <- c("macos-highsierra-release-cran", "fedora-clang-devel",
              "linux-x86_64-rocker-gcc-san")
devtools::check_rhub(interactive = FALSE,
                     platforms = my_platf,
                     env_vars = c(`_R_CHECK_FORCE_SUGGESTS_` = "false", ## rgl not working on Ubuntu Linux 16.04 LTS, R-release, GCC
                                  `_R_CHECK_CRAN_INCOMING_USE_ASPELL_` = "true"))


devtools::spell_check(use_wordlist = TRUE)

cat(devtools::spell_check(use_wordlist = TRUE)$word, sep = "\n")

#
# usethis::use_cran_comments(open = FALSE)



## reverse dependency? not working well as of Nov 2020
# usethis::use_revdep()
# revdepcheck::revdep_check(num_workers = 1)

# Error in if (worker$process$get_exit_status()) { : 
#     argument is not interpretable as logical
#   In addition: Warning message:
#     call dbDisconnect() when finished working with a connection 

## Release
devtools::release_checks()
devtools::release()
# direct: devtools::submit_cran()