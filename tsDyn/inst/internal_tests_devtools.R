
## Online Checks
devtools::check_rhub(interactive = FALSE)
devtools::check_win_devel()
devtools::check_win_release()

devtools::spell_check(use_wordlist = TRUE)

#
# usethis::use_cran_comments(open = FALSE)



## reverse dependency? not working well as of Nov 2020
# usethis::use_revdep()
# revdepcheck::revdep_check(num_workers = 4)

## Release
devtools::release_checks()
devtools::release()
