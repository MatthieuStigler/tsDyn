
## Online Checks
# _R_CHECK_FORCE_SUGGESTS_
devtools::check_rhub(interactive = FALSE,
                     env_vars = c(`_R_CHECK_FORCE_SUGGESTS_` = "false", ## rgl not working on Ubuntu Linux 16.04 LTS, R-release, GCC
                                  `_R_CHECK_CRAN_INCOMING_USE_ASPELL_` = "true"))
devtools::check_win_devel()
devtools::check_win_release()
devtools::check_win_oldrelease()

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
