# Prepare for CRAN ----

# Update dependencies in DESCRIPTION
attachment::att_amend_desc()

# Run tests and examples
devtools::run_examples()

# Check package as CRAN
rcmdcheck::rcmdcheck(
  args = c("--no-manual", "--as-cran")
)

# Check content
# remotes::install_github("ThinkR-open/checkhelper")
checkhelper::find_missing_tags()

# Check spelling
# usethis::use_spell_check()
spelling::spell_check_package()

# Check URL are correct
# remotes::install_github("r-lib/urlchecker")
urlchecker::url_check()
urlchecker::url_update()

# check on other distributions
platforms <- c('debian-gcc-devel', 'ubuntu-gcc-devel', 'windows-x86_64-devel-ucrt', 'macos-m1-bigsur-release')

devtools::check_rhub(
  platforms = platforms, 
  email = 'stefan.filges@gu.se'
)

# _win devel
devtools::check_win_devel()

# Verify you're ready for release, and release
devtools::release()
