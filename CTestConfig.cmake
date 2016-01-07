# This file should be placed in the root directory of your project.
set(CTEST_PROJECT_NAME "dynac")
set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "abp-cdash.web.cern.ch/abp-cdash/")
set(CTEST_DROP_LOCATION "submit.php?project=dynac")
set(CTEST_DROP_SITE_CDASH TRUE)
set(UPDATE_TYPE "git")
