Dear cran members.

This submission correspond to version 0.6 of package greed. The check results obtained are similar to those of the previous version. I added to my tests a check on fedora-clang-rdevel since the current cran check result page for greed give a warning with it. The results obtained with R-hub and this image do not raise a warning anymore. The main changes are listed in the NEWS.md file and the check result details follows.

## Test environments

* ubuntu 20.04 (on github-action), R 4.1.3 
* ubuntu 20.04 (on github-action), R-devel 
* windows server 2019 (on github-action), R 4.1.3
* macOS-latest (on github-action), R 4.1.3
* solaris 10 (on rhub), R 4.1.3
* fedora-clang (on rhub), R-devel

## R CMD check results

I still have one note on installed package size that comes from my usage of RcppArmadillo and on UTF-8 usage on SOLARIS.

### fedora-clang-devel
── R CMD check results ──────────────────────────────────────── greed 0.6 ────

❯ checking installed package size ... NOTE
    installed size is  24.4Mb
    sub-directories of 1Mb or more:
      libs   22.1Mb
❯ checking data for non-ASCII characters ... NOTE
  Note: found 989 marked UTF-8 strings
    
0 errors ✔ | 0 warnings ✔ | 2 note ✖

### solaris-release
── R CMD check results ──────────────────────────────────────── greed 0.6 ────
      
❯ checking data for non-ASCII characters ... NOTE
  Note: found 989 marked UTF-8 strings
      
0 errors ✔ | 0 warnings ✔ | 1 note ✖

### macOS-release
── R CMD check results ──────────────────────────────────────── greed 0.6.0 ────

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

### windows-release
── R CMD check results ──────────────────────────────────────── greed 0.5.1 ────

❯ checking installed package size ... NOTE
    installed size is  7.2Mb
    sub-directories of 1Mb or more:
      libs   4.6Mb
      
0 errors ✔ | 0 warnings ✔ | 1 note ✖

### Ubuntu 20.04-release
── R CMD check results ──────────────────────────────────────── greed 0.6.0 ────

❯ checking installed package size ... NOTE
    installed size is 45.5Mb
    sub-directories of 1Mb or more:
      libs  42.8Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

### Ubuntu 20.04-devel
── R CMD check results ──────────────────────────────────────── greed 0.6 ────

❯ checking installed package size ... NOTE
    installed size is 45.5Mb
    sub-directories of 1Mb or more:
      libs  42.8Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖


