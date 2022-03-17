Dear cran members.

This submission correspond to version 0.6 of package greed. The check results obtained are similar to those of the previous version.

## Test environments

* ubuntu 20.04 (on github-action), R 4.1.3 
* ubuntu 20.04 (on github-action), R-devel 
* windows (on github-action), R 4.1.3
* macOS-11 (on github-action), R 4.1.3
* solaris 10 (on rhub), R 4.1.3

## R CMD check results

I still have one note on installed package size that comes from my usage of RcppArmadillo and on UTF-8 usage on SOLARIS.

### solaris-release

── R CMD check results ──────────────────────────────────────── greed 0.5.1 ────
Duration: 6m 58.3s

❯ checking data for non-ASCII characters ... NOTE
  Note: found 6693 marked UTF-8 strings
    
0 errors ✔ | 0 warnings ✔ | 1 note ✖

### macOS-release

── R CMD check results ──────────────────────────────────────── greed 0.5.1 ────
Duration: 7m 9.9s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

### windows-release
── R CMD check results ──────────────────────────────────────── greed 0.5.1 ────
Duration: 13m 24.6s

❯ checking installed package size ... NOTE
    installed size is 5.5Mb
    sub-directories of 1Mb or more:
      libs  4.6Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

### Ubuntu 20.04-release

── R CMD check results ──────────────────────────────────────── greed 0.6.0 ────
Duration: 6m 23.9s

❯ checking installed package size ... NOTE
    installed size is 45.5Mb
    sub-directories of 1Mb or more:
      libs  42.8Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

### Ubuntu 20.04-devel


── R CMD check results ──────────────────────────────────────── greed 0.5.1 ────
Duration: 6m 58.3s

❯ checking installed package size ... NOTE
    installed size is 43.5Mb
    sub-directories of 1Mb or more:
      libs  42.5Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖


