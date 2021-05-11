Dear cran members,

I corrected the problems found with valgrind. All the provided tests did not raise issues anymore with valgrind. I have performed the same tests as previously with similar results (same notes; 1 on UTF8 on Solaris, and 1 on package size on Windows and Ubuntu).  

## Test environments

* ubuntu 20.04 (on github-action), R 4.0.5 
* ubuntu 20.04 (on github-action), R-devel (2021-04-06 r80146)
* windows (on github-action), R 4.0.5
* macOS darwin 17 (on github-action), R 4.0.5 
* solaris 10 (on rhub), R 4.0.5

## R CMD check results

I still have one note on installed package size that comes from my usage of RcppArmadillo and one on UTF-8 usage on SOLARIS.

### solaris-release

── R CMD check results ──────────────────────────────────────── greed 0.5.1 ────
Duration: 6m 58.3s

❯ checking data for non-ASCII characters ... NOTE
  Note: found 6693 marked UTF-8 strings
    
0 errors ✔ | 0 warnings ✔ | 1 note ✖

### macOS-release

── R CMD check results ──────────────────────────────────────── greed 0.5.1 ────
Duration: 8m 5.1s

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


── R CMD check results ──────────────────────────────────────── greed 0.5.1 ────
Duration: 9m 33.5s

❯ checking installed package size ... NOTE
    installed size is 48.1Mb
    sub-directories of 1Mb or more:
      libs  45.8Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

### Ubuntu 20.04-devel

 ── R CMD check results ──────────────────────────────────────── greed 0.5.1 ────
Duration: 6m 59.9s

❯ checking installed package size ... NOTE
    installed size is 48.1Mb
    sub-directories of 1Mb or more:
      libs  45.8Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖


