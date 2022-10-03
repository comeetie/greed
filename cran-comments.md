Dear cran members.

This submission correspond to version 0.6.1 of package greed. This is minor modification to get rid of deprecated << initialization in Armadillo following RcppArmadillo advices (see  https://github.com/RcppCore/RcppArmadillo/issues/391). The check results obtained are similar to those of the previous version, note on package size and non-ascii characters. 

My previous submission failed the pre-tests because of URL problems (a pnas link was not responding). I rechecked all the URL with urlchecker and all URLs are correct locally. In the pre-test a one note about detritus in temp directory also appears with  windows-devel. I thus have run a check with rhub on the same platform (windows-devel) which did not raise errors.


## Test environments

* ubuntu 20.04 (on github-action), R 4.2.1
* ubuntu 20.04 (on github-action), R-devel 
* windows server 2019 (on github-action), R 4.2.1
* macOS-latest (on github-action), R 4.2.1
* windows-devel (on github-action), R-devel
* fedora-clang (on rhub), R-devel
* windows-devel (on winbuilder), R-devel

## R CMD check results

I still have one note on installed package size that comes from my usage of RcppArmadillo and on UTF-8 usage on some system.

### fedora-clang-devel
── R CMD check results ──────────────────────────────────────── greed 0.6.1 ────

❯ checking installed package size ... NOTE
    installed size is  24.8Mb
    sub-directories of 1Mb or more:
      libs   22.4Mb
❯ checking data for non-ASCII characters ... NOTE
  Note: found 989 marked UTF-8 strings
    
0 errors ✔ | 0 warnings ✔ | 2 note ✖



### macOS-release
 ── R CMD check results ──────────────────────────────────────── greed 0.6.1 ────
Duration: 8m 30.6s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

### windows-release


── R CMD check results ──────────────────────────────────────── greed 0.6.1 ────
Duration: 7m 4.8s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

### windows-devel

── R CMD check results ──────────────────────────────────────── greed 0.6.1 ────
Duration: 7m 53.3s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

### Ubuntu 20.04-release

 ── R CMD check results ──────────────────────────────────────── greed 0.6.1 ────
Duration: 8m 43.3s

❯ checking installed package size ... NOTE
    installed size is 45.4Mb
    sub-directories of 1Mb or more:
      libs  42.9Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

### Ubuntu 20.04-devel

 ── R CMD check results ──────────────────────────────────────── greed 0.6.1 ────
Duration: 6m 6.6s

❯ checking installed package size ... NOTE
    installed size is 45.6Mb
    sub-directories of 1Mb or more:
      libs  42.9Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

