Dear cran members.

I corrected the two remarks raised by my previous submission. The modification are described bellow. Thanks in advance for your help.

Remark : Please write the reference in the description text in form author(s)
(year) <arXiv:2002:11577>
Response : The reference was update in Readme.md to the published version of the paper and in the correct format.

Remark : Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\value{No return value, called for side effects} or similar)
Missing Rd-tags:
      print-icl_path-method.Rd: \value
Response : The value tag was added to the print-icl_path-method Rd file with the following description (None (invisible NULL). No return value, called for side effects.).

## Test environments

* ubuntu 20.04 (on github-action), R 4.0.5 
* ubuntu 20.04 (on github-action), R-devel (2021-04-06 r80146)
* windows (on github-action), R 4.0.5
* macOS darwin 17 (on github-action), R 4.0.5 

## R CMD check results

I had one note on installed package size that comes from my usage of RcppArmadillo.

### macOS-release

── R CMD check results ──────────────────────────────────────── greed 0.5.0 ────
Duration: 7m 9.9s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

### windows-release
── R CMD check results ──────────────────────────────────────── greed 0.5.0 ────
Duration: 13m 24.6s

❯ checking installed package size ... NOTE
    installed size is 5.5Mb
    sub-directories of 1Mb or more:
      libs  4.6Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

### Ubuntu 20.04-release

── R CMD check results ──────────────────────────────────────── greed 0.5.0 ────
Duration: 6m 23.5s

❯ checking installed package size ... NOTE
    installed size is 43.5Mb
    sub-directories of 1Mb or more:
      libs  42.5Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

### Ubuntu 20.04-devel


── R CMD check results ──────────────────────────────────────── greed 0.5.0 ────
Duration: 6m 58.3s

❯ checking installed package size ... NOTE
    installed size is 43.5Mb
    sub-directories of 1Mb or more:
      libs  42.5Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖


