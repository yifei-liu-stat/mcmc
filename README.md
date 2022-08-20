# mcmc
Bayesian analysis of sports league data from 2019–2020 NFL Season;
Implementation of [STAT8056 doing assignment 4](https://www.stat.umn.edu/geyer/8054/hw/nfl.html) by Yifei Liu

* `.Rhistory` and `.Rprofile` are used for convenience only for my PC.
* `mcmchw.R` contains all R codes I used for this assignment.
* `mout.RData`, `mout_Q1.RData` and `mout_Q2.RData` contain R objects `mout`, `mout_1` and `mout_2` returned from `mcmc::metrop()`,
  I saved them so that I can save compiling time of `8054.Rmd` simply by loading them into R workspace when they are needed, with
  * `mout` for the chain of $\boldsymbol \beta$
  * `mout_1` for the chain in Question 1 (replay the final game), and
  * `mout_2` for the chain in Question 2 (replay the whole playoffs)
* `8054hw.Rmd` is used for generating the report `8054hw.html`, one can reproduce the work via:
  ```R
  library(rmarkdown)
  render("8054hw.Rmd", html_document(number_sections = TRUE, toc = TRUE))
  ```
  and make sure you have aforementioned three `.RData` files in the same work directory before the compiling.
 
