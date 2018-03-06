[![Build Status](https://travis-ci.org/smwindecker/decay.svg?branch=master)](https://travis-ci.org/smwindecker/decay)
[![codecov](https://codecov.io/gh/smwindecker/decay/branch/master/graph/badge.svg)](https://codecov.io/gh/smwindecker/decay)

# decay

R package to run decay models with stan. 

Currently supports the negative exponential and weibull decay functions, the option to include random effects, fixed effects, and perform cross-validation by specified folds. 

```{r}
if (!require("devtools")) {
  install.packages("devtools")
  library("devtools")
}
install_github("smwindecker/decay")
```
