# TSAutils
The package implements some useful functions in time series analysis, mainly for evaluation of time series models. Examples include rolling origin evaluation and bootstrap methods for time series.

## Installation
You can install the package with
```
library(remotes)
remotes::install_github("rwigren/TSAutils")
```
or with
```
library(devtools)
devtools::install_github("rwigren/TSAutils")
```

## To be added:
* Multivariate rolling origin evaluation
* Easy evaulation of multiple models on the same series
* Exhaustive search of SARIMA model order by information criteria
    + E.g.: Get k model orders with smalles AIC 
* Possibly some simulation functions for:
    + Time-series regression
    + Multivariate linear processes (VARMA)



## List of functions 
