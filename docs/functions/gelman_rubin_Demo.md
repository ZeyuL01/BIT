# `gelman_rubin_Demo`

gelman-rubin demo


## Description

Gelman-Rubin convergence diagnostic for Demo dataset, recommend for N >= 5000


## Usage

```r
gelman_rubin_Demo(
  N = 5000,
  dat = c("CTCF", "ZBTB7A"),
  nchains = parallel::detectCores() - 1
)
```


## Arguments

Argument      |Description
------------- |----------------
`N`     |     number of rounds
`dat`     |     pre-processed Demo data, including CTCF, ZBTB7A. For details please check the manual for the data.
`nchains`     |     number of chains.


## Value

A list object contains gelman-rubin diagnostic results.


