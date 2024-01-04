# `BIMTR_multi`

BIMTR parallel computation.


## Description

the multi-cores parallel version of BIMTR to acquire parallel MCMC chains for gelman-rubin diagnostic.


## Usage

```r
BIMTR_multi(
  file,
  format = c("bed", "narrowPeak", "broadPeak", "bigNarrowPeak"),
  N = 5000,
  bin_width = c(100, 500, 1000),
  option = c("ALL", "CREs", "PLS", "ELS"),
  numCores
)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     file path to the user-input.
`format`     |     format can be .bed, .narrowPeak, .broadPeak and .bigNarrowPeak.
`N`     |     number of iterations for gibbs sampler, recommended for > 5000.
`bin_width`     |     desired width of bin, should be in 100/500/1000.
`option`     |     option to filter peaks with candidate cis-regulatory elements from ENCODE.
`numCores`     |     number of cores used in parallel computation.


## Value

A list object of list objects that contains results of Gibbs sampler.


