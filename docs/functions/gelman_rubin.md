# `gelman_rubin`

gelman-rubin convergence diagnostic


## Description

gelman-rubin convergence diagnostic based on parallel computation, used chains can be varied but must less or equal to the usable cores.


## Usage

```r
gelman_rubin(
  file,
  format = c("bed", "narrowPeak", "broadPeak", "bigNarrowPeak"),
  N = 2000,
  bin_width = 1000,
  nchains = parallel::detectCores() - 1
)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     file path to the user-input.
`format`     |     format can be .bed, .narrowPeak, .broadPeak and .bigNarrowPeak.
`N`     |     number of iterations for gibbs sampler, recommended for > 5000.
`bin_width`     |     desired width of bin, should be in 100/500/1000.
`nchains`     |     number of chains used in diagnostic.


## Value

A list object contains multiple mcmc chains.


