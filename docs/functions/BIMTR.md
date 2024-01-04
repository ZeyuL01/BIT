# `BIMTR`

BIMTR


## Description

Main interface to run BIMTR method, please set the input file path, input file format, number of iterations and bin width.
 users can also change default parameters used in Gibbs sampler.


## Usage

```r
BIMTR(
  file,
  format = c("bed", "narrowPeak", "broadPeak", "bigNarrowPeak", "csv"),
  N = 5000,
  bin_width = c(100, 500, 1000),
  option = c("ALL", "CREs", "PLS", "ELS")
)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     file path to the user-input.
`format`     |     format can be .bed, .narrowPeak, .broadPeak and .bigNarrowPeak.
`N`     |     number of iterations for gibbs sampler, recommended for >= 5000.
`bin_width`     |     desired width of bin, should be in 100/500/1000.
`option`     |     option to filter peaks with candidate cis-regulatory elements from ENCODE.


## Value

A list object contains results of Gibbs sampler.


