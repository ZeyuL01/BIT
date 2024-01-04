# `filter_peaks`

filter peaks based on CREs, Promoter, Enhancer


## Description

filter peaks based on CREs, Promoter, Enhancer


## Usage

```r
filter_peaks(
  input_vec,
  option = c("ALL", "CREs", "CLS", "ELS"),
  bin_width = c(100, 500, 1000)
)
```


## Arguments

Argument      |Description
------------- |----------------
`input_vec`     |     input vector that contains indices from import_peaks
`option`     |     option, ALL: Non-filtering, CREs: All cis-regulatory elements, PLS: Promoters, ELS: Enhancers.
`bin_width`     |     bin width


## Value

a vector contains filtered bin indices.


