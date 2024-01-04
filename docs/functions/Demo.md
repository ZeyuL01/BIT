# `Demo`

Demo to show BIMTR, which can be used without loading the ChIP-seq data.


## Description

Demo to show BIMTR, which can be used without loading the ChIP-seq data.


## Usage

```r
Demo(N, dat = c("CTCF", "ZBTB7A"))
```


## Arguments

Argument      |Description
------------- |----------------
`N`     |     number of iterations for gibbs sampler, recommended for >= 5000.
`dat`     |     pre-processed Demo data, including CTCF, ZBTB7A. For details please check the manual for the data.


## Value

A list object contains results of Gibbs sampler.


