# `Show_Results`

To show the ranking table by inspecting the results of Gibbs sampler.


## Description

To show the ranking table by inspecting the results of Gibbs sampler.


## Usage

```r
Show_Results(dat, burnin)
```


## Arguments

Argument      |Description
------------- |----------------
`dat`     |     results from Main_sampling.
`burnin`     |     number of samples used for burn-in.


## Value

a list object contains two data.frame. one by ranking the theta_ij, the other one by ranking theta_i.


