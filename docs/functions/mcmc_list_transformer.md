# `mcmc_list_transformer`

mcmc_list_transformer


## Description

function to transform parallel chains into a mcmc list.


## Usage

```r
mcmc_list_transformer(dat, label, nchains)
```


## Arguments

Argument      |Description
------------- |----------------
`dat`     |     input results from parallel chains.
`label`     |     label of TFs.
`nchains`     |     number of chains.


## Value

A list object contains multiple mcmc chains.


