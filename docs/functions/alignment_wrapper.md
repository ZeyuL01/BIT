# `alignment_wrapper`

Count the 'good' and 'total' informative cases between input peak set with reference database.


## Description

Count the 'good' and 'total' informative cases between input peak set with reference database.


## Usage

```r
alignment_wrapper(input_vec, bin_width, option)
```


## Arguments

Argument      |Description
------------- |----------------
`input_vec`     |     A input vector contains index of peaks transformed by applying import_peaks.
`bin_width`     |     width of bin, should be in 100/500/1000.


## Value

a data frame has three columns, TF labels, number of 'good' windows, number of 'total' informative cases.


