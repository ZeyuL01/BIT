# `import_peaks`

Transform the input file to vector


## Description

Transform the input file to vector


## Usage

```r
import_peaks(
  file,
  format = c("bed", "narrowPeak", "broadPeak", "bigNarrowPeak", "csv"),
  bin_width = 1000
)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     file path to the user-input.
`format`     |     format can be .bed, .narrowPeak, .broadPeak and .bigNarrowPeak.
`bin_width`     |     desired width of bin, should be in 100/500/1000.


## Value

A numeric vector contains the index of peaks with pre-specified number of bins in each chromosome.


