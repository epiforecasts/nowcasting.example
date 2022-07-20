
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Nowcasting example

## Installation

You can install the utility functions contained in this repo with

``` r
# install.packages("devtools")
devtools::install_github("epiforecasts/nowcasting.example")
```

## Render report

To render the report, run

``` r
rmarkdown::render("inst/reports/est_trunc_epinow2.Rmd")
```

Afterwards the report can be viewed in
`inst/reports/est_trunc_epinow2.html`. A rendered version of the report
is also available in
[inst/reports/est_trunc_epinow2.md](the%20repository).

## Create mock dataset

A mock dataset can be created using

``` r
source("inst/scripts/create_mock_dataset.r")
```
