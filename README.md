
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Nowcasting examples

This repository contains two example approaches to nowcasting, one using
the [EpiNow2](https://epiforecasts.io/EpiNow2/) package, and one using
the [epinowcast](https://epiforecasts.io/epinowcast/) package which is
currently being developed to eventually replace `EpiNow2` for real-time
applications.

## Installation

You can install the utility functions contained in this repo with

``` r
# install.packages("devtools")
devtools::install_github("epiforecasts/nowcasting.example")
```

## Render the `EpiNow2` report

To render the report, run

``` r
rmarkdown::render("inst/reports/epinow2.Rmd")
```

Afterwards the report can be viewed in `inst/reports/epinow2.html`. A
[rendered version](inst/reports/epinow2.md) of the report is also
available in the repository.

## Render the `epinowcast` report

To render the report, run

``` r
rmarkdown::render("inst/reports/epinowcast.Rmd")
```

Afterwards the report can be viewed in `inst/reports/epinowcast.html`. A
[rendered version](inst/reports/epinowcast.md) of the report is also
available in the repository.

## Create mock dataset

A mock dataset can be created using

``` r
source("inst/scripts/create_mock_dataset.r")
```
