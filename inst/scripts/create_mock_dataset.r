## Script to create a mock line list

library("EpiNow2")
library("dplyr")
library("here")

## length of time series
length <- 60
## proportion zeroes
zero <- 0.4
## proportion negative
neg <- 0.02
## proportion onsets missing
miss <- 0.1


## delay distribution parameters
meanlog <- 2
sdlog <- 0.5

## downsample example
onsets <- example_confirmed |>
  slice(seq_len(length)) |>
  mutate(confirm = rbinom(
    n = n(),
    size = confirm,
    p = 0.002
  ))

max_delay <- 26

## sampling function for truncated discretised log normal
rtdislnorm <- function(n, meanlog, sdlog, trunc = Inf) {
  res <- rep(trunc, n)
  while(any(res >= trunc)) {
    res[res >= trunc] <- floor(rlnorm(sum(res >= trunc), meanlog, sdlog))
  }
  return(res + 1L)
}

linelist <- tibble(date_onset = rep(onsets$date, times = onsets$confirm)) |>
  ## create report dates
  mutate(
    report_date = date_onset + rtdislnorm(n(), meanlog, sdlog, max_delay),
    id = 1:n()
  ) |>
  ## create zeroes
  mutate(report_date = if_else(
    id %in% sample(1:n(), round(zero * n())),
    date_onset,
    report_date
  )) |>
  ## create negative
  mutate(report_date = if_else(
    id %in% sample(1:n(), round(neg * n())),
    date_onset - (report_date - date_onset),
    report_date
  )) |>
  ## create missing
  mutate(date_onset = if_else(
    id %in% sample(1:n(), round(miss * n())),
    as.Date(NA_character_),
    date_onset
  )) |>
  dplyr::select(-id) |>
  ## truncate
  filter(report_date <= max(date_onset, na.rm = TRUE))

saveRDS(linelist, here::here("data", "linelist.rds"))
