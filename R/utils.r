##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title n/a
##' @return a data frame
##' @importFrom here here
##' @author Sebastian Funk
load_data <- function() {
  df <- readRDS(here::here("data", "linelist_sf.rds"))
  return(df)
}

##' Create data snapshots from a reporting date
##'
##' .. content for \details{} ..
##' @title n/a
##' @param x data frame containing `report_date` and `date_onset` columns
##' @param max_delay maximum delay
##' @return a list of data frames to feed into `estimate_truncation`
##' @importFrom purrr map
##' @importFrom dplyr arrange count
##' @author Sebastian Funk
create_snapshots <- function(x, max_delay) {
  ## remove NA and sort
  x <- x |>
    filter(!is.na(date_onset), !is.na(report_date),
           report_date >= date_onset) |>
    dplyr::arrange(report_date, date_onset)
  min_report_date <- min(x$date_onset) + max_delay
  report_dates <- x |>
    filter(report_date >= min_report_date) |>
    pull(report_date) |>
    unique()
  ## generate snapshots
  snapshots <- purrr::map(
    report_dates,
    ~ x |>
      filter(report_date <= .x) |>
      dplyr::count(date_onset, name = "confirm") |>
      tidyr::complete(
        date_onset = seq(min(date_onset), max(date_onset), by = "day"),
        fill =  list(confirm = 0)
      ) |>
      rename(date = date_onset)
    )
  return(snapshots)
}
