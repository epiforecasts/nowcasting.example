##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title n/a
##' @return a data frame
##' @importFrom here here
##' @importFrom dplyr mutate
##' @author Sebastian Funk
load_data <- function() {
  df <- readRDS(here::here("data", "linelist_sf.rds")) |>
    dplyr::mutate(delay = as.integer(report_date - date_onset))
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
##' @importFrom dplyr arrange count pull rename filter
##' @importFrom tidyr complete
##' @author Sebastian Funk
create_snapshots <- function(x, max_delay) {
  ## remove NA and sort
  x <- x |>
    dplyr::arrange(report_date, date_onset)
  min_report_date <- min(x$date_onset) + max_delay
  report_dates <- x |>
    dplyr::filter(report_date >= min_report_date) |>
    dplyr::pull(report_date) |>
    unique()
  ## generate snapshots
  snapshots <- purrr::map(
    report_dates,
    ~ x |>
      dplyr::filter(report_date <= .x) |>
      dplyr::count(date_onset, name = "confirm") |>
      tidyr::complete(
        date_onset = seq(
          min(date_onset),
          max(.x, min(date_onset) + max_delay),
          by = "day"),
        fill =  list(confirm = 0)
      ) |>
      dplyr::rename(date = date_onset)
    )
  return(snapshots)
}
