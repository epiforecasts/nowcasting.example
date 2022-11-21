##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title n/a
##' @return a data frame
##' @importFrom here here
##' @importFrom dplyr mutate
##' @author Sebastian Funk
##' @export
load_data <- function() {
  df <- readRDS(here::here("data", "linelist.rds")) |>
    dplyr::mutate(delay = as.integer(report_date - date_onset))
  return(df)
}

##' Create data snapshots from a reporting date
##'
##' .. content for \details{} ..
##' @title n/a
##' @param x data frame containing `report_date` and `date_onset` columns
##' @param max_delay maximum delay
##' @param second name of the secondary (reporting) column
##' @param first name of the primary (outcome) column
##' @return a list of data frames to feed into `estimate_truncation`
##' @importFrom purrr map
##' @importFrom dplyr arrange count pull rename filter
##' @importFrom tidyr complete
##' @importFrom rlang sym
##' @author Sebastian Funk
##' @export
create_snapshots <- function(x, max_delay, second = "report_date",
                             first = "date_onset") {
  ## remove NA and sort
  x <- x |>
    dplyr::rename(..first = sym(!!first),
                  ..second = sym(!!second)) |>
    dplyr::arrange(..second, ..first)
  min_report_date <- min(x$`..first`) + max_delay
  report_dates <- x |>
    dplyr::filter(..second >= min_report_date) |>
    dplyr::pull(..second) |>
    unique()
  ## generate snapshots
  snapshots <- purrr::map(
    report_dates,
    ~ x |>
      dplyr::filter(..second <= .x) |>
      dplyr::count(..first, name = "confirm") |>
      tidyr::complete(
        `..first` = seq(
          min(..first),
          max(.x, min(..first) + max_delay - 1),
          by = "day"),
        fill = list(confirm = 0)
      ) |>
      dplyr::rename(date = `..first`)
    )
  return(snapshots)
}
