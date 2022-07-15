library("EpiNow2")
library("mpx.nowcasting")
library("dplyr")
library("ggplot2")

df <- load_data()
delays <- df |>
  mutate(diff = as.integer(report_date - date_onset)) |>
  filter(!is.na(diff)) |>
  pull(diff)
max_delay <- max(delays)
snapshots <- create_snapshots(df, max_delay)

est <- estimate_truncation(
  snapshots,
  max_truncation = max_delay,
  chains = 2, iter = 2000
)
