elapsed_months <- function(end_date, start_date) {
  ed <- as.Date(end_date)
  sd <- as.Date(start_date)

  (12 * (as.numeric(substr(ed, 1,4)) - as.numeric(substr(sd, 1,4)))) +
      (as.numeric(substr(ed, 6,7)) - as.numeric(substr(sd, 6,7)))
}
