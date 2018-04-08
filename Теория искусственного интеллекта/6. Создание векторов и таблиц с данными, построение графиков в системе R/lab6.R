handle_vector <- function (x, k = length(x) %% 10) {
  sorted_x <- sort(x)
  intervals <- seq(min(x), max(x), length = k + 1)
  for(interval_index in seq_along(intervals)) {
    if (interval_index < length(intervals)) {
      print(intervals[interval_index])
    }
  }
}