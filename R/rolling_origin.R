#' A test function
#'
#' This function allows you to express your love of cats.
#' @param a integer wallah
#' @param b integer brah
#' @keywords add
#' @export
#' @examples
#' test(1, 2)
test <- function(a, b, ...){ return(a +  b)}

#' Rolling origin forecast evaluation, a.k.a. time-series cross validation.
#'
#' This function does something
#' @param series series used for ...
#' @param start starting origin for the rolling origin forecast
#' @param forecast_fn forecasting function of the series
#' @param h forecasting horizon. Will evaluate all ... from 1 to h
#' @keywords rolling origin evaluation
#' @export
#' @examples
#' test(1, 2)
rolling_origin_forecast <- function(series, start, forecast_fn, h=1, ...){
  time_points <- time(window(series, start=start, end=end(series)))
  forecasts <- NULL
  for (t in head(time_points, -1)){
    H <- (frequency(series) * (tail(time_points, n=1) - t)) %>%
          round() %>%
          min(h)

    forecasts <- fn(window(series, end = t), h=H, ...) %>%
      c(rep(NA, h-H), .) %>%
      rbind(forecasts, .)
  }

  forecasts <- ts(forecasts, end=end(series), frequency = frequency(series))

  for (i in 1:h){
    forecasts[, i] <-
      lag(forecasts[, i], k=-i+1) %>% # Shifts forecast to correct row
      window(end = end(series)) %>%   # Cuts off unwanted values
      c(rep(NA, i-1), .)              # Fills non-forecasted entries with NA
  }
  colnames(forecasts) <- paste('h=', 1:h, sep = '')
  errors <- forecasts - series
  colnames(errors) <- colnames(forecasts)
  return(list(forecasts=forecasts, errors=errors))
}
