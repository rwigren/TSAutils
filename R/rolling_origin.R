#' Rolling origin forecast evaluation
#'
#' Rolling origin forecast evaluation, a.k.a. time-series cross validation, of a model or method.
#' Computes errors and prediction of a forecast function applied to a time series according to the rolling origin scheme.
#'
#' This method implements the rolling origin forecast evaluation (see e.g. Hyndman and Athanasopoulos, 2018).
#' Returns a list of two matrices (multivariate time series), one containing the errors and one the predictions,
#' where columns represent forecast horizon and rows represent time points.
#'
#' The method starts with a subset X[1:t] and forecasts X[(t+1):(t+h)] based on this subset.
#' Then forecasts X[(t+2):(t+h+1)] based on X[1:(t+1)] and so on.
#'
#' @param series Univariate time series used for fitting and computing forecast errors.
#' @param start Time point used as the starting point for the rolling origin forecast.
#' @param forecast_fn Function which returns forecasts. Takes a time series as its first argument and has an argument
#' h representing the forecast horizon. Returns a vector of length equal to the given forecast horizon h.
#' @param h Forecast horizon. This function will evaluate all horizons 1 to h.
#' @param ... Additional arguments passed to forecast_fn
#' @keywords rolling origin evaluation
#' @examples
#' series <- arima.sim(model=list(ar=c(0.6, -0.4)), n=150)
#' arma_forecast <- function(x, h) {return(
#'    forecast(arima(x, order=c(2,0,0)), h=h)$mean
#' )}
#' roe.res <- rolling_origin_eval(series, 120, arma_forecast, h=5)
#' @references Hyndman, R.J., & Athanasopoulos, G. (2018)
#' \emph{Forecasting: principles and practice, 2nd edition}, OTexts: Melbourne, Australia.
#' OTexts.com/fpp2. Accessed on 2019-10-18
#' @export
rolling_origin_eval <- function(series, start, forecast_fn, h=1, ...){
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
