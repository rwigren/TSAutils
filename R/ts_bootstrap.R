#' Block-wise bootstrap sample
#'
#' Gets a single sample according to the block-wise bootstrap scheme. Used in the blockwise_bootstrap function.
#'
#' @param series A stationary time series from which the block-wise sampling will be done. Can be univariate or multivariate.
#' @param l Length of the blocks
#' @param sim.length Length of the new, sampled series
#' @examples
#' series <- arima.sim(model=list(ar=c(0.6, -0.4)), n=100)
#' block_BS_sample(series, 5)
#' @export
block_BS_sample <- function(series, l, sim.length=NROW(series)){
  k <- ceiling(sim.length / l) # Amount of blocks to be sampled

  samples <- sample(0:(sim.length-l), k, replace=T) # Sample block start indeces
  samples <- samples %x% rep(1, l) # Repeat each l times through kronecker product
  samples <- samples + rep(1, k) %x% 1:(l) # Add increasing stuff, to fix indices.
  samples <- samples[1:sim.length]

  return(as.matrix(series)[samples, ])
}

#' Block-wise bootstrap
#'
#' Procedure for evaluating a statistic on a stationary time series.
#'
#' Returns a bootstrap distribution of the chosen statistic.
#'
#' The implementation follows the block-wise bootstrap method as presented by Künsch (1989).
#' Resamples blocks of the time series and uses these to create new time series, which the
#' statistic is evaluated on. Can for example be used to estimate the variance or confidence
#' intervals of a statistic (e.g. mean of the series).
#'
#' Note that the implemented method is not the 'block of blocks'-bootstrap but the more general
#' 'naive' block-wise bootstrap method.
#'
#' @param series A stationary time series from which the block-wise sampling will be done.
#' Can be univariate or multivariate.
#' @param statistic Function: the statistic for which a bootstrap-distribution will be created.
#' @param B Amount of bootstrap samples. Choice depends on statistic and purpose
#' @param l Length of blocks in the procedure. Good values depend on the series and the statistic,
#' however, good default values ... n^(1/3) (ref Künsch) or 2n^(1/3) ref (Bühlmann).
#' @param sim.length Length of the sampled series. Will usually be the same length as the original series,
#' but can in some cases (such as residual bootstrapping) be useful to be set to another value.
#' @param ... Additional values passed on to the statistic function.
#' @keywords bootstrap blockwise block
#' @examples
#' series <- arima.sim(model=list(ar=c(0.6, -0.4)), n=100)
#' bw.res <- blockwise_bootstrap(series, mean, B=500, l=5)
#' boxplot(bw.res)
#' @references Kunsch, Hans R. (1989) "The Jackknife and the Bootstrap for General Stationary Observations".
#' \emph{The Annals of Statistics}, \bold{17}(3), 1217--1241.
#'
#' @seealso \code{\link{ARsieve_bootstrap}}
#' @export
blockwise_bootstrap <- function(series, statistic, B, l, sim.length=NROW(series), ...){
  sample_stats <- NULL
  for (i in 1:B){
    sample.data <- block_BS_sample(series, l, sim.length)
    sample_stats <- rbind(sample_stats, statistic(sample.data, ...))
  }
  return(sample_stats)
}


AR_innov_sample <- function(ar.fit){
  innovations <- ar.fit$resid
  p <- ar.fit$order
  n_innov <- length(innovations)

  sample_fn <- function(n){
    samples <- sample((p+1):n_innov, n, replace = T)
    return(innovations[samples])
  }
  return(sample_fn)
}

#' AR-sieve sample
#'
#' Gets a single sample according to the AR-sieve bootstrap scheme.
#'
#' @param mean Mean of the series
#' @param ar.coef Coefficients of the AR process
#' @param sim.length Length of the new, sampled series
#' @examples
#' series <- arima.sim(model=list(ar=c(0.6, -0.4)), n=100)
#' block_BS_sample(series, 5)
#' @export
ARsieve_BS_sample <- function(mean, ar.coef, innov.fn, sim.length){
  return(arima.sim(model=list(ar=ar.coef), rand.gen=innov.fn, n=sim.length) + mean)
}

#' AR-sieve bootstrap
#'
#' Bootstrapping method for linear stationary univariate time series. Returns a bootstrap distribution of the chosen statistic.
#'
#' The implementation follows the AR-sieve bootstrap method as presented by Bühlmann (1997).
#' Fits an AR process based on AIC minimization and creates samples based on fitted coefficients
#' and empirical residuals. The statistic is then evaluated on these new samples. Can for example
#' be used to estimate the variance or confidence intervals of a statistic (e.g. auto-regressive
#' coefficients of the series).
#'
#' @param series A univariate linear stationary time series from which the AR-sieve sampling will be done.
#' @param statistic Function: the statistic for which a bootstrap-distribution will be created.
#' @param B Amount of bootstrap samples.
#' @param sim.length Length of the sampled series. Will usually be the same length as the original series,
#' but can in some cases (such as residual bootstrapping) be useful to be set to another value.
#' @param max.p Maximum order the fitted AR model can take.
#' @param ar.method Method used to fit the auto-regressive process. ()
#' @param ... Extra values passed on to the statistic function.
#' @keywords bootstrap AR-sieve sieve ARsieve
#' @examples
#' series <- arima.sim(model=list(ar=c(0.6, -0.4)), n=100)
#' ars.res <- ARsieve_bootstrap(series, mean, B=500)
#' boxplot(ars.res)
#' @references Bühlmann, Peter (1997) "Sieve bootstrap for time series".
#' \emph{Bernoulli}, \bold{3}(2), 123--148.
#' @seealso \code{\link{blockwise_bootstrap}}, \code{\link{stats::ar}}
#' @export
ARsieve_bootstrap <- function(series, statistic, B, sim.length=NROW(series),
                              max.p=10*log10(NROW(series)), ar.method='yw', ...){
  ar.fit <- ar(series, order.max = max.p, method=ar.method)
  innov.fn <- AR_innov_sample(ar.fit)

  sample_stats <- NULL
  for (i in 1:B){
    sample.data <- ARsieve_BS_sample(mean=ar.fit$x.mean, ar.coef=ar.fit$ar, innov.fn=innov.fn, sim.length=sim.length)
    sample_stats <- rbind(sample_stats, statistic(sample.data, ...))
  }
  return(sample_stats)
}

#' Monte-Carlo simulation
#'
#' Method to evaluate the distribution of a statistic, dependent on the distribution of the data.
#' Returns B samples from the distribution of the statistic.
#'
#' Equal to drawing samples from the theoretical distribution of the statistic given the distribution of the data.
#' This is useful when the theoretical distribution of the statistic, given the data distribution, does not have
#' an analytical solution or to save time by computing an empirical distribution rather than calculating the theoretical.
#'
#' Can for example be used to evaluate other bootstrap methods accuracy in a simulation study.
#'
#' @param sim.fn Simulation function (distribution) from which the time series is drawn from. Takes only one argument
#' representing the length of the drawn sample
#' @param statistic Function: the statistic for which a sample-distribution will be created.
#' @param B Amount of samples.
#' @param sim.length Length of the sampled series.
#' @param ... Additional values passed on to the statistic function.
#' @keywords monte carlo
#' @export
#' @examples
#' arima_sim <- function(n){arima.sim(model=list(ar=c(0.6, -0.4)), n=n)}
#' mcs.res <- monte_carlo_sim(series, mean, 5)
#'
#' series <- arima_sim(100)
#' ars.res <- ARsieve_bootstrap(series, mean, B=500)
#' boxplot(cbind(mcs.res, ars.res))
monte_carlo_sim <- function(sim.fn, statistic, B, sim.length, ...){
  sample_stats <- NULL
  for (i in 1:B){
    sample.data <- sim.fn(sim.length)
    sample_stats <- rbind(sample_stats, statistic(sample.data, ...))
  }
  return(sample_stats)
}
