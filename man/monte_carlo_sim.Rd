% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts_bootstrap.R
\name{monte_carlo_sim}
\alias{monte_carlo_sim}
\title{Monte-Carlo simulation}
\usage{
monte_carlo_sim(sim.fn, statistic, B, sim.length, ...)
}
\arguments{
\item{sim.fn}{Simulation function (distribution) from which the time series is drawn from. Takes only one argument
representing the length of the drawn sample}

\item{statistic}{Function: the statistic for which a sample-distribution will be created.}

\item{B}{Amount of samples.}

\item{sim.length}{Length of the sampled series.}

\item{...}{Additional values passed on to the statistic function.}
}
\description{
Method to evaluate the distribution of a statistic, dependent on the distribution of the data.
Returns B samples from the distribution of the statistic.
}
\details{
Equal to drawing samples from the theoretical distribution of the statistic given the distribution of the data.
This is useful when the theoretical distribution of the statistic, given the data distribution, does not have
an analytical solution or to save time by computing an empirical distribution rather than calculating the theoretical.

Can for example be used to evaluate other bootstrap methods accuracy in a simulation study.
}
\examples{
arima_sim <- function(n){arima.sim(model=list(ar=c(0.6, -0.4)), n=n)}
mcs.res <- monte_carlo_sim(series, mean, 5)

series <- arima_sim(100)
ars.res <- ARsieve_bootstrap(series, mean, B=500)
boxplot(cbind(mcs.res, ars.res))
}
\keyword{carlo}
\keyword{monte}
