% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predicted_y_quantiles_onechain.R
\name{predicted_y_quantiles_onechain}
\alias{predicted_y_quantiles_onechain}
\title{Get quantiles for y (given x), mean quantiles for all chains}
\usage{
predicted_y_quantiles_onechain(
  x,
  mcmc,
  chain_no,
  quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)
)
}
\arguments{
\item{x}{x value}

\item{mcmc}{MCMC object}

\item{chain_no}{Chain number}

\item{quantiles}{Quantiles. Default is c(0.025, 0.25, 0.5, 0.75, 0.975)}
}
\value{
Matrix of one column and (with default quantiles) 5 rows, containing quantiles for the given x value
}
\description{
Get quantiles for y (given x), mean quantiles for all chains
}
\examples{
sim <- lc_simulate()
result <- lc_linear(sim$data)
predicted_y_quantiles_onechain(2000, result$model, 1)

}
