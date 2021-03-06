% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lc_simulate.R
\name{lc_simulate}
\alias{lc_simulate}
\title{Simulate censored data}
\usage{
lc_simulate(
  n = 100,
  intercept = 30,
  slope = -3,
  sigma = 4,
  threshold_1 = 20,
  threshold_2 = 10,
  threshold_change = 4,
  plot = TRUE,
  seed = NULL
)
}
\arguments{
\item{n}{Number of simulated observations}

\item{intercept}{Intercept of the linear relationship between y and x}

\item{slope}{Slope of the linear relationship between y and x}

\item{sigma}{Standard deviance around the linear relationship (the s.d. of y for a given x)}

\item{threshold_1}{threshold when x <= threshold_change}

\item{threshold_2}{threshold when x > threshold_change}

\item{threshold_change}{The level of x when threshold changes from threshold_1 to threshold_2}

\item{plot}{If TRUE, makes a plot of the data}

\item{seed}{If you want the simulation to be reproducible, set seed to some integer}
}
\value{
A data frame with 5 columns is returned invisibly. The variables are 1) x; 2) y (the real, unobserved value);
3) y_uncens (the observed y values, which equals y for values above threshold, and is NA for values under threshold);
4) uncensored (equals 0 for values below threshold and 1 for values above threshold); and 5) threshold (the threshold
value for each observation).
}
\description{
Simulates data where the real value y depends linearly on x, with random variation added.
The values are censored, so we only observe the actual value above some threshold.
The threshold has two different levels depending on x. This is a typical situation if y is chemical
measurements and x is time in years, as chemical methods tends to improve (decrease the limit of detection) over the years.
}
\examples{

# Default parameters
lc_simulate()

# Modified parameters
sim <- lc_simulate(slope = -1.5, threshold_1 = 27, threshold_2 = 20)

# The data object of the output can be directly used in lm_leftcensored 
result <- lm_linear(sim$data)

}
\seealso{
\code{\link{lm_leftcensored}}
}
