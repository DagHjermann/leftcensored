% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/leftcensored_clean.R
\name{lc_flag2}
\alias{lc_flag2}
\title{Flag data before running trend analyses, using Rob's rule 2}
\usage{
lc_flag2(data, show_result = TRUE)
}
\arguments{
\item{data}{}
}
\value{
The same data set, with added TRUE/FALSE variable "Rule1"
}
\description{
Rule 2. If a linear/smooth trend is fitted, the first year must be non-censored
}
\examples{

# Simulate some data
testdata <- data.frame(
  x = rep(2009:2020, each = 3),
  threshold = NA,
  uncensored = 1)

# Set the  data for the first 7 years to be left-censored   
sel <- testdata$x \%in\% c(2010:2016)
testdata$threshold[sel] <- testdata$y[sel]
testdata$y[sel] <- NA
testdata$uncensored[sel] <- 0

testdata_flagged <- lc_flag2(testdata)
}
