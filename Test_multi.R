
library(devtools)
load_all()

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Thin plate splines ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

#
# . simulate data ----
#
set.seed(2) ## simulate some data... 
n <- 50
dat <- mgcv::gamSim(1,n=n,dist="normal",scale=1)  

# we will use only x2 and y, and x2 is renamed 'x'
dat <- dat[c("x2", "y")]
names(dat)[1] <- "x"

ggplot(dat, aes(x, y)) +
  geom_point()

# Here: fixed threshold (but that is not necessary)
thresh <- 4
dat_cens <- dat[c("x","y")]
dat_cens$y_orig <- dat_cens$y       # original (will not be used)
sel_uncens <- dat_cens$y > thresh
dat_cens$y[!sel_uncens] <- NA
dat_cens$cut <- thresh
dat_cens$cut[sel_uncens] <- NA
dat_cens$uncensored <- 0
dat_cens$uncensored[sel_uncens] <- 1

# Data
ggplot() +
  geom_point(data = dat_cens[sel_uncens,], aes(x = x, y = y)) +
  geom_point(data = dat_cens[!sel_uncens,], aes(x = x, y = cut), shape = 6)


#
# . quick tests ----
#

# Quick test that JAGS runs
# debugonce(lc_fixedsplines_tp)
test <- lc_fixedsplines_tp(data = dat_cens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                           normalize = TRUE, k = 3, initialize_only = TRUE)

#
# Back to basic jagam ----
#

## the following illustrates a typical workflow. To run the 
## 'Not run' code you need rjags (and JAGS) to be installed.
require(mgcv)

set.seed(2) ## simulate some data... 
n <- 50
dat <- gamSim(1,n=n,dist="normal",scale=2)
## regular gam fit for comparison...
b0 <- gam(y ~ s(x2), data=dat, method="REML")

plot(b0, res = T, pch = 19, cex = 1)


## Set directory and file name for file containing jags code.
## In real use you would *never* use tempdir() for this. It is
## only done here to keep CRAN happy, and avoid any chance of
## an accidental overwrite. Instead you would use
## setwd() to set an appropriate working directory in which
## to write the file, and just set the file name to what you
## want to call it (e.g. "test.jags" here). 

jags.file <- paste(tempdir(),"/test.jags",sep="") 

## Set up JAGS code and data. In this one might want to diagonalize
## to use conjugate samplers. Usually call 'setwd' first, to set
## directory in which model file ("test.jags") will be written.



jd <- jagam(y ~ s(x2), data=dat, file=jags.file,
            sp.prior="gamma", diagonalize=TRUE)

## In normal use the model in "test.jags" would now be edited to add 
## the non-standard stochastic elements that require use of JAGS....

## Not run: 
require(rjags)
load.module("glm") ## improved samplers for GLMs often worth loading
jm <-jags.model(jags.file,data=jd$jags.data,inits=jd$jags.ini,n.chains=1)
# list.samplers(jm)
sam <- jags.samples(jm, c("b","rho","scale"), n.iter=10000, thin=10)
jam <- sim2jam(sam,jd$pregam)
plot(jam,pages=1)
jam
pd <- data.frame(x0=c(.5,.6),x1=c(.4,.2),x2=c(.8,.4),x3=c(.1,.1))
fv <- predict(jam,newdata=pd)
## and some minimal checking...
require(coda)
effectiveSize(as.mcmc.list(sam$b))

## End(Not run)

#
# Code directly from 'jags.file' (file name 'test.jags')
#
jagscode_txt <- '
model {
  mu <- X %*% b ## expected response
  for (i in 1:n) { 
    y[i] ~ dnorm(mu[i],tau)  ## response
  }
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/88^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00013) }
  ## prior for s(x2)... 
  for (i in c(2:9)) { b[i] ~ dnorm(0, lambda[1]) }
  for (i in c(10)) { b[i] ~ dnorm(0, lambda[2]) }
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

str(jd$jags.data, 1)
jm <-jags.model(textConnection(jagscode_txt), data = jd$jags.data, inits=jd$jags.ini, n.chains=1)


#
# Code directly from 'jags.file' (file name 'test.jags')
#
jagscode_txt <- '
model {
  mu <- X %*% b ## expected response
  for (i in 1:n) { 
    y[i] ~ dnorm(mu[i],tau)  ## response
  }
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/88^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00013) }
  ## prior for s(x2)... 
  for (i in c(2:9)) { b[i] ~ dnorm(0, lambda[1]) }
  for (i in c(10)) { b[i] ~ dnorm(0, lambda[2]) }
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

str(jd$jags.data, 1)
jm <-jags.model(textConnection(jagscode_txt), data = jd$jags.data, inits=jd$jags.ini, n.chains=1)


#
# Linear model
#
jagscode_txt <- '
model {
  mu <- X %*% b ## expected response
  for (i in 1:n) { 
    y[i] ~ dnorm(mu[i],tau)  ## response
  }
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/88^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00013) }
  ## prior for s(x2)... 
  for (i in c(2:9)) { b[i] ~ dnorm(0, lambda[1]) }
  for (i in c(10)) { b[i] ~ dnorm(0, lambda[2]) }
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

dat$x2

str(jd$jags.data, 1)
jm <-jags.model(textConnection(jagscode_txt), data = jd$jags.data, inits=jd$jags.ini, n.chains=1)



#
# Andrew C Parnell multiple regression example ----
#

# . Header ------------------------------------------------------------------

# Fitting a multiple response regression in JAGS
# Andrew Parnell

# In this file we fit a Bayesian muliple response regression model

# Some boiler plate code to clear the workspace, and load in required packages
library(MASS) # For the multivariate normal distribution

# . Maths -------------------------------------------------------------------

# Description of the Bayesian model fitted in this file
# Notation:
# y_{i} = J-vector of response variables for observation i, i = 1, ..., N
# x_{i} = K-vector of explanatory variables for observation i
# Together Y is an N by J matrix of response variables (each observation is of dimenstion J) and X is an N by K matrix of explanatory variables (each observation has K response variables)

# Parameters
# A = intercept vector of length J
# B = slope matrix of dimension K by J, i.e. this is the matrix of slopes for explanatory variable k on dimensino j
# Sigma = residual variance matrix of dimension J by J

# Likelihood
# y[i,] ~ N(A + B*x[i,], Sigma)

# Priors - all vague
# A[j] ~ normal(0, 100)
# B[j,k] ~ normal(0,100)
# Smarter priors might tie these values, especially the slopes together
# This borrowing of strength might take place over dimensions (j) or covariates (k)

# . Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model
N <- 100
J <- 3
K <- 5
set.seed(123)
X <- matrix(rnorm(N * K), nrow = N, ncol = K)

# Simulate parameters
Sigma <- rWishart(1, df = J + 1, Sigma = diag(J))[, , 1]
A <- rnorm(J)
B <- matrix(rnorm(J * K) * 5, ncol = J, nrow = K)

# Get the means and simulate data
mean <- y <- matrix(NA, ncol = J, nrow = N)
for (i in 1:N) {
  mean[i, ] <- A + X[i, ] %*% B
  y[i, ] <- mvrnorm(1, mean[i, ], Sigma)
}

# Very hard to visualise!
pairs(data.frame(
  y = y,
  X = X
))

# . Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data
model_code <- "
model
{
  # Likelihood
  for (i in 1:N) {
    y[i, ] ~ dmnorm(mu[i, ], Sigma.Inv)
    mu[i, 1:J] <- A + X[i, ]%*%B
  }
  Sigma.Inv ~ dwish(I, J+1)
  Sigma <- inverse(Sigma.Inv)
  # Priors
  for(j in 1:J) {
    A[j] ~ dnorm(0, 100^-2)
    for(k in 1:K) {
      B[k,j] ~ dnorm(0, 100^-2)
    }
  }
}
"

# Wishhart distribution:
# install.packages("MCMCpack")
# ?MCMCpack::dwish

# Set up the data
model_data <- list(
  N = N, y = y, X = X, K = ncol(X), J = ncol(y),
  I = diag(J)
)

# Choose the parameters to watch
model_parameters <- c("A", "B", "Sigma")

# Run the model
# Change to 'runjags' - original was with R2jags 
model_run <- runjags::run.jags(
  data = model_data,
  monitor = model_parameters,
  model = model_code
)

# . Simulated results -------------------------------------------------------

# model_result
model_mcmc <- coda::as.mcmc(model_run)

summary <- summary(model_mcmc)
quants <- summary$quantiles
stats <- summary$statistics



#
# Andrew C Parnell linear regression ----
#

# . Header ------------------------------------------------------------------

# Fitting a linear regression in JAGS
# Andrew Parnell

# In this code we generate some data from a simple linear regression model and fit is using jags. We then interpret the output.


# . Maths -------------------------------------------------------------------

# Description of the Bayesian model fitted in this file
# Notation:
# y_i = repsonse variable for observation t=i,..,N
# x_i = explanatory variable for obs i
# alpha, beta = intercept and slope parameters to be estimated
# sigma = residual standard deviation

# Likelihood:
# y[i] ~ N(alpha + beta * x[i], sigma^2)
# Prior
# alpha ~ N(0,100) - vague priors
# beta ~ N(0,100)
# sigma ~ U(0,10)

# . Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model
n <- 100
alpha <- 2
beta <- 3
sigma <- 1
# Set the seed so this is repeatable
set.seed(123)
x <- sort(runif(n, 0, 10)) # Sort as it makes the plotted lines neater
y <- rnorm(n, mean = alpha + beta * x, sd = sigma)

# Also creat a plot
plot(x, y)
lines(x, alpha + beta * x)

# . Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data

model_code <- "
model
{
  # Likelihood
  for (i in 1:n) {
    y[i] ~ dnorm(alpha + beta * x[i], sigma^-2)
  }
  # Priors
  alpha ~ dnorm(0, 100^-2)
  beta ~ dnorm(0, 100^-2)
  sigma ~ dunif(0, 10)
}
"

# Set up the data
model_data <- list(n = n, y = y, x = x)

# Choose the parameters to watch
model_parameters <- c("alpha", "beta", "sigma")

# Run the model
model_run <- runjags::run.jags(
  data = model_data,
  monitor = model_parameters,
  model = model_code,
  n.chains = 4, # Number of different starting positions
  sample = 1000, # Number of iterations
  burnin = 200, # Number of iterations to remove at start
  thin = 2
) # Amount of thinning


# . Simulated results -------------------------------------------------------

# model_result
model_mcmc <- coda::as.mcmc(model_run)

summary <- summary(model_mcmc)
quants <- summary$quantiles
stats <- summary$statistics



#
# ACP - linear regr. for sum (underlying variable) ----
#
# What we observe is three variables (y) that together
#   is y_latent (the sum) 
# Proportion 
#

# . Header ------------------------------------------------------------------

# Fitting a linear regression in JAGS
# Andrew Parnell

# In this code we generate some data from a simple linear regression model and fit is using jags. We then interpret the output.


# . Maths -------------------------------------------------------------------

# Description of the Bayesian model fitted in this file
# Notation:
# y_i = repsonse variable for observation t=i,..,N
# x_i = explanatory variable for obs i
# alpha, beta = intercept and slope parameters to be estimated
# sigma = residual standard deviation

# Likelihood:
# y[i] ~ N(alpha + beta * x[i], sigma^2)
# Prior
# alpha ~ N(0,100) - vague priors
# beta ~ N(0,100)
# sigma ~ U(0,10)

# . Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model
n <- 100
alpha <- 20
beta <- 1
sigma_lat <- 3
# Set the seed so this is repeatable
set.seed(123)
x <- sort(runif(n, 0, 10)) # Sort as it makes the plotted lines neater
y_latent <- rnorm(n, mean = alpha + beta * x, sd = sigma_lat)

# Add components  
y <- matrix(nrow = n, ncol = 3)
p_comp <- c(0.2, 0.3, 0.5)
sigma_comp <- c(3,2,1)
for (i in 1:3){
  y[,i] <- rnorm(n, mean = y_latent*p_comp[i], sd = sigma_comp[i])
}

# Data frame of true parameter values    
true_df <- data.frame(
  true = c(alpha, beta, sigma_lat,
           sigma_comp, p_comp))


# Also create a plot
ylim <- c(0,35)
par(mfrow = c(3,2), mar = c(3,4,2,1))
plot(x, y_latent, ylim = ylim, main = "latent (unobserved)")
lines(x, alpha + beta * x)
for (i in 1:3){
  plot(x, y[,i], ylim = ylim, main = paste("component", i))
}
plot(x, apply(y, 1, sum), ylim = ylim, main = "sum of components")

pairs(data.frame(
  y = y,
  x = x
))


# . Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data

model_code <- "
model
{
  # Likelihood
  for (i in 1:n) {
    y_latent[i] ~ dnorm(alpha + beta * x[i], sigma_lat^-2)
    for (j in 1:3) {
      y[i,j] ~ dnorm(prop[j] * y_latent[i], sigma_comp[j]^-2)
    }
  }
  prop[3] <- 1- prop[1] - prop[2]
  # Priors
  alpha ~ dnorm(0, 100^-2)
  beta ~ dnorm(0, 100^-2)
  sigma_lat ~ dunif(0, 10)
  for (k in 1:3) {
    sigma_comp[k] ~ dunif(0, 10)
  }
  for (m in 1:2) {
    prop[m] ~ dunif(0, 1)
  }
}
"

# Set up the data
model_data <- list(n = n, y = y, x = x)

# Choose the parameters to watch
model_parameters <- c("alpha", "beta", "sigma_lat", "sigma_comp", "prop")

# Run the model
model_run <- runjags::run.jags(
  data = model_data,
  monitor = model_parameters,
  model = model_code,
  n.chains = 4, # Number of different starting positions
  sample = 1000, # Number of iterations
  burnin = 200, # Number of iterations to remove at start
  thin = 2
) # Amount of thinning


# . Simulated results -------------------------------------------------------

# model_result
model_mcmc <- coda::as.mcmc(model_run)

summary <- summary(model_mcmc)
quants <- summary$quantiles
stats <- summary$statistics
quants_df <- bind_cols(
  data.frame(parameter = rownames(quants)),
  as.data.frame(quants),
  true_df
)
stats_df <- bind_cols(
  data.frame(parameter = rownames(stats)),
  as.data.frame(stats),
  true_df
)

# Plot estimates (black) and true values (red)  
ggplot(quants_df %>% filter(parameter != "alpha"), aes(parameter, `50%`)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  geom_point(aes(y = true), color = "red", shape = 1, size = 3) +
  coord_flip() +
  ylab("Estimates (black) and true values (red)")
  


