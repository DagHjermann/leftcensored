#' Linear regression for interval-censored data
#' 
#' This function runs linear regression when the dependent (y) variable is left-censored. That is,
#' for some or all measurements, we only know that the actual value is below some upper bound (y_up) 
#' and above some lower bound (y_lo). The model does not take measurement error into account. The model
#' formulation of Qi et al. (2022) is used.
#' 
#' @param data The data must have (at least) four columns: the predictor variable, the (uncensored) response variable, 
#' a variable which is 1 for every uncensored observation, and a variable with the threhold of censoring for every censored observation
#' @param x Variable name for the predictor (independent) variable        
#' @param y_lo Variable name for the lower bound of the response (dependent) variable    
#' @param y_up Variable name for the lower bound of the response (dependent) variable    
#' @param uncensored Variable name for a variable which is 1 for uncensored values (for chemical measurements, above LOQ) 
#' and 0 for censored values (for chemical measurements: less than LOQ, but above some other value). 
#' values.     
#' @param resolution The number of points along the x axis used to describe the spline. 
#' @param n.chains The number of MCMC chains (replicates) to run. The default is 4. Using more than 1 chain enables us to say whether 
#' @param n.iter The number of iterations for each MCMC chains. The default is 5000, which is usually sufficient for this application.
#' @param n.burnin The number of iterations to remove at start of each MCMC chain, before results are collected for statistics. The default is 1000.
#' If n.burnin is too small, plots of the trace (see examples) will show whether the chains are homogeneous along the chain (there should be from 
#' no decreasing or increasing trend in the traceplot). One can also use R2jags::traceplot to assess whether the chains behave differently, 
#' depending on their different starting point. If they do, n.burnin should be increased.
#' @param n.thin The number of MCMC iterations that are kept for statistics
#' 
#' @keywords Statistics, regression, censored
#' 
#' @references 
#' Qi X, Zhou S and Plummer M. 2022. On Bayesian modeling of censored data in JAGS. BMC Bioinformatics 23: 102
#' (doi: 10.1186/s12859-021-04496-8). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8944154/  
#' 
#' Plummer N (2017). JAGS Version 4.3.0 user manual. https://sourceforge.net/projects/mcmc-jags/files/  
#' 
#' @return The function returns a list with two parts: \code{summary} and \code{model}. \code{summary} shows a summary of the result, i.e., estimates
#' of the parameters of the linear regression: intercept, slope and sigma (the estimated standard deviation of the data around the regression line).
#' It is common to use the quantiles for parameter estimates, i.e., using the  
#' 50% quantile as the "best estimate" of the parameters, and using the 2.5% and 97.5% quantiles as endpoints of a 95% confidence interval. 
#' \code{model} is the output of the jags() command, which is what \code{lc_linear()} runs under the hood for estimation. 
#' It is an MCMC object, which has methods for functions such as plot (see ?mcmc). If you have some knowledge of the MCMC technique for     
#' state-space models, this can be used for diagnostic plots of the model. See examples.
#'   
#' @examples
#' # Simulate data (using function made for left-censored data)
#' set.seed(11)
#' sim <- lc_simulate(n = 30, plot = FALSE)
#' 
#' # The upper bound is equal to 'y_plot' 
#' sim$data$y_up <- sim$data$y_plot
#' # Initially, we set lower bound to the same
#' sim$data$y_lo <- sim$data$y_up
#' # For censored data, let us say that the lower bound is 2 less than 
#' # the upper bound:  
#' sel <- sim$data$uncensored == 0
#' sim$data$y_lo[sel] <- sim$data$y_up[sel] + log(0.1)
#' 
#' # Plot data
#' plot(y_up~x, data = sim$data, ylim = range(sim$data$y_lo, sim$data$y_up), type = "n")
#' segments(x0 = sim$data$x, y0 = sim$data$y_lo, y1 = sim$data$y_up)
#' points(y_lo~x, data = sim$data, pch = 16)
#' points(y_up~x, data = sim$data, pch = 16)
#' result <- lc_linear(sim$data)
#' 
#' Estimate regression parameters:
#' result <- int_linear_qi(sim$data
#' 
#' # Get best estimates and plot its regression line on top of the plot 
#' a <- result$intercept["50%"]
#' b <- result$slope["50%"]
#' abline(a, b, col = "green4")
#' Add confidence interval  
#' lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green4")
#' lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green4")
#' 
#' @export
int_linear_qi <- function(data,
                          x = "x", 
                          y_lo = "y_lo", 
                          y_up = "y_up", 
                          uncensored = "uncensored",
                          resolution = 50,
                          n.chains = 4, 
                          n.iter = 5000, 
                          n.burnin = 1000, 
                          n.thin = 2,
                          mean_y = NULL,
                          sd_y = NULL,
                          plot_input = FALSE,
                          plot_norm = FALSE,
                          detailed = FALSE,
                          model_parameters_for_convergence = c("intercept", "slope")){
  
  # Censoring vs truncation:
  # https://stats.stackexchange.com/a/144047/13380 
  # Also see here for an alternative way to set up the model:
  # https://stats.stackexchange.com/questions/6870/censoring-truncation-in-jags
  # NOTE: CHECK THIS for multi-level statistical models:
  # https://stats.stackexchange.com/questions/185254/multi-level-bayesian-hierarchical-regression-using-rjags
  
  if (sum(data[[y_up]] < data[[y_lo]]) > 0)
    stop("Values of ", sQuote(y_up), " must be higher than the values of ", sQuote(y_lo))
  
  # Set 'uncensored' and 'y'  
  data$uncensored <- ifelse(data[[y_up]]-data[[y_lo]] < 0.000001, 1, 0)
  data$y <- data[[y_up]]
  
  # Set all censored data to NA (if not already done)
  # Important! Otherwise all LOQ stuff is ignored
  data$y[data$uncensored %in% 0] <- NA
  
  if (plot_input){
    plot(data$x, data$y, ylim = range(data$y, data$y_up, na.rm = TRUE))
    points(data$x[data$uncensored == 0], data$y_hi[data$uncensored == 0], pch = 20, col = 2)
    points(data$x[data$uncensored == 0], data$y_lo[data$uncensored == 0], pch = 20, col = 4)
  }
  
  # For making predicted lines (with SE) 
  xmin <- min(data[[x]], na.rm = TRUE)
  xmax <- max(data[[x]], na.rm = TRUE)
  x.out <- seq(xmin, xmax, length = resolution)
  
  # Jags code to fit the model to the simulated data
  
  if (sum(data$uncensored == 1) > 0){
    
    model_code = '
model
{
  # Uncensored observations 
  for (o in 1:O) {
    y.uncens[o] ~ dnorm(intercept + slope * x[o], sigma^-2)
  }
  # Censored observations 
  for (c in 1:C) {
    Z1[c] ~ dbern(p[c])
    # make sure p is never approximatly zero (below 0.01):
    p[c] <- max(p_below_hi[c]*p_above_lo[c], 0.01)     
    p_below_hi[c] <- pnorm(cut_up[c], intercept + slope * x[O+c], sigma^-2)
    p_above_lo[c] <- 1 - pnorm(cut_lo[c], intercept + slope * x[O+c], sigma^-2)
  }
  for (i in 1:resolution) {
    y.hat.out.norm[i] <- intercept + slope * x.out[i]
    y.hat.out[i] <- y.hat.out.norm[i]*sd_y + mean_y
  }
  # Priors
  intercept ~ dnorm(0, 100^-2)
  slope ~ dnorm(0, 100^-2)
  sigma ~ dunif(0, 10)
}
'
    
  }  else {
    
    model_code = '
model
{
  # No uncensored observations 
  # Censored observations 
  for (c in 1:C) {
    Z1[c] ~ dbern(p[c])
    # make sure p is never approximatly zero (below 0.01):
    p[c] <- max(p_below_hi[c]*p_above_lo[c], 0.01)      
    p_below_hi[c] <- pnorm(cut_up[c], intercept + slope * x[O+c], sigma^-2)
    p_above_lo[c] <- 1 - pnorm(cut_lo[c], intercept + slope * x[O+c], sigma^-2)
  }
  for (i in 1:resolution) {
    y.hat.out.norm[i] <- intercept + slope * x.out[i]
    y.hat.out[i] <- y.hat.out.norm[i]*sd_y + mean_y
  }
  # Priors
  intercept ~ dnorm(0, 100^-2)
  slope ~ dnorm(0, 100^-2)
  sigma ~ dunif(0, 10)
}
'
    
    
  }
  
  ### Set up data and parameters
  
  # Normalize y
  # Achieves mean = 0
  if (is.null(mean_y)){
    mean_y1 <- mean(data[[y_lo]], na.rm = TRUE)
    mean_y2 <- mean(data[[y_up]], na.rm = TRUE)
    mean_y <- mean(mean_y1, mean_y2)
  }
  if (is.null(sd_y)){
    sd_y1 <- sd(data[[y_lo]], na.rm = TRUE)
    sd_y2 <- sd(data[[y_up]], na.rm = TRUE)
    sd_y <- mean(sd_y1, sd_y2)
  }
  # Functions for normalization and "un-normalization" (back-transformation)
  norm_y <- function(x) (x-mean_y)/sd_y
  unnorm_y <- function(x) x*sd_y + mean_y
  
  # Normalize (or more correctly centralize) x
  # Achieves mean = 0
  mean_x <- mean(data[[x]], na.rm = TRUE)
  # Functions for normalization and "un-normalization" (back-transformation)
  norm_x <- function(x) x-mean_x
  
  # Set up the data for Qi's method
  # Data has normalized y values and centralized x values
  
  # Split the data into uncensored and censored parts
  data_obs <- data[data[[uncensored]] %in% 1,]
  data_cen <- data[data[[uncensored]] %in% 0,]
  data_all <- rbind(data_obs, data_cen)
  
  if (sum(data$uncensored == 1) > 0){
    
    model_data <- list(x = norm_x(data[[x]]),
                       y.uncens = norm_y(data_obs[["y"]]),
                       O = nrow(data_obs),
                       Z1 = rep(1, nrow(data_cen)),  # because all are left-censored, see text below 'Model 2' in Qi' et al. 2022's paper
                       cut_up = norm_y(data_cen[[y_up]]),
                       cut_lo = norm_y(data_cen[[y_lo]]),
                       C = nrow(data_cen),
                       x.out = norm_x(x.out),
                       resolution = resolution,
                       mean_y = mean_y,
                       sd_y = sd_y)
    
  } else {
    
    model_data <- list(x = norm_x(data[[x]]),
                       O = nrow(data_obs),
                       Z1 = rep(1, nrow(data_cen)),  # because all are left-censored, see text below 'Model 2' in Qi' et al. 2022's paper
                       cut_up = norm_y(data_cen[[y_up]]),
                       cut_lo = norm_y(data_cen[[y_lo]]),
                       C = nrow(data_cen),
                       x.out = norm_x(x.out),
                       resolution = resolution,
                       mean_y = mean_y,
                       sd_y = sd_y)
    
  }
  
  if (plot_norm){
    # plot(model_data$x, model_data$y.uncens, ylim = range(model_data$y.uncens, model_data$cut_up, na.rm = TRUE))
    # points(model_data$x[model_data$uncensored == 0], model_data$cut_up, pch = 20, col = 2)
    # points(model_data$x[model_data$uncensored == 0], model_data$cut_lo, pch = 20, col = 4)
    # Un-normalized data:
    # plot(data_obs$x, data_obs$y, ylim = range(data_obs$y, data_cen$cut_up, na.rm = TRUE))
    # points(data_cen$x, data_cen$cut_up, pch = 20, col = 2)
    # points(data_cen$x, data_cen$cut_lo, pch = 20, col = 4)
  }
  
  
  
  # Choose the parameters to watch
  if (detailed){
    model_parameters <-  c("intercept", "slope", "sigma", 
                           "y.hat.out.norm", "y.hat.out", "y.uncens")
  } else {
    model_parameters <-  c("intercept", "slope", "sigma", "y.hat.out")
  }
  
  # Initial values  
  if (sum(data$uncensored == 1) > 0){
    init_model_df <- data.frame(
      x = model_data$x, y = c(model_data$y.uncens, model_data$cut_up))
  } else {
    init_model_df <- data.frame(
      x = model_data$x, y = model_data$cut_up)
  }
  init_model <- lm(y ~ x, data = init_model_df)
  init_summ <- summary(init_model)$coef
  jags.inits <- function(){
    list("intercept" = rnorm(1, mean = init_summ[1,1], sd = init_summ[1,2]), 
         "slope" =  rnorm(1, mean = init_summ[2,1], sd = init_summ[2,2]),
         "sigma" = runif(1))
  }
  
  #
  ### Run model
  #
  
  # Alt. 1. Run the model using R2jags::jags
  
  # model_run <- R2jags::jags(
  #   data = model_data,
  #   init = jags.inits,
  #   parameters.to.save = model_parameters,
  #   model.file=textConnection(model_code),
  #   n.chains=n.chains,   # Number of different starting positions
  #   n.iter = n.iter,     # Number of iterations
  #   n.burnin = n.burnin, # Number of iterations to remove at start
  #   n.thin = n.thin)     # Amount of thinning
  
  # R2jags::jags workflow continues with:  
  #   model_mcmc <- coda::as.mcmc(model_run)
  #   summary <- summary(model_mcmc)
  
  # Alt. 2: Run the model using rjags::jags.model
  #   As used in 'Binomial Data' in
  #   https://github.com/xinyue-qi/Censored-Data-in-JAGS/blob/main/R_program.md
  # parameters.to.save - specified in coda.samples
  # n.burnin           - specified in coda.samples 
  # n.thin             - specified in coda.samples
  
  # Choose the parameters to watch
  model_parameters <-  c('intercept', 'slope', 'sigma', 'y.hat.out')
  
  ### Run model
  # Initial run, using just sigma and dic 
  model_converged <- runjags::autorun.jags(
    data = model_data,
    monitor = model_parameters_for_convergence,     
    inits = jags.inits,
    model = model_code,
    n.chains = n.chains,    # Number of different starting positions
    startsample = 4000,     # Number of iterations
    startburnin = n.burnin, # Number of iterations to remove at start
    thin = n.thin)          # Amount of thinning
  
  # Add all model parameters and get samples for them
  model_result <- runjags::extend.jags(model_converged, 
                                       add.monitor = model_parameters,
                                       sample = n.iter)
  
  # model_result
  model_mcmc <- coda::as.mcmc(model_result)
  
  summary <- summary(model_mcmc)
  
  #
  # DIC
  #
  dic.pd <- rjags::dic.samples(model = runjags::as.jags(model_result), n.iter=1000, type="pD")
  
  # Not used now:
  # dic.popt <- dic.samples(model=model_run, n.iter=30000, type="popt"); dic.popt
  
  # Select the observations for which we got penalties
  dic.sel.pd <- !is.nan(dic.pd$penalty )
  
  # Get penalties and deviances for those
  pd <- dic.pd$penalty[dic.sel.pd]
  deviance <- dic.pd$deviance[dic.sel.pd]
  
  # Calculate DIC
  dic <- sum(deviance) + sum(pd)
  
  #
  # Get predicted line 
  #
  quants <- summary$quantiles
  length.out <- length(x.out)
  pick_rownames <- sprintf("y.hat.out[%i]", 1:length.out)
  # y and lower and upper CI  values are back-transformed (un-normalized) using unnorm:
  plot_data <- data.frame(
    x = x.out, 
    y = quants[pick_rownames,"50%"],
    y_lo = quants[pick_rownames,"2.5%"],
    y_hi = quants[pick_rownames,"97.5%"]  )
  
  # Get regression coefficients, originals:
  intercept.norm <- summary$quantiles["intercept",]
  slope.norm <- summary$quantiles["slope",]
  #
  # Get regression coefficients, back-transformed
  intercept <- intercept.norm*sd_y - slope.norm*sd_y*mean_x + mean_y
  slope <- slope.norm*sd_y
  #
  # Basis for formulae above:
  #   y' = (y - mean_y)/sd_y   (1)
  #   x' = x - mean_x          (2)
  # Slope formula used for normalized data:
  #   y' = a' + b'x'           (3)
  #   - where a' and b' are the intercept and slope found for normalized data
  # To get the formulae used above, substitute (1) and (2) into (3):
  #   (y - mean_y)/sd_y = a' + b'(x - mean_x)               (4)
  # And solve for y on the left side. This results in
  #   y = [a'*sd_y - b'*sd_y*mean_x + mean_y] + [b*sd_y]*x  (5)
  # Where the two parentheses are the back-transformed intercept and slope, respectively
  #
  
  list(summary = summary,
       plot_data = plot_data,
       model = model_mcmc,
       model_from_jags = model_result,
       intercept = intercept,
       slope = slope,
       mean_y = mean_y,
       sd_y = sd_y,
       norm_y = norm_y,
       dic_all = dic.pd,
       dic = dic)  
  
}
