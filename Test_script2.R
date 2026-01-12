
#
# When starting testing ----
#
if (FALSE){
  library(devtools)
  load_all()
  
  packageVersion("mgcv")
  # was '1.8.38'
  # after update 12.9.2022: '1.8.40' (see Appendix 2)
  
  # not needed for just running code below
  install()

  }

#
# Linear, no measurement error ----
#


# Simulate data and estimate regression 
set.seed(11)
sim <- lc_simulate(n = 30)   # also plots the data

# Perform estimation
result <- lc_linear(sim$data)

# Get best estimate fitted line       
a <- result$intercept["50%"]
b <- result$slope["50%"]
# Add regression line to the plot  
abline(a, b, col = "green3")
# Add confidence interval  
lines(y_q2.5 ~ x, data = result$plot_data, lty = "dashed", col = "green3")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green3")

# DIC
result$dic


#
# . No rows with censored data ----
#

# Simulate data and estimate regression 
set.seed(11)
# Change the thresholds so no data are censored  
sim <- lc_simulate(n = 30, threshold_1 = -5, threshold_2 = -5)   # also plots the data

# Perform estimation
result <- lc_linear(sim$data)

# Get best estimate fitted line       
a <- result$intercept["50%"]
b <- result$slope["50%"]
# Add regression line to the plot  
abline(a, b, col = "green3")
# Add confidence interval  
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green3")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green3")

# DIC
result$dic

# Using good old lm()  
# A bit different result (because data are not normalized before estimation?)  
result_lm <- lm(y~x, data = sim$data)
result$plot_data_lm <- result$plot_data
result_lm_pred <- predict(result_lm, result$plot_data_lm, se.fit = TRUE) 
result$plot_data_lm$y <- result_lm_pred$fit
result$plot_data_lm$y_lo <- result_lm_pred$fit - 1.96*result_lm_pred$se.fit
result$plot_data_lm$y_hi <- result_lm_pred$fit + 1.96*result_lm_pred$se.fit
## Add regression line to the plot  
abline(coef = coef(result_lm), col = "red", lty = "dotted")
# Add confidence interval  
lines(y_lo ~ x, data = result$plot_data_lm, lty = "dotted", col = "red")
lines(y_hi ~ x, data = result$plot_data_lm, lty = "dotted", col = "red")



#
# Linear, measurement error ----
#

#
# . Simplest case ----
#

# Simulate data and estimate regression 
set.seed(11)
sim <- lc_simulate(n = 30)   # also plots the data
sim$data$meas_error <- 5
result_me <- lc_linear(sim$data, measurement_error = "meas_error")



# Get best estimate fitted line       
a <- result_me$intercept["50%"]
b <- result_me$slope["50%"]
# Add regression line to the plot  
abline(a, b, col = "purple")
# Add confidence interval  
lines(y_lo ~ x, data = result_me$plot_data, lty = "dashed", col = "purple")
lines(y_hi ~ x, data = result_me$plot_data, lty = "dashed", col = "purple")

# DIC
result_me$dic

#
# . Additive measurement errors ----
#

error_size <- c(0.1, 2, 5, 10)

# Function for setting measurement error and then estimate regression
get_result <- function(error_size){
  sim$data$meas_error <- error_size
  lc_linear(sim$data, measurement_error = "meas_error")
}

# Estimate regression for all error sizes  
result_me_list <- purrr::map(error_size, get_result)
names(result_me_list) <- paste("Error =", error_size)

# Plot results
lc_plot(sim$data, results = result_me_list)
lc_plot(sim$data, results = result_me_list, facet = "cols")


#
# . Measurement error as percentage, log-transformed data ----
#

# Data are log-transformed -> proportional error becomes additive
# Assume that error sd is proportional, e.g. 20% of the value
#   (errorfraction = 0.2)
# Log-transforming Y to Y' = log(Y) means that 
#   standard error becomes additive with sd = exp(errorfraction)-1
# See "Theory" below!


#
# .. Test/demonstration ----
#

# Simulate data and estimate regression 
set.seed(11)
sim <- lc_simulate(n = 30, 
                   intercept = 3, slope = -0.2, sigma = 1, 
                   threshold_1 = 1.5, threshold_2 = 0.8, threshold_change = 4)
sim$data$y_orig <- exp(sim$data$y)
# plot(y_orig ~ x, data = sim$data)


# Assume 40% measurement error (fraction = 0.4)
sim$data <- sim$data %>%
  mutate(
    y_orig = exp(y),
    error_orig = y_orig*0.4,
    error_log = exp(0.4)-1
  )

# Plot on original scale
ggplot(sim$data, aes(x, y_orig)) +
  geom_pointrange(aes(ymin = y_orig - error_orig, 
                      ymax = y_orig + error_orig))

# Plot on log scale
ggplot(sim$data, aes(x, y)) +
  geom_pointrange(aes(ymin = y - error_log, 
                      ymax = y + error_log))

result_me <- lc_linear(sim$data, measurement_error = "error_log")
lc_plot(sim$data, results = result_me)

# Simulate data and estimate regression 
set.seed(11)
sim <- lc_simulate(n = 30, 
                   intercept = 3, slope = -0.2, sigma = 1, 
                   threshold_1 = 1.5, threshold_2 = 0.8, threshold_change = 4)
sim$data$y_orig <- exp(sim$data$y)
# plot(y_orig ~ x, data = sim$data)

# Assume 20% measurement error (fraction = 0.2)
sim$data <- sim$data %>%
  mutate(
    y_orig = exp(y),
    error_orig = y_orig*0.2,
    error_log = exp(0.2)-1
  )

# Plot on original scale
ggplot(sim$data, aes(x, y_orig)) +
  geom_pointrange(aes(ymin = y_orig - error_orig, 
                      ymax = y_orig + error_orig))

# Plot on log scale
ggplot(sim$data, aes(x, y)) +
  geom_pointrange(aes(ymin = y - error_log, 
                      ymax = y + error_log))

result_me <- lc_linear(sim$data, measurement_error = "error_log")

lc_plot(sim$data, results = result_me)


#
# .. same data, different percentages ----
#

# Simulate data (without measurement error given)
set.seed(11)
sim <- lc_simulate(n = 30, 
                   intercept = 3, slope = -0.2, sigma = 1, 
                   threshold_1 = 1.5, threshold_2 = 0.8, threshold_change = 4)

get_data <- function(error_percent){
  sim$data <- sim$data %>%
    mutate(
      y_orig = exp(y),
      error_orig = y_orig*error_percent/100,
      error_log = exp(error_percent/100)-1
    )
}

library(purrr)
results_me <- data.frame(
  error_percent = c(0.1, 15, 45)) %>%
  mutate(
    data = map(error_percent, get_data),
  )

# Estimate regression for all error percentages  
result_me_list <- purrr::map(results_me$data, lc_linear, measurement_error = "error_log")
names(result_me_list) <- paste("Error =", results_me$error_percent)

# Plot results
# Narrower confidence intervals with higher error_percent??
lc_plot(sim$data, results = result_me_list, facet = "wrap")

# The data look as expected  
results_me %>%
  select(error_percent, data) %>%
  tidyr::unnest(cols = c(error_percent, data)) %>%
  ggplot(aes(x, y)) +
  geom_pointrange(aes(ymin = y - error_log, 
                      ymax = y + error_log)) +
  facet_wrap(vars(error_percent))

# Get normalized data used as input to the Jags modelling
model_data <- map(result_me_list, "model_data")

# Get "observed" (uncensored) values only (number 1:O) 
jags_input <- purrr::map_dfr(model_data, ~data.frame(x = .$x[1:.$O], 
                                                     y = .$y.uncens.error[1:.$O], 
                                                     error = .$meas_error), .id = "error_percent")
# - plot
ggplot(jags_input, aes(x, y)) +
  geom_pointrange(aes(ymin = y - error, 
                      ymax = y + error)) +
  facet_wrap(vars(error_percent))

#
# . No rows with censored data ----
#

# Simulate data and estimate regression 
set.seed(11)
sim <- lc_simulate(n = 30, threshold_1 = -5, threshold_2 = -5)   # also plots the data
sim$data$meas_error <- 5
result_me <- lc_linear(sim$data, measurement_error = "meas_error")

# Get best estimate fitted line       
a <- result_me$intercept["50%"]
b <- result_me$slope["50%"]
# Add regression line to the plot  
abline(a, b, col = "purple")
# Add confidence interval  
lines(y_lo ~ x, data = result_me$plot_data, lty = "dashed", col = "purple")
lines(y_hi ~ x, data = result_me$plot_data, lty = "dashed", col = "purple")

# DIC
result_me$dic

#
# . No rows with censored data, different errors ----
#

error_size <- c(0.1, 2, 5, 10)

# Function for setting measurement error and then estimate regression
get_result <- function(error_size){
  sim$data$meas_error <- error_size
  lc_linear(sim$data, measurement_error = "meas_error")
}

# Estimate regression for all error sizes  
result_me_list <- purrr::map(error_size, get_result)
names(result_me_list) <- paste("Error =", error_size)

# Estimate regression without error
ordinary_lm <- predict.lm(lm(y~x, data = sim$data), se.fit = TRUE) 
ordinary_lm_plot_data <- data.frame(
  x = sim$data$x,
  y = ordinary_lm$fit,
  y_lo = ordinary_lm$fit - 1.96*ordinary_lm$se.fit,
  y_hi = ordinary_lm$fit + 1.96*ordinary_lm$se.fit
)

# Add to 'ordinary_lm_list' (creating 'ordinary_lm_list2')
ordinary_lm_list <- list(list(plot_data = ordinary_lm_plot_data))
names(ordinary_lm_list) <- "No meas. error"
result_me_list2 <- append(ordinary_lm_list, result_me_list)

# Plot results
lc_plot(sim$data, results = result_me_list2, facet = "wrap")


#
# Linear, actual data with proportional error (e.g. 20%) ----
#

# Get one station
data_test_orig <- subset(polybrom, station %in% "23B")

# Prepare data
# debugonce(lc_prepare)
data_test_prep <- lc_prepare(data_test_orig, 
                             x = "year",
                             y = "concentration", 
                             censored = "LOQ_flag",
                             log = TRUE,
                             keep_original_columns = TRUE)

# Plot
lc_plot(data_test_prep)



#
# . Try different percentages of measurement error ----
#

get_data <- function(error_percent, data = data_test_prep){
  data$meas_error <- exp(error_percent/100) - 1
  data
}
# get_data(15, data_test_prep)
# get_data(45, data_test_prep)

# X <- get_data(15, data_test_prep)
# ggplot(X, aes(x, y)) + 
#   geom_pointrange(aes(ymin = y-meas_error, ymax = y+meas_error))
# X <- get_data(45, data_test_prep)
# ggplot(X, aes(x, y)) + 
#   geom_pointrange(aes(ymin = y-meas_error, ymax = y+meas_error))

# Function for setting measurement error and then estimate regression
get_result <- function(data_for_analysis){
  lc_linear(data_for_analysis, measurement_error = "meas_error")
}

results_error <- data.frame(
  error_percent = c(0.1, 15, 45)) %>%
  mutate(
    data = map(error_percent, get_data),
    lc_result = map(data, lc_linear, measurement_error = "meas_error")
  )


# Estimate regression for all error sizes  
result_me_list <- purrr::map(error_percent, get_result)
names(result_me_list) <- paste("Error =", error_percent)

str(result_me_list, 1)
# Plot results
lc_plot(data_test_prep, results = result_me_list, facet = "wrap")

map_dfr(result_me_list, "slope")


result_ord <- lc_linear(data_test_prep)   
result_me <- lc_linear(data_test_prep, measurement_error = 0.2)   

lc_plot(data_test_prep)

lc_plot(data_test_prep, 
        results = list(Ordinary = result_ord,
                       Meas.error.20perc = result_me))



# Get best estimate fitted line       
a <- result_me$intercept["50%"]
b <- result_me$slope["50%"]
# Add regression line to the plot  
abline(a, b, col = "purple")
# Add confidence interval  
lines(y_lo ~ x, data = result_me$plot_data, lty = "dashed", col = "purple")
lines(y_hi ~ x, data = result_me$plot_data, lty = "dashed", col = "purple")

# DIC
result_me$dic


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

load_all()

#
# . quick tests ----
#

# Quick test that JAGS runs
# debugonce(lc_fixedsplines_tp)
test <- lc_fixedsplines_tp(data = dat_cens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                           normalize = TRUE, k = 3, initialize_only = TRUE)

# Quick test that JAGS runs with measurement error
dat_cens$error <- 1

ggplot() +
  geom_pointrange(data = dat_cens[sel_uncens,], aes(x = x, y = y, ymin = y-error, ymax = y+error)) +
  geom_point(data = dat_cens[!sel_uncens,], aes(x = x, y = cut), shape = 6)


# debugonce(lc_fixedsplines_tp)
test <- lc_fixedsplines_tp(data = dat_cens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                           measurement_error = "error", 
                           normalize = TRUE, k = 3, initialize_only = TRUE)

#
# . full tests, no measurement error ----
#

# Test without measurement error

k_values <- 1:7
results <- purrr::map(k_values, 
                      ~lc_fixedsplines_tp(data = dat_cens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                                          normalize = TRUE, k = .x, raftery = TRUE)
)
results <- purrr::map(k_values, 
                      ~lc_fixedsplines_tp(data = dat_cens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                                          normalize = TRUE, k = .x, raftery = FALSE)
)

names(results) <- paste("k =", k_values)


# Plot
lc_plot(dat_cens, threshold = "cut", results = results, facet = "wrap")

# DIC
dic1 <- purrr::map_dbl(results, "dic")
plot(k_values, dic1, type = "b")

#
# . full tests, with measurement error ----
#


results_me <- purrr::map(k_values, 
                         ~lc_fixedsplines_tp(data = dat_cens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                                             normalize = TRUE, k = .x, raftery = FALSE, measurement_error = "error")
)
names(results_me) <- paste("k =", k_values)

# Plot
lc_plot(dat_cens, threshold = "cut", results = results_me, facet = "wrap")

dic2 <- purrr::map_dbl(results_me, "dic")
plot(k_values, dic2, type = "b")

purrr::map_dbl(results_me, ~sum(.x$pd))
purrr::map_dbl(results_me, ~sum(.x$deviance))

# Compare DIC with and without measurement error  
plot(k_values, dic1, type = "b", ylim = range(dic1, dic2))
lines(k_values, dic2, type = "b", pch = 19)



#
# . actual data with proportional error (e.g. 20%) ----
#

# Get one station
data_test_orig <- subset(polybrom, station %in% "23B")

# Prepare data
# debugonce(lc_prepare)
data_test_prep <- lc_prepare(data_test_orig, 
                             x = "year",
                             y = "concentration", 
                             censored = "LOQ_flag",
                             log = TRUE,
                             keep_original_columns = TRUE)

# Plot
lc_plot(data_test_prep)

data_test_prep$meas_error <- exp(0.3) - 1

#
# . - test 'reference_x'  ----
#
# debugonce(lc_fixedsplines_tp)
test <- lc_fixedsplines_tp(data = data_test_prep,
                           normalize = FALSE, k = 3, raftery = FALSE, measurement_error = "meas_error",
                           predict_x = seq(min(data_test_prep$x), max(data_test_prep$x)),
                           reference_x = 2010)


# debugonce(get_jagam_object)
test <- lc_fixedsplines_tp(data = data_test_prep,
                           normalize = FALSE, k = 5, raftery = FALSE, measurement_error = "meas_error",
                           predict_x = seq(min(data_test_prep$x), max(data_test_prep$x)),
                           reference_x = 2010, make_data_only = T)
str(test$jagam_object$jags.data, 1)

ggplot(test$diff_data, aes(x, y)) +
  geom_ribbon(aes(ymin = y_q2.5, ymax = y_q97.5), fill = "grey70") +
  geom_line() +
  geom_point(aes(color = p_category))

ggplot(test$diff_data, aes(x, p)) +
  geom_line() + geom_point() +
  scale_y_log10() +
  geom_hline(yintercept = c(0.05, 0.01, 0.001), linetype = "dashed")


#
# . - test all k's ----  
#

k_values <- 1:5

# Note: we hd to set normalize = FALSE
results_me <- purrr::map(k_values, 
                         ~lc_fixedsplines_tp(data = data_test_prep, 
                                             normalize = FALSE, k = .x, raftery = FALSE, measurement_error = "meas_error", 
                                             predict_x = seq(min(data_test_prep$x), max((data_test_prep$x))))
                         
)
names(results_me) <- paste("k =", k_values)

# Plot
lc_plot(data_test_prep,  results = results_me, facet = "wrap")

# DIC
dic <- purrr::map_dbl(results_me, "dic")
plot(k_values, dic, type = "b")

# debugonce(lc_fixedsplines_tp)
test <- lc_fixedsplines_tp(data = data_test_prep, 
                           normalize = FALSE, k = 2, raftery = FALSE, measurement_error = "meas_error", 
                           predict_x = seq(min(data_test_prep$x), max((data_test_prep$x))), compare_with_last = TRUE)
lc_plot(data_test_prep,  results = test, facet = "wrap")
ggplot(test$diff_data, aes(x, y)) +
  geom_ribbon(aes(ymin = y_q2.5, ymax = y_q97.5), fill = "grey70") +
  geom_line()


#
# . - test all k's, no censored data ----
#

data_test_prep2 <- data_test_prep %>% filter(uncensored==1)

# Note: we hd to set normalize = FALSE
results_me2 <- purrr::map(k_values, 
                         ~lc_fixedsplines_tp(data = data_test_prep2, 
                                             normalize = FALSE, k = .x, raftery = FALSE, measurement_error = "meas_error", 
                                             predict_x = seq(min(data_test_prep$x), max((data_test_prep$x))))
                         
)
names(results_me) <- paste("k =", k_values)

# Plot
lc_plot(data_test_prep2,  results = results_me, facet = "wrap")



#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Difficult case  ----
#
# Thin plate splines, censoring
#
# This is a series that resulted in error when JAGS initialized
# The reason was that 'prob[j]' in the model became appx. = 1,
#   result in error about problem with "parent of Z"
# The solution was to wrap the 'prob[j]' definition in max(..., 0.99) 
#   I.e. replace
#   max(pnorm(cut[j], mu[n+j], tau), 0.01)
#   with
#   min(max(pnorm(cut[j], mu[n+j], tau), 0.01),0.99)
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Get one station
data_all <-  readRDS("../leftcensored_testdata/125_results_2021_04_input/125_dat_all_prep3.rds")
data_test_prep <- data_all %>%
  filter(STATION_CODE == "30B" & PARAM == "CB101")
data_test_prep$meas_error <- exp(0.3) - 1

# data_test_prep <- data_all %>%
#   filter(STATION_CODE == "19N" & PARAM == "CB118" & TISSUE_NAME %in% "Egg")
# data_test_prep$meas_error <- exp(0.3) - 1

# Plot
lc_plot(data_test_prep)

data_test_prep$meas_error <- exp(0.3) - 1

# debugonce(lc_fixedsplines_tp)
test <- lc_fixedsplines_tp(data = data_test_prep,
                           normalize = FALSE, k = 3, raftery = FALSE, measurement_error = "meas_error",
                           predict_x = seq(min(data_test_prep$x), max(data_test_prep$x)),
                           reference_x = 2010)

# runjags::failed.jags('output')
# runjags::failed.jags('data')

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# HG WWa 36B ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

data_all <-  readRDS("../../seksjon 212/Milkys2_pc/Files_from_Jupyterhub_2021/Raw_data/109_adjusted_data_2022-09-23.rds")

data_all %>%
  filter(PARAM == "HG" & TISSUE_NAME %in% "Muskel" & MYEAR == 2021) %>%
  count(STATION_CODE, VALUE_WWa_lacking = is.na(VALUE_WWa)) %>% View("2021")

data_all %>%
  filter(PARAM == "HG" & TISSUE_NAME %in% "Muskel") %>%
  group_by(STATION_CODE, VALUE_WWa_lacking = is.na(VALUE_WWa)) %>%
  summarize(n = length(unique(MYEAR))) %>%
  View("")

data_all %>%
  filter(PARAM == "HG" & TISSUE_NAME %in% "Muskel") %>%
  distinct(STATION_CODE, MYEAR, VALUE_WWa_lacking = is.na(VALUE_WWa)) %>%
  count(STATION_CODE, VALUE_WWa_lacking) %>% View()

data_test_orig <- data_all %>%
  filter(STATION_CODE == "36B" & PARAM == "HG" & TISSUE_NAME %in% "Muskel")

ggplot(data_test_orig, aes(MYEAR, VALUE_WWa))+
  geom_point()

# Prepare data
# debugonce(lc_prepare)
data_test_prep <- lc_prepare(data_test_orig, 
                             x = "MYEAR",
                             y = "VALUE_WWa", 
                             censored = "FLAG1",
                             log = TRUE,
                             keep_original_columns = TRUE)

lc_plot(data_test_prep)

# 20% measurement error
data_test_prep$meas_error <- exp(0.2) - 1

table(data_test_prep$uncensored)

# Delete the 2 values under LOQ
data_test_prep2 <- data_test_prep %>%
  filter(uncensored == 1)

# debugonce(lc_fixedsplines_tp)
model2 <- lc_fixedsplines_tp(data = data_test_prep2,
                             normalize = FALSE, k = 2, raftery = FALSE, measurement_error = "meas_error",
                             predict_x = seq(min(data_test_prep$x), max(data_test_prep$x)),
                             reference_x = 2021)
model3 <- lc_fixedsplines_tp(data = data_test_prep2,
                           normalize = FALSE, k = 3, raftery = FALSE, measurement_error = "meas_error",
                           predict_x = seq(min(data_test_prep$x), max(data_test_prep$x)),
                           reference_x = 2021)
lc_plot(data_test_prep, list(k2 = model2, k3 = model3))
lc_plot(data_test_prep, results = list(k2 = model2, k3 = model3))  

str(model3, 1)
str(model3$diff_data, 1)
View(model3$diff_data)


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Thin plate splines, no censoring ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

#
# . simulate data ----
#
set.seed(2) ## simulate some data... 
n <- 50
dat <- mgcv::gamSim(1,n=n,dist="normal",scale=1)  

# we will use only x2 and y, and x2 is renamed 'x'
dat_uncens <- dat[c("x2", "y")]
names(dat_uncens)[1] <- "x"

dat_uncens$cut <- NA
dat_uncens$uncensored <- 1

ggplot(dat_uncens, aes(x, y)) +
  geom_point()

load_all()

#
# . quick tests ----
#

# Quick test for making data
# debugonce(lc_fixedsplines_tp)
test <- lc_fixedsplines_tp(data = dat_uncens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                           normalize = TRUE, k = 3, make_data_only = TRUE)

str(test, 1)
str(test$jagam_object$jags.data, 1)

# Quick test that JAGS runs
# debugonce(lc_fixedsplines_tp)
test <- lc_fixedsplines_tp(data = dat_uncens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                           normalize = TRUE, k = 3, initialize_only = TRUE)

# Quick test that JAGS works
# debugonce(lc_fixedsplines_tp)
test <- lc_fixedsplines_tp(data = dat_uncens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                           normalize = TRUE, k = 3, raftery = FALSE)
# Plot
lc_plot(dat_uncens, threshold = "cut", results = test, facet = "wrap")

# Quick test that JAGS runs with measurement error
dat_uncens$error <- 1

ggplot() +
  geom_pointrange(data = subset(dat_uncens, uncensored == 1), aes(x = x, y = y, ymin = y-error, ymax = y+error))


test <- lc_fixedsplines_tp(data = dat_uncens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                           measurement_error = "error", 
                           normalize = TRUE, k = 5, initialize_only = TRUE)


#
# . full tests, no measurement error ----
#

# Test without measurement error

k_values <- 2:7
results <- purrr::map(2:7, 
                      ~lc_fixedsplines_tp(data = dat_uncens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                                          normalize = TRUE, k = .x, raftery = FALSE)
)

names(results) <- paste("k =", k_values)


# Plot
lc_plot(dat_uncens, threshold = "cut", results = results, facet = "wrap")

# DIC
purrr::map_dbl(results, "dic")
plot(k_values, purrr::map_dbl(results, "dic"), type = "b")

#
# . full tests, with measurement error ----
#


results_me <- purrr::map(k_values, 
                         ~lc_fixedsplines_tp(data = dat_uncens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                                             normalize = TRUE, k = .x, raftery = FALSE, measurement_error = "error")
)
names(results_me) <- paste("k =", k_values)

# Plot
lc_plot(dat_uncens, threshold = "cut", results = results_me, facet = "wrap")

dic2 <- purrr::map_dbl(results_me, "dic")
plot(k_values, dic2, type = "b")

purrr::map_dbl(results_me, ~sum(.x$pd))
purrr::map_dbl(results_me, ~sum(.x$deviance))

# Compare DIC with and without measurement error  
dic1 <- purrr::map_dbl(results, "dic")
dic2 <- purrr::map_dbl(results_me, "dic")
plot(k_values, dic1, type = "b", ylim = range(dic1, dic2))
lines(k_values, dic2, type = "b", pch = 19)





#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Parallel processing using furrr ----
#
# Using thin-plate splines
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# install.packages("furrr")
library(furrr)

# Check number cores
future::availableCores()

# Set a "plan" for how many cores to use:
plan(multisession, workers = 6)


results_me <- future_map(k_values, 
                         ~leftcensored::lc_fixedsplines_tp(data = dat_cens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                                             normalize = TRUE, k = .x, raftery = TRUE, measurement_error = "error")
)

## Explicitly close multisession workers by switching plan
plan(sequential)




library(doParallel)
cl <- makeCluster(6)
registerDoParallel(cores = 6)

my_data <- iris[which(iris[,5] != "setosa"), c(1,5)]
trials <- 3

calculate_coefficients <- function(){
  ind <- sample(100, 100, replace=TRUE)
  result1 <- glm(my_data[ind,2] ~ my_data[ind,1], family=binomial(logit))
  coefficients(result1)
}
calculate_coefficients()

result <- foreach(icount(trials), 
                  .combine=cbind, 
                  .export = c("my_data","calculate_coefficients")) %dopar%
  calculate_coefficients()

result


calculate_coefficients <- function(i){
  ind <- sample(100, 100, replace=TRUE)
  result1 <- glm(my_data[ind,2] ~ my_data[ind,1], family=binomial(logit))
  coefficients(result1)
}
calculate_coefficients()

result <- foreach(icount(trials), 
                  .combine=cbind, 
                  .export = c("my_data","calculate_coefficients")) %dopar%
  calculate_coefficients()

result

foreach(i=1:5, j=11:15) %dopar% (i+j)



result <- foreach(icount(trials), 
                  .export = c("my_data","calculate_coefficients")) %dopar%
  calculate_coefficients()

result


result <- foreach(icount(trials), 
                  .export = c("my_data")) %dopar%
  calculate_coefficients()

result


#
# .. get DIC only ----
#
get_dic <- function(k){
  result_list <- lc_fixedsplines_tp(data = dat_cens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                                    normalize = TRUE, k = k, raftery = FALSE)
  result_list$dic
}
get_dic(4)


# WORKS
test <- foreach(k_number = 2:3, 
                .export = c("dat_cens","lc_fixedsplines_tp", "get_dic"),
                .packages = c("dplyr", "mgcv", "runjags", "rjags", "splines", "stats", "leftcensored")) %dopar%
  get_dic(k_number)


#
# .. get full result only ----
#
get_full_result <- function(k){
  lc_fixedsplines_tp(data = dat_cens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                                    normalize = TRUE, k = k, raftery = FALSE)
}

# WORKS!
test <- foreach(k = 2:3, 
                .export = c("dat_cens","lc_fixedsplines_tp", "get_dic"),
                .packages = c("dplyr", "mgcv", "runjags", "rjags", "splines", "stats", "leftcensored")) %dopar%
  get_full_result(k)

purrr::map_dbl(test, "dic")


#
# .. get full result from specified data ----
#







#
# . test effect of dropping Raftery ----
# (and normalizaton)
#

# Test with measurement error, no normalization, and no Raftery check   
# Much faster and seems ok enough  
res_k5_me2 <- lc_fixedsplines_tp(data = dat_cens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                                 normalize = FALSE, k = 5, raftery = FALSE, measurement_error = "error")

# Test again res_k5 (no measurement error) but with no normlization and Raftery
# Slightly wider confidence interval at the far right side
res_k5_noraft <- lc_fixedsplines_tp(data = dat_cens, x = "x", y = "y", uncensored = "uncensored", threshold = "cut",
                                    normalize = FALSE, k = 5, raftery = FALSE)

# Big effect on DIC though....
res_k3$dic
res_k4$dic
res_k5$dic
res_k5_noraft$dic

sum(res_k3$deviance)
sum(res_k4$deviance)
sum(res_k5$deviance)
sum(res_k5_noraft$deviance)

ggplot() +
  geom_ribbon(data = res_k5_me2$plot_data, aes(x= x, ymin = y_q2.5, ymax = y_q97.5), fill = "blue", alpha = 0.5) +
  geom_ribbon(data = res_k5$plot_data, aes(x= x, ymin = y_q2.5, ymax = y_q97.5), fill = "red", alpha = 0.5) +
  geom_ribbon(data = res_k5_noraft$plot_data, aes(x= x, ymin = y_q2.5, ymax = y_q97.5), fill = "green", alpha = 0.5) +
  geom_point(data = dat_cens[sel_uncens,], aes(x = x, y = y)) +
  geom_point(data = dat_cens[!sel_uncens,], aes(x = x, y = cut), shape = 6)


#
# . real data ----
#

# Get one station
data_test_orig <- subset(polybrom, station %in% "23B")

# Prepare data
# debugonce(lc_prepare)
data_test_prep <- lc_prepare(data_test_orig, 
                             x = "year",
                             y = "concentration", 
                             censored = "LOQ_flag",
                             log = TRUE,
                             keep_original_columns = TRUE)

# Plot
lc_plot(data_test_prep)  


# Assume 30% measurement error  
data_test_prep$error <- log(1.3)

debugonce(lc_fixedsplines_tp)
# normalize = TRUE gives problems in the JAGS process
results_me <- purrr::map(3:7, 
                         ~lc_fixedsplines_tp(data = data_test_prep, 
                                             normalize = FALSE, k = .x, raftery = TRUE, measurement_error = "error")
)
names(results_me) <- paste("k =", 3:7)

saveRDS(results_me, "../../temp/leftcensored_Test_script2_resut_me.rds")

# result_me_lin <- lc_linear(data_test_prep, measurement_error = "meas_error")


# Plot
lc_plot(data_test_prep, results = results_me, facet = "wrap")

# DIC
purrr::map_dbl(results_me, "dic")



#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Show that normalizing Y means that standard error is divided by sd(Y)
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Random variable with mean 100 and sd 10
n <- 5000
y <- rep(100, n)
err <- rnorm(n, 0, 10)
# Observed y (ym):
ym <- y + err 
sd(ym)

# Normalize 'ym'
ym_norm <- (ym-100)/10
mean(ym_norm)
sd(ym_norm)

# Recalculate normalized ym from y and err
ym_norm2 <- (y-100)/10 + err/10

# ym_norm and ym_norm2 seems exactly the same
mean(ym_norm2)
sd(ym_norm2)

# More 'proof':
plot(ym_norm[1:50])
points(ym_norm2[1:50], pch = 18, col = 2)

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Theory: assume that error sd is proportional, e.g. 20% of the value
#   (errorfraction = 0.2)
# Show that log-transforming Y to Y' = log(Y) means that 
#   standard error becomes additive with sd = exp(errorfraction)-1
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Random variable with mean 100 and sd 10
n <- 10000
y_mean <- rep(100, n)
sigma <- rnorm(n, 0, 10)
# Observed y (ym):
y <- y_mean + sigma 
sd(y)

logy <- log(y)
sd(logy)
log(sd(y))

# Random variable with mean 100 and sd 10
y_actual <- runif(1000, 20, 200)
errorfraction = 0.2
sd <- y_actual*errorfraction
error <- rnorm(1000, 0, sd)
y_obs <- y_actual + error

# Plot SD relative to value
plot(y_actual, sd)
plot(y_obs, sd)

# Plot some points withiut and with error
plot(head(y_actual, 50))
points(head(y_obs, 50), pch = 18, col = 2)

# Plot all 
plot(y_actual, y_obs)

# Check some large numbers
sel <- y_obs > 150
sd(y_obs[sel] - y_actual[sel])
175*0.2  # the "expected", quite close

# Check some smaller numbers
sel <- y_obs < 40
sd(y_obs[sel] - y_actual[sel])
30*0.2  # the "expected", quite close

# Log transform  
log_y_actual <- log(y_actual)
log_y_obs <- log(y_obs)

# Plot some points withiut and with error 
plot(head(log_y_actual, 50))
points(head(log_y_obs, 50), pch = 18, col = 2)

# Plot all 
plot(log_y_actual, log_y_obs)

# Expected additive approximate errors
error_logscale <- rnorm(1000, 0, exp(errorfraction) - 1)
log_y_obs2 <- log_y_actual + error_logscale

# Plot all 
points(log_y_actual, log_y_obs2, pch = 18, col = 2)

# Check sigma from variation around line - very close  
mod1 <- lm(log_y_obs ~ log_y_actual)
summary(mod1)$sigma
mod2 <- lm(log_y_obs2 ~ log_y_actual)
summary(mod2)$sigma
# [1] 0.2172967
# [1] 0.2218082



#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Preprocessing ----
#
# - Following Rob's rules 
# - See "Refinements" on https://dome.ices.dk/ohat/trDocuments/2021/help_methods_less_thans.html  
#
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o


#
# . - test rule 1+2 ----
#

#
# Using test data
#

# Make 'base' test data to use for the testthat tests
test_data_base <- data.frame(
  x = rep(2009:2020, each = 3),
  threshold = NA,
  uncensored = 1
)
test_data_base$y <- 3 + (test_data_base$x-2009)*0.5 + rnorm(12*3)

# saveRDS(test_data_base, "tests/testthat/fixtures/test_data_clean.rds")

lc_plot(test_data_base)  
load_all()

# debugonce(lc_clean)
test <- lc_clean1(test_data_base)  

# set_x_to_censored() is defined in tests/testthat/fixtures

test_data_1 <- test_data_base %>% set_x_to_censored(c(2010:2011, 2013, 2016))
test_data_2 <- test_data_base %>% set_x_to_censored(c(2010:2011, 2013, 2016:2019))
test_data_3 <- test_data_base %>% set_x_to_censored(c(2010:2016))
test_data_4 <- test_data_base %>% set_x_to_censored(c(2009, 2012:2016))
test_data_5 <- test_data_base %>% set_x_to_censored(c(2009:2013))
lc_plot(test_data_1)  

#
# Test lc_clean1
#

test <- lc_clean1(test_data_1)  
# No changes made to the data (lc_clean1)

test <- lc_clean1(test_data_2)  
# Deleted data with x < 2019 (30 rows of data)
# Final data has 2 unique x values (1 with uncensored data), and 6 rows of data.

test <- lc_flag1(test_data_2)  

test <- lc_clean1(test_data_3)  
# Deleted data with x < 2013 (12 rows of data)
# Final data has 8 unique x values (4 with uncensored data), and 24 rows of data.

# Deleted data with x < 2013 (12 rows of data)

# c(2010:2011, 2013, 2016:2019)
# Deleted data with x < 2013 (12 rows of data)


#
# Test lc_clean2
#

test <- lc_clean2(test_data_1)  
test <- lc_clean2(test_data_2)  
test <- lc_clean2(test_data_3)  
# No changes made to the data (lc_clean1)

test <- lc_clean2(test_data_4)  
test <- lc_clean2(test_data_5)  


#
# . - test rule 3: setting the last years to constant ----
#
# With real data
#

dat_all <- readRDS("../../seksjon 212/Milkys2_pc//Files_from_Jupyterhub_2021/Raw_data/109_adjusted_data_2022-06-04.rds")

dat <- dat_all %>%
  filter(PARAM %in% "HBCDB", 
         STATION_CODE %in% "13B")


# Prepare data
# debugonce(lc_prepare)
data_test_prep <- lc_prepare(dat, 
                             x = "MYEAR",
                             y = "VALUE_WW", 
                             censored = "FLAG1",
                             log = TRUE,
                             keep_original_columns = TRUE)
data_test_prep$meas_error <- exp(0.3) - 1

# Plot
lc_plot(data_test_prep, sampleinfo = TRUE)  

last_year_with_uncens <- max(subset(data_test_prep, uncensored %in% 1)$x)

# debugonce(lc_fixedsplines_tp)
test <- lc_fixedsplines_tp(data = data_test_prep,
                           normalize = FALSE, k = 2, raftery = FALSE, measurement_error = "meas_error",
                           predict_x = seq(min(data_test_prep$x), max(data_test_prep$x)),
                           reference_x = last_year_with_uncens, set_last_equal_x = last_year_with_uncens)

test <- lc_fixedsplines_tp(data = data_test_prep,
                           normalize = FALSE, k = 1, raftery = FALSE, measurement_error = "meas_error",
                           predict_x = seq(min(data_test_prep$x), max(data_test_prep$x)),
                           reference_x = last_year_with_uncens, set_last_equal_x = last_year_with_uncens)


ggplot(test$plot_data, aes(x, y)) +
  geom_ribbon(aes(ymin = y_q2.5, ymax = y_q97.5), fill = "grey70") +
  geom_line()

ggplot(test$diff_data, aes(x, y)) +
  geom_ribbon(aes(ymin = y_q2.5, ymax = y_q97.5), fill = "grey70") +
  geom_line() +
  geom_point(aes(color = p_category))





#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Cubic splines, no measurement error  ----   
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Cubic splines + Qi method version with simulated data

#
# Simulate strongly non-linear data
#

# Function for simulation
lc_simulate_nonlin <- function(n = 50, n_knots = 4, linear_part = 0.2,
                               constant = 2,
                               sigma = 0.15, 
                               censored_fraction = 0.25){
  x <- seq(0, 10, length = n) # Sort as it makes the plotted lines neater
  knots <- seq(0, 10, length = n_knots)
  B <- t(splines::bs(x, knots=knots, degree=3, intercept = TRUE)) # creating the B-splines (dim 7 x 81)
  num_data <- length(X)   # 81
  num_basis <- nrow(B)    # 7
  a0 <- linear_part # 'intercept'
  a <- rnorm(num_basis, 0, 1) # coefficients of B-splines
  Y_true <- as.vector(a0*x + a%*%B) + constant # generating the output
  Y <- Y_true + rnorm(length(x), 0, sigma) # adding noise
  dat_sim <- data.frame(x = x, y_uncensored = Y, y_true = Y_true)
  censoring_limit <- quantile(Y, censored_fraction)
  dat_sim %>%
    mutate(
      y = case_when(
        y_uncensored >= censoring_limit ~ y_uncensored,
        TRUE ~ as.numeric(NA)),
      threshold = case_when(
        y_uncensored < censoring_limit ~ y_uncensored,
        TRUE ~ as.numeric(NA)),
      uncensored = case_when(
        y_uncensored >= censoring_limit ~ 1,
        TRUE ~ 0)
    )
}
# debugonce(lc_simulate_nonlin)

# Simulate data
set.seed(13)
dat_sim <- lc_simulate_nonlin(linear_part = 0)

# Plot simulated data
lc_plot(dat_sim, y_true = "y_true")

# Estimate model
result_5knots <- lc_fixedsplines(dat_sim, knots = 5)

# Plot 
lc_plot(dat_sim, 
        y_true = "y_true", 
        results = list(Nonlin_5knots = result_5knots))




#
# Cubic splines, measurement error  ----   
#

dat_sim$meas_error <- 0.4

# debugonce(lc_fixedsplines_qi_measerror)
result_5knots_me <- lc_fixedsplines(dat_sim, knots = 5, measurement_error = "meas_error")

lc_plot(dat_sim, 
        y_true = "y_true", 
        results = list(No_error = result_5knots,
                       'Error = 0.4' = result_5knots_me), facet = "wrap")



#
# . no rows with censored data  ----   
#

# Simulate data
set.seed(13)
dat_sim_uncens <- lc_simulate_nonlin(linear_part = 0, censored_fraction = 0)
dat_sim_uncens$meas_error <- 0.4
# lc_plot(dat_sim)

# debugonce(lc_fixedsplines_qi_measerror)
result_5knots_me_uncens <- lc_fixedsplines(dat_sim_uncens, knots = 5, measurement_error = "meas_error")

lc_plot(dat_sim, 
        y_true = "y_true", 
        results = result_5knots_me_uncens, facet = "wrap")



#
# Cubic splines, test DIC values ----
#
# Does work as expected
#
# debugonce(lc_fixedsplines_qi)
result_linear <- lc_linear(dat_sim)
result_2knots <- lc_fixedsplines(dat_sim, knots = 2)
result_3knots <- lc_fixedsplines(dat_sim, knots = 3)
result_4knots <- lc_fixedsplines(dat_sim, knots = 4)
result_5knots <- lc_fixedsplines(dat_sim, knots = 5)
dic_values <- data.frame(
  Model = c("Linear", "2 knots", "3 knots", "4 knots", "5 knots"),
  DIC <- c(result_linear_qi$dic, result_2knots_qi$dic, result_3knots_qi$dic,
           result_4knots_qi$dic, result_5knots_qi$dic)
)
barplot(dic_values$DIC, names.arg = dic_values$Model)

#
# . plot data ---- 
#
#   true model (if existing), and fitted line of model(s)
#

lc_plot(dat_sim, y_true = "y_true", results = result_3knots_qi)
lc_plot(dat_sim, y_true = "y_true", results = list(result_3knots_qi))
lc_plot(dat_sim, 
        y_true = "y_true", 
        results = list(Linear = result_linear_qi,
                       Nonin_2knots = result_2knots_qi,
                       Nonlin_3knots = result_3knots_qi,
                       Nonlin_4knots = result_4knots_qi,
                       Nonlin_5knots = result_5knots_qi))




#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Appendix 1: Test log-likelihood ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

test <- dat_cens[15:20,]

# y[i] ~ dnorm(mu[i], total_var[i]^-1)       ## response

#   Z[j] ~ dbern(prob[j])
#   prob[j] <- max(pnorm(cut[j], mu[n+j], tau), 0.01)

fitted <- c(9,4,5,7,3.5,4)  

# uncensored part, test
dnorm(x = 8.73, mean = 9, sd = 2)

# uncensored part, all
i <- c(1,2,4)
log(dnorm(x = test$y[i], mean = fitted[i], sd = 2))
sum(log(dnorm(x = test$y[i], mean = fitted[i], sd = 2)))

# sensored part, test
prob <- pnorm(q = 4, mean = 1, sd = 2)
mc2d::dbern(1, prob)

# sensored part, all
i <- c(3,5,6)
prob <- pnorm(q = test$cut[i], mean = fitted[i], sd = 2)
prob
mc2d::dbern(1, prob)
log(mc2d::dbern(1, prob))


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Appendix 2: Model and data output of mgcv::jagam in mgcv_1.8-38 (Windows) vs. 1.8-39 (Linux)  ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# The jagam reults differ only for k = 3, not k = 4
# 
# Results of jagam for 
# y_comb ~ s(x, bs = "tp", k = 3)

#
#  mgcv_1.8-38 (Windows)
#

# model {                                                        
#   mu <- X %*% b ## expected response                           
#   for (i in 1:n) { y[i] ~ dnorm(mu[i],tau) } ## response      
#   scale <- 1/tau ## convert tau to standard GLM scale          
#   tau ~ dgamma(.05,.005) ## precision parameter prior          
#   ## Parametric effect priors CHECK tau=1/21^2 is appropriate!
#   for (i in 1:1) { b[i] ~ dnorm(0,0.0022) }                    
#   ## prior for s(x)...                                         
#   K1 <- S1[1:2,1:2] * lambda[1]  + S1[1:2,3:4] * lambda[2]    
#   b[2:3] ~ dmnorm(zero[2:3],K1)                                
#   ## smoothing parameter priors CHECK...                       
#   for (i in 1:2) {                                            
#     lambda[i] ~ dgamma(.05,.005)                               
#     rho[i] <- log(lambda[i])                                   
#   }                                                           
# }  

# Data:
#   List of 5
# $ y   : num [1:316(1d)] 0.262 -0.844 0.916 -0.713 -1.309 ...
# $ n   : int 316
# $ X   : num [1:316, 1:3] 1 1 1 1 1 1 1 1 1 1 ...
# ..- attr(*, "dimnames")=List of 2
# $ S1  : num [1:2, 1:4] 1.73e+01 1.01e-13 1.01e-13 5.85e-28 0.00 ...
# $ zero: num [1:3] 0 0 0


#
#  mgcv_1.8-39 (Linux)
#

# model {                                                        
#   mu <- X %*% b ## expected response                          
#   for (i in 1:n) { y[i] ~ dnorm(mu[i],tau) } ## response       
#   scale <- 1/tau ## convert tau to standard GLM scale         
#   tau ~ dgamma(.05,.005) ## precision parameter prior          
#   ## Parametric effect priors CHECK tau=1/12^2 is appropriate!
#   for (i in 1:1) { b[i] ~ dnorm(0,0.0068) }                    
#   ## prior for s(x)...                                        
#   for (i in c(2)) { b[i] ~ dnorm(0, lambda[1]) }               
#   for (i in c(3)) { b[i] ~ dnorm(0, lambda[2]) }              
#   ## smoothing parameter priors CHECK...                       
#   for (i in 1:2) {                                            
#     lambda[i] ~ dgamma(.05,.005)                               
#     rho[i] <- log(lambda[i])                                  
#   }   
  
  # Data: 
  #   List of 3
  # $ y: num [1:80(1d)] 0.24 0.305 0.249 0.98 3.31 0.685 0.604 0.302 0.341 0.284 ...
  # $ n: int 80
  # $ X: num [1:80, 1:3] 1 1 1 1 1 1 1 1 1 1 ...
  
  
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Appendix 3: Rob's rule 2 using dplyr  ----
#
# https://dplyr.tidyverse.org/reference/cumall.html
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

library(dplyr)

test <- data.frame(
  a = rep(1:2, each = 10),
  i = rep(1:10, 2),
  x = c(1,0,0,0,0,1,0,1,0,1,
        0,1,1,0,0,0,1,0,1,0))

test %>%
  group_by(a) %>%
  arrange(a, desc(i)) %>%
  mutate(cummean = cummean(x)) %>%
  arrange(a, i) %>%
  mutate(Test2 = cumany(cummean >= 0.5))
  


