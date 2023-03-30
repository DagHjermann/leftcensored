
#
# Start up ----
#

#
# Note 5.8.2022: For some reason, after updating R and rjags,
#   rjags will no longer start in RStudio. 
#   However it starts in base R.

# Error message:
#   library(rjags)
#   Error: package or namespace load failed for ‘rjags’ in get(method, envir = envir):
#    lazy-load database 'C:/R/Library/rjags/R/rjags.rdb' is corrupt
#   In addition: Warning message:
#   In get(method, envir = envir) : internal error -3 in R_decompress1


# For Base R use (outside RStudio):
# setwd("C:/Data/R_packages/leftcensored")

# Also see "zeroes trick", e.g. here:
# https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/17f237b2/

library(devtools)
# library(ggplot2)
# library(mgcv)
# library(splines)
load_all()

#
# Checking ..... -----
#
# adding packages / fixing functions / fixing documentation between checks  
#

check()
# use_package('rjags')
# use_package('dplyr')
# use_package('ggplot2', 'suggests')
# use_mit_license()
# check()
# use_package('ggplot2')
# check()
# use_package('splines')
# use_package('graphics')
# use_package('stats')

# Update documentation (including examples)  
document()

# Update functions
load_all()

# Final check
check()

# Installing package  
install()

# help(package = "leftcensored")
# help(lc_linear)

#
# Linear regression (LR) ----
# 
# Works ok
#

# Simulate data and estimate regression 
set.seed(11)
sim <- lc_simulate(n = 30)   # also plots the data
# debugonce(lc_linear)
result <- lc_linear(sim$data)


# Get best estimate fitted line       
a <- result$intercept["50%"]
b <- result$slope["50%"]
# Add regression line to the plot  
abline(a, b, col = "green3")
# Add confidence interval  
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green3")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green3")

# The result is a list of 9 things:  
str(result, 1)
str(result, 2)
str(result$summary, 1)
head(result$summary$statistics)
head(result$summary$quantiles)
rownames(result$summary$quantiles)

pick_rownames <- sprintf("y.hat.out[%i]", 1:length.out)

#' # Make a standard MCMC plot: the trace and the density for each estimated parameter  
#' par(mar = c(2,4,3,1))
#' plot(result$model)


#
# . Qi version ----
#

set.seed(11)
sim <- lc_simulate(n = 30)

# debugonce(lc_linear_qi)
result <- lc_linear_qi(sim$data)

# Get best estimates and plot its regression line on top of the plot  
a <- result$intercept["50%"]
b <- result$slope["50%"]
abline(a, b, col = "green4")
# Add confidence interval  
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green4")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green4")

# DIC
result$dic



#
# . Qi version, interval ----
#

set.seed(11)
sim <- lc_simulate(n = 30)

sim$data$y_up <- sim$data$y
sim$data$y_up[sim$data$uncensored == 0] <- sim$data$threshold[sim$data$uncensored == 0]
sim$data$y_lo <- sim$data$y
sim$data$y_lo[sim$data$uncensored == 0] <- 0
sim$data

# debugonce(int_linear_qi)
result <- int_linear_qi(sim$data)

# Get best estimates and plot its regression line on top of the plot  
a <- result$intercept["50%"]
b <- result$slope["50%"]
abline(a, b, col = "green4")
# Add confidence interval  
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green4")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green4")

# DIC
result$dic

#
# LR with measurement error ----
#
# Slope is significantly biased downwards
# Or: y is underestimated for high x
# - actual y-mean is within the CI of estimate y-mean, though  
#

set.seed(11)
sim <- lc_simulate(n = 30)  

result <- lc_linear_measerror(sim$data, measurement_error = 0.1)  
# result <- lc_linear_measerror(sim$data, measurement_error = 0.1, detailed = TRUE)  

# Get best estimates and plot its regression line on top of the plot  
a <- result$intercept["50%"]
b <- result$slope["50%"]
abline(a, b, col = "green2")
# Add confidence interval  
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")

#
# With measurement error and a minimum Y ----
#
# This increased the bias, not reduced it....
#

set.seed(11)
sim <- lc_simulate(n = 30)
result <- lc_linear_measerror_min(sim$data, measurement_error = 0.1, minimum_y = -2)  
# result <- lc_linear_measerror(sim$data, measurement_error = 0.1, detailed = TRUE)  

# Get best estimates and plot its regression line on top of the plot  
a <- result$intercept["50%"]
b <- result$slope["50%"]
abline(a, b, col = "green2")
# Add confidence interval  
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")


#
# . - Bias? Test for many simulations     ----
#

#
# 1. base example using 'lc_linear_measerror'    
#
sim_slope <- function(i){
  sim <- lc_simulate(n = 30, plot = FALSE)
  result <- lc_linear_measerror(sim$data, measurement_error = 0.1)  
  c(i=i, result$slope)
}

if (FALSE){
  # Main test - around 2-3 minutes run  
  set.seed(11)
  sim_result <- 1:20 %>% purrr::map_dfr(sim_slope)
}

ggplot(sim_result, aes(i, `50%`)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  geom_hline(yintercept = -3, color = "red")



#
# 2. give the actual (normally unknown) mean/sd of y  
#
# Didn't help. So hat's not the problem
#
sim_slope_givenmean <- function(i){
  sim <- lc_simulate(n = 30, plot = FALSE)
  result <- lc_linear_measerror(sim$data, measurement_error = 0.1, 
                                mean_y = mean(sim$data$y_real),
                                sd_y = sd(sim$data$y_real))  
  c(i=i, result$slope)
}

if (FALSE){
  # Main test
  set.seed(11)
  sim_result_givenmean <- 1:20 %>% purrr::map_dfr(sim_slope_givenmean)
}

ggplot(sim_result_givenmean, aes(i, `50%`)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  geom_hline(yintercept = -3, color = "red")

#
# 3. try with minimum value
#

# Did not work ell at all, on the contrary  

sim_slope <- function(i){
  sim <- lc_simulate(n = 30, plot = FALSE)
  result <- lc_linear_measerror_min(sim$data, measurement_error = 0.1, minimum_y = -2)  
  c(i=i, result$slope)
}

if (FALSE){
  # Main test
  set.seed(11)
  sim_result_miny <- 1:20 %>% purrr::map_dfr(sim_slope)
}

ggplot(sim_result_miny, aes(i, `50%`)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  geom_hline(yintercept = -3, color = "red")



#
# 4. try using 'lc_linear' (not 'lc_linear_measerror'), but instead add random noise to data  
#    
#

# Add random error (sd = 0.10) for uncensored data
# Mimicking measurement error
create_one_random <- function(i, simulation, sd_fraction = 0.1){
  sel <- simulation$data$uncensored == 1
  sd <- simulation$data$y[sel]*sd_fraction 
  simulation$data$y[sel] <- simulation$data$y[sel] + rnorm(sum(sel), 0, sd)
  # Make sure that no data are under the threshold
  sel2 <- sel & (simulation$data$y < simulation$data$threshold)
  simulation$data$y[sel2] <- simulation$data$threshold[sel2]*1.01
  simulation
}
create_one_random(1, sim1) %>% str(1)

sim_list <- lapply(1:2, create_one_random, simulation = sim1)

# Test that 'lc_linear' works
result <- lc_linear(sim1$data)
result <- lc_linear(sim_list[[1]]$data)

# Check how to get slope "raw" numbers (so they lter can be concayenated for each random data)
str(result$model)
str(result$model[[1]])
dimnames(result$model[[1]])
result$model[[1]][,"slope"] %>% head(50)  # note: normalized numbers

# View(sim1$data)
# View(sim_list[[1]]$data)

# This gets you only the quantiles som not so useful  
get_slope <- function(simulation){
  result <- lc_linear(simulation$data)
  result$slope
}
#  sim_list %>% purrr::map_dfr(get_slope)

# This gets you the raw data (not only the quantiles)    
# NOTE: still the untransformed data  
get_slope_raw <- function(i, simulation){
  result <- lc_linear(simulation$data)
  data.frame(
    i = i,
    slope = c(result$model[[1]][,"slope"],
      result$model[[2]][,"slope"],
      result$model[[3]][,"slope"],
      result$model[[4]][,"slope"]
    )
  )
}
test <- get_slope_raw(42, sim_list[[1]])
str(test)

result_rand <- seq_along(sim_list) %>% purrr::map_dfr(~get_slope_raw(.x, sim_list[[.x]]))
str(result_rand, 1)

#
# GOTTEN THIS FAR
#




# create_one_random(sim1)
sim_slopes_rand <- function(simulation, n){
  sim_result_slopes <- 1:n %>% purrr::map_dfr(create_one_random(simulation))
  sim_result_slopes
}
sim_slopes_rand(sim1, 2)



  
rnorm(10, mean = 1:10, sd = seq(100,0.1,length=10))


sim_slope_rand <- function(i){
  sim1 <- lc_simulate(n = 30, plot = FALSE)
  result1 <- lc_linear(sim$data, measurement_error = 0.1)  
  result2 <- lc_linear_measerror(sim$data, measurement_error = 0.1)  
  c(i=i, result$slope)
}

if (FALSE){
  # Main test
  set.seed(11)
  sim_result_miny <- 1:20 %>% purrr::map_dfr(sim_slope)
}

ggplot(sim_result_miny, aes(i, `50%`)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  geom_hline(yintercept = -3, color = "red")





# Test 
sim_slope(1)






# Truncae distribution  
y_min <- 0
(y_min - result$mean_y)/result$sd_y

# Test plot
plot(result$model_data$x, result$model_data$y.uncens.error, 
     ylim = range(result$model_data$y.uncens.error, result$model_data$threshold, na.rm = TRUE))
sel <- is.na(result$model_data$y.uncens.error)
points(result$model_data$x[sel], result$model_data$threshold[sel], pch = 18, col = "red")
  

set.seed(11)
sim <- lc_simulate(n = )
head(sim$data)

# Normalized y, centralized x - no data under threshold   
set.seed(22)
sim <- lc_simulate(x = seq(-4,4,length=50), intercept = 0, slope = -0.25, sigma = 0.75, 
                   threshold_1 = -3, threshold_2 = -3, threshold_change = 0)
mean(sim$data$y_real, na.rm = TRUE)
sd(sim$data$y_real, na.rm = TRUE)
# - analysis   
result <- lc_linear_measerror(sim$data, measurement_error = 0.1)  
a <- result$intercept["50%"]
b <- result$slope["50%"]
abline(a, b, col = "green2")
# Add confidence interval  
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")

sim_est <- function(..., mean_y = NULL, sd_y = NULL, seed = 22){
  set.seed(seed)
  sim <- lc_simulate(...)  
  mean(sim$data$x, na.rm = TRUE) %>% cat("mean x:", ., "\n")
  mean(sim$data$y_real, na.rm = TRUE) %>% cat("mean y:", ., "\n")
  sd(sim$data$y_real, na.rm = TRUE) %>% cat("sd y:", ., "\n")
  # - analysis   
  result <- lc_linear_measerror(sim$data, measurement_error = 0.1, mean_y = mean_y, sd_y = sd_y)  
  a <- result$intercept["50%"]
  b <- result$slope["50%"]
  abline(a, b, col = "green2")
  # Add confidence interval  
  lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
  lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")
  invisible(result)
}

# 0a. Normalized y, centralized x - no data under threshold   
X <- sim_est(x = seq(-4,4,length=50), intercept = 0, slope = -0.25, sigma = 0.75, 
             threshold_1 = -1000, threshold_2 = -1000, threshold_change = 0)

# 0b. As 0a, lower N     
X <- sim_est(x = seq(-4,4,length=30), intercept = 0, slope = -0.25, sigma = 0.75, 
             threshold_1 = -1000, threshold_2 = -1000, threshold_change = 0)

# run for k samples with n = 30
if (FALSE){
  seeds <- round(runif(n = 2, min = 1, max = 10000))
  X_list <- seeds %>%
    purrr::map(~sim_est(x = seq(-4,4,length=30), intercept = 0, slope = -0.25, sigma = 0.75, 
                 threshold_1 = -1000, threshold_2 = -1000, threshold_change = 0, seed = .)
    )
}

# 1. normalized y, centralized x - some data under threshold   
X <- sim_est(x = seq(-4,4,length=50), intercept = 0, slope = -0.25, sigma = 0.75, 
             threshold_1 = -0.1, threshold_2 = -0.5, threshold_change = 0)

# 2a. as 1 but lower N  
X <- sim_est(x = seq(-4,4,length=30), intercept = 0, slope = -0.25, sigma = 0.75, 
             threshold_1 = -0.1, threshold_2 = -0.5, threshold_change = 0, seed = 40)

# 2b. as 2a (lower N) but mean_y and sd_y are given    
X <- sim_est(x = seq(-4,4,length=30), intercept = 0, slope = -0.25, sigma = 0.75, 
             threshold_1 = -0.1, threshold_2 = -0.5, threshold_change = 0, mean_y = 0, sd_y = 0)

# 3. As 1 but "un-centralize" x (adding 10)       
# - threshold parameters also adjusted
# - works fine  
X <- sim_est(x = seq(-4,4,length=50) + 10, intercept = 0, slope = -0.25, sigma = 0.75, 
             threshold_1 = -2.7, threshold_2 = -3.2, threshold_change = 10)

# 4. As 3 but also "un-centralize" y (adding 10)       
# - threshold parameters also adjusted
# - still works fine  
X <- sim_est(x = seq(-4,4,length=50) + 10, intercept = 10, slope = -0.25, sigma = 0.75, 
             threshold_1 = -2.7 + 10, threshold_2 = -3.2 + 10, threshold_change = 10)

# 5. As 4 but also "un-normalize" y (adjusting slope)       
# - threshold parameters also adjusted
# - still works fine  
X <- sim_est(x = seq(-4,4,length=50) + 10, intercept = 10, slope = -0.55, sigma = 0.75, 
             threshold_1 = 4.7, threshold_2 = 3.2, threshold_change = 10)

# 6. As 5 but lower sample size      
# - threshold parameters also adjusted
# - still works fine  
X <- sim_est(x = seq(-4,4,length=30) + 10, intercept = 10, slope = -0.55, sigma = 0.75, 
             threshold_1 = 4.7, threshold_2 = 3.2, threshold_change = 10)

# 7. As 6 but "un-centralize" y       
X <- sim_est(x = seq(-4,4,length=50), intercept = 0, slope = -0.25, sigma = 0.75, 
             threshold_1 = -0.1, threshold_2 = -0.5, threshold_change = 0)


result <- lc_linear_measerror(sim$data, measurement_error = 0.1, detailed = TRUE)  


rownames(result$summary$quantiles)

#
# With measurement error, less censoring ----
#

set.seed(11)
sim <- lc_simulate(n = 30, threshold_1 = 15, threshold_2 = 5)

# Without measurement error (green)
result <- lc_linear(sim$data)
a <- result$intercept["50%"]
b <- result$slope["50%"]
abline(a, b, col = "green2")
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")

# With measurement error (blue)
result <- lc_linear_measerror(sim$data, measurement_error = 0.1)
a <- result$intercept["50%"]
b <- result$slope["50%"]
abline(a, b, col = "blue3")
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "blue3")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "blue3")


#
# With measurement error, test: without error vs. with different error sizes ----
#

# Simulate data  
set.seed(6)
sim0 <- sim1 <- sim2 <- sim3 <- lc_simulate(n = 30)

# sd_values   
sd_values <- c(0, 0.10, 0.25, 0.50)

# Make list of simulation objects
sim_list <- list(sim0, sim1, sim2, sim3)

# Estimate models
result_list <- list()
result_list[[1]] <- lc_linear(sim_list[[1]]$data)
for (i in 2:4){
  result_list[[i]] <- lc_linear_measerror(sim_list[[i]]$data, 
                                                measurement_error = sd_values[i],
                                                detailed = TRUE)
}

names(result_list) <- sprintf("SD = %.2f", sd_values)

# Extract quantiles of intercept and slope  
intercepts <- purrr::map_dfr(result_list, c("intercept"), .id = "Model")
slopes <- purrr::map_dfr(result_list, c("slope"), .id = "Model")

# Plot slope estimates
ggplot(slopes, aes(Model, `50%`)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`))

# Plot graphs for 
for (i in 1:4){
  gg <- ggplot(sim_list[[i]]$data, aes(x, y_plot, color = factor(uncensored))) +
    geom_point() +
    geom_abline(
      intercept = intercepts$`50%`[i],
      slope = slopes$`50%`[i]
    ) +
    labs(title = names(result_list)[i])
  print(gg)
}


#
# Add example data ----
#

# ?lc_prepare

# Get actual data  
fn <- "../../seksjon 212/Milkys2_pc/Files_from_Jupyterhub_2020/Raw_data/109_adjusted_data_2021-09-15.rds"                         # data FROM Milkys2 on PC 
dat_all <- readRDS(fn)

# Tables for less-than's
dat_all %>%
  filter(MYEAR > 2012) %>%
  xtabs(~PARAM + addNA(FLAG1), .)

# as above, BDEs only 
dat_all %>%
  filter(MYEAR > 2012 & substr(PARAM,1,3) == "BDE") %>%
  xtabs(~PARAM + addNA(FLAG1), .)

# One BDE only, stations 
param <- "BDE99"
tab <- dat_all %>%
  filter(MYEAR > 2012 & PARAM == param) %>%
  xtabs(~STATION_CODE + addNA(FLAG1), .)
# Stations with total sample size >= 80 and at least 20 less-thans: 
tab2 <- tab[apply(tab, 1, sum) >= 50 & tab[,1] >= 15,]
tab2

# One BDE only, data for plot
dat1 <- dat_all %>%
  filter(PARAM == param & TISSUE_NAME == "Lever",
         STATION_CODE %in% rownames(tab2) & !is.na(VALUE_WW)) %>%
  mutate(Over_LOQ = is.na(FLAG1)) %>%
  arrange(desc(Over_LOQ)) 

# One BDE only, plot
ggplot(dat1, aes(MYEAR, VALUE_WW, shape = Over_LOQ, color = Over_LOQ)) +
  geom_jitter(width = 0.1) +
  scale_y_log10() +
  facet_wrap(vars(STATION_CODE))+
  labs(title = param)

if (FALSE){
  
  # Adding example data 'polybrom'  
  
  # Make example dataset
  polybrom <- dat1 %>% 
    filter(STATION_CODE %in% c("13B", "23B", "53B")) %>%
    rename(
      station = STATION_CODE,
      year = MYEAR,
      concentration = VALUE_WW,
      LOQ_flag = FLAG1,
    ) %>%
    arrange(station, year) %>%
    select(station, year, concentration, LOQ_flag)
  
  # Add to package  
  usethis::use_data(polybrom)
  
  # File R/data.R manually edited (as far as I remember?)
  
  # Check
  ggplot(polybrom, aes(year, concentration, 
                       shape = (LOQ_flag %in% "<"), 
                       colour = (LOQ_flag %in% "<"))) +
    geom_jitter(width = 0.1) +
    scale_shape_manual(values = c(16,6)) +
    scale_y_log10() +
    facet_wrap(vars(station))

}



#
# lc_prepare (data from station 23B) ----
#


xtabs(~station + addNA(LOQ_flag), polybrom)

# Prepare data
data_test_orig <- subset(polybrom, station %in% "23B")

# debugonce(lc_prepare)
data_test_prep <- lc_prepare(data_test_orig, 
                        x = "year",
                        y = "concentration", 
                        censored = "LOQ_flag",
                        log = TRUE)

debugonce(lc_plot)
lc_plot(data_test_prep)

# . analysis using lc_linear ----  
result <- lc_linear(data_test_prep, plot_input = TRUE, plot_norm = TRUE)  

# Qi version of the same (2 minutes)
# debugonce(lc_linear_qi)
result <- lc_linear_qi(data_test_prep)  

# Check
str(result, 1)
str(result$model_data, 1)

# . plot data and lc_linear result (base plot) ----

sel_uncens <- !data_test_orig$LOQ_flag %in% "<" 
plot(log(concentration) ~ year, 
     data = data_test_orig[sel_uncens,], 
     ylim = range(log(data_test_orig$concentration), na.rm = TRUE),
     pch = 16, col = 4)
points(log(concentration) ~ year, 
     data = data_test_orig[!sel_uncens,], 
     pch = 6, col = 2)
# Get best estimate fitted line       
a <- result$intercept["50%"]
b <- result$slope["50%"]
# Add regression line to the plot  
abline(a, b, col = "green2")
# Add confidence interval  
lines(y_lo ~ x, data = result$plot_data, lty = "dashed", col = "green2")
lines(y_hi ~ x, data = result$plot_data, lty = "dashed", col = "green2")


# . plot data and lc_linear result (ggplot) ----

ggplot(data_test_orig, aes(x = year)) +
  geom_ribbon(data = result$plot_data, aes(x=x, ymin = y_lo, ymax = y_hi),
              fill = "lightgreen") +
  geom_abline(intercept = a, slope = b, color = "green3") +
  geom_point(aes(y = log(concentration), color = is.na(LOQ_flag)))


#
# Linear MCMC (station 23B), comparing with ordinary LM ---- 
#

# debugonce(lc_linear)
# debugonce(R2jags::jags)
result <- lc_linear(data_test_prep, plot_input = TRUE, plot_norm = TRUE)  

# Plot with regression line  
a <- result$intercept["50%"]
b <- result$slope["50%"]
ggplot(data_test_prep, aes(x, y)) +
  geom_point(aes(y = threshold), shape = 1, size = 3, color = "red") +
  geom_point(aes(color = factor(uncensored))) +
  geom_abline(intercept = result$intercept["50%"], slope = result$slope["50%"])

# Plot with regression line and confidence interval   
a <- result$intercept["50%"]
b <- result$slope["50%"]
ggplot(data_test_prep, aes(x = x)) +
  geom_ribbon(data = result$plot_data, aes(ymin = y_lo, ymax = y_hi), fill = "grey80") + 
  geom_point(aes(y = threshold), shape = 1, size = 3, color = "red") +
  geom_point(aes(y = y, color = factor(uncensored))) +
  geom_abline(intercept = result$intercept["50%"], slope = result$slope["50%"])

# Slope from lc_linear
result$slope

# Slope from linear regression, ignoring "<"
mod <- lm(concentration ~ year, data = subset(data_test_orig, station == "23B"))
ols <- summary(mod)$coef["year",]

# Compare slopes  
df_slopes <- bind_rows(
  data.frame(
    Analysis = "lc_linear",
    slope_min = result$slope["2.5%"],
    slope = result$slope["50%"],
    slope_max = result$slope["97.5%"]),
  data.frame(
    Analysis = "lm",
    slope_min = ols["Estimate"] - 2*ols["Std. Error"],
    slope = ols["Estimate"],
    slope_max = ols["Estimate"] + 2*ols["Std. Error"])
)
ggplot(df_slopes, aes(x = Analysis)) +
  geom_pointrange(aes(y = slope, ymin = slope_min, ymax = slope_max))


#
# Spline MCMC (station 23B) ---- 
#

result_2knots <- lc_fixedsplines(data_test_prep, 
                                 x = "x", y = "y", uncensored = "uncensored", threshold = "threshold",
                                 knots = 2)
result_3knots <- lc_fixedsplines(data_test_prep, 
                                 x = "x", y = "y", uncensored = "uncensored", threshold = "threshold",
                                 knots = 3)
result_4knots <- lc_fixedsplines(data_test_prep, 
                                 x = "x", y = "y", uncensored = "uncensored", threshold = "threshold",
                                 knots = 4)
result_5knots <- lc_fixedsplines(data_test_prep, 
                                 x = "x", y = "y", uncensored = "uncensored", threshold = "threshold",
                                 knots = 5)
result_6knots <- lc_fixedsplines(data_test_prep, 
                                 x = "x", y = "y", uncensored = "uncensored", threshold = "threshold",
                                 knots = 6)

#
# Qi version ----
#
# debugonce(lc_fixedsplines_qi)
result_linear_qi <- lc_linear_qi(data_test_prep, 
                                 x = "x", y = "y", uncensored = "uncensored", threshold = "threshold")
result_2knots_qi <- lc_fixedsplines_qi(data_test_prep, 
                                       x = "x", y = "y", uncensored = "uncensored", threshold = "threshold",
                                       knots = 2)
result_3knots_qi <- lc_fixedsplines_qi(data_test_prep, 
                                    x = "x", y = "y", uncensored = "uncensored", threshold = "threshold",
                                    knots = 3)
result_4knots_qi <- lc_fixedsplines_qi(data_test_prep, 
                                       x = "x", y = "y", uncensored = "uncensored", threshold = "threshold",
                                       knots = 4)
result_linear_qi$dic
result_2knots_qi$dic
result_3knots_qi$dic
result_4knots_qi$dic

plot_prediction <- function(data, jagsresult, title){
  sel_uncens <- is.na(data$LOQ_flag)
  plot(log(concentration) ~ year, 
       data = data[sel_uncens,], 
       ylim = range(log(data$concentration), na.rm = TRUE),
       pch = 16, col = 4, main = title)
  points(log(concentration) ~ year, 
         data = data[!sel_uncens,], 
         pch = 6, col = 2)
  lines(y ~ x, data = jagsresult$plot_data, col = "red")
  lines(y_lo ~ x, data = jagsresult$plot_data, lty = "dashed", col = "red")
  lines(y_hi ~ x, data = jagsresult$plot_data, lty = "dashed", col = "red")  
}

plot_prediction(data_test_orig, result_2knots, "Station 23B, 2 knots")
plot_prediction(data_test_orig, result_3knots, "Station 23B, 3 knots")

plot_prediction(data_test_orig, result_linear_qi, "Station 23B, linear, Qi's method")
plot_prediction(data_test_orig, result_2knots_qi, "Station 23B, 2 knots, Qi's method")
plot_prediction(data_test_orig, result_3knots_qi, "Station 23B, 3 knots, Qi's method")
plot_prediction(data_test_orig, result_4knots_qi, "Station 23B, 4 knots, Qi's method")

plot_prediction(result_4knots, "Station 23B, 4 knots")
plot_prediction(result_5knots, "Station 23B, 5 knots")
plot_prediction(result_6knots, "Station 23B, 6 knots")

#
# . DIC calcultion 1 ----
#
# - DIC values are very similar for 2,3 and 4 knots, and for 5 and 6 knots.
# - Also gives the following warning:
#
# Warning message:
#   In rjags::dic.samples(result_2knots$model_from_jags$model, n.iter = 1000,  :
#                           Failed to set mean monitor for pD
#                         Support of observed nodes is not fixed
#
# The reason for this is that dinterval() function has a limitation in deviance calculation.
# Qi et al. (2022) writes (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8944154/):
#   In the presence of censored outcomes, even though the dinterval() function can 
#   generate the proper posterior distribution of the parameters in JAGS, the likelihood 
#   function is misspecified with the wrong focus of inference on the censored outcome 
#   variable [22]. Instead, a constant value of 1 for the likelihood function, 
#   or equivalently, a constant value of 0 for the deviance function, is misspecified 
#   for the censored outcomes in the deviance monitor.
# 
# Also see Qi's Github:
# https://github.com/xinyue-qi/Censored-Data-in-JAGS  

rjags::dic.samples(result_2knots$model_from_jags$model, n.iter = 1000, type = "pD")
rjags::dic.samples(result_3knots$model_from_jags$model, n.iter = 1000, type = "pD")
rjags::dic.samples(result_4knots$model_from_jags$model, n.iter = 1000, type = "pD")

#
# . DIC calcultion 2  ----
#
# Has the same limitation as 

dic <- c(
  AICcmodavg::DIC(result_2knots$model_from_jags),
  AICcmodavg::DIC(result_3knots$model_from_jags),
  AICcmodavg::DIC(result_4knots$model_from_jags),
  AICcmodavg::DIC(result_5knots$model_from_jags),
  AICcmodavg::DIC(result_6knots$model_from_jags))
plot(2:6, dic)


#
# Test Qi version with simulated data ----   
#
# Make strongly non-linear data
#

X <- seq(from=-1, to=1, by=.025) # generating inputs
B <- t(splines::bs(X, knots=seq(-1,1,1), degree=3, intercept = TRUE)) # creating the B-splines
num_data <- length(X); num_basis <- nrow(B)
a0 <- 0.2 # intercept

set.seed(991)
# num_basis <- 6
a <- rnorm(num_basis, 0, 1) # coefficients of B-splines
n_param <- length(a)

Y_true <- as.vector(a0*X + a%*%B) # generating the output
Y <- Y_true + rnorm(length(X),0,.1) # adding noise

dat_sim <- data.frame(x = X, y_uncensored = Y, y_true = Y_true)
# ggplot(dat_sim, aes(x, y_uncensored)) +
#   geom_point() +
#   geom_line(aes(y = y_true), color = "blue")

# Add censoring 
dat_sim$y <- dat_sim$y_uncensored
dat_sim$uncensored <- 1
threshold_fixed <- -0.3
sel <- dat_sim$y_uncensored < threshold_fixed
dat_sim$y[sel] <- NA  
dat_sim$uncensored[sel] <- 0  
dat_sim$threshold <- threshold_fixed

# Plot
lc_plot(dat_sim)


#
# . fit models and test DIC values ----
#
# Does work as expected
#
# debugonce(lc_fixedsplines_qi)
result_linear_qi <- lc_linear_qi(dat_sim)
result_2knots_qi <- lc_fixedsplines_qi(dat_sim, knots = 2)
result_3knots_qi <- lc_fixedsplines_qi(dat_sim, knots = 3)
result_4knots_qi <- lc_fixedsplines_qi(dat_sim, knots = 4)
result_5knots_qi <- lc_fixedsplines_qi(dat_sim, knots = 5)

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


#
# APPENDIX? - Test truncation ---- 
#
# Based on
# http://www.johnmyleswhite.com/notebook/2010/08/20/using-jags-in-r-with-the-rjags-package/
# Linear regression example

library(rjags)

N <- 100
x <- runif(N, 0, 10)
epsilon <- rnorm(N, , 1)
y <- x + epsilon
# y[y<0] <- 0
  
dat <- data.frame(X = x, Y = y, Epsilon = epsilon)
plot(Y ~ X, dat)

write.table(
  dat,
  file = "C:/data/temp/example2.data",
  row.names = FALSE,
  col.names = TRUE
)

writeLines("model {
    for (i in 1:N) {
        y[i] ~ dnorm(y.hat[i], tau) T(0,)
        y.hat[i] <- a + b * x[i]
    }
    a ~ dnorm(0, 0.0001)
    b ~ dnorm(0, 0.0001)
    tau <- pow(sigma, -2)
    sigma ~ dunif(0, 100)
}", "C:/data/temp/example2.bug")

jags <- jags.model(
  'C:/data/temp/example2.bug',
  data = list(
    'x' = x,
    'y' = y,
    'N' = N
  ),
  n.chains = 4,
  n.adapt = 100
)

update(jags, 1000)

X <- jags.samples(jags,
             c('a', 'b'),
             1000)
hist(X$b[1, ,1])


