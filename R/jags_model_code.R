
get_jags_model_code <- function(bs = "tp",
                                k_code = 5,
                                type = "leftcensored"){
  
  #
  # Censored without measurement error ----
  #
  
  if (bs == "tp" & k_code == 1 & type == "leftcensored"){
    
    # 
    code <- '
model {
  b_new <- c(b[1], b2, 0)                           # Added line
  mu <- X %*% b_new ## expected response
  for (i in 1:n) { y[i] ~ dnorm(mu[i], tau) } ## response 
  for (j in 1:m) {
    Z[j] ~ dbern(prob[j])
    prob[j] <- min(max(pnorm(cut[j], mu[n+j], tau), 0.01),0.99)
  }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  # K1 <- S1[1,1] * lambda[1]  + S1[1,3] * lambda[2] # Changed from the code for k_code = 3
  # b[2] ~ dmnorm(zero[3], K1)                       # Changed from the code for k_code = 3
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'
  
  } else if (bs == "tp" & k_code == 2 & type == "leftcensored"){
    
    # 
    code <- '
model {
  b_new <- c(b[1], b2, b[2])                       # Added line
  mu <- X %*% b_new ## expected response
  # Difference mu - reference mu (observed only if reference_x is set): 
  for (t in 1:(t2-t1+1)) { 
    dmu[t] <- mu[t1+t-1]-mu[t_ref]
  }
  for (i in 1:n) { y[i] ~ dnorm(mu[i], tau) } ## response 
  for (j in 1:m) {
    Z[j] ~ dbern(prob[j])
    prob[j] <- min(max(pnorm(cut[j], mu[n+j], tau), 0.01),0.99)
  }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  K1 <- S1[1,1] * lambda[1]  + S1[1,3] * lambda[2] # Changed from the code for k_code = 3
  b[2] ~ dmnorm(zero[3], K1)                       # Changed from the code for k_code = 3
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'
  
  } else if (bs == "tp" & k_code >= 3 & type == "leftcensored"){
    
    code <- '
model {
  mu <- X %*% b ## expected response
  # Difference mu - reference mu (observed only if reference_x is set): 
  for (t in 1:(t2-t1+1)) { 
    dmu[t] <- mu[t1+t-1]-mu[t_ref]
  }
  for (i in 1:n) { y[i] ~ dnorm(mu[i], tau) } ## response 
  for (j in 1:m) {
    Z[j] ~ dbern(prob[j])
    prob[j] <- min(max(pnorm(cut[j], mu[n+j], tau), 0.01),0.99)
  }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  for (i in c(2:(k-1))) { b[i] ~ dnorm(0, lambda[1]) }
  for (i in c(k)) { b[i] ~ dnorm(0, lambda[2]) }
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

#
# Censored with measurement error ----
#

  } else if (bs == "tp" & k_code == 1 & type == "leftcensored_measerror"){
    
    code <- '
model {
  b_new <- c(b[1], b2, 0)                       # Added line
  mu <- X %*% b_new ## expected response
  for (i in 1:n) { 
    y[i] ~ dnorm(mu[i], total_var[i]^-1)       ## response
    total_var[i] <- scale^2 + meas_error[i]^2
    }  
  for (j in 1:m) {
    Z[j] ~ dbern(prob[j])
    prob[j] <- min(max(pnorm(cut[j], mu[n+j], tau), 0.01),0.99)
  }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior   
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  # K1 <- S1[1,1] * lambda[1]  + S1[1,3] * lambda[2] # Changed from the code for k_code = 3
  # b[2] ~ dmnorm(zero[3], K1)                       # Changed from the code for k_code = 3
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

  } else if (bs == "tp" & k_code == 2 & type == "leftcensored_measerror"){
  
  code <- '
model {
  b_new <- c(b[1], b2, b[2])                       # Added line
  mu <- X %*% b_new ## expected response
  # Difference mu - reference mu (observed only if reference_x is set): 
  for (t in 1:(t2-t1+1)) { 
    dmu[t] <- mu[t1+t-1]-mu[t_ref]
  }
  for (i in 1:n) { 
    y[i] ~ dnorm(mu[i], total_var[i]^-1)       ## response
    total_var[i] <- scale^2 + meas_error[i]^2
    }  
  for (j in 1:m) {
    Z[j] ~ dbern(prob[j])
    prob[j] <- min(max(pnorm(cut[j], mu[n+j], tau), 0.01),0.99)
  }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  K1 <- S1[1,1] * lambda[1]  + S1[1,3] * lambda[2] # Changed from the code for k_code = 3
  b[2] ~ dmnorm(zero[3], K1)                       # Changed from the code for k_code = 3
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

} else if (bs == "tp" & k_code >= 3 & type == "leftcensored_measerror"){
  
  code <- '
model {
  mu <- X %*% b ## expected response
  # Difference mu - reference mu (observed only if reference_x is set): 
  for (t in 1:(t2-t1+1)) { 
    dmu[t] <- mu[t1+t-1]-mu[t_ref]
  }
  for (i in 1:n) { 
    y[i] ~ dnorm(mu[i], total_var[i]^-1)       ## response
    total_var[i] <- scale^2 + meas_error[i]^2
    }  
  for (j in 1:m) {
    Z[j] ~ dbern(prob[j])
    prob[j] <- min(max(pnorm(cut[j], mu[n+j], tau), 0.01),0.99)
  }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  for (i in c(2:(k-1))) { b[i] ~ dnorm(0, lambda[1]) }
  for (i in c(k)) { b[i] ~ dnorm(0, lambda[2]) }
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

#
# Uncensored with measurement error ----
#
# Code: as censored, just deleting the 4 lines:
# for (j in 1:m) {
#   Z[j] ~ dbern(prob[j])
#   prob[j] <- min(max(pnorm(cut[j], mu[n+j], tau), 0.01),0.99)
# }  
#
# and not supplying "cut", "Z" and "m" in the data, see get_jagam_object
#   for the case of m == 0 


} else if (bs == "tp" & k_code == 1 & type == "uncensored"){
  
  # 
  code <- '
model {
  b_new <- c(b[1], b2, 0)                       # Added line
  mu <- X %*% b_new ## expected response
  for (i in 1:n) { y[i] ~ dnorm(mu[i], tau) } ## response 
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  # K1 <- S1[1,1] * lambda[1]  + S1[1,3] * lambda[2] # Changed from the code for k_code = 3
  # b[2] ~ dmnorm(zero[3], K1)                       # Changed from the code for k_code = 3
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'
} else if (bs == "tp" & k_code == 2 & type == "uncensored"){
  
  # 
  code <- '
model {
  b_new <- c(b[1], b2, b[2])                       # Added line
  mu <- X %*% b_new ## expected response
  # Difference mu - reference mu (observed only if reference_x is set): 
  for (t in 1:(t2-t1+1)) { 
    dmu[t] <- mu[t1+t-1]-mu[t_ref]
  }
  for (i in 1:n) { y[i] ~ dnorm(mu[i], tau) } ## response 
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  K1 <- S1[1,1] * lambda[1]  + S1[1,3] * lambda[2] # Changed from the code for k_code = 3
  b[2] ~ dmnorm(zero[3], K1)                       # Changed from the code for k_code = 3
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'
  } else if (bs == "tp" & k_code >= 3 & type == "uncensored"){
    
    code <- '
model {
  mu <- X %*% b ## expected response
  # Difference mu - reference mu (observed only if reference_x is set): 
  for (t in 1:(t2-t1+1)) { 
    dmu[t] <- mu[t1+t-1]-mu[t_ref]
  }
  for (i in 1:n) { y[i] ~ dnorm(mu[i], tau) } ## response 
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  for (i in c(2:(k-1))) { b[i] ~ dnorm(0, lambda[1]) }
  for (i in c(k)) { b[i] ~ dnorm(0, lambda[2]) }
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

#
# Uncensored with measurement error ----
#

  } else if (bs == "tp" & k_code == 1 & type == "uncensored_measerror"){
    
    code <- '
model {
  b_new <- c(b[1], b2, 0)                       # Added line
  mu <- X %*% b_new ## expected response
  for (i in 1:n) { 
    y[i] ~ dnorm(mu[i], total_var[i]^-1)       ## response
    total_var[i] <- scale^2 + meas_error[i]^2
    }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  # K1 <- S1[1,1] * lambda[1]  + S1[1,3] * lambda[2] # Changed from the code for k_code = 3
  # b[2] ~ dmnorm(zero[3], K1)                       # Changed from the code for k_code = 3
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

  } else if (bs == "tp" & k_code == 2 & type == "uncensored_measerror"){
    
    code <- '
model {
  b_new <- c(b[1], b2, b[2])                       # Added line
  mu <- X %*% b_new ## expected response
  # Difference mu - reference mu (observed only if reference_x is set): 
  for (t in 1:(t2-t1+1)) { 
    dmu[t] <- mu[t1+t-1]-mu[t_ref]
  }
  for (i in 1:n) { 
    y[i] ~ dnorm(mu[i], total_var[i]^-1)       ## response
    total_var[i] <- scale^2 + meas_error[i]^2
    }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  K1 <- S1[1,1] * lambda[1]  + S1[1,3] * lambda[2] # Changed from the code for k_code = 3
  b[2] ~ dmnorm(zero[3], K1)                       # Changed from the code for k_code = 3
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

  } else if (bs == "tp" & k_code >= 3 & type == "uncensored_measerror"){
    
    code <- '
model {
  mu <- X %*% b ## expected response
  # Difference mu - reference mu (observed only if reference_x is set): 
  for (t in 1:(t2-t1+1)) { 
    dmu[t] <- mu[t1+t-1]-mu[t_ref]
  }
  for (i in 1:n) { 
    y[i] ~ dnorm(mu[i], total_var[i]^-1)       ## response
    total_var[i] <- scale^2 + meas_error[i]^2
    }  
  scale <- 1/tau ## convert tau to standard GLM scale
  tau ~ dgamma(.05,.005) ## precision parameter prior 
  ## Parametric effect priors CHECK tau=1/79^2 is appropriate!
  for (i in 1:1) { b[i] ~ dnorm(0,0.00016) }
  ## prior for s(x2)... 
  for (i in c(2:(k-1))) { b[i] ~ dnorm(0, lambda[1]) }
  for (i in c(k)) { b[i] ~ dnorm(0, lambda[2]) }
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
}
'

  }  else {
    
    stop("The given combination of bs = ", sQuote(bs), ", k_code = ", k_code, ", and type = ", sQuote(type),
         " has not been implemented.")
    
  }

code

}


