############ 00. Clean workspace

rm(list=ls())    

############ 0. Load packages and/or functions
library(quantmod)
getwd()

############ 1. Load data

sep500_data <- read.delim("SeP500.txt")
y_price <- as.numeric(sep500_data[,])
y <- diff(log(y_price))

sstd.lik <- function(par, y){
  
  n <- length(y)
  logl <- numeric(length = n)
  lambda <- numeric(length = n)
  #lambda_star <- numeric(length = n+1)
  u <- numeric(length = n)
  eps <- numeric(length =n)
  #u_score <- numeric(length =n)
  sigma <- numeric(length = n)
  
  #set of parameters
  nu <- par[1]
  phi <- par[2]
  kappa <- par[3]
  omega <- par[4]
  gamma <- par[5]
  
  lambda <-rep(0,n+1)
  lambda[1] <- omega / (1-phi)
  
  #lambda_star <- rep(0,n+1)
 
  u <- rep(0,n)
  #u_score <- rep(0,n)
  sigma <- rep(0,n)
  
  for(t in 1:n){
    
    eps[t] <- (y[t])/exp(lambda[t])
    
    sigma[t] <- (gamma + 1/gamma)*exp(lambda[t])*sqrt(nu/4*(nu-2))
    lambda = log(sigma)
    
    logl[t] <- log(2) - log(gamma + (1/gamma)) + lgamma((nu+1)/2) - lgamma(nu/2) - (1/2)*log(pi*nu) -
      lambda[t] - (nu+1)/2 * log(1+ ((y[t])^2 / (xi^(2*sign(y[t]))*nu*exp(2*lambda[t]))))
    
    #omega <- mu/(1-phi)
    #score function (equation 2)
    #u[t] <- ((nu + 1)*(y[t]-mu)^2)/(nu*exp(2*lambda[t])+(y[t]-mu)^2) -1
    

    b_plus <- ((y[t])^2/(nu*gamma^2*exp(2*lambda[t]))) /
      (1+((y[t])^2/(nu*gamma^2 * exp(2*lambda[t]))))
    
    b_minus <- ((y[t])^2/(nu*gamma^-2*exp(2*lambda[t]))) /
      (1+((y[t])^2/(nu*gamma^-2 * exp(2*lambda[t]))))
    
    #conditional score
    u_score <- (nu+1)*b_plus + (nu+1)*b_minus  - 1
    #u_score[t] <- (nu+1)*b - 1
    #u_plus <- (nu+1)*b_plus - 1
    #u_min <- (nu+1)*b_minus - 1
    
    #indicator variable 
    if((y[t]) >= 0){
      #b = b_plus
      b_minus <- 0
      #u_plus <- u_score
      if(eps[t] >= 0)
      I_plus <- 1
      I_min <- 0
      #u_plus <- u_score
      #u_min <- 0
      if(eps[t] < 0)
      I_plus <- 0
      I_min <- 1
    }
    else if(y[t]< 0){
      #b = b_minus
      b_plus <- 0
      #u_min <- u_score
      if(eps[t] >= 0)
      I_min <- 1
      I_plus <- 0
      #u_plus <- 0
      if(eps[t] < 0)
      I_min <- 0
      I_plus <- 1
    }
  
    u[t] <- u_score * I_plus * y[t]  +  u_score * I_min * y[t] 
    
    # without leverage, dynamic equation: 
    lambda[t+1] <- omega*(1-phi) + phi*lambda[t] + kappa* u[t]
    #lambda_star[t+1] <- phi* lambda_star[t] + kappa* u[t]
    #lambda[t+1] <- omega + lambda_star[t+1]
    }

  ## Calculate Log Likelihood Values
  #construct sequence of log-likelihood contributions 

  LogLik <- sum(logl)     # obtain the sum log-likelihood
  
  if(is.na(LogLik) | !is.finite(LogLik)) LogLik <- -1e10
  if(LogLik != LogLik) LogLik <- -1e10
  
  return(-LogLik)   
  
}


# Initial parameter values for the Newton Raphson optimization
nu <- 10
phi <- 0.9
kappa <- 0.04
omega <- 0.4
gamma <- 0.03

# Transform initial values using the inverse of the link functions
par_ini <- c(nu, phi, kappa, omega, gamma)

# Optimize log likelihood function

lower <- c(2.1, -1e10, -0.999999, -1e10, -1e10)
upper <- c(100, 1e10, 0.999999, 1e10, 1e10 )

est <- optimx(par = par_ini, fn = sstd.lik, 
              y = y, method = "nlminb",
              control = list(trace = 1, maxit = 1000), hessian = TRUE,
              lower = lower, upper = upper)

#est <- optim(par=par_ini,fn=function(par) sstd.lik(par,y), method = "BFGS")

#est <- optimx(par = par_ini, fn = sstd.lik, 
             #y = y, method = "nlminb",
             #control = list(trace = 1, maxit = 1000), hessian = TRUE,
             #lower = lower, upper = upper)


# Obtain parameter estimates
nu_hat <- est$par[1]
phi_hat <-  est$par[2]
kappa_hat <- est$par[3]
omega_hat <- est$par[4]
gamma_hat <- est$par[5]

theta_hat <- c(nu_hat, phi_hat, kappa_hat, omega_hat, gamma_hat)

cat("The parameter estimates are:")
round(theta_hat,5)

cat("The log-likelihood value is:")
-est$value*length(y)


# display the exit flag to see if convergence of the algorithm has been attained

cat("Exit flag:")
est$convergence # zero indicates succesfull optimization