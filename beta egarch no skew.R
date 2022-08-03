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
  u <- numeric(lenght = n)
  z <- numeric(lenght =n)
  
  #set parameter values from the input par using link functions for restrictions  
  mu <- par[1] 
  nu <- exp(par[2])
  phi <- (2/pi)*tan(par[3])
  kappa <- exp(par[4])
  omega <- par[5]
  
  
  #do we need to add the restriction for a Beta-t-Egarch (not skewed) for which
  # -1<= u[t] <= nu , nu >0 ??
  u[1] <- 0
  
  lambda <-rep(0,n+1)
  
  lambda[1] <- mu / (1-phi)
  
  lambda_star <- rep(0,n+1)
  u <- rep(0,n)
  z <- rep(0,n)
  
  for(t in 1:n){
    
    omega <- mu/(1-phi)
    
    b[t] <- ((y[t]-mu)^2/(nu*exp(2*lambda[t]))) /
      (1+((y[t]-mu)^2/(nu*exp(2*lambda[t]))))
    
    #score function
    u[t] <- (nu +1)*b -1
  

    
    #lambda[t] <- omega + lambda_plus[t]
    #lamda_plus[t+1] <- phi*lamda_plus[t] + kappa*u[t]
    
    sigma_u^2 <- (2*nu)/(nu+3)
    psi <- rbind(kappa_tilde, phi_tilde, omega_tilde)
    a <- phi - (kappa*2*nu)/(nu+3)
    b <- phi^2 - (phi*kappa*4*nu)/(nu+3) + (kappa^2*12*nu*(nu+1)*(nu+2))/((nu+7)*(nu+5)*(nu+3))
    c <- (kappa*4*nu*(1-nu)) / ((nu+5)*(nu+3))
    
    A <- sigma_u^2
    B <- (kappa^2 * sigma_u^2 * (1 + a*phi))/((1-phi^2)(1-a*phi))
    C <- ((1-phi)^2*(1+a))/(1-a)
    D <- (a*kappa*sigma_u^2) / (1-a*phi)
    E <- C*(1-phi)/(1-a)
    F <- a * C * kappa* (1-phi)/ (1-a)*(1-a*phi)
    
    D(psi) = (1/(1-b)) * matrix(c(A, D, E, D, B, F, E, F, C),
                                nrow = 3,
                                ncol = 3,
                                byrow = TRUE)
    
    I(psi) <-sigma_u^2 * D(psi)
    
    
    
    u[t] <- u[t] * I_positive *(y[t] -mu) + u[t] * I_negative * (y[t]-mu)
    
    # without leverage, dynamic equation:  (mistake in paper; w/(1-phi))
    lambda_star[t+1] <- phi* lambda_star[[t] + kappa* u[t]]
    lambda[t] <- omega + lambda_star[t]
    
    
    ## Calculate Log Likelihood Values
    #construct sequence of log-likelihood contributions    
    logl[t] <- log(mu) + log(e[t])*(lambda[t])
  }
  
  LogLik <- sum(logl)     # obtain the sum log-likelihood
  
  if(is.na(LogLik) | !is.finite(LogLik)) LogLik <- -1e10
  if(LogLik != LogLik) LogLik <- -1e10
  
  return(-LogLik)         # return the average sum log-likelihood as output
}

# Initial parameter values for the Newton Raphson optimization
mu <- 0 # initial value for mu
nu <- 0.05
phi <- 0.02
kappa <- 0.02
omega <- 0.05
# Transform initial values using the inverse of the link functions
par_ini <- c(mu, log(nu), pi/2*atan(phi), log(kappa), omega)

# Optimize log likelihood function
#est <- optim(par=par_ini,fn=function(par)sstd.lik(par,y), method = "L-BFGS-B")
lower <- c(-10, 1e-3, 1e-3, 2.01)
upper <- c(10, 100, 5, 200)

est <- optimx(par = par_ini, fn = sstd.lik, 
              y = y, method = "nlminb",
              control = list(trace = 1, maxit = 1000), hessian = TRUE,
              lower = lower, upper = upper)


# Obtain parameter estimates
mu_hat <- est$par[1]
nu_hat <- exp(est$par[2])
phi_hat <-  tan(est$par[3])*(2/pi) 
kappa_hat <- exp(est$par[4])
psi_hat <- est$par[5]


theta_hat <- c(mu_hat, nu_hat, phi_hat, kappa_hat, psi_hat)

cat("The parameter estimates are:")
round(theta_hat,5)

cat("The log-likelihood value is:")
-est$value*length(y)


# display the exit flag to see if convergence of the algorithm has been attained

cat("Exit flag:")
est$convergence # zero indicates succesfull optimization