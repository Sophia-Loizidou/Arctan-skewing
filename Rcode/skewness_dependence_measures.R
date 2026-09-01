#' Arctan Skewing: Skewness and Dependence Measures
#'
#' This file contains functions to calculate skewness and dependence measures for
#' bivariate distributions with arctan skewing transformation. 
#'
#' @author Sophia Loizidou
#' @license CC-BY-4.0 (Attribution 4.0 International)

## Code for densities: dcosine, dsine, dbwc
source('https://raw.githubusercontent.com/Sophia-Loizidou/Symmetry_test_on_hypertorus/main/Rcode/Generating_data.R')


trig_moments_beta <- function(p1, p2, dens, lambda1, lambda2, params, c = 1){
  require(pracma)
  
  int1 <- function(theta1, theta2){
    sin(p1*(theta1 - params$mu[1]) + p2*(theta2 - params$mu[2])) * atan(lambda1*sin(theta1 - params$mu[1]) / (1 + sin(theta1 - params$mu[1])^2)) * 
      eval(parse(text=paste("d",dens,"(theta1,theta2,params$mu[1],params$mu[2],params$kappa[1],params$kappa[2],params$rho)",sep="")))
  }
  
  int2 <- function(theta1, theta2){
    sin(p1*(theta1 - params$mu[1]) + p2*(theta2 - params$mu[2])) * atan(lambda2*sin(theta2 - params$mu[2]) / (1 + sin(theta2 - params$mu[2])^2)) * 
      eval(parse(text=paste("d",dens,"(theta1,theta2,params$mu[1],params$mu[2],params$kappa[1],params$kappa[2],params$rho)",sep="")))
  }
  
  beta_p <- 1/(pi*c) * (integral2(int1, xmin = -pi, xmax = pi, ymin = -pi, ymax = pi)[[1]] + 
                          integral2(int2, xmin = -pi, xmax = pi, ymin = -pi, ymax = pi)[[1]])
  
  return(beta_p)
}

trig_moments_alpha <- function(p1, p2, dens, params){
  require(pracma)
  
  int <- function(theta1, theta2){
    cos(p1*(theta1 - params$mu[1]) + p2*(theta2 - params$mu[2])) *
      eval(parse(text=paste("d",dens,"(theta1,theta2,params$mu[1],params$mu[2],params$kappa[1],params$kappa[2],params$rho)",sep="")))
  }
  
  alpha_p <- integral2(int, xmin = -pi, xmax = pi, ymin = -pi, ymax = pi)[[1]]
  
  return(alpha_p)
}

skewness_dependence_arctan_skewing <- function(dens, lambda, params){
  
  alpha1_0 <- trig_moments_alpha(1,0,dens,params)
  alpha2_0 <- trig_moments_alpha(2,0,dens,params)
  alpha0_1 <- trig_moments_alpha(0,1,dens,params)
  alpha0_2 <- trig_moments_alpha(0,2,dens,params)
  alpha1_1 <- trig_moments_alpha(1,1,dens,params)
  alpha1_minus1 <- trig_moments_alpha(1,-1,dens,params)
  alphaminus1_1 <- trig_moments_alpha(1,-1,dens,params)
  
  beta1_0 <- trig_moments_beta(1,0,dens,lambda[1],lambda[2],params)
  beta2_0 <- trig_moments_beta(2,0,dens,lambda[1],lambda[2],params)
  beta0_1 <- trig_moments_beta(0,1,dens,lambda[1],lambda[2],params)
  beta0_2 <- trig_moments_beta(0,2,dens,lambda[1],lambda[2],params)
  beta1_1 <- trig_moments_beta(1,1,dens,lambda[1],lambda[2],params)
  beta1_minus1 <- trig_moments_beta(1,-1,dens,lambda[1],lambda[2],params)
  betaminus1_1 <- trig_moments_beta(-1,1,dens,lambda[1],lambda[2],params)
  
  
  ## skewness
  rho_sq1 <- alpha1_0^2 + beta1_0^2
  rho_sq2 <- alpha0_1^2 + beta0_1^2
  
  beta_centered1 <- ((alpha1_0^2 - beta1_0^2) * beta2_0 - 2 * alpha1_0 * beta1_0 * alpha2_0) / rho_sq1
  beta_centered2 <- ((alpha0_1^2 - beta0_1^2) * beta0_2 - 2 * alpha0_1 * beta0_1 * alpha0_2) / rho_sq2
  
  sk1 <- beta_centered1 / (1- sqrt(rho_sq1))^{3/2}
  sk2 <- beta_centered2 / (1- sqrt(rho_sq2))^{3/2}
  
  rho_fl_num <- alpha1_minus1^2 - alpha1_1^2 - (beta1_1 + beta1_minus1)^2
  rho_fl_den <- sqrt((1 - alpha2_0^2 - beta2_0^2) * (1 - alpha0_2^2 - beta0_2^2))
    
  return(list(skewness1 = sk1, skewness2 = sk2, rho_fl = rho_fl_num/rho_fl_den))
}

## Example
# skewness_dependence_arctan_skewing(model_best, c(-1,2), params=list(mu = c(1, 2), kappa = c(55,10), rho = -0.3))


