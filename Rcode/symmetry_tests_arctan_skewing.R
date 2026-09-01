#' Arctan Skewing: Symmetry Tests
#'
#' This file implements symmetry tests for circular and toroidal data with arctan
#' skewing transformation. Includes optimal tests for specified and unspecified medians.
#'
#' @author Sophia Loizidou
#' @license CC-BY-4.0 (Attribution 4.0 International)

source("~/Library/CloudStorage/OneDrive-UniversityofLuxembourg/PhD/Projects/arctan_skewing/R_code/necessary_functions/basic_functions.R")


####################### Optimal tests #######################

# Test for specified median --------------------------------------------------------
specified_median_test_arctan <- function(theta, mu){
  if(length(mu) != length(theta[1,])){stop('mu must be a vector of same length as the number of random variables')}
  n <- length(theta[1,])
  d <- length(mu)
  
  gamma <- matrix(rep(NA), nrow = d, ncol = d)
  for(i in 1:d){
    for(j in i:d){
      gamma[i,j] <- (2/(pi * d))^2 * (1/n)*sum(sin(theta[,i]-mu[i])*sin(theta[,j]-mu[j]) / ((1 + sin(theta[,i]-mu[i])^2) * (1 + sin(theta[,j]-mu[j])^2)))
      gamma[j,i] <- gamma[i,j]
    }
  }
  
  gamma_inv <- solve(gamma)
  delta <- 2/(pi * d) * matrix(apply(sin(t(t(theta) - mu)) / (1 + sin(t(t(theta) - mu))^2), 2, sum), nrow = d)
  
  statistic <- (1/n)*t(delta)%*%gamma_inv%*%delta
  
  # the test has asymptotic chi-square distribution
  p_value <- as.numeric(1-pchisq(statistic, d))
  names(statistic) = "Specified median optimal test"
  lamb0 = 0
  names(lamb0) = "value of the vector lambda"
  totval = list(statistic = statistic, p.value = p_value, null.value = lamb0, alternative = "two.sided",
                method = "Symmetry Optimal Test for Specified Median", sample.size = n, 
                data.name=deparse(substitute(x)))
  class(totval) <- "htest"
  
  return(totval)
}

# Test for unspecified median --------------------------------------------------------
## Quantities required for f0 = sine model
phi1_sine_model <- function(x, y, mu1, mu2, kappa1, kappa2, rho){
  kappa1*sin(x-mu1) - rho*cos(x-mu1)*sin(y-mu2)
}

phi2_sine_model <- function(x, y, mu1, mu2, kappa1, kappa2, rho){
  kappa2*sin(y-mu2) - rho*sin(x-mu1)*cos(y-mu2)
}

phi1_dtheta1_sine_model <- function(x, y, mu1, mu2, kappa1, kappa2, rho){
  kappa1*cos(x-mu1) + rho*sin(x-mu1)*sin(y-mu2)
}

phi1_dtheta2_sine_model <- function(x, y, mu1, mu2, kappa1, kappa2, rho){
  - rho*cos(x-mu1)*cos(y-mu2)
}

phi2_dtheta2_sine_model <- function(x, y, mu1, mu2, kappa1, kappa2, rho){
  kappa2*cos(y-mu2) + rho*sin(x-mu1)*sin(y-mu2)
}

## Quantities required for f0 = cosine model
phi1_cosine_model <- function(x, y, mu1, mu2, kappa1, kappa2, rho){
  kappa1*sin(x-mu1) - rho*sin(x-mu1-y+mu2)
}

phi2_cosine_model <- function(x, y, mu1, mu2, kappa1, kappa2, rho){
  kappa2*sin(y-mu2) + rho*sin(x-mu1-y+mu2)
}

phi1_dtheta1_cosine_model <- function(x, y, mu1, mu2, kappa1, kappa2, rho){
  kappa1*cos(x-mu1) - rho*cos(x-mu1-y+mu2)
}

phi1_dtheta2_cosine_model <- function(x, y, mu1, mu2, kappa1, kappa2, rho){
  rho*cos(x-mu1-y+mu2)
}

phi2_dtheta2_cosine_model <- function(x, y, mu1, mu2, kappa1, kappa2, rho){
  kappa2*cos(y-mu2) - rho*cos(x-mu1-y+mu2)
}

## Quantities required for f0 = bivariate wrapped cauchy model
phi1_bwc_model <- function(x, y, mu1, mu2, kappa1, kappa2, rho){
  c0= (1+rho^2)*(1+kappa1^2)*(1+kappa2^2)-8*abs(rho)*kappa1*kappa2
  c1= 2*(1+rho^2)*(1+kappa2^2)*kappa1-4*abs(rho)*(1+kappa1^2)*kappa2
  c2= 2*(1+rho^2)*(1+kappa1^2)*kappa2-4*abs(rho)*(1+kappa2^2)*kappa1
  c3= -4*(1+rho^2)*kappa1*kappa2+2*abs(rho)*(1+kappa1^2)*(1+kappa2^2)
  c4= 2*rho*(1-kappa1^2)*(1-kappa2^2)
  denom <- c0-c1*cos(x-mu1)-c2*cos(y-mu2)-c3*cos(x-mu1)*cos(y-mu2)-c4*sin(x-mu1)*sin(y-mu2)
  num <- -c1*sin(x-mu1)-c3*sin(x-mu1)*cos(y-mu2)+c4*cos(x-mu1)*sin(y-mu2)
  return( - num / denom)
}

phi2_bwc_model <- function(x, y, mu1, mu2, kappa1, kappa2, rho){
  c0= (1+rho^2)*(1+kappa1^2)*(1+kappa2^2)-8*abs(rho)*kappa1*kappa2
  c1= 2*(1+rho^2)*(1+kappa2^2)*kappa1-4*abs(rho)*(1+kappa1^2)*kappa2
  c2= 2*(1+rho^2)*(1+kappa1^2)*kappa2-4*abs(rho)*(1+kappa2^2)*kappa1
  c3= -4*(1+rho^2)*kappa1*kappa2+2*abs(rho)*(1+kappa1^2)*(1+kappa2^2)
  c4= 2*rho*(1-kappa1^2)*(1-kappa2^2)
  denom <- c0-c1*cos(x-mu1)-c2*cos(y-mu2)-c3*cos(x-mu1)*cos(y-mu2)-c4*sin(x-mu1)*sin(y-mu2)
  num <- -c2*sin(y-mu2)-c3*cos(x-mu1)*sin(y-mu2)+c4*sin(x-mu1)*cos(y-mu2)
  return( - num / denom)
}

phi1_dtheta1_bwc_model <- function(x, y, mu1, mu2, kappa1, kappa2, rho){
  c0= (1+rho^2)*(1+kappa1^2)*(1+kappa2^2)-8*abs(rho)*kappa1*kappa2
  c1= 2*(1+rho^2)*(1+kappa2^2)*kappa1-4*abs(rho)*(1+kappa1^2)*kappa2
  c2= 2*(1+rho^2)*(1+kappa1^2)*kappa2-4*abs(rho)*(1+kappa2^2)*kappa1
  c3= -4*(1+rho^2)*kappa1*kappa2+2*abs(rho)*(1+kappa1^2)*(1+kappa2^2)
  c4= 2*rho*(1-kappa1^2)*(1-kappa2^2)
  M <- c0-c1*cos(x-mu1)-c2*cos(y-mu2)-c3*cos(x-mu1)*cos(y-mu2)-c4*sin(x-mu1)*sin(y-mu2)
  num <- (c1*cos(x-mu1)+c3*cos(x-mu1)*cos(y-mu2)+c4*sin(x-mu1)*sin(y-mu2))*M - (c1*sin(x-mu1)+c3*sin(x-mu1)*cos(y-mu2)-c4*cos(x-mu1)*sin(y-mu2))^2
  return(num / M^2)
}

phi1_dtheta2_bwc_model <- function(x, y, mu1, mu2, kappa1, kappa2, rho){
  c0= (1+rho^2)*(1+kappa1^2)*(1+kappa2^2)-8*abs(rho)*kappa1*kappa2
  c1= 2*(1+rho^2)*(1+kappa2^2)*kappa1-4*abs(rho)*(1+kappa1^2)*kappa2
  c2= 2*(1+rho^2)*(1+kappa1^2)*kappa2-4*abs(rho)*(1+kappa2^2)*kappa1
  c3= -4*(1+rho^2)*kappa1*kappa2+2*abs(rho)*(1+kappa1^2)*(1+kappa2^2)
  c4= 2*rho*(1-kappa1^2)*(1-kappa2^2)
  M <- c0-c1*cos(x-mu1)-c2*cos(y-mu2)-c3*cos(x-mu1)*cos(y-mu2)-c4*sin(x-mu1)*sin(y-mu2)
  num <- (-c3*sin(x-mu1)*sin(y-mu2)-c4*cos(x-mu1)*cos(y-mu2))*M - (c1*sin(x-mu1)+c3*sin(x-mu1)*cos(y-mu2)-c4*cos(x-mu1)*sin(y-mu2))*(c2*sin(y-mu2)+c3*cos(x-mu1)*sin(y-mu2)-c4*sin(x-mu1)*cos(y-mu2))
  return(num / M^2)
}

phi2_dtheta2_bwc_model <- function(x, y, mu1, mu2, kappa1, kappa2, rho){
  c0= (1+rho^2)*(1+kappa1^2)*(1+kappa2^2)-8*abs(rho)*kappa1*kappa2
  c1= 2*(1+rho^2)*(1+kappa2^2)*kappa1-4*abs(rho)*(1+kappa1^2)*kappa2
  c2= 2*(1+rho^2)*(1+kappa1^2)*kappa2-4*abs(rho)*(1+kappa2^2)*kappa1
  c3= -4*(1+rho^2)*kappa1*kappa2+2*abs(rho)*(1+kappa1^2)*(1+kappa2^2)
  c4= 2*rho*(1-kappa1^2)*(1-kappa2^2)
  M <- c0-c1*cos(x-mu1)-c2*cos(y-mu2)-c3*cos(x-mu1)*cos(y-mu2)-c4*sin(x-mu1)*sin(y-mu2)
  num <- (c2*cos(y-mu2)+c3*cos(x-mu1)*cos(y-mu2)+c4*sin(x-mu1)*sin(y-mu2))*M-(c2*sin(y-mu2)+c3*cos(x-mu1)*sin(y-mu2)-c4*sin(x-mu1)*cos(y-mu2))^2
  return(num / M^2)
}

## Quantities required for wrapped Cauchy distribution
cdf_wc <- function(theta, rho, mu){
  if((rho >= 1) | (rho < 0)){stop('rho should be in the interval [0,1)')}
  
  theta <- theta %% (2*pi); mu <- mu %% (2*pi)
  
  pdf <- function(x){dwrpcauchy(x, mu, rho)}; cdf <- numeric(length(theta))
  
  for(i in 1:length(theta)){cdf[i] <- integrate(pdf, 0, theta[i])$value}
  return(cdf)
}

wc_dtheta <- function(theta, rho, mu){
  -(1 - rho^2) * 2 * rho * sin(theta-mu) / (2 * pi * (1 + rho^2 - 2*rho*cos(theta-mu))^2)
}

wc_dtheta_dtheta <- function(theta, rho, mu){
  -(1 - rho^2) * (2*rho*cos(theta-mu)*(1 + rho^2 - 2*rho*cos(theta-mu)) - 8*rho^2*sin(theta-mu)^2) /(2*pi*(1 + rho^2 - 2*rho*cos(theta-mu))^3)
}

phi_wc_model <- function(theta, rho, mu){
  rho*sin(theta - mu) / (1 + rho^2 - 2*rho*cos(theta - mu))
}

phi_dtheta_wc_model <- function(theta, rho, mu){
  (rho*cos(theta - mu)*(1 + rho^2 - 2*rho*cos(theta - mu)) - rho^2*cos(theta - mu)^2) / (1 + rho^2 - 2*rho*cos(theta - mu))^2
}

## Quantities required for cardioid distribution
phi_card_model <- function(theta, rho, mu){
  -sin(theta-mu)/(1-cos(theta-mu))
}

phi_dtheta_card_model <- function(theta, rho, mu){
  1/(1-cos(theta-mu))
}

## Quantities required for f0 = trivariate wrapped Cauchy distribution
phi_twd_model_wc_marginals <- function(M, M_dtheta, f, f_dtheta){
  -(f_dtheta*M - f*M_dtheta)/(f*M)
}

phi1_dtheta1_model_wc_marginals <- function(M, M_dtheta, M_dtheta2, f, f_dtheta, f_dtheta2){
  -((f_dtheta2*M - f*M_dtheta2)*f*M - (f_dtheta*M - f*M_dtheta)*(f_dtheta*M + f*M_dtheta)) / (f*M)^2
}

phi1_dtheta2_model_wc_marginals <- function(M, M_dtheta1, M_dtheta2, M_dtheta3, f, f_dtheta){
  -((f_dtheta*M_dtheta2 - f*M_dtheta3)*M - (f_dtheta*M - f*M_dtheta1)*M_dtheta2) / (f*(M^2))
}

## Test for unspecified median
unspecified_median_test_arctan <- function(x, params, 
              model = c("vonmises", "wrappedcauchy", "cardioid", "wrappednormal", "katojones", 
                        "sine", "cosine", "bwc", 
                        "twcc")){
  require(circular)
  require(pracma)
  
  model <- match.arg(model)
  
  n <- dim(x)[1]; d <- dim(x)[2]
  
  if(d > 3){stop("The test has been developed for d = 1, 2 or 3.")}
  if(!is.matrix(x)) {stop("Argument 'x' must be a numeric matrix.")}
  if(model %in% c("vonmises", "wrappedcauchy", "cardioid", "wrappednormal", "katojones")){
    if(d != 1){stop(paste0("For model = ", model, " argument 'x' must be a matrix with one column."))}
    if(model %in% c("vonmises", "wrappednormal", "katojones")){stop(paste0("Optimal tests for unspecified median have not been employed for model = ", model))}
  } else if (model  %in% c("sine", "cosine", "bwc")){
    if(d != 2){stop(paste0("For model = ", model, " argument 'x' must be a matrix with two columns."))}
  } else {
    if(d != 3){stop(paste0("For model = ", model, " argument 'x' must be a matrix with three columns."))}
  }
  
  #estimate mu using circular median
  mu_hat <- rep(NA, d)
  for(i in 1:d){
    mu_hat[i] <- as.numeric(median.circular(circular(x[, i])))
  }
  
  if(d == 1) {
    if(model == 'wrappedcauchy'){
      
      phi1 = phi_wc_model(x[, i], params$rho, mu_hat)
      phi1_dtheta1 = phi_dtheta_wc_model(x[, i], params$rho, mu_hat)
      
    } else if (model == 'cardioid') {
      x = x %% (2*pi)
      
      ## remove points that are 0 or 2pi
      prob = c()
      for(i in 1:n){
        if(any(abs(x[i,] - c(mu_hat[1], mu_hat[2])) < 1e-5) | any(abs(x[i,] - c(mu_hat[1], mu_hat[2]) - 2*pi) < 1e-5)){
          prob = c(prob, i)
        }
      }
      if(!is.null(prob)){
        x = x[-prob,]; n = dim(x)[1]
      }
      
      phi1 = phi_card_model(x[, i], params$rho, mu_hat)
      phi1_dtheta1 = phi_dtheta_card_model(x[, i], params$rho, mu_hat)
    }
    
  } else if (d == 2){
    phi1 <- eval(parse(text=paste("phi1_",model,"_model(x[,1],x[,2],mu_hat[1],mu_hat[2],params$kappa[1],params$kappa[2],params$rho)",sep="")))
    phi2 <- eval(parse(text=paste("phi2_",model,"_model(x[,1],x[,2],mu_hat[1],mu_hat[2],params$kappa[1],params$kappa[2],params$rho)",sep="")))
    phi1_dtheta1 <- eval(parse(text=paste("phi1_dtheta1_",model,"_model(x[,1],x[,2],mu_hat[1],mu_hat[2],params$kappa[1],params$kappa[2],params$rho)",sep="")))
    phi1_dtheta2 <- eval(parse(text=paste("phi1_dtheta2_",model,"_model(x[,1],x[,2],mu_hat[1],mu_hat[2],params$kappa[1],params$kappa[2],params$rho)",sep="")))
    phi2_dtheta2 <- eval(parse(text=paste("phi2_dtheta2_",model,"_model(x[,1],x[,2],mu_hat[1],mu_hat[2],params$kappa[1],params$kappa[2],params$rho)",sep="")))
    
  } else {
    require(CircStats)
    for(i in 1:d){
      eval(parse(text=paste("u",i,"= 2*pi*cdf_wc(x[,",i,"], params$params", i, "$rho, mu_hat[i])",sep="")))
      eval(parse(text=paste("f",i,"= dwrpcauchy(x[,",i,"], params$params", i, "$rho, mu_hat[i])",sep="")))
      eval(parse(text=paste("f",i,"_dtheta", i, "= wc_dtheta(x[,",i,"], params$params", i, "$rho, mu_hat[i])",sep="")))
      eval(parse(text=paste("f",i,"_dtheta", i, i, "= wc_dtheta_dtheta(x[,",i,"], params$params", i, "$rho, mu_hat[i])",sep="")))
    }
    
    rho12 = params$rho12; rho23 = params$rho23; rho13 = params$rho13
    c1 = rho12*rho13/rho23 + rho12*rho23/rho13 + rho13*rho23/rho12
    M = c1 + 2*rho12*cos(u1 - u2) + 2*rho13*cos(u1 - u3) + 2*rho23*cos(u2 - u3)
    M_dtheta1 = 4*pi*f1*(-rho12*sin(u1 - u2) - rho13*sin(u1 - u3))
    M_dtheta2 = 4*pi*f2*(rho12*sin(u1 - u2) - rho23*sin(u2 - u3))
    M_dtheta3 = 4*pi*f3*(rho13*sin(u1 - u3) + rho23*sin(u2 - u3))
    M_dtheta11 = 4*pi*f1_dtheta1*(-rho12*sin(u1 - u2) - rho13*sin(u1 - u3)) + 8*pi^2*f1^2*(-rho12*cos(u1 - u2) - rho13*cos(u1 - u3))
    M_dtheta22 = 4*pi*f2_dtheta2*(rho12*sin(u1 - u2) - rho23*sin(u2 - u3)) + 8*pi^2*f2^2*(-rho12*cos(u1 - u2) - rho23*cos(u2 - u3))
    M_dtheta33 = 4*pi*f3_dtheta3*(rho13*sin(u1 - u3) + rho23*sin(u2 - u3)) + 8*pi^2*f3^2*(-rho13*cos(u1 - u3) - rho23*cos(u2 - u3))
    M_dtheta12 = 8*pi^2*f1*f2*rho12*cos(u1 - u2)
    M_dtheta13 = 8*pi^2*f1*f3*rho13*cos(u1 - u3)
    M_dtheta23 = 8*pi^2*f2*f3*rho23*cos(u2 - u3)
    
    for(i in 1:d){
      eval(parse(text=paste("phi",i,"<- phi_twd_model_wc_marginals(M, M_dtheta", i, ", f", i, ", f", i, "_dtheta", i, ")",sep="")))
      
      for(j in i:d){
        if(i == j){
          eval(parse(text=paste("phi",i,"_dtheta",i,"<-phi1_dtheta1_model_wc_marginals(M, M_dtheta", i, ", M_dtheta", i,i,", f", i, ", f", i, "_dtheta", i, ", f", i, "_dtheta", i,i,")",sep="")))
        } else {
          eval(parse(text=paste("phi",i,"_dtheta",j,"<-phi1_dtheta2_model_wc_marginals(M, M_dtheta", i, ", M_dtheta", j, ", M_dtheta", i,j,", f", i, ", f", i, "_dtheta", i, ")",sep="")))
        }
      }
    }
  }
  
  delta_lambda <- rep(NA, d)
  delta_mu <- rep(NA, d)
  c_mu_lambda <- matrix(0, nrow = d, ncol = d)
  c_mu_mu <- matrix(0, nrow = d, ncol = d)
  for(i in 1:d){
    delta_lambda[i] <- eval(parse(text=paste("2 /(pi * d) * sum(sin(x[,",i,"] - mu_hat[i])/ (1 + sin(x[,",i,"] - mu_hat[i])^2)) / sqrt(n)",sep="")))
    delta_mu[i] <- eval(parse(text=paste("sum(phi",i,") / sqrt(n)",sep="")))
    
    c_mu_lambda[i,i] = eval(parse(text=paste("sum(log(sec((x[,", i, "] - mu_hat[i])/2)^2 - 2*sqrt(2) + 2) - log((sec(x[,", i, "] - mu_hat[i])/2)^2 + 2*sqrt(2) + 2)) / (n * pi * d * sqrt(2))", sep="")))
    
    for(j in i:d){
      c_mu_mu[i,j] = eval(parse(text=paste("sum(phi",i,"_dtheta",j,") / n",sep="")))
      c_mu_mu[j,i] = c_mu_mu[i,j]
    }
  }
  
  c_mu_mu_inv = solve(c_mu_mu)
  
  central_seq <- delta_lambda - c_mu_lambda %*% c_mu_mu_inv %*% delta_mu
  var <- matrix(0, nrow = d, ncol = d)
  for(i in 1:n){
    phi <- c()
    for(j in 1:d){
      eval(parse(text=paste("phi <- c(phi, phi",j,"[", i, "])",sep="")))
    }
    delta <- as.numeric(2 /(pi * d) * sin(x[i,] - mu_hat) / (1 + sin(x[i,] - mu_hat)^2) - c_mu_lambda %*% c_mu_mu_inv %*% phi)
    var <- var + delta %*% t(delta)
  }
  var <- var / n
  
  statistic <- t(central_seq) %*% solve(var) %*% central_seq
  
  p_value <- 1 - pchisq(statistic, d)
  
  names(statistic) = "Unspecified median optimal test"
  lamb0 = 0
  names(lamb0) = "value of the vector lambda"
  totval = list(statistic = statistic, p.value = p_value, null.value = lamb0, alternative = "two.sided",
                method = paste("Symmetry Optimal Test for Unspecified Median with", model, "model"), sample.size = n, 
                data.name = deparse(substitute(x)))
  class(totval) <- "htest"
  
  return(totval)
}



# ## Examples
# set.seed(1)
# x <- r_arctan(n = 5000, d = 1, lambda = c(0), model = 'wrappedcauchy', params = list(mu = 0, rho = 0.1))
# specified_median_test_arctan(x, mu = 0)
# unspecified_median_test_arctan(x, model = 'wrappedcauchy', params = list(mu = 0, rho = 0.1))
# 
# set.seed(1)
# x <- r_arctan(n = 5000, d = 2, lambda = c(0.6, 0.1), model = 'sine', params = list(mu = c(0,0), kappa = c(0.1, 0.1), rho = 0.2))
# specified_median_test_arctan(x, mu = c(0,0))
# unspecified_median_test_arctan(x, model = 'sine', params = list(kappa = c(0.1, 0.1), rho = 0.2))
# unspecified_median_test_arctan(x, model = 'bwc', params = list(kappa = c(0.1, 0.1), rho = 0.2))

# set.seed(1)
# x <- r_arctan(n = 1000, d = 3, lambda = c(10,20,30), model = 'twcc',
#               params = list(rho12 = 0.25, rho13 = 1, rho23 = 4, 
#                             marginals = rep('wrapped cauchy', 3),
#                             params1 = list(mu = 0, rho = 0.1),
#                             params2 = list(mu = 0, rho = 0.2),
#                             params3 = list(mu = 0, rho = 0.1)))
# specified_median_test_arctan(x, mu = c(0,0,0))
# unspecified_median_test_arctan(x, model = 'twcc', params = list(rho12 = 0.25, rho13 = 1, rho23 = 4, 
#                                                                 marginals = rep('wrapped cauchy', 3),
#                                                                 params1 = list(mu = 0, rho = 0.1),
#                                                                 params2 = list(mu = 0, rho = 0.1),
#                                                                 params3 = list(mu = 0, rho = 0.1)))

