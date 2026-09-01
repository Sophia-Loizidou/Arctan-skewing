#' Arctan Skewing: Basic Functions
#'
#' This file contains core functions for the arctan skewing transformation on circular
#' and toroidal data. It provides density estimation, random sampling generation, and
#' visualization tools for various circular distributions with arctan-based skewing.
#'
#' @author Sophia Loizidou
#' @license CC-BY-4.0 (Attribution 4.0 International)
#' @references 
#'   - Symmetry test on hypertorus: https://github.com/Sophia-Loizidou/Symmetry_test_on_hypertorus
#'   - Trivariate wrapped Cauchy copula: https://github.com/Sophia-Loizidou/Trivariate-wrapped-Cauchy-copula

## Read data generating code from Symmetry test and TWCC repositories
source('https://raw.githubusercontent.com/Sophia-Loizidou/Symmetry_test_on_hypertorus/main/Rcode/Generating_data.R')
source('https://raw.githubusercontent.com/Sophia-Loizidou/Trivariate-wrapped-Cauchy-copula/main/R%20code/functions%20for%20trivariate%20wrapped%20Cauchy%20copula.R')

## Parameter checks
param_checks <- function(d, lambda, model, params){
  ## checks for d
  if(all(d != c(1,2,3))){stop("The dimension 'd' must be 1, 2 or 3.")}
  
  ## checks for lambda
  for(i in 1:d){
    if (!is.numeric(lambda[i])) {stop("Argument 'lambda' must be a vector of real numbers.")}  
  }
  if (length(lambda) != d) {stop("Argument 'lambda' must be a d-dimensional vector.")}  
  
  if(d == 1){
    ## checks for model
    if (is.null(model)){
      warning("Argument 'model' is missing. By default, the von Mises distribution is employed, i.e., 'model='vonmises''.")
      model="vonmises"
    }
    if(all(model != c('vonmises', 'wrappedcauchy', 'cardioid', 'wrappednormal', 'katojones'))){
      stop("For d=1, argument 'model' must be one of 'vonmises', 'wrappedcauchy', 'cardioid', 'wrappednormal' or 'katojones'.")
    }
    
    ## checks for mu
    if(is.null(params$mu)){stop(paste0("The list 'params' should include a circular parameter 'mu'."))}
    if(!is.numeric(params$mu)){stop(paste0("The list 'params' should include a circular parameter 'mu'."))}
    if(length(params$mu) != d){stop("The list 'params' should include a circular parameter 'mu'.")}
    
    if(model == 'vonmises'){
      
      if(is.null(params$kappa)){stop(paste0("When model = '", model, "', the list 'params' should include a real number 'kappa'."))}
      if(!is.numeric(params$kappa)){stop(paste0("When model = '", model, "', the list 'params' should include a real number 'kappa'."))}
      if(length(params$kappa) != 1){stop(paste0("When model = '", model, "', the list 'params' should include a real number 'kappa'."))}
      
    } else if(model %in% c('cardioid', 'wrappedcauchy')){
      
      if(is.null(params$rho)){stop(paste0("When model = '", model, "', the list 'params' should include a real number 'rho'."))}
      if(!is.numeric(params$rho)){stop(paste0("When model = '", model, "', the list 'params' should include a real number 'rho'."))}
      if(length(params$rho) != 1){stop(paste0("When model = '", model, "', the list 'params' should include a real number 'rho'."))}
      
    } else if(model == 'wrappednormal'){
      
      if((is.null(params$rho) & is.null(params$sd)) | (!is.null(params$rho) & !is.null(params$sd))){
        stop(paste0("When model = '", model, "', the list 'params' should include either a parameter 'rho' in the interval [0,1] or a positive number 'sd'."))
      }
      
      if(!is.null(params$rho)){
        if(!is.numeric(params$rho)){stop(paste0("When model = '", model, "', the list 'params' should include a parameter 'rho' in the interval [0,1]."))}
        if(length(params$rho) != 1){stop(paste0("When model = '", model, "', the list 'params' should include a parameter 'rho' in the interval [0,1]."))}
        if((params$rho < 0) | (params$rho >= 1)){stop(paste0("When model = '", model, "', the list 'params' should include a parameter 'rho' in the interval [0,1]."))}
      } else if(!is.null(params$sd)){
        if(!is.numeric(params$sd)){stop(paste0("When model = '", model, "', the list 'params' should include a positive number 'sd'."))}
        if(length(params$sd) != 1){stop(paste0("When model = '", model, "', the list 'params' should include a positive number 'sd'."))}
        if(params$sd < 0){stop(paste0("When model = '", model, "', the list 'params' should include a positive number 'sd'."))}
      }
      
    } else if (model == 'katojones'){
      
      if(is.null(params$kappa)){stop(paste0("When model = '", model, "', the list 'params' should include a positive number 'kappa'."))}
      if(!is.numeric(params$kappa)){stop(paste0("When model = '", model, "', the list 'params' should include a positive number 'kappa'."))}
      if(length(params$kappa) != 1){stop(paste0("When model = '", model, "', the list 'params' should include a positive number 'kappa'."))}
      if(params$kappa < 0){stop(paste0("When model = '", model, "', the list 'params' should include a positive number 'kappa'."))}
      
      if(is.null(params$nu)){stop(paste0("When model = '", model, "', the list 'params' should include a circular parameter 'nu'."))}
      if(!is.numeric(params$nu)){stop(paste0("When model = '", model, "', the list 'params' should include a circular parameter 'nu'."))}
      if(length(params$nu) != 1){stop(paste0("When model = '", model, "', the list 'params' should include a circular parameter 'nu'."))}
      
      if(is.null(params$r)){stop(paste0("When model = '", model, "', the list 'params' should include a parameter 'r' in the interval [0,1)."))}
      if(!is.numeric(params$r)){stop(paste0("When model = '", model, "', the list 'params' should include a parameter 'r' in the interval [0,1)."))}
      if(length(params$r) != 1){stop(paste0("When model = '", model, "', the list 'params' should include a parameter 'r' in the interval [0,1)."))}
      if((params$r < 0) | (params$r >= 1)){stop(paste0("When model = '", model, "', the list 'params' should include a parameter 'r' in the interval [0,1)."))}
      
    }
  } else if (d == 2){
    ## checks for model
    if (is.null(model)){
      warning("Argument 'model' is missing. By default, the Sine distribution is employed, i.e., 'model='sine''.")
      model="sine"
    }
    if(all(model != c('sine', 'cosine', 'bwc'))){stop("For d=2, argument 'model' must be one of 'sine', 'cosine' or 'bwc'.")}
    
    ## checks for mu
    if(is.null(params$mu)){stop(paste0("The list 'params' should include a 2-dimensional vector 'mu'."))}
    if(!is.numeric(params$mu)){stop(paste0("The list 'params' should include a 2-dimensional vector 'mu'."))}
    if(length(params$mu) != d){stop("The list 'params' should include a 2-dimensional vector 'mu'.")}
    
    ## checks for kappa
    if(is.null(params$kappa)){stop(paste0("When model = '", model, "', the list 'params' should include a 2-dimensional vector 'kappa'."))}
    if(!is.numeric(params$kappa)){stop(paste0("When model = '", model, "', the list 'params' should include a 2-dimensional vector 'kappa'."))}
    if(length(params$kappa) != 2){stop(paste0("When model = '", model, "', the list 'params' should include a 2-dimensional vector 'kappa'."))}
    
    ## checks for rho
    if(is.null(params$rho)){stop(paste0("When model = '", model, "', the list 'params' should include a real number 'rho'."))}
    if(!is.numeric(params$rho)){stop(paste0("When model = '", model, "', the list 'params' should include a real number 'rho'."))}
    if(length(params$rho) != 1){stop(paste0("When model = '", model, "', the list 'params' should include a real number 'rho'."))}
    
  } else if (d == 3){
    
    ## checks for model
    if (is.null(model)){
      warning("Argument 'model' is missing. By default, the TWCC distribution is employed, i.e., 'model='twcc''.")
      model="twcc"
    }
    if(model != 'twcc'){stop("For d=3, argument 'model' must be 'twcc'")}
    
    marginals = params$marginals
    if(is.null(params$marginals)){
      warning("Argument 'params$marginals' is missing. By default, uniform marginals are employed, i.e., params$marginals = rep('uniform', 3).")
      marginals = rep('uniform', 3)
    }
    
    ## parameter checks are performed within the functions
    
    if(is.null(params$rho12) | is.null(params$rho13) | is.null(params$rho23)){
      stop(paste0("When model = '", model, "', the list 'params' should include 'rho12', 'rho13' and 'rho23'."))
    }
  }
  
  if(d != 3){marginals = NULL}
  
  return(list(model = model, marginals = marginals))
}

## Density estimation

d_arctan <- function(theta, lambda, model = NULL, params = list(mu = c(0,0), kappa = c(0.1, 0.1), rho = 0.2), param_check = T){
  if (!is.numeric(theta)) {stop("'theta', 'params$mu', and 'lambda' must be numeric")}
  if(is.null(dim(theta))){theta <- t(as.matrix(theta))}
  
  n <- dim(theta)[1]
  d <- dim(theta)[2]
  
  if(param_check){
    temp <- param_checks(d = d, lambda = lambda, model = model, params = params)
    model = temp[[1]]; params$marginals = temp[[2]]
  }
  
  require(circular)
  
  skewing_transf = function(theta, mu, lambda){
    temp_sum <- 0
    for(i in 1:d){
      temp_sum <- temp_sum + atan((lambda[i] * sin(theta[,i] - mu[i])) / (1 + sin(theta[,i] - mu[i])^2))
    }
    return(1 + (2 / (pi * d)) * temp_sum)
  }
  
  if(d == 1){
    
    if(model == 'vonmises'){
      devalmodel = dvonmises(theta, circular(params$mu), params$kappa)
    } else if(model %in% c('cardioid', 'wrappedcauchy')){
      devalmodel = eval(parse(text=paste("d",model,"(theta, circular(params$mu), params$rho)",sep="")))
    } else if(model == 'wrappednormal'){
      devalmodel = dwrappednormal(theta, circular(params$mu), params$rho, params$sd)
    }  else if(model == 'katojones'){
      devalmodel = dkatojones(theta, circular(params$mu), circular(params$nu), params$r, params$kappa)
    }
    
  } else if(d == 2){
    
    devalmodel = eval(parse(text=paste("d",model,"(theta[,1], theta[,2], params$mu[1], params$mu[2], params$kappa[1], params$kappa[2], params$rho)",sep="")))
    
  } else if(d == 3){
    
    devalmodel = dtri(theta, params$rho12, params$rho13, params$rho23, params$marginals, params$params1, params$params2, params$params3)
    
    params$mu = ifelse(is.null(params$params1$mu), 0, params$params1$mu)
    params$mu = c(params$mu, ifelse(is.null(params$params2$mu), 0, params$params2$mu))
    params$mu = c(params$mu, ifelse(is.null(params$params3$mu), 0, params$params3$mu))
  }
  
  return(devalmodel*skewing_transf(theta, params$mu, lambda))
}

# Examples
# require(circular)
# theta <- cbind(rvonmises(n = 200, mu = circular(0), kappa = 3), rvonmises(n = 200, mu = circular(0), kappa = 3), rvonmises(n = 200, mu = circular(0), kappa = 3))
# d_arctan(theta, lambda = c(1,10,1), model = "twcc", params = list(rho12 = 4, rho13=-0.25, rho23=-1))


## Generating mechanism
r_arctan <- function(n, d, lambda, model = NULL, 
                params = list(mu = c(0,0), kappa = c(0.1, 0.1), rho = 0.2)){ 
  
  require(circular)
  temp <- param_checks(d = d, lambda = lambda, model = model, params = params)
  model = temp[[1]]; params$marginals = temp[[2]]
  
  if(d == 1){
    
    if(model == 'vonmises'){
      newsample = rvonmises(n, circular(params$mu), params$kappa)
    } else if(model %in% c('cardioid', 'wrappedcauchy')){
      newsample = eval(parse(text=paste("r",model,"(n, circular(params$mu), params$rho)",sep="")))
    } else if(model == 'wrappednormal'){
      newsample = rwrappednormal(n, circular(params$mu), params$rho, params$sd)
    }  else if(model == 'katojones'){
      newsample = rkatojones(n, circular(params$mu), circular(params$nu), params$r, params$kappa)
    }
    newsample = as.matrix(newsample, ncol = 1)
    
  } else if(d == 2){
    
    newsample=eval(parse(text=paste("r",model,"(n, params$mu[1], params$mu[2], params$kappa[1], params$kappa[2], params$rho)",sep="")))
    
  } else if(d == 3){
    
    newsample = rtri(n, params$rho12, params$rho13, params$rho23, params$marginals, 
                     params$params1, params$params2, params$params3)
  }
  
  ssnewsample = newsample
  un = runif(n)
  sum_temp <- 0
  for(i in 1:d){
    sum_temp <- sum_temp + atan(lambda[i]*sin(newsample[,i] - params$mu[i]) / (1 + sin(newsample[,i] - params$mu[i])^2))
  }
  cons = un > (1 + 2/(pi*d)*sum_temp)/2
  
  for(i in 1:d){
    ssnewsample[cons,i]=-ssnewsample[cons,i]+2*params$mu[i]
  }
  
  return(ssnewsample)
  
}

# Examples
# r_arctan(n = 100, d = 1, lambda = c(1), model = 'katojones', 
#     params = list(mu = c(0), rho = c(0.2), kappa = 0.2, nu = 3))
# r_arctan(n = 100, d = 3, lambda = c(1,1,2), model = 'twcc',
#     params = list(rho12 = 1, rho13 = 0.25, rho23 = 4,
#                   params1 = list(mu = 0, rho = 0.2), params2 = list(mu = 0, rho = 0.2),
#                   params3 = list(mu = 0, rho = 0.2), marginals = c('wrapped Cauchy', 'von mises', 'von mises')))



## Contour plots
plot_arctan_contour = function(lambda = NULL, params, model = NULL, points = NULL, 
                               m = NULL, nn = 12, get_data = FALSE, cex = 1, add = FALSE){
  if(is.null(m)){m=50}
  zz=matrix(1:m^2,m)
  aa=NULL
  xx=seq(-pi,pi,len=m)
  yy=seq(-pi,pi,len=m)
  
  a=0.2 # 0.22
  ca=1.25
  
  for(j in 1:m){
    for(k in 1:m){
      zz[j,k] = d_arctan(c(xx[j],yy[k]), lambda, model, params)
    }
  }
  contour(xx,yy,zz,xlim=c(a-pi,pi-a),ylim=c(a-pi,pi-a),xlab=expression(x),
          ylab='',nlevels=nn,xaxt="n",yaxt="n", cex.lab = cex,
          main = bquote(paste(f = .(model), ", ", lambda, ' = ', f = .(lambda[1]), ', ', f = .(lambda[2]), ", ",
                              kappa, ' = ', f = .(params$kappa[1]), ', ', f = .(params$kappa[2]), ", ",
                              rho, ' = ', f = .(params$rho))),
          cex.main = cex, add = add)
  axis(1, at=(-4:4)*2*pi/8,labels=c(expression(paste(-pi)),"",expression(paste(-pi/2)),"","0","",expression(paste(pi,"/2")),"",expression(paste(pi))), las=1,col.axis="black", cex.axis = cex)
  axis(2, at=(-4:4)*2*pi/8,labels=c(expression(paste(-pi)),"",expression(paste(-pi/2)),"","0","",expression(paste(pi,"/2")),"",expression(paste(pi))), las=1,col.axis="black", cex.axis = cex)
  
  if(!is.null(points)){
    points(points, col = rgb(0, 0, 1, alpha = 0.1))
  }
  
  if(get_data){
    return(list(xx, yy, zz))
  } else {
    return(c(max(zz), min(zz)))
  }
}

