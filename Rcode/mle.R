#' Arctan Skewing: Maximum Likelihood Estimation
#'
#' This file implements maximum likelihood estimation (MLE) functions for the arctan
#' skewing transformation on circular and toroidal data. Uses the Rsolnp library for
#' optimization with support for symmetric and asymmetric submodels across dimensions.
#'
#' @author Sophia Loizidou
#' @license CC-BY-4.0 (Attribution 4.0 International)

## Basic functions for arctan skewing
source('~/basic_functions.R')

# MLE functions for TWCM
source('https://raw.githubusercontent.com/Sophia-Loizidou/Trivariate-wrapped-Cauchy-copula/main/R%20code/MLE.R')

'
Maximum likelihood estimation
Obtained with the Rsolnp library 
For the symmetric submodels set symmetric=TRUE 
ninipar is the number of initial starting points employed in the optimization algorithm
vecmax and vecmin are contain the largest and smallest values employed for the random initial points 
lowercons and uppercons are the lower and upper parameter constrains
'
mle_d1 <- function(x, symmetric, ninipar, model, vecmax, vecmin, lowercons, uppercons,
                   lambdamax, lambdamin, verbose){
  if (!is.null(vecmax)) {
    warning("The parameter 'vecmax' is only used for d=2, with model = c('sine', 'cosine', 'bwc')")
  }
  if (!is.null(vecmin)) {
    warning("The parameter 'vecmin' is only used for d=2, with model = c('sine', 'cosine', 'bwc')")
  }
  
  if(model == 'katojones'){
    p <- if (symmetric) 4L else 5L
  } else {
    p <- if (symmetric) 2L else 3L
  }
  
  if (!is.null(lowercons)) {
    warning("The parameter 'lowercons' is only used for d=2, with model = c('sine', 'cosine', 'bwc')")
  }
  if(model == 'vonmises') {
    lowercons <- c(-pi, 0)
  } else if (model == 'cardioid'){
    lowercons <- c(-pi, -0.5)
  } else if (model %in% c('wrappedcauchy', 'wrappednormal')){
    lowercons <- c(-pi, 0)
  } else if (model == 'katojones'){
    lowercons <- c(-pi, -pi, 0, 0)
  }
  
  if (!symmetric) lowercons <- c(lowercons, lambdamin)
  
  if (!is.null(uppercons)) {
    warning("The parameter 'uppercons' is only used for d=2, with model = c('sine', 'cosine', 'bwc')")
  }
  if(model == 'vonmises') {
    uppercons <- c(pi, Inf)
  } else if (model == 'cardioid'){
    uppercons <- c(pi, 0.5)
  } else if (model %in% c('wrappedcauchy', 'wrappednormal')){
    uppercons <- c(pi, 1)
  } else if (model == 'katojones'){
    uppercons <- c(pi, pi, 1, Inf)
  }
  
  if (!symmetric) uppercons <- c(uppercons, lambdamax)
  
  
  if (!is.numeric(lowercons) || length(lowercons) != p || anyNA(lowercons)) {
    stop(sprintf("`lowercons` must be a numeric vector of length %d.", p))
  }
  if (!is.numeric(uppercons) || length(uppercons) != p || anyNA(uppercons)) {
    stop(sprintf("`uppercons` must be a numeric vector of length %d.", p))
  }
  if (any(lowercons >= uppercons)) {
    stop("Each element of `lowercons` must be strictly smaller than `uppercons`.")
  }
  
  if(model == 'vonmises'){
    negloglik <- function(par) {
      lambda <- if (symmetric) c(0) else c(par[3])
      dens <- d_arctan(x, lambda = lambda, model = model, 
                       params = list(mu = par[1], kappa = par[2]),
                       param_check = F
      )
      
      if (any(!is.finite(dens)) || any(dens <= 0)) return(Inf)
      -sum(log(dens))
    }
  } else if (model == 'katojones'){
    negloglik <- function(par) {
      lambda <- if (symmetric) c(0) else c(par[5])
      dens <- d_arctan(x, lambda = lambda, model = model, 
                       params = list(mu = par[1], nu = par[2], r = par[3], kappa = par[4]),
                       param_check = F
      )
      
      if (any(!is.finite(dens)) || any(dens <= 0)) return(Inf)
      -sum(log(dens))
    }
  } else {
    negloglik <- function(par) {
      lambda <- if (symmetric) c(0) else c(par[3])
      dens <- d_arctan(x, lambda = lambda, model = model, 
                       params = list(mu = par[1], rho = par[2]),
                       param_check = F
      )
      
      if (any(!is.finite(dens)) || any(dens <= 0)) return(Inf)
      -sum(log(dens))
    }
  }
  
  if(model == 'katojones'){
    draw_init <- function() {
      init <- numeric(p)
      init[1] <- runif(1, -pi, pi)
      init[2] <- runif(1, lowercons[2], uppercons[2])
      init[3] <- runif(1, lowercons[3], uppercons[3])
      init[4] <- runif(1, lowercons[4], min(100, uppercons[4]))
      if (!symmetric) {
        init[5] <- runif(1, max(-10,lambdamin), min(10,lambdamax))
      }
      init
    }
  } else {
    draw_init <- function() {
      init <- numeric(p)
      init[1] <- runif(1, -pi, pi)
      init[2] <- runif(1, lowercons[2], min(100, uppercons[2]))
      if (!symmetric) {
        init[3] <- runif(1, max(-10,lambdamin), min(10,lambdamax))
      }
      init
    }
  }
  
  best_value <- Inf
  best_par <- rep(NA_real_, p)
  best_fit <- NULL
  
  for (i in seq_len(ninipar)) {
    if(verbose){if (i %% 20 == 0){cat(i/ninipar*100, '% \n')}}
    
    init <- draw_init()
    
    fit <- tryCatch(
      R.utils::withTimeout(
        Rsolnp::solnp(
          pars = init,
          fun = negloglik,
          LB = lowercons,
          UB = uppercons,
          control = list(trace = 0)
        ),
        timeout = 10,
        onTimeout = "error"
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    value <- fit$values[length(fit$values)]
    if (is.finite(value) && value < best_value) {
      best_value <- value
      best_par <- fit$pars
      best_fit <- fit
    }
  }
  
  if (is.null(best_fit)) {
    stop("All optimization attempts failed.")
  }
  
  out <- list(
    model = model,
    mu = best_par[1]
  )
  
  if(model == 'katojones'){
    out$nu = best_par[2]
    out$r = best_par[3]
    out$kappa = best_par[4]
    
    if (!symmetric) {
      out$lambda <- best_par[5]
    }
  } else {
    if(model == 'vonmises'){
      out$kappa = best_par[2]
    } else {
      out$rho = best_par[2]
    }
    if (!symmetric) {
      out$lambda <- best_par[3]
    }
  }
  
  out$LL = -best_value
  out$convergence = best_fit$convergence
  out$pars = best_par
  out$fit = best_fit
  
  out
}

mle_d2 <- function(x, symmetric, ninipar, model, vecmax, vecmin, lowercons, uppercons,
                   lambdamax, lambdamin, verbose){
  
  if (is.null(vecmax)) {
    valmax <- if (model %in% c("sine", "cosine")) 10 else 1
    vecmax <- rep(valmax, 3L)
  }
  if (is.null(vecmin)) {
    valmin <- if (model %in% c("sine", "cosine")) -10 else -1
    vecmin <- c(0, 0, valmin)
  }
  
  if (!is.numeric(vecmax) || length(vecmax) != 3L || anyNA(vecmax)) {
    stop("`vecmax` must be a numeric vector of length 3.")
  }
  if (!is.numeric(vecmin) || length(vecmin) != 3L || anyNA(vecmin)) {
    stop("`vecmin` must be a numeric vector of length 3.")
  }
  
  p <- if (symmetric) 5L else 7L
  
  if (is.null(lowercons)) {
    lowercons <- c(-pi, -pi, vecmin)
    if (model %in% c("sine")) lowercons[5] <- -Inf
    if (model %in% c("cosine")) lowercons[5] <- -15 ## numerical issues with smaller numbers
    if (!symmetric) lowercons <- c(lowercons, lambdamin, lambdamin)
  }
  if (is.null(uppercons)) {
    uppercons <- c(pi, pi, vecmax)
    if (model %in% c("sine", "cosine")) uppercons[3:5] <- Inf
    if (!symmetric) uppercons <- c(uppercons, lambdamax, lambdamax)
  }
  
  if (!is.numeric(lowercons) || length(lowercons) != p || anyNA(lowercons)) {
    stop(sprintf("`lowercons` must be a numeric vector of length %d.", p))
  }
  if (!is.numeric(uppercons) || length(uppercons) != p || anyNA(uppercons)) {
    stop(sprintf("`uppercons` must be a numeric vector of length %d.", p))
  }
  if (any(lowercons >= uppercons)) {
    stop("Each element of `lowercons` must be strictly smaller than `uppercons`.")
  }
  
  negloglik <- function(par) {
    lambda <- if (symmetric) c(0, 0) else c(par[6], par[7])
    dens <- d_arctan(x, lambda = lambda, model = model, 
                     params = list(mu = c(par[1], par[2]),kappa = c(par[3], par[4]), rho = par[5]),
                     param_check = F
    )
    
    if (any(!is.finite(dens)) || any(dens <= 0)) return(Inf)
    -sum(log(dens))
  }
  
  draw_init <- function() {
    init <- numeric(p)
    init[1] <- runif(1, -pi, pi)
    init[2] <- runif(1, -pi, pi)
    init[3] <- runif(1, vecmin[1], vecmax[1])
    init[4] <- runif(1, vecmin[2], vecmax[2])
    init[5] <- runif(1, vecmin[3], vecmax[3])
    if (!symmetric) {
      init[6] <- runif(1, max(-10, lambdamin), min(10,lambdamax))
      init[7] <- runif(1, max(-10, lambdamin), min(10,lambdamax))
    }
    init
  }
  
  best_value <- Inf
  best_par <- rep(NA_real_, p)
  best_fit <- NULL
  
  for (i in seq_len(ninipar)) {
    if(verbose){if (i %% 20 == 0){cat(i/ninipar*100, '% \n')}}
    
    init <- draw_init()
    
    fit <- tryCatch(
      R.utils::withTimeout(
        Rsolnp::solnp(
          pars = init,
          fun = negloglik,
          LB = lowercons,
          UB = uppercons,
          control = list(trace = 0)
        ),
        timeout = 10,
        onTimeout = "error"
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    value <- fit$values[length(fit$values)]
    if (is.finite(value) && value < best_value) {
      best_value <- value
      best_par <- fit$pars
      best_fit <- fit
    }
  }
  
  if (is.null(best_fit)) {
    stop("All optimization attempts failed.")
  }
  
  out <- list(
    model = model,
    mu1 = best_par[1],
    mu2 = best_par[2],
    kappa1 = best_par[3],
    kappa2 = best_par[4],
    rho = best_par[5]
  )
  
  if (!symmetric) {
    out$lambda1 <- best_par[6]
    out$lambda2 <- best_par[7]
  }
  
  out$LL = -best_value
  out$convergence = best_fit$convergence
  out$pars = best_par
  out$fit = best_fit
  
  out
}


mle.copula_arctan_skewing = function(x, symmetric, ninipar, marginals, vecmax, vecmin, lambdamax, lambdamin, verbose){ 
  
  ## Checks for the inputs
  
  if(!is.matrix(x)) { stop("Argument 'x' must be a numeric matrix")}
  if(dim(x)[2]!=3) { stop("Argument 'x' must be a matrix with three columns")}
  
  if (!is.numeric(ninipar)) {
    warning("Argument 'ninipar' must be a positive integer number. Default value of 'ninipar' was used")
    ninipar=100
  }
  if ((length(ninipar) != 1) | (ninipar%%1 != 0) | (ninipar <= 0)){
    warning("Argument 'ninipar' must be a positive integer number Default value of 'ninipar' was used")
    ninipar=100
  }
  
  if(is.null(vecmax)) vecmax = 10
  if (!is.numeric(vecmax)) {stop("Argument 'vecmax' must be a real number")}
  if (length(vecmax) != 1) {stop("Argument 'vecmax' must be a real number")} 
  
  if(is.null(vecmin)) vecmin = -10
  if (!is.numeric(vecmin)) {stop("Argument 'vecmin' must be a real number")}
  if (length(vecmin) != 1) {stop("Argument 'vecmin' must be a real number")}
  
  marginals <- tolower(marginals)
  if(length(marginals) != 3){stop("'marginals' should be a vector of length 3 with elements 'uniform', 'wrapped cauchy', 'cardioid', 
         'vonMises' or 'katojones' for the circular part and 'weibull' for the linear part")}
  if(any(!(marginals %in% c('uniform', 'wrapped cauchy', 'cardioid', 'vonmises', 'katojones', 'weibull')))){
    stop("'marginals' should be a vector of length 3 with elements 'uniform', 'wrapped cauchy', 'cardioid', 
         'vonMises' or 'katojones' for the circular part and 'weibull' for the linear part")
  }
  
  require(Rsolnp)
  require(CircStats)
  require(Directional)
  require(VGAM)
  require(circular)
  
  ## First calculate parameters for the marginals
  
  pars1 <- list(); pars2 <- list(); pars3 <- list()
  if (any(marginals == 'uniform')){
    index <- which(marginals == 'uniform')
    for(i in 1:length(index)){
      eval(parse(text = paste('uniform_x', index[i], ' = x[,', index[i], ']', sep ='')))
      eval(parse(text = paste('dens', index[i], ' = 1', sep ='')))
    }
  } 
  if (any(marginals == 'wrapped cauchy')){
    index <- which(marginals == 'wrapped cauchy')
    for(i in 1:length(index)){
      eval(parse(text = paste('pars', index[i], ' = mle.marginals.wrpcauchy(x[,', index[i], '], ninipar = ninipar, verbose = verbose)', sep ='')))
      eval(parse(text = paste('uniform_x', index[i], ' = vector_to_uniform(x[,', index[i], "], marginals = 'wrapped cauchy', params = pars", index[i], ')', sep ='')))
      eval(parse(text = paste('dens', index[i], ' = (2*pi) * dwrpcauchy(x[,', index[i],'],pars', index[i],'$mu,pars', index[i], '$rho)', sep ='')))
    }
  } 
  if (any(marginals == 'cardioid')){
    index <- which(marginals == 'cardioid')
    for(i in 1:length(index)){
      eval(parse(text = paste('pars', index[i], ' = cardio.mle(x[,', index[i], '], rads = TRUE)', sep ='')))
      eval(parse(text = paste('pars', index[i], " = list('mu' = as.numeric(pars", index[i], "$param[1]), 'rho' = as.numeric(pars", index[i], "$param[2]))", sep ='')))
      eval(parse(text = paste('uniform_x', index[i], ' = vector_to_uniform(x[,', index[i], "], marginals = 'cardioid', params = pars", index[i], ')', sep ='')))
      eval(parse(text = paste('dens', index[i], ' = (2*pi) * dcard(x[,', index[i],'],pars', index[i],'$mu,pars', index[i], '$rho)', sep ='')))
    }
  } 
  if (any(marginals == 'vonmises')){
    index <- which(marginals == 'vonmises')
    for(i in 1:length(index)){
      suppressWarnings({
        eval(parse(text = paste('pars', index[i], ' = mle.vonmises(circular(x[,', index[i], ']))', sep ='')))
        eval(parse(text = paste('pars', index[i], " = list('mu' = as.numeric(pars", index[i], "$mu), 'kappa' = pars", index[i], "$kappa)", sep ='')))
        eval(parse(text = paste('uniform_x', index[i], ' = vector_to_uniform(x[,', index[i], "], marginals = 'vonmises', params = pars", index[i], ')', sep ='')))
        eval(parse(text = paste('dens', index[i], ' = (2*pi) * dvonmises(circular(x[,', index[i],']),circular(pars', index[i],'$mu),pars', index[i], '$kappa)', sep ='')))
      })
    }
  } 
  if(any(marginals == 'katojones')){
    index <- which(marginals == 'katojones')
    for(i in 1:length(index)){
      eval(parse(text = paste('pars', index[i], ' = mle.marginals.katojones(x[,', index[i], '], ninipar = ninipar, verbose = verbose)', sep ='')))
      eval(parse(text = paste('uniform_x', index[i], ' = vector_to_uniform(x[,', index[i], "], marginals = 'katojones', params = pars", index[i], ')', sep ='')))
      eval(parse(text = paste('dens', index[i], ' = (2*pi) * dkatojones(x[,', index[i],'],pars', index[i],'$mu,pars', index[i], '$gamma,pars', index[i], '$rho,pars', index[i], '$lambda)', sep ='')))
    }
  }
  if(any(marginals == 'weibull')){
    require(EnvStats)
    index <- which(marginals == 'weibull')
    for(i in 1:length(index)){
      eval(parse(text = paste('pars', index[i], ' = eweibull(x[,', index[i], '], method = "mle")', sep ='')))
      eval(parse(text = paste('pars', index[i], " = list('shape' = as.numeric(pars", index[i], "$parameters[1]), 'scale' = as.numeric(pars", index[i], "$parameters[2]))", sep ='')))
      eval(parse(text = paste('uniform_x', index[i], ' = vector_to_uniform(x[,', index[i], "], marginals = 'weibull', params = pars", index[i], ')', sep ='')))
      eval(parse(text = paste('dens', index[i], ' = dweibull(x[,', index[i],'],pars', index[i],'$shape,pars', index[i], '$scale)', sep ='')))
    }
  }
  if(marginals[1] != 'uniform') names(pars1) <- paste(names(pars1), 1, sep = '')
  if(marginals[2] != 'uniform') names(pars2) <- paste(names(pars2), 2, sep = '')
  if(marginals[3] != 'uniform') names(pars3) <- paste(names(pars3), 3, sep = '')
  
  dens_marginals <- dens1*dens2*dens3
  
  uniform_x <- cbind(uniform_x1, uniform_x2, uniform_x3)
  
  copula_params <- mle.copula.params_arctan_skewing(uniform_x, symmetric, ninipar, vecmax, vecmin, lambdamax, lambdamin, verbose)
  
  out <- list(
    marginals = marginals,
    rho12 = copula_params$final_rho12,
    rho13 = copula_params$final_rho13,
    rho23 = copula_params$final_rho23
  )
  
  if (!symmetric) {
    out$lambda1 <- copula_params$final_lambda1
    out$lambda2 <- copula_params$final_lambda2
    out$lambda3 <- copula_params$final_lambda3
  }
  
  out$params1 = pars1
  out$params2 = pars2
  out$params3 = pars3
  out$LL = copula_params$final_LL + sum(log(dens_marginals))
  out$convergence = copula_params$final_fit$convergence
  # out$pars = copula_params$final_fit$pars
  out$fit = copula_params$final_fit
  
  out
}

mle.copula.params_arctan_skewing <- function(uniform_x, symmetric, ninipar, vecmax, vecmin, lambdamax, lambdamin, verbose){
  p <- ifelse(symmetric, 2, 5)
  
  LB = c(-1, -Inf); UB = c(1, Inf)
  if(!symmetric){
    LB <- c(LB, rep(lambdamin, 3)); UB <- c(UB, rep(lambdamax, 3))
  }
  
  draw_init <- function() {
    init <- numeric(p)
    init[1] = runif(1,-1,1)
    init[2] = runif(1,vecmin,vecmax)
    if (!symmetric) {
      init[3] <- runif(1, max(-10, lambdamin), min(10,lambdamax))
      init[4] <- runif(1, max(-10, lambdamin), min(10,lambdamax))
      init[5] <- runif(1, max(-10, lambdamin), min(10,lambdamax))
    }
    init
  }
  
  ## (k,l) = (1,2)
  
  ## - Log-likelihood because the function used finds the minimum of the input function
  llpar12=function(par){
    rho23 = par[2]
    rho13 = (1 + sqrt(1 + 4*abs(rho23)^3)) / (2*abs(rho23)^2) * 1 / par[1]
    rho12 = 1/(rho23*rho13)
    
    lambda <- if (symmetric) c(0, 0, 0) else c(par[3], par[4], par[5])
    dens <- d_arctan(uniform_x, lambda = lambda, model = 'twcc', 
                     params = list(rho12 = rho12, rho13 = rho13, rho23 = rho23, marginals = rep('uniform', 3)),
                     param_check = F
    )
    
    if (any(!is.finite(dens)) || any(dens <= 0)) return(Inf)
    -sum(log(dens))
  }
  
  valllfin12=Inf
  paramfin12=numeric()
  fit12 = NULL
  
  for(iterip in 1:ninipar){
    if(verbose){print(iterip)}
    
    inipar = draw_init()
    
    fit <- tryCatch(
      R.utils::withTimeout(
        Rsolnp::solnp(
          pars = inipar,
          fun = llpar12,
          LB = LB,
          UB = UB,
          control = list(trace = 0)
        ),
        timeout = 10,
        onTimeout = "error"
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    valll12 <- fit$values[length(fit$values)]
    if (is.finite(valll12) && valll12 < valllfin12) {
      valllfin12 <- valll12
      paramfin12 <- fit$pars
      fit12 <- fit
    }
  }
  
  
  
  ## (k,l) = (2,3)
  
  ## - Log-likelihood because the function used finds the minimum of the input function
  llpar23=function(par){
    rho13 = par[2]
    rho12 = (1 + sqrt(1 + 4*abs(rho13)^3)) / (2*abs(rho13)^2) * 1 / par[1]
    rho23 = 1/(rho12*rho13)
    
    lambda <- if (symmetric) c(0, 0, 0) else c(par[3], par[4], par[5])
    dens <- d_arctan(uniform_x, lambda = lambda, model = 'twcc', 
                     params = list(rho12 = rho12, rho13 = rho13, rho23 = rho23, marginals = rep('uniform', 3)),
                     param_check = F
    )
    
    if (any(!is.finite(dens)) || any(dens <= 0)) return(Inf)
    -sum(log(dens))
  }
  
  valllfin23=Inf
  paramfin23=numeric()
  fit23 = NULL
  
  for(iterip in 1:ninipar){
    if(verbose){print(iterip)}
    
    inipar = draw_init()
    
    fit <- tryCatch(
      R.utils::withTimeout(
        Rsolnp::solnp(
          pars = inipar,
          fun = llpar23,
          LB = LB,
          UB = UB,
          control = list(trace = 0)
        ),
        timeout = 10,
        onTimeout = "error"
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    valll23 <- fit$values[length(fit$values)]
    if (is.finite(valll23) && valll23 < valllfin23) {
      valllfin23 <- valll23
      paramfin23 <- fit$pars
      fit23 <- fit
    }
  }
  
  
  ## (k,l) = (3,1)
  
  ## - Log-likelihood because the function used finds the minimum of the input function
  llpar13=function(par){
    rho12 = par[2]
    rho23 = (1 + sqrt(1 + 4*abs(rho12)^3)) / (2*abs(rho12)^2) * 1 / par[1]
    rho13 = 1/(rho23*rho12)
    
    lambda <- if (symmetric) c(0, 0, 0) else c(par[3], par[4], par[5])
    dens <- d_arctan(uniform_x, lambda = lambda, model = 'twcc', 
                     params = list(rho12 = rho12, rho13 = rho13, rho23 = rho23, marginals = rep('uniform', 3)),
                     param_check = F
    )
    
    if (any(!is.finite(dens)) || any(dens <= 0)) return(Inf)
    -sum(log(dens))
  }
  
  valllfin13=Inf
  paramfin13=numeric()
  fit13 = NULL
  
  for(iterip in 1:ninipar){
    if(verbose) {print(iterip)}
    
    inipar = draw_init()
    
    fit <- tryCatch(
      R.utils::withTimeout(
        Rsolnp::solnp(
          pars = inipar,
          fun = llpar13,
          LB = LB,
          UB = UB,
          control = list(trace = 0)
        ),
        timeout = 10,
        onTimeout = "error"
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    valll13 <- fit$values[length(fit$values)]
    if (is.finite(valll13) && valll13 < valllfin13) {
      valllfin13 <- valll13
      paramfin13 <- fit$pars
      fit13 <- fit
    }
  }
  
  if((valllfin12 < valllfin13) && (valllfin12 < valllfin23)){
    paramfin = paramfin12
    final_rho23 = paramfin[2]
    final_rho13 = (1 + sqrt(1 + 4*abs(paramfin[2])^3)) / (2*abs(paramfin[2])^2) * 1/paramfin[1]
    final_rho12 = 1/(final_rho23*final_rho13)
    
    final_LL=-valllfin12 ## Estimated LL value
    final_fit = fit12
    
  } else if((valllfin23 < valllfin12) && (valllfin23 < valllfin13)){
    paramfin = paramfin23
    final_rho13 = paramfin[2]
    final_rho12 = (1 + sqrt(1 + 4*abs(paramfin[2])^3)) / (2*abs(paramfin[2])^2) * 1/paramfin[1]
    final_rho23 = 1/(final_rho12*final_rho13)
    
    final_LL=-valllfin23 ## Estimated LL value
    final_fit = fit23
    
  } else {
    paramfin = paramfin13
    final_rho12 = paramfin[2]
    final_rho23 = (1 + sqrt(1 + 4*abs(paramfin[2])^3)) / (2*abs(paramfin[2])^2) * 1/paramfin[1]
    final_rho13 = 1/(final_rho23*final_rho12)
    
    final_LL=-valllfin13 ## Estimated LL value
    final_fit = fit13
    
  }
  if(symmetric){
    final_lambda1 = 0
    final_lambda2 = 0
    final_lambda3 = 0
  } else {
    final_lambda1 = paramfin[3]
    final_lambda2 = paramfin[4]
    final_lambda3 = paramfin[5]
  }
  
  return(list('final_rho12' = final_rho12, 'final_rho23' = final_rho23, 'final_rho13' = final_rho13, 
              'final_lambda1' = final_lambda1, 'final_lambda2' = final_lambda2, 'final_lambda3' = final_lambda3, 
              'final_LL' = final_LL, 'final_fit' = final_fit))
}

mle_arctan <- function(x, symmetric = FALSE, ninipar = NULL, 
                       model = c("vonmises", "wrappedcauchy", "cardioid", "wrappednormal", "katojones", 
                                 "sine", "cosine", "bwc",
                                 "twcc"), 
                       marginals = NULL,
                       vecmax = NULL, vecmin = NULL, lowercons = NULL, uppercons = NULL,
                       lambdamax = Inf, lambdamin = -Inf, seed = NULL, verbose = F) {
  
  model <- match.arg(model)
  
  if(model %in% c("vonmises", "wrappedcauchy", "cardioid", "wrappednormal", "katojones")){
    d = 1
  } else if (model %in% c("twcc")){
    d = 3
  } else {
    d = 2
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix.")
  }
  if (ncol(x) != d) {
    stop(paste("For model =", model, ", `x` must have exactly", d, "columns."))
  }
  if (!is.logical(symmetric) || length(symmetric) != 1L || is.na(symmetric)) {
    stop("`symmetric` must be TRUE or FALSE.")
  }
  
  default_ninipar <- if (model %in% c("cosine", "sine")) 1000L else 100L
  if (is.null(ninipar)) ninipar <- default_ninipar
  if (!is.numeric(ninipar) || length(ninipar) != 1L || ninipar <= 0 || ninipar %% 1 != 0) {
    stop("`ninipar` must be a positive integer.")
  }
  ninipar <- as.integer(ninipar)
  
  if (!is.numeric(lambdamax) || length(lambdamax) != 1L || is.na(lambdamax)) {
    stop("`lambdamax` must be a numeric scalar.")
  }
  if (!is.numeric(lambdamin) || length(lambdamin) != 1L || is.na(lambdamin)) {
    stop("`lambdamin` must be a numeric scalar.")
  }
  if (lambdamin >= lambdamax) {
    stop("Require `lambdamin < lambdamax`.")
  }
  
  if(d == 1){
    mle <- mle_d1(x, symmetric, ninipar, model, vecmax, vecmin, lowercons, uppercons,
                  lambdamax, lambdamin, verbose)
  } else if (d == 2){
    mle <- mle_d2(x, symmetric, ninipar, model, vecmax, vecmin, lowercons, uppercons,
                  lambdamax, lambdamin, verbose)
  } else {
    if(is.null(marginals)){
      warning("Argument 'marginals' is missing. By default, circular uniform distribution is emplyed, i.e. 'marginals = rep('uniform', 3)'.")
      marginals = rep('uniform', 3)
    }
    mle <- mle.copula_arctan_skewing(x, symmetric, ninipar, marginals, vecmax, vecmin, 
                                     lambdamax, lambdamin, verbose)
  }
  
  return(mle)
}

# Examples
# set.seed(1)
# data <- r_arctan(n = 5000, d = 2, lambda = c(1, 2), model = 'sine', params = list(mu = c(0,0), kappa = c(0.1, 0.1), rho = 0.2))
# mle_arctan(data, seed = 1, model = 'sine')
# x <- r_arctan(n = 5000, d = 1, lambda = c(1), model = 'vonmises', params = list(mu = 0, kappa = 0.1))
# mle_arctan(x, seed = 1, model = 'vonmises', symmetric=F)
# 
# set.seed(1)
# data <- r_arctan(n = 5000, d = 3, lambda = c(10,10,10), model = 'twcc',
#                  params = list(rho12 = 1, rho13 = 0.25, rho23 = 4,
#                                params1 = list(mu = 0, rho = 0.2), params2 = list(mu = 0, rho = 0.2),
#                                params3 = list(mu = 0, rho = 0.2),
#                                marginals = rep('wrapped cauchy', 3)))
# mle_arctan(data, seed = 1, model = 'twcc', marginals = rep('wrapped Cauchy', 3))
# mle_arctan(data, seed = 1, model = 'twcc', marginals = c('wrapped Cauchy', 'vonMises', 'vonMises'), symmetric = T)
