
## loo.ivbart <-
##     function(x,
##              ...,
##              r_eff = NULL,
##              save_psis = FALSE,
##              cores = getOption("mc.cores", 1),
##              moment_match=FALSE,
##              k_threshold=0.7,
##              max_iters=25L)
## {  
##     if (is.null(r_eff))
##     {
##       r_eff=loo::relative_eff(exp(x$dLL),cores=cores,chain_id=x$chain)  
##     } 
##     loo <- loo.matrix(x$dLL,r_eff,save_psis,cores)
##     if(moment_match)
##     {
##         loo <- tryCatch(expr={loo_moment_match(x=x,loo=loo,cores=cores,
##                                                    max_iters=max_iters,r_eff=r_eff,
##                                                    k_threshold=k_threshold,
##                                                    post_draws=post_draws_ivbart,
##                                                    log_lik_i=log_lik_i_ivbart,
##                                                    unconstrain_pars=unconstrain_pars_ivbart,
##                                                    log_prob_upars=log_prob_upars_ivbart,
##                                                    log_lik_i_upars=log_lik_i_upars_ivbart)},
##                             error=function(e){
##                                 print(e)
##                                 return(loo)
##                             })
##     }
##     return(loo)
## }

## loo.hetivbart <-
##     function(x,
##              r_eff = NULL,
##              save_psis = FALSE,
##              cores = getOption("mc.cores", 1),
##              moment_match=FALSE,
##              k_threshold=0.7,
##              max_iters=25L,
##              ...)
## {  
##     if (is.null(r_eff))
##     {
##       r_eff=loo::relative_eff(exp(x$dLL),cores=cores,chain_id=x$chain)  
##     } 
##     loo.out <- loo.matrix(x$dLL,r_eff,save_psis,cores)
##     if(moment_match)
##     {
##         loo.out <- loo_moment_match(x=x,loo=loo.out,cores=cores,
##                                     max_iters=max_iters,r_eff=r_eff,
##                                     k_threshold=k_threshold,
##                                     post_draws=post_draws_hetivbart,
##                                     log_lik_i=log_lik_i_hetivbart,
##                                     unconstrain_pars=unconstrain_pars_hetivbart,
##                                     log_prob_upars=log_prob_upars_hetivbart,
##                                     log_lik_i_upars=log_lik_i_upars_hetivbart)
##     }
##     return(loo.out)
## }

## #setMethod("loo", "hetivbart", loo.hetivbart)

## loo.npivbart <-
##     function(x,
##              r_eff = NULL,
##              save_psis = FALSE,
##              cores = getOption("mc.cores", 1),
##              moment_match=FALSE,
##              k_threshold=0.7,
##              max_iters=25L,
##              ...)
## {  
##     if (is.null(r_eff))
##     {
##       r_eff=loo::relative_eff(exp(x$dLL),cores=cores,chain_id=x$chain)  
##     } 
##     loo.out <- loo.matrix(x$dLL,r_eff,save_psis,cores)
##     if(moment_match)
##     {
##         loo.out <- loo_moment_match(x=x,loo=loo.out,cores=cores,
##                                     max_iters=max_iters,r_eff=r_eff,
##                                     k_threshold=k_threshold,
##                                     post_draws=post_draws_npivbart,
##                                     log_lik_i=log_lik_i_npivbart,
##                                     unconstrain_pars=unconstrain_pars_npivbart,
##                                     log_prob_upars=log_prob_upars_npivbart,
##                                     log_lik_i_upars=log_lik_i_upars_npivbart)
##     }
##     return(loo.out)
## }

## loo.npivbart <-
##     function(x,
##              r_eff = NULL,
##              save_psis = FALSE,
##              cores = getOption("mc.cores", 1),
##              moment_match=FALSE,
##              k_threshold=0.7,
##              max_iters=25L,
##              ...)
## {  
##     if (is.null(r_eff))
##     {
##       r_eff=loo::relative_eff(exp(x$dLL),cores=cores,chain_id=x$chain)  
##     } 
##     loo.out <- loo.matrix(x$dLL,r_eff,save_psis,cores)
##     if(moment_match)
##     {
##         loo.out <- loo_moment_match(x=x,loo=loo.out,cores=cores,
##                                     max_iters=max_iters,r_eff=r_eff,
##                                     k_threshold=k_threshold,
##                                     post_draws=post_draws_npivbart,
##                                     log_lik_i=log_lik_i_npivbart,
##                                     unconstrain_pars=unconstrain_pars_npivbart,
##                                     log_prob_upars=log_prob_upars_npivbart,
##                                     log_lik_i_upars=log_lik_i_upars_npivbart)
##     }
##     return(loo.out)
## }

## loo.hetnpivbart <-
##     function(x,
##              r_eff = NULL,
##              save_psis = FALSE,
##              cores = getOption("mc.cores", 1),
##              moment_match=FALSE,
##              k_threshold=0.7,
##              max_iters=25L,
##              ...)
## {
##     if (is.null(r_eff))
##     {
##       r_eff=loo::relative_eff(exp(x$dLL),cores=cores,chain_id=x$chain)  
##     } 
##     loo.out <- loo.matrix(x$dLL,r_eff,save_psis,cores)
##     if(moment_match)
##     {
##         loo.out <- loo_moment_match(x=x,loo=loo.out,cores=cores,
##                                     max_iters=max_iters,r_eff=r_eff,
##                                     k_threshold=k_threshold,
##                                     post_draws=post_draws_hetnpivbart,
##                                     log_lik_i=log_lik_i_hetnpivbart,
##                                     unconstrain_pars=unconstrain_pars_hetnpivbart,
##                                     log_prob_upars=log_prob_upars_hetnpivbart,
##                                     log_lik_i_upars=log_lik_i_upars_hetnpivbart)
##     }
##     return(loo.out)
## }

## post_draws_ivbart <- function(x, ...)
## {
##     return(cbind(x$df,x$yhat))
## }
## log_lik_i_ivbart <- function(x, i, ...)
## {
##     return(x$dLL[,i])
## }
## unconstrain_pars_ivbart <- function(x, pars, ...)
## {
##     return(pars)
## }
## log_prob_upars_ivbart <- function(x, upars, ...)
## {
##     n <- 0.5*ncol(upars)
##     nd <- nrow(upars)
##     Tm <- matrix(rep(x$data[[2]],each=nd),ncol=n)
##     Ym <- matrix(rep(x$data[[1]],each=nd),ncol=n)
##     tf2 <- x$tau[[1]]^2
##     th12 <- x$tau[[2]]^2
##     sT2 <- x$dsigma1^2
##     sY2 <- x$dsigma2^2
##     vB <- 1/x$iv.prior[[2]]
##     m2 <- x$m2
##     lvarT <- log(x$m1)+log(tf2)+log(sT2)-log(x$m1*tf2+sT2)
##     muT <- x$m1*tf2*Tm/(x$m1*tf2+sT2)
##     varY <- sY2*(Tm*Tm*vB+th12*m2)/(Tm*Tm*vB+th12*m2+sY2)
##     muY <- (Tm*Tm*vB+m2*th12)*Ym/(Tm*Tm*vB+th12*m2+sY2)
##     out <- -0.5*log(2*pi)-0.5*lvarT-1/(2*exp(lvarT))*(upars[,1:n]-muT)^2+
##         -0.5*log(2*pi)-0.5*log(varY)-1/(2*varY)*(upars[,1:n+n]-muY)^2
##     return(rowSums(out))
## }
## log_lik_i_upars_ivbart <- function(x,upars,i,...)
## {
##     n <- 0.5*ncol(upars)
##     nd <- nrow(upars)
##     out <- rep(0,nd)
##     for(mcmc in 1:nd){
##         val <- matrix(c(x$data$T[i],x$data$Y[i]),ncol=1)
##         mu <- matrix(c(upars[mcmc,i],upars[mcmc,i+n]),ncol=1)
##         chl <- matrix(c(x$dsigma1[mcmc,i],x$dgamma[mcmc,i],0,x$dsY[mcmc,i]),ncol=2)
##         sig <- chl%*%t(chl)
##         sig <- round(sig,8)
##         out[mcmc] <- -log(2*pi)-0.5*log(det(sig))-0.5*t(val-mu)%*%solve(sig)%*%(val-mu)
##     }
##     return(out)
## }
## post_draws_hetivbart <- function(x, ...)
## {
##     return(cbind(x$df,x$yhat))
## }
## log_lik_i_hetivbart <- function(x, i, ...)
## {
##     return(x$dLL[,i])
## }
## unconstrain_pars_hetivbart <- function(x, pars, ...)
## {
##     return(pars)
## }
## log_prob_upars_hetivbart <- function(x, upars, ...)
## {
##     n <- 0.5*ncol(upars)
##     nd <- nrow(upars)
##     Tm <- matrix(rep(x$data[[2]],each=nd),ncol=n)
##     Ym <- matrix(rep(x$data[[1]],each=nd),ncol=n)
##     tf2 <- x$tau[[1]]^2
##     th12 <- x$tau[[2]]^2
##     sT2 <- x$dsigma1^2
##     sY2 <- x$dsigma2^2
##   #  vB <- 1/x$iv.prior[[2]]
##     m2 <- x$m2
##     lvarT <- log(x$m1)+log(tf2)+log(sT2)-log(x$m1*tf2+sT2)
##     muT <- x$m1*tf2*Tm/(x$m1*tf2+sT2)
##     varY <- sY2*m2*th12*(1+Tm*Tm)/(m2*th12*(1+Tm*Tm)+sY2)
##     muY <- m2*th12*(1+Tm*Tm)*Ym/(m2*th12*(1+Tm*Tm)+sY2)
##     out <- -0.5*log(2*pi)-0.5*lvarT-1/(2*exp(lvarT))*(upars[,1:n]-muT)^2+
##         -0.5*log(2*pi)-0.5*log(varY)-1/(2*varY)*(upars[,1:n+n]-muY)^2
##     return(rowSums(out))
## }
## log_lik_i_upars_hetivbart <- function(x,upars,i,...)
## {
##     n <- 0.5*ncol(upars)
##     nd <- nrow(upars)
##     out <- rep(0,nd)
##     for(mcmc in 1:nd){
##         val <- matrix(c(x$data$T[i],x$data$Y[i]),ncol=1)
##         mu <- matrix(c(upars[mcmc,i],upars[mcmc,i+n]),ncol=1)
##         chl <- matrix(c(x$dsigma1[mcmc,i],x$dgamma[mcmc,i],0,x$dsY[mcmc,i]),ncol=2)
##         sig <- chl%*%t(chl)
##         sig <- round(sig,8)
##         out[mcmc] <- -log(2*pi)-0.5*log(det(sig))-0.5*t(val-mu)%*%solve(sig)%*%(val-mu)
##     }
##     return(out)
## }
## post_draws_npivbart <- function(x, ...)
## {
##     return(cbind(x$df,x$yhat))
## }
## log_lik_i_npivbart <- function(x, i, ...)
## {
##     return(x$dLL[,i])
## }
## unconstrain_pars_npivbart <- function(x, pars, ...)
## {
##     return(pars)
## }
## log_prob_upars_npivbart <- function(x, upars, ...)
## {
##     n <- 0.5*ncol(upars)
##     nd <- nrow(upars)
##     Tm <- matrix(rep(x$data[[2]],each=nd),ncol=n)
##     Ym <- matrix(rep(x$data[[1]],each=nd),ncol=n)
##     tf2 <- x$tau[[1]]^2
##     th12 <- x$tau[[2]]^2
##     th22 <- x$tau[[3]]^2
##     sT2 <- x$dsigma1^2
##     sY2 <- x$dsigma2^2
##   #  vB <- 1/x$iv.prior[[2]]
##     m2 <- x$m2
##     lvarT <- log(x$m1)+log(tf2)+log(sT2)-log(x$m1*tf2+sT2)
##     muT <- x$m1*tf2*Tm/(x$m1*tf2+sT2)
##     varY <- sY2*m2*(th22+th12)/(th22*m2+th12*m2+sY2)
##     muY <- (th22*m2+m2*th12)*Ym/(th22*m2+th12*m2+sY2)
##     out <- -0.5*log(2*pi)-0.5*lvarT-1/(2*exp(lvarT))*(upars[,1:n]-muT)^2+
##         -0.5*log(2*pi)-0.5*log(varY)-1/(2*varY)*(upars[,1:n+n]-muY)^2
##     return(rowSums(out))
## }
## log_lik_i_upars_npivbart <- function(x,upars,i,...)
## {
##     n <- 0.5*ncol(upars)
##     nd <- nrow(upars)
##     out <- rep(0,nd)
##     for(mcmc in 1:nd){
##         val <- matrix(c(x$data$T[i],x$data$Y[i]),ncol=1)
##         mu <- matrix(c(upars[mcmc,i],upars[mcmc,i+n]),ncol=1)
##         chl <- matrix(c(x$dsigma1[mcmc,i],x$dgamma[mcmc,i],0,x$dsY[mcmc,i]),ncol=2)
##         sig <- chl%*%t(chl)
##         sig <- round(sig,8)
##         out[mcmc] <- -log(2*pi)-0.5*log(det(sig))-0.5*t(val-mu)%*%solve(sig)%*%(val-mu)
##     }
##     return(out)
## }
## post_draws_hetnpivbart <- function(x, ...)
## {
##     return(cbind(x$df,x$yhat))
## }
## log_lik_i_hetnpivbart <- function(x, i, ...)
## {
##     return(x$dLL[,i])
## }
## unconstrain_pars_hetnpivbart <- function(x, pars, ...)
## {
##     return(pars)
## }
## log_prob_upars_hetnpivbart <- function(x, upars, ...)
## {
##     n <- 0.5*ncol(upars)
##     nd <- nrow(upars)
##     Tm <- matrix(rep(x$data[[2]],each=nd),ncol=n)
##     Ym <- matrix(rep(x$data[[1]],each=nd),ncol=n)
##     tf2 <- x$tau[[1]]^2
##     th12 <- x$tau[[2]]^2
##     sT2 <- x$dsigma1^2
##     sY2 <- x$dsigma2^2
##   #  vB <- 1/x$iv.prior[[2]]
##     m2 <- x$m2
##     lvarT <- log(x$m1)+log(tf2)+log(sT2)-log(x$m1*tf2+sT2)
##     muT <- x$m1*tf2*Tm/(x$m1*tf2+sT2)
##     varY <- sY2*m2*th12/(th12*m2+sY2)
##     muY <- (m2*th12)*Ym/(th12*m2+sY2)
##     out <- -0.5*log(2*pi)-0.5*lvarT-1/(2*exp(lvarT))*(upars[,1:n]-muT)^2+
##         -0.5*log(2*pi)-0.5*log(varY)-1/(2*varY)*(upars[,1:n+n]-muY)^2
##     return(rowSums(out))
## }
## log_lik_i_upars_hetnpivbart <- function(x,upars,i,...)
## {
##     n <- 0.5*ncol(upars)
##     nd <- nrow(upars)
##     out <- rep(0,nd)
##     for(mcmc in 1:nd){
##         val <- matrix(c(x$data$T[i],x$data$Y[i]),ncol=1)
##         mu <- matrix(c(upars[mcmc,i],upars[mcmc,i+n]),ncol=1)
##         chl <- matrix(c(x$dsigma1[mcmc,i],x$dgamma[mcmc,i],0,x$dsY[mcmc,i]),ncol=2)
##         sig <- chl%*%t(chl)
##         sig <- round(sig,8)
##         out[mcmc] <- -log(2*pi)-0.5*log(det(sig))-0.5*t(val-mu)%*%solve(sig)%*%(val-mu)
##     }
##     return(out)
## }


## throw_loo_r_eff_warning <- function() {
##   warning(
##     "Relative effective sample sizes ('r_eff' argument) not specified.\n",
##     "For models fit with MCMC, the reported PSIS effective sample sizes and \n",
##     "MCSE estimates will be over-optimistic.",
##     call. = FALSE
##   )
## }

## prepare_psis_r_eff <- function(r_eff, len) {
##   if (isTRUE(is.null(r_eff) || all(is.na(r_eff)))) {
##     if (is.null(r_eff)) {
##       throw_psis_r_eff_warning()
##     }
##     r_eff <- rep(1, len)
##   } else if (length(r_eff) != len) {
##     stop("'r_eff' must have one value per observation.", call. = FALSE)
##   } else if (anyNA(r_eff)) {
##     stop("Can't mix NA and not NA values in 'r_eff'.", call. = FALSE)
##   }
##   return(r_eff)
## }
