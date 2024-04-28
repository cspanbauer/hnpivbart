## ivbart: Instrumental Variable BART with Diricihlet Process Mixtures
## Copyright (C) 2020 Robert McCulloch and Rodney Sparapani

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/GPL-2

ivbart = function(Z, T, Y, X=NULL,
                  Z.test=NULL,
                  T.test=NULL, 
                  X.test=NULL,
                  type1=1,type2=1,
                  burn=1000, nd=2000, burnf1=1000, burnf22=1000,
                  keepevery=10,
                  betabar=0, Abeta=NA,
                  m1=200, m2=200, nc=100,
                  power=2, base=0.95,
                  k1=2, k2=2,
                  sigmaf1=NA, sigmaf2=NA,
                  doDP=TRUE, do.tsls=TRUE,
                  v=0.17, nu=2.004, a=0.016,
                  sparse1=0,sparse2=0,
                  Imin=1, Imax=floor(0.1*n)+1, psi=0.5, gs=100,
                  centermeans=TRUE, in_sample=TRUE, aggregate_sigma=FALSE,
                  n=length(Y),
                  f1s = NULL, f22s = NULL,
                  mTs=0, mYs=0,
                  betas=NA, sTs=NA, sYs=NA, gammas=NA,
                  null.f1=FALSE, null.f21=FALSE,
                  ## betas=0, mTs=0, mYs=0, sTs=1, sYs=1, gammas=0,
                  printevery=100,quiet=TRUE
                  )
{

    ## Data checking
    X <- check_build_training_data(Z,T,X,n,"ivbart")
    ZX.ls <- check_build_test_data(Z.test,T.test,X.test,"ivbart")
    Z.test <- ZX.ls[[1]]; X.test <- ZX.ls[[2]]
    rm(ZX.ls)

    ## Getting tau values (prior varince of terminal node value)
    if(type2==1) {
        tau.ls <- get_tau(sigmaf1,sigmaf2,k1,k2,m1,m2,max(T)-min(T),max(Y)-min(Y))
        tauf1 <- tau.ls[1]; tauf2 <- tau.ls[2]; Abeta <- tauf2/(max(T)-min(T))
    }
    else {
        tau.ls <- get_tau(sigmaf1,sigmaf2,k1,k2,m1,m2,max(T)-min(T),6)
        tauf1 <- tau.ls[1]; tauf2 <- tau.ls[2]; Abeta <- tauf2/6
    }
    rm(tau.ls)

    ## Get DPM priors
    priold = getaprior(length(Y),Imin,Imax,psi,gs)
    priag = priold$p
    ag = priold$a
    rm(priold)

    ## Do TSLS for coefficients which can be used for starting values
    if(do.tsls){
        tsls.ls <- get_tsls_coef(T,Z,X,Y)
        tsls2.coef <- tsls.ls[[1]]
        L <- tsls.ls[[2]]
        rm(tsls.ls)
    }
    else{
        tsls2.coef <- 0
        L <- matrix(c(1,rep(stats::rnorm(1,0,.1),2),1),nrow=2)
    }

    
    ## Offset
    if(type2==1)
        offset=mean(Y)
    else if(type2==2)
        offset=qnorm(mean(Y))

    ## Starting values
    starts.ls <- get_starting_values(betas,sTs,sYs,gammas,L,null.f21)
    betas=starts.ls[1]; sTs=starts.ls[2]; sYs=starts.ls[3]; gammas=starts.ls[4]
    if(length(f1s)==0) f1s=rep(offset,n)
    if(length(f22s)==0) f22s=rep(offset,n)
    rm(starts.ls)

    

    include_output <- 0
    
    res = .Call("cbiv",
                t(Z),
                t(Z.test),
                t(X),
                t(X.test),
                T,
                Y,
                type1,
                type2,
                burn,
                nd,
                keepevery,
                burnf1,
                burnf22,
                m1, m2, nc,
                power, base,
                tauf1,
                tauf2,
                betabar, Abeta,
                doDP,
                v, nu, a, #base prior
                sparse1, sparse2,
                ag, priag, #alpha prior
                centermeans,
                offset, in_sample, aggregate_sigma,
                f1s, f22s,
                betas, mTs, mYs, sTs, gammas, sYs,
                null.f1, null.f21,
                printevery, quiet
                )
    
    res$check = NULL
    res$dnpart = res$dnpart
    res$dalpha = res$dalpha
    res$dbeta = res$dbeta
    if(in_sample){
        res$df1 = res$df1
        res$df22 = res$df22
        res$df2 = res$df22+matrix(rep(res$dbeta,n),nrow=nd)*matrix(rep(T,each=nd),nrow=nd)
    }
    res$dsigma1 = res$dsigma1
    res$dsigma2 = res$dsigma2
    res$dcov = res$dcov
    res$df1.test = res$df1.test
    res$df22.test = res$df22.test
    res$df2.test = matrix(rep(T.test,each=nd),nrow=nd)*matrix(rep(res$dbeta,length(T.test)),nrow=nd)+res$df22.test
    res$df1.varcount = res$df1.varcount
    res$df2.varcount = res$df2.varcount
    class(res) <- 'ivbart'
    return(res)
}
