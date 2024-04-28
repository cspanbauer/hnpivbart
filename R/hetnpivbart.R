
## ivbart: Instrumental Variable BART with Dirichlet Process Mixtures
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

hetnpivbart = function(Z, T, Y, X=NULL,
                       Z.test=NULL, 
                       T.test=NULL,
                       X.test=NULL,
                       type1=1,type2=1,
                       burn=1000, nd=2000, burnf1=1000, burnf2=1000,
                       keepevery=10,
                       m1=200, m2=200, nc=100,
                       power=2, base=0.95,
                       k1=2, k2=2,
                       sigmaf1=NA, sigmaf2=NA,
                       doDP=TRUE, do.tsls=TRUE,
                       v=0.17, nu=2.004, a=0.016,
                       sparse1=0, sparse2=0,
                       Imin=1, Imax=floor(0.1*n)+1, psi=0.5, gs=100,
                       centermeans=TRUE, in_sample=TRUE, aggregate_sigma=FALSE,
                       n = length(Y),
                       f1s = rep(0, n), f2s = rep(0, n),
                       mTs=0, mYs=0, sTs=NA, sYs=NA, gammas=NA,
                       null.f1=FALSE, null.f2=FALSE,
                                        #   betas=rep(0,2), mTs=0, mYs=0, sTs=1, sYs=1, gammas=0,
                       printevery=100, quiet=TRUE
                       )
{


    
    ## Data checking
    TX <- check_build_training_data(Z,T,X,n,"hetnpivbart")
    if(null.f2) TX <- X

    ZXT.ls <- check_build_test_data(Z.test,T.test,X.test,"hetnpivbart")
    Z.test <- ZXT.ls[[1]];
    if(null.f2) {
        TX.test <- X.test
    }
    else TX.test <- ZXT.ls[[2]] 
    rm(ZXT.ls)

    ## Getting tau values (prior variance of terminal node value)
    if(type2==1) {
        tau.ls <- get_tau(sigmaf1,sigmaf2,k1,k2,m1,m2,max(T)-min(T),max(Y)-min(Y))
        tauf1 <- tau.ls[1]; tauf2 <- tau.ls[2]
    }
    else {
        tau.ls <- get_tau(sigmaf1,sigmaf2,k1,k2,m1,m2,max(T)-min(T),6)
        tauf1 <- tau.ls[1]; tauf2 <- tau.ls[2]
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
        L <- matrix(c(1,rep(0.01,2),1),nrow=2)
    }

    ## Starting values
    starts.ls <- get_starting_values(0,sTs,sYs,gammas,L,null.f2)
    sTs=starts.ls[2]; sYs=starts.ls[3]; gammas=starts.ls[4]
    rm(starts.ls)

    include_output <- 0

    ## Offset
    if(type2==1) {
        offset=mean(Y)
    } else {
        offset=qnorm(mean(Y))
    }

    res = .Call("cbhetnpiv",
                t(Z),
                t(Z.test),
                T,
                t(TX),
                t(TX.test),
                Y,
                1,
                type2,
                burn,
                nd,
                keepevery,
                burnf1,
                burnf2,
                m1, m2, nc,
                power, base,
                tauf1,
                tauf2,
                doDP,
                v, nu, a, #base prior
                sparse1,
                sparse2,
                ag, priag, #alpha prior
                centermeans,
                offset,
                in_sample,
                aggregate_sigma,
                f1s, f2s,
                mTs, mYs, sTs, gammas, sYs,
                null.f1, null.f2,
                printevery, quiet
                )

    res$check = NULL
    res$ag = ag
    res$priag = priag
    res$dnpart = res$dnpart
    res$dalpha = res$dalpha
    if(in_sample){
        res$df1 = res$df1
        res$df2 = res$df2
    }
    res$dsigma1 = res$dsigma1
    res$dsigma2 = res$dsigma2
    res$dcov = res$dcov
    res$df1.test = res$df1.test
    res$df2.test = res$df2.test
    res$df1.varcount = res$df1.varcount
    res$df2.varcount = res$df2.varcount
    class(res) <- 'hetnpivbart'
    return(res)
}
