
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

hetivbart = function(Z, T, Y, X=NULL,
                     Z.test=NULL,
                     T.test=NULL, 
                     X.test=NULL,
                     burn=1000, nd=2000, burnf1=1000, burnf2=1000,
                     keepevery=10,
                     m1=200, m2=200, nc=100,
                     power=2, base=0.95,
                     k=2, sigmaf1=NA, sigmaf2=NA,
                     doDP=TRUE, do.tsls=TRUE,
                     v=0.17, nu=2.004, a=0.016,
                     Imin=1, Imax=floor(0.1*n)+1, psi=0.5, gs=100,
                     centermeans=TRUE,
                     n=length(Y),
                     f1s = rep(0, n), f2s = rep(0, n),
                     mTs=0, mYs=0, sTs=NA, sYs=NA, gammas=NA,
                     null.f1=FALSE, null.f2=FALSE,
                     printevery=100, quiet=TRUE
                     )
{


    ## Data checking
    X <- check_build_training_data(Z,T,X,n,"ivbart")
    ZX.ls <- check_build_test_data(Z.test,T.test,X.test,"ivbart")
    Z.test <- ZX.ls[[1]]; X.test <- ZX.ls[[2]]
    rm(ZX.ls)

    ## Getting tau values (prior variance of terminal node value)
    tau.ls <- get_tau(sigmaf1,sigmaf2,k,m1,m2,max(T)-min(T),max(Y)-min(Y))
    tauf1 <- tau.ls[1]; tauf2 <- tau.ls[2]
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

    ## Starting values
    starts.ls <- get_starting_values(0,sTs,sYs,gammas,L,null.f2)
    sTs=starts.ls[2]; sYs=starts.ls[3]; gammas=starts.ls[4]
    rm(starts.ls)

    include_output <- 0
    
    res = .Call("cbhetiv",
                t(Z),
                t(Z.test),
                t(X),
                t(X.test),
                T,
                Y,
                burn,
                nd*keepevery,
                burnf1,
                burnf2,
                m1, m2, nc,
                power, base,
                tauf1,
                tauf2,
                doDP,
                v, nu, a, #base prior
                ag, priag, #alpha prior
                centermeans,
                include_output,
                f1s, f2s,
                mTs, mYs, sTs, gammas, sYs,
                null.f1, null.f2,
                printevery, quiet
                )

    res$check = NULL
    thin <- seq(1,nd*keepevery,keepevery)
    res$dnpart = res$dnpart[thin]
    res$dalpha = res$dalpha[thin]
    if(include_output==1){
        res$dsigma1 = res$dsigma1[thin, ]
        res$dsigma2 = res$dsigma2[thin, ]
        res$dcov = res$dcov[thin, ]
        res$df1 = res$df1[thin, ]
        res$df21 = res$dbeta0[thin,]
        res$df22 = res$dbeta1[thin,]
        res$df2 = res$df21+matrix(rep(T,each=nd),nrow=nd)*res$df22
        res$dLL = res$dLL[thin, ]
    }
    if(NROW(Z.test)!=0){
        res$df1.test = res$df1.test[thin, ]
    }
    if(NROW(X.test)!=0){
        res$df21.test = res$dbeta0.test[thin,]
        res$df22.test = res$dbeta1.test[thin,]
        if(length(T.test)!=0)
            res$df2.test = res$df21.test+matrix(rep(T.test,each=nd),nrow=nd)*res$df22.test
    }
    class(res) <- 'hetivbart'
    return(res)
}
