# ivbart: Instrumental Variable BART with Diricihlet Process Mixtures
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

ivbart2s = function(Z1, Z2, X1, X2, T1, Y, X.test=NULL, Z.test=NULL,
   burn=1000, nd=2000, burnf=1000, burnh=1000,
   keepevery=10,
   betabar=0, Abeta=NA,
   m=200, nc=100,
   power=2, base=0.95,
   k=2, sigmaf=NA, sigmah=NA,
   v=0.17, nu=2.004, a=0.016,
   n1 = length(T1), n2 = length(Y),
   Imin1=1, Imax1=floor(0.1*n1)+1, psi=0.5, gs=100,
   Imin2=1, Imax2=floor(0.1*n2)+1,
   centermeans=TRUE,
   fs = rep(0, n2), hs = rep(0, n2),
   mTs=0, mYs=0,
   betas=NA, sTs=NA, sYs=NA, gammas=NA, T2s=rep(mean(T1),nrow(Z2)),
## betas=0, mTs=0, mYs=0, sTs=1, sYs=1, gammas=0,
   printevery=100
)
{
### data
    if(length(T1)!=nrow(Z1)) stop('The length of T1 must equal the number of rows of Z1')
    if(length(X1)==1) X = rep(X1, n1)
    if(length(X2)==1) X = rep(X2, n2)
    X1 = cbind(X1)
    X2 = cbind(X2)    
    Z1 = cbind(Z1)
    Z2 = cbind(Z2)
    ZX1 = cbind(Z1, X1)
    ZX2 = cbind(Z2, X2)
    ZX12 = rbind(ZX1,ZX2)
if(nrow(Z2)!=n2) stop('The rows of Z2 must equal the length of Y')
    if(nrow(X2)!=n2) stop('Either X is a constant (no measured confounders)\n or the rows of X must equal the length of Y')
    if(ncol(Z1)!=ncol(Z2)) stop('The number of columns of Z1 must be the same as Z2')
    if(ncol(X1)!=ncol(X2)) stop('The number of columns of X1 must be the same as X2')
    if(length(X.test)==0) X.test <- matrix(0,nrow=0,ncol=0)
    if(length(Z.test)==0) Z.test <- matrix(0,nrow=0,ncol=0)
    ZX.test <- cbind(Z.test,X.test)
    
## tau
if(is.na(sigmaf)) {
   tauf = (max(T)-min(T))/(2*k*sqrt(m));
} else {
   tauf = sigmaf/sqrt(m)
}
if(is.na(sigmah)) {
   tauh = (max(Y)-min(Y))/(2*k*sqrt(m));
} else {
   tauh = sigmah/sqrt(m)
}

### prior for alpha
priold = getaprior(length(T1),Imin1,Imax1,psi,gs)
priag = priold$p
ag = priold$a

    prioldB = getaprior(length(Y),Imin2,Imax2,psi,gs)
    priagB = prioldB$p
    agB = prioldB$a
    
### unless provided starting values from TSLS
## first stage manually
    tsls1 = lm(T~., data.frame(T = T1, Z1, X1))
    print(dim(cbind(1,ZX2)))
    print(matrix(tsls1$coefficients,ncol=1))
    Ts <- cbind(1,ZX2)%*%matrix(tsls1$coefficients,ncol=1)
## second stage manually
    tsls2 = lm(Y~T+., data.frame(Y = Y, T = Ts, X2))
L = t(chol(var(cbind(tsls1$resid, tsls2$resid))))
    if(is.na(betas)) betas=tsls2$coeff[2]
    
if(is.na(sTs)) sTs=L[1, 1]
if(is.na(sYs)) sYs=L[2, 2]
if(is.na(gammas)) gammas=L[2, 1]
if(is.na(Abeta)) Abeta=1/(sqrt(m)*tauh/(max(T)-min(T)))
    
res = .Call("cbiv2s",
            t(ZX1),
            t(ZX2),
            t(ZX12),
            t(ZX.test),
            t(X1),
            t(X2),
   t(X.test),
   T,
   Y,
   burn,
   nd*keepevery,
   burnf,
   burnh,
   m, nc,
   power, base,
   tauf,
   tauh,
   betabar, Abeta,
   v, nu, a, #base prior
   ag, priag, #alpha prior
   agB, priagB,
   centermeans,
   fs, hs,
   betas, mTs, mYs, sTs, gammas, sYs, T2s,
   printevery
)
res$check = NULL
#res$ag = ag
#res$priag = priag
thin = seq(1,nd*keepevery,keepevery)
res$dnpart = res$dnpart[thin]
res$dalpha = res$dalpha[thin]
res$dmu1 = res$dmu1[thin, ]
res$dsigma1 = res$dsigma1[thin, ]
res$dmu2 = res$dmu2[thin, ]
    res$dsigma2 = res$dsigma2[thin, ]
    res$dgamma = res$dgamma[thin, ]
res$dcov = res$dcov[thin, ]
    res$dcor = res$dcor[thin, ]
    res$dsY = res$dsY[thin, ]
res$dbeta = res$dbeta[thin]
res$dh = res$dh[thin, ]
    res$df = res$df[thin, ]
    res$df.test = res$df.test[thin, ]
res$dh.test = res$dh.test[thin,]
res$betas = betas
res$sTs = sTs
res$sYs = sYs
res$gammas = gammas
res$tsls1 = tsls1
    res$tsls2 = tsls2
    res$tau <- list(tauf,tauh)
    res$prior = list(res$dnu,res$da,res$dv)
return(res)
}
