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

## Checks if input training data has correct dimensions.
## if X is NULL, then TX is equal to T for npivBART-h
## if X is NULL, then it is translated to a 0x0 matrix to pass via Rcpp for other methods
check_build_training_data <- function(Z,T,X,n,method){
    if(method=="hetnpivbart"){        
        if(length(T)!=n) stop('The length of T must equal the length of Y')
        if(NROW(Z)!=n) stop('The number of rows of Z must equal the length of Y')
        if(length(X)==0)
        {
            X = matrix(0,nrow=0,ncol=0)
            TX = T
        }
        else if(length(X)!=0)
        {
            if(NROW(X)!=n) stop('The number of rows of X must equal the length of Y')
            else TX = cbind(T, X)
        }
        return(TX)
    }
    else{
        if(length(T)!=n) stop('The length of T must equal the length of Y')
        if(NROW(Z)!=n) stop('The number of rows of Z must equal the length of Y')
        if(length(X)==0)
        {
            X = matrix(0,nrow=0,ncol=0)
        }
        else if(length(X)!=0)
        {
            if(nrow(X)!=n) stop('The number of rows of X must equal the length of Y')
        }    
        return(X)
    }
}
