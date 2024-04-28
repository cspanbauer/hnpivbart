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

# Checks if input test data has correct dimensions.

check_build_test_data <- function(Z.test,T.test,X.test,method){
    if(method=="hetnpivbart"){
        if(length(Z.test)==0) Z.test = matrix(0,nrow=0,ncol=0)
        if(length(X.test)==0&length(T.test)==0)
        {
            X.test = matrix(0,nrow=0,ncol=0)
            TX.test = X.test
        }
        else if(length(X.test)==0&length(T.test)!=0)
        {
            X.test = matrix(0,nrow=0,ncol=0)
            TX.test = T.test
        }
        else
        {
            if(nrow(X.test)!=length(T.test)) stop('The number of rows of X.test must equal the number of rows of T.test')
            else {
                TX.test = cbind(T.test,X.test)
            }
        }        
        return(list(Z.test,TX.test))
    }
    else{
        if(length(X.test)==0) X.test = matrix(0,nrow=0,ncol=0)
        if(length(Z.test)==0) Z.test = matrix(0,nrow=0,ncol=0)
        if(length(T.test)==0) T.test = matrix(0,nrow=0,ncol=0)
        return(list(Z.test,X.test,T.test))
    }
}
