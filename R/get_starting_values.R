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

## Gets 2SLS coefficients with data

get_starting_values <- function(betas,sTs,sYs,gammas,L,null2){
    starts <- rep(0,4)
    if(null2){
        starts[1] <- 0
    }
    if(is.na(sTs)) starts[2] <- L[1,1]
    if(is.na(sYs)) starts[3] <- L[2,2]
    if(is.na(gammas)) starts[4] <- L[2,1]
    return(starts)
}
