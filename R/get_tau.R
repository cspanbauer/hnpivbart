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

## Gets tau values

get_tau <- function(sigmaf1,sigmaf2,k1,k2,m1,m2,rT,rY){
    if(is.na(sigmaf1)) {
        tauf1 = rT/(2*k1*sqrt(m1));
    } else {
        tauf1 = sigmaf1/sqrt(m1)
    }
    if(is.na(sigmaf2)) {
        tauf2 = rY/(2*k2*sqrt(m2));
    } else {
        tauf2 = sigmaf2/sqrt(m2)
    }
    return(c(tauf1,tauf2))
}
