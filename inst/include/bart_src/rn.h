/*
 *  hetnpivbart: Heterogenous Nonparametric IV Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Charles Spanbauer
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#ifndef RN_H
#define RN_H

#define RTPI 2.506628274631000502415765284811

class rn
{
 public:
 rn():df(1.),shape(1.),scale(1.) {}
  ~rn() {}
  double normal() {return R::rnorm(0.,1.);}
  double uniform() {return R::runif(0.,1.);}
  double chi_square() {return R::rchisq(df);}
  double gamma() {return R::rgamma(shape,scale);}
  double exp() {return R::rgamma(1.,scale);}
  void set_df(double _df) {this->df=_df;}
  void set_alpha(double _shape) {this->shape=_shape;}
  void set_scale(double _scale) {this->scale=_scale;}
  
 private:
  double df, shape, scale;
};
    
#endif  
