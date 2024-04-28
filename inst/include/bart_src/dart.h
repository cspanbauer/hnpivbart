/*
 *  hetnpivbart: Heterogenous Nonparametric IV Bayesian Additive Regression Trees
 *  Copyright (C) 2023 Charles Spanbauer
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

double draw_concentration(double *lpv, size_t p)
{
  double theta=1./log((double)p);
  return theta;
}

void draw_log_dirichlet(size_t *vc, double *lpv, double theta, size_t p, rn gen)
{
  gen.set_alpha(theta+(double)vc[0]);
  lpv[0]=log(gen.gamma());
  double mx_lpv=lpv[0];
  for(size_t j=1;j<p;j++){
    gen.set_alpha(theta+(double)vc[j]);
    lpv[j]=log(gen.gamma());
    if(lpv[j]>mx_lpv) mx_lpv=lpv[j];
  }
  // first get sum of gamma
  double log_sum_gamma=0.;
  for(size_t j=0;j<p;j++){
    log_sum_gamma+=exp(lpv[j]-mx_lpv);
  }
  // log-sum-exp trick
  log_sum_gamma=mx_lpv+log(log_sum_gamma);
  for(size_t j=0;j<p;j++){
    lpv[j]=lpv[j]-log_sum_gamma;
  }
}
