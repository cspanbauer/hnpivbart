
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>

using std::endl;
using std::cout;

#include <Rcpp.h>

#include <bart_src/rn.h>
#include <bart_src/rtnorm.h>
#include <bart_src/tree.h>
#include <bart_src/treefuns.h>
#include <bart_src/info.h>
#include <bart_src/bartfuns.h>
#include <bart_src/bd.h>
#include <bart_src/bart.h>
#include <bart_src/heterbart.h>
#include <bart_src/dart.h>
//#include <bart_src/semibart.h>
//#include <bart_src/semibd.h>
//#include <bart_src/semibartfuns.h>
//#include <bart_src/hetersemibart.h>


#include <bart_src/Dp.h>
#include <bart_src/DpMuSigma.h>
#include <bart_src/DpMuTau.h>
//#include <bart_src/DpMuSigmaLKJ.h>

#define LTPI 1.837877066409345483560659472811
#define PI   3.141592653589793238462643383279

double lreg(int n, double *x, double *y, double sigma, double betabar, double Abeta, rn& gen);

double lreg(int n, double *x, double *y, double sigma, double betabar, double Abeta, rn& gen)
{
   double sxx=0.0,sxy=0.0;
   for(int i=0;i<n;i++) {
      sxx += x[i]*x[i];
      sxy += x[i]*y[i];
   }
   double sig2 = sigma*sigma;
   double v = 1.0/((sxx/sig2) + Abeta);
   double m = v*((sxy/sig2) + Abeta*betabar);

   return (m + sqrt(v)*gen.normal());
}

double bvnorm(double x, double y, double mux, double muy, double sigx, double sigy, double sigxy);

double bvnorm(double x, double y, double mux, double muy, double sigx, double sigy, double sigxy)
{
  double rho = sigxy/(sigx*sigy);
  double T1MR = 2.0*(1.0-rho*rho);
  return -log(2.)-log(PI)-log(sigx)-log(sigy)-0.5*log(1-rho*rho)-(x-mux)*(x-mux)/(T1MR*sigx*sigx)+2*rho*(x-mux)*(y-muy)/(T1MR*sigx*sigy)-(y-muy)*(y-muy)/(T1MR*sigy*sigy);
}

double log_rho_posterior(double r, double sumY1sq, double sumY2sq, double sumY1Y2, double s, double a, double b, size_t n){
  return -0.5*((double)n)*(log(s*s)+log(1.-r*r))-0.5*(sumY1sq/(s*s)-2.*r*sumY1Y2/s+sumY2sq)/(1.-r*r) + (a-1)*log(1.+r) + (b-1)*log(1.-r);
}

double log_sig_posterior(double s, double sumY1sq, double sumY2sq, double sumY1Y2, double r, double m, double v, size_t n){
  return -0.5*((double)n)*(log(s*s)) - 0.5*(sumY1sq/(s*s)-2*r*sumY1Y2/s)/(1.-r*r) - 0.5*(v+1.)*log(1+s*s/v*m);
}


