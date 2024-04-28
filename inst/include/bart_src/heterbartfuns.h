#ifndef GUARD_heterbartfuns_h
#define GUARD_heterbartfuns_h

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "rn.h"
#include <map>


//--------------------------------------------------
//heterlh, replacement for lil that only depends on sum y.
double heterlh(double b, double M, double tau);
//--------------------------------------------------
//compute b and M  for left and right give bot and v,c
void hetergetsuff(tree<double>& x, typename tree<double>::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& srl, double& sryl, size_t& nr, double& srr, double& sryr, double *sigma);
//--------------------------------------------------
//compute b and M for left and right bots
void hetergetsuff(tree<double>& x, typename tree<double>::tree_p left, typename tree<double>::tree_p right, xinfo& xi, dinfo& di, double& srl, double& sryl, double& srr, double& sryr, double *sigma);
//compute sufficient statistics for left and right give bot and v,c
void hetergetsuff(tree<std::vector<double>>& x, typename tree<std::vector<double>>::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& srl, double& sryl, double& srwl, double& srwyl, double& srwwl, size_t& nr,  double& srr, double& sryr, double& srwr, double& srwyr, double& srwwr, double *sigma);
//compute sufficient statistics for left and right bots
void hetergetsuff(tree<std::vector<double>>& x, typename tree<std::vector<double>>::tree_p left, typename tree<std::vector<double>>::tree_p right, xinfo& xi, dinfo& di, double& srl, double& sryl, double& srwl, double& srwyl, double& srwwl, double& srr, double& sryr, double& srwr, double& srwyr, double& srwwr, double *sigma);
//--------------------------------------------------
//draw one mu from post
double heterdrawnodetheta(double sr, double sry, double tau, rn& gen);
//--------------------------------------------------
// draw two mu from post (Linear)
std::vector<double> heterdrawnodetheta(double sr, double sry, double srw, double srwy, double srww, double tau, double tau2, rn& gen);
//--------------------------------------------------
//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void heterallsuff(tree<double>& x, xinfo& xi, dinfo& di, tree<double>::npv& bnv, std::vector<double>& srv, std::vector<double>& sryv,double *sigma);
//--------------------------------------------------
//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void heterallsuff(tree<std::vector<double>>& x, xinfo& xi, dinfo& di, tree<std::vector<double>>::npv& bnv, std::vector<double>& srv, std::vector<double>& sryv, std::vector<double>& srwv, std::vector<double>& srwyv, std::vector<double>& srwwv, double *sigma);
//--------------------------------------------------
//heter version of drtheta, need b and M instead of n and sy
void heterdrtheta(tree<double>& t, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, rn& gen);
//--------------------------------------------------
//heter version of drtheta, need b and M instead of n and sy
void heterdrtheta(tree<std::vector<double>>& t, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, rn& gen);

//#include "heterbartfuns.h"

//--------------------------------------------------
//heterlh, replacement for lil that only depends on sum y.
double heterlh(double b, double M, double tau) {
   double t2 = tau*tau;
   double k = b*t2+1;
   return -.5*log(k)+.5*M*M*t2/k;
}
//--------------------------------------------------
//compute b and M  for left and right give bot and v,c
void hetergetsuff(tree<double>& x, tree<double>::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& srl, double& sryl, size_t& nr,  double& srr, double& sryr, double *sigma)
{
   double *xx;//current x
   srl=0; srr=0.0; sryl=0; sryr=0.0; nl=0; nr=0;
   double r;
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node
         r= 1.0/(sigma[i]*sigma[i]);
         if(xx[v] < xi[v][c]) {
               nl+=1;
               srl+=r;
               sryl+=r*di.y[i];
          } else {
               nr+=1;
               srr+=r;
               sryr+=r*di.y[i];
	 }
      }
   }
}
//--------------------------------------------------
//compute b and M for left and right bots
void hetergetsuff(tree<double>& x, tree<double>::tree_p left, tree<double>::tree_p right, xinfo& xi, dinfo& di, double& srl, double& sryl, double& srr, double& sryr, double *sigma)
{

   double *xx;//current x
   srl=0; sryl=0.0; srr=0; sryr=0.0;
   double r;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tree<double>::tree_cp bn = x.bn(xx,xi);
      if(bn==left) {
         r = 1.0/(sigma[i]*sigma[i]);
         srl+=r;
         sryl += r*di.y[i];
      }
      if(bn==right) {
         r = 1.0/(sigma[i]*sigma[i]);
         srr+=r;
         sryr += r*di.y[i];
      }
   }
}
//compute sufficient statistics for left and right give bot and v,c
void hetergetsuff(tree<std::vector<double>>& x, tree<std::vector<double>>::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& srl, double& sryl, double& srwl, double& srwyl, double& srwwl, size_t& nr,  double& srr, double& sryr, double& srwr, double& srwyr, double& srwwr, double *sigma)
{
   double *xx;//current x
   nl=0;srl=0.0;sryl=0.0;srwl=0.0;srwyl=0.0;srwwl=0.0;
   nr=0;srr=0.0;sryr=0.0;srwr=0.0;srwyr=0.0;srwwr=0.0;
   
   double r;
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node
         r = 1.0/(sigma[i]*sigma[i]);
         if(xx[v] < xi[v][c]) {
               nl+=1;
               srl+=r;
               sryl += r*di.y[i];
               srwl += r*di.w[i];
               srwyl += r*di.w[i]*di.y[i];
               srwwl += r*di.w[i]*di.w[i];
         }
         else {
               nr+=1;
               srr+=r; 
               sryr += r*di.y[i];
               srwr += r*di.w[i];
               srwyr += r*di.w[i]*di.y[i];
               srwwr += r*di.w[i]*di.w[i];
         }
      }
   }
}
//compute sufficient statistics for left and right bots
void hetergetsuff(tree<std::vector<double>>& x, typename tree<std::vector<double>>::tree_p left, typename tree<std::vector<double>>::tree_p right, xinfo& xi, dinfo& di, double& srl, double& sryl, double& srwl, double& srwyl, double& srwwl, double& srr, double& sryr, double& srwr, double& srwyr, double& srwwr, double *sigma)
{
  
  double *xx;//current x
   srl=0.0; sryl=0.0; srwl=0.0; srwyl=0.0; srwwl=0.0;
   srr=0.0; sryr=0.0; srwr=0.0; srwyr=0.0; srwwr=0.0;
   double r;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      typename tree<std::vector<double>>::tree_cp bn = x.bn(xx,xi);
      if(bn==left) {
        r = 1.0/(sigma[i]*sigma[i]);
        srl+=r;
        sryl += r*di.y[i];
        srwl += r*di.w[i];
        srwyl += r*di.w[i]*di.y[i];
        srwwl += r*di.w[i]*di.w[i];
      }
      if(bn==right) {
        r = 1.0/(sigma[i]*sigma[i]);
        srr+=r;
        sryr += r*di.y[i];
        srwr += r*di.w[i];
        srwyr += r*di.w[i]*di.y[i];
        srwwr += r*di.w[i]*di.w[i];
      }
   }
}
//--------------------------------------------------
//draw one mu from post
double heterdrawnodetheta(double sr, double sry, double tau, rn& gen)
{
   double muhat = sry/sr;
   double a = 1.0/(tau*tau);
   return (sr*muhat)/(a+sr) + gen.normal()/sqrt(a+sr);
}
//--------------------------------------------------
// draw two mu from post (Linear)
std::vector<double> heterdrawnodetheta(double sr, double sry, double srw, double srwy, double srww, double tau, double tau2, rn& gen)
{
  double t2 = tau*tau;
  double inv_V00 = sr + 1./t2;
  double inv_V01 = srw;
  double inv_V11 = srww + 1./(tau2*tau2);
  double det_invV = inv_V00*inv_V11-inv_V01*inv_V01;
  double V00 = (1./det_invV)*inv_V11;
  double V11 = (1./det_invV)*inv_V00;
  double V01 = -(1./det_invV)*inv_V01;
  double wty0 = sry;
  double wty1 = srwy;
  double m0 = wty0*V00+wty1*V01;
  double m1 = wty0*V01+wty1*V11;
  double chol_V00 = sqrt(V00);
  double chol_V01 = V01/chol_V00;
  double chol_V11 = sqrt(V11-chol_V01*chol_V01);
  std::vector<double> out(2,0.);
  double Z0 = gen.normal();
  double Z1 = gen.normal();
  // cout << "M0: " << m0 << '\n';
  // cout << "M1: " << m1 << '\n';
  // cout << "V00: " << V00 << '\n';
  // cout << "V01: " << V01 << '\n';
  // cout << "V11: " << V11 << '\n';
  out[0] = chol_V00*Z0 + m0;
  out[1] = chol_V01*Z0 + chol_V11*Z1 + m1;
  return out;
}
//--------------------------------------------------
//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void heterallsuff(tree<double>& x, xinfo& xi, dinfo& di, tree<double>::npv& bnv, std::vector<double>& srv, std::vector<double>& sryv, double *sigma)
{
   tree<double>::tree_cp tbn; //the pointer to the bottom node for the current observations
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x

   bnv.clear();
   x.getbots(bnv);

   typedef tree<double>::npv::size_type bvsz;
   bvsz nb = bnv.size();
   srv.resize(nb);
   sryv.resize(nb);

   std::map<tree<double>::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) {bnmap[bnv[i]]=i;srv[i]=0;sryv[i]=0.0;}

   double r;
   for(size_t i=0;i<di.n;i++) {
      r = 1.0/(sigma[i]*sigma[i]);
      xx = di.x + i*di.p;
      tbn = x.bn(xx,xi);
      ni = bnmap[tbn];

      srv[ni] += r;
      sryv[ni] += r*di.y[i];
   }
}
//--------------------------------------------------
//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void heterallsuff(tree<std::vector<double>>& x, xinfo& xi, dinfo& di, tree<std::vector<double>>::npv& bnv, std::vector<double>& srv, std::vector<double>& sryv, std::vector<double>& srwv, std::vector<double>& srwyv, std::vector<double>& srwwv, double *sigma)
{
  tree<std::vector<double>>::tree_cp tbn; //the pointer to the bottom node for the current observations
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x

   bnv.clear();
   x.getbots(bnv);

   typedef tree<std::vector<double>>::npv::size_type bvsz;
   bvsz nb = bnv.size();
   srv.resize(nb);
   sryv.resize(nb);
   srwv.resize(nb);
   srwyv.resize(nb);
   srwwv.resize(nb);

   std::map<tree<std::vector<double>>::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) {
     bnmap[bnv[i]]=i;
     srv[i]=0.0;
     sryv[i]=0.0;
     srwv[i]=0.0;
     srwyv[i]=0.0;
     srwwv[i]=0.0;
   }

   double r;
   for(size_t i=0;i<di.n;i++) {
      r = 1.0/(sigma[i]*sigma[i]);
      xx = di.x + i*di.p;
      tbn = x.bn(xx,xi);
      ni = bnmap[tbn];

      srv[ni] += r;
      sryv[ni] += r*di.y[i];
      srwv[ni] += r*di.w[i];
      srwyv[ni] += r*di.w[i]*di.y[i];
      srwwv[ni] += r*di.w[i]*di.w[i];
   }
}
//--------------------------------------------------
//heter version of drtheta, need b and M instead of n and sy
void heterdrtheta(tree<double>& t, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, rn& gen)
{
   tree<double>::npv bnv;
   std::vector<double> srv;
   std::vector<double> sryv;
   heterallsuff(t,xi,di,bnv,srv,sryv,sigma);
   for(tree<double>::npv::size_type i=0;i!=bnv.size();i++)
     bnv[i]->settheta(heterdrawnodetheta(srv[i],sryv[i],pi.tau,gen));
}
//--------------------------------------------------
//heter version of drtheta, need b and M instead of n and sy
void heterdrtheta(tree<std::vector<double>>& t, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, rn& gen)
{
  tree<std::vector<double>>::npv bnv;
  std::vector<double> srv;
  std::vector<double> sryv;
  std::vector<double> srwv;
  std::vector<double> srwyv;
  std::vector<double> srwwv;
  heterallsuff(t,xi,di,bnv,srv,sryv,srwv,srwyv,srwwv,sigma);
  for(tree<std::vector<double>>::npv::size_type i=0;i!=bnv.size();i++)
    bnv[i]->settheta(heterdrawnodetheta(srv[i],sryv[i],srwv[i],srwyv[i],srwwv[i],pi.tau,pi.tau2,gen));
}

#endif
