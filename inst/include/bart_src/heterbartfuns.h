#ifndef GUARD_heterbartfuns_h
#define GUARD_heterbartfuns_h

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "rn.h"
#include <map>


//--------------------------------------------------
//heterlh, replacement for lil that only depends on sum y.
double heterlh(double sr, double sry, double tau);
//--------------------------------------------------
//compute sufficient statistics for left and right give bot and v,c
void hetergetsuff(tree& x, tree::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& srl, double& sryl, double& srwl, double& srwyl, double& srwwl, size_t& nr,  double& srr, double& sryr, double& srwr, double& srwyr, double& srwwr, double *sigma);
//--------------------------------------------------
//compute sufficient statistics for left and right bots
void hetergetsuff(tree& x, tree::tree_p left, tree::tree_p right, xinfo& xi, dinfo& di, double& srl, double& sryl, double& srwl, double& srwyl, double& srwwl, double& srr, double& sryr, double& srwr, double& srwyr, double& srwwr, double *sigma);
//--------------------------------------------------
//draw one mu from post (Constant)
double heterdrawnodemu0(double sr, double sry, double tau, rn& gen);
//draw two mu from post (Linear)
std::vector<double> heterdrawnodemu1(double sr, double sry, double srw, double srwy, double srww, double tau, rn& gen);
//--------------------------------------------------
//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void heterallsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<double>& srv, std::vector<double>& sryv, std::vector<double>& srwv, std::vector<double>& srwyv, std::vector<double>& srwwv, double *sigma);
//--------------------------------------------------
//heter version of drmu, need b and M instead of n and sy
void heterdrmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, rn& gen);

//#include "heterbartfuns.h"

//--------------------------------------------------
//heterlh, replacement for lil that only depends on sum y.
double heterlh(double sr, double sry, double tau) {
   double t2 =tau*tau;
   double k = sr*t2+1;
   return -.5*log(k)+.5*sry*sry*t2/k;
}
//--------------------------------------------------
//compute b and M  for left and right give bot and v,c
void hetergetsuff(tree& x, tree::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& srl, double& sryl, double& srwl, double& srwyl, double& srwwl, size_t& nr,  double& srr, double& sryr, double& srwr, double& srwyr, double& srwwr, double *sigma)
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
//--------------------------------------------------
//compute b and M for left and right bots
void hetergetsuff(tree& x, tree::tree_p left, tree::tree_p right, xinfo& xi, dinfo& di, double& srl, double& sryl, double& srwl, double& srwyl, double& srwwl, double& srr, double& sryr, double& srwr, double& srwyr, double& srwwr, double *sigma)
{

  double *xx;//current x
   srl=0.0; sryl=0.0; srwl=0.0; srwyl=0.0; srwwl=0.0;
   srr=0.0; sryr=0.0; srwr=0.0; srwyr=0.0; srwwr=0.0;
   double r;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tree::tree_cp bn = x.bn(xx,xi);
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
//draw one mu from post (Constant)
double heterdrawnodemu0(double sr, double sry, double tau, rn& gen)
{
  double muhat = sry/sr;
  double a = 1.0/(tau*tau);
  double v = 1.0/(a+sr);
  double m = v*sry;
  return(m + sqrt(v)*gen.normal());
}
//--------------------------------------------------
// draw two mu from post (Linear)
std::vector<double> heterdrawnodemu1(double sr, double sry, double srw, double srwy, double srww, double tau, rn& gen)
{
  std::vector<double> _out (2,0.0);
  double t2 = tau*tau;
  double vinv11=sr+1.0/t2; double vinv12=srw; double vinv22=srww+1.0/t2;
  double dtrmntinv = 1.0/(vinv11*vinv22-vinv12*vinv12);
  double vbeta11=dtrmntinv*vinv22;
  double vbeta12=-dtrmntinv*vinv12;
  double vbeta22=dtrmntinv*vinv11;
  double zscore1=gen.normal();
  double zscore2=gen.normal();
  double sdbeta0=sqrt(vbeta11);
  double cholfactor=vbeta12/sdbeta0;
  double sdbeta1=sqrt(vbeta22-cholfactor*cholfactor);
  _out[0] = sdbeta0*zscore1;
  _out[1] = zscore1*cholfactor+sdbeta1*zscore2;
  _out[0] = _out[0]+vbeta11*sry+vbeta12*srwy;
  _out[1] = _out[1]+vbeta12*sry+vbeta22*srwy;
  return _out;
}
//--------------------------------------------------
//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void heterallsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<double>& srv, std::vector<double>& sryv, std::vector<double>& srwv, std::vector<double>& srwyv, std::vector<double>& srwwv, double *sigma)
{
   tree::tree_cp tbn; //the pointer to the bottom node for the current observations
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x

   bnv.clear();
   x.getbots(bnv);

   typedef tree::npv::size_type bvsz;
   bvsz nb = bnv.size();
   srv.resize(nb);
   sryv.resize(nb);
   srwv.resize(nb);
   srwyv.resize(nb);
   srwwv.resize(nb);

   std::map<tree::tree_cp,size_t> bnmap;
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
//heter version of drmu, need b and M instead of n and sy
void heterdrmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, rn& gen)
{
   tree::npv bnv;
   std::vector<double> srv;
   std::vector<double> sryv;
   std::vector<double> srwv;
   std::vector<double> srwyv;
   std::vector<double> srwwv;
   heterallsuff(t,xi,di,bnv,srv,sryv,srwv,srwyv,srwwv,sigma);

   if(di.size_beta==1){
     for(tree::npv::size_type i=0;i!=bnv.size();i++){
       bnv[i]->setbeta(0,heterdrawnodemu0(srv[i],sryv[i],pi.tau,gen));
     }
   }
   else if(di.size_beta==2){
     for(tree::npv::size_type i=0;i!=bnv.size();i++){
       std::vector<double> tmpout;
       tmpout=heterdrawnodemu1(srv[i],sryv[i],srwv[i],srwyv[i],srwwv[i],pi.tau,gen);
       bnv[i]->setbeta(0,tmpout[0]);
       bnv[i]->setbeta(1,tmpout[1]);
     }
   }
   else {}
}

     

#endif
