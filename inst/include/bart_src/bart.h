#ifndef GUARD_bart_h
#define GUARD_bart_h

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <string>

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"

class bart {
public:
   //------------------------------
   //friends
   friend bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen);
   //------------------------------
   //constructor/destructor
   bart();
  bart(size_t m,size_t size_beta);
   bart(const bart&);
   ~bart();
   //------------------------------
   //operators
   bart& operator=(const bart&);
   //------------------------------
   //get,set
   size_t getm() {return m;}
   void setm(size_t m);
  void setdata(size_t p,size_t n, double *x, double *y, double *w, size_t nc=100);
   void setpi(pinfo& pi) {this->pi = pi;}
   void setprior(double alpha, double beta, double tau) 
      {pi.alpha=alpha; pi.mybeta = beta; pi.tau=tau;} 
   void settau(double tau) {pi.tau=tau;}
  double getbeta0(size_t i) {return beta[0][i];}
  double getbeta1(size_t i) {return beta[1][i];}
   tree& gettree(size_t i ) { return t[i];}
   xinfo& getxinfo() {return xi;}
  std::vector<size_t>& getnv() {return nv;}
  std::vector<double>& getpv() {return pv;}
  void setpv(std::vector<double>& varprob) {
    for(size_t j=0;j<p;j++) pv[j]=varprob[j];
  }
   //------------------------------
   //public methods
  void birth(size_t i, size_t nid,size_t v, size_t c, std::vector<double> betal, std::vector<double> betar)
         {t[i].birth(nid,v,c,betal,betar);}
  void death(size_t i,size_t nid, std::vector<double> beta)
         {t[i].death(nid,beta);}
   void pr();
   void tonull() {for(size_t i=0;i!=t.size();i++) t[i].tonull();}
  void predict(size_t p, size_t n, double *x, double *w, double *fp);
  void predictb(size_t p, size_t n, double *xp, double *b0p, double *b1p);
   void draw(double sigma, rn& gen);
   double f(size_t i) {return allfit[i];}
protected:
  size_t sparse;
   size_t m;  //number of trees
  std::vector<tree> t; //the trees 
   pinfo pi; //prior and mcmc info
   //data
   size_t p,n; //x has dim p, n obserations
  double *x,*y, *w;  //x is column stack, pxn
   xinfo xi; //cutpoint info
   //working
   double *allfit; //if the data is set, should be f(x)
   double *r;
   double *ftemp;
  double lpp;
   dinfo di;
  size_t size_beta;
  std::vector < std::vector < double > > beta;
  std::vector<size_t> nv;
  std::vector<double> pv, lpv;
};

#include <iostream>
//#include "bart.h"

//--------------------------------------------------
//constructor
bart::bart():m(200),t(m),pi(),p(0),n(0),x(0),y(0),w(0),xi(),allfit(0),r(0),ftemp(0),lpp(0.0),di() {}
bart::bart(size_t im, size_t size_beta):m(im),t(m),pi(),p(0),n(0),x(0),y(0),w(0),xi(),allfit(0),r(0),ftemp(0),lpp(0.0),di(),size_beta(size_beta) {
  for(size_t h=0;h<m;h++){
    t[h].setbetasize(size_beta);
  }
}
bart::bart(const bart& ib):m(ib.m),t(m),pi(ib.pi),p(0),n(0),x(0),y(0),w(0),xi(),allfit(0),r(0),ftemp(0),di(),lpp(0.0)
{
   this->t = ib.t;
}
bart::~bart() 
{
   if(allfit) delete[] allfit;
   if(r) delete[] r;
   if(ftemp) delete[] ftemp;
}

//--------------------------------------------------
//operators
bart& bart::operator=(const bart& rhs)
{
   if(&rhs != this) {

      this->t = rhs.t;
      this->m = t.size();

      this->pi = rhs.pi;

      p=0;n=0;x=0;y=0;w=0;
      xi.clear();

      if(allfit) {delete[] allfit; allfit=0;}
      if(r) {delete[] r; r=0;}
      if(ftemp) {delete[] ftemp; ftemp=0;}

   }
   return *this;
}
//--------------------------------------------------
//get,set
void bart::setm(size_t m)
{
   t.resize(m);
   this->m = t.size();

   if(allfit && (xi.size()==p)) predict(p,n,x,w,allfit);
}
//--------------------------------------------------
void bart::setdata(size_t p,size_t n, double *x, double *y, double *w, size_t nc)
{
  this->p=p; this->n=n; this->x=x; this->y=y; this->w=w; 
  if(xi.size()==0) makexinfo(p,n,&x[0],xi,nc);
  if(allfit) delete[] allfit;
  allfit = new double[n];
  predict(p,n,x,w,allfit);
  if(r) delete[] r;
  r = new double[n];
  
  if(ftemp) delete[] ftemp;
  ftemp = new double[n*size_beta];

  this->beta.resize(size_beta);
  for(size_t i=0;i<size_beta;i++) {
    this->beta[i].resize(n);
    std::fill(this->beta[i].begin(),this->beta[i].end(),0.0);
  }
  
  di.n=n; di.p=p; di.x = &x[0]; di.y=r; di.w=w; di.size_beta=this->size_beta;
  for(size_t j=0;j<p;j++){
    nv.push_back(0);
    pv.push_back(1./(double)p);
    lpv.push_back(-log(p));
  }
}
//--------------------------------------------------
void bart::predict(size_t p, size_t n, double *x, double *w, double *fp)
//uses: m,t,xi
{
  double *fptemp = new double[n*size_beta];
   for(size_t j=0;j<n;j++) fp[j]=0.0;
   for(size_t j=0;j<m;j++) {
     fit(t[j],xi,p,n,x,size_beta,fptemp);
     if(size_beta==1) {
       for(size_t k=0;k<n;k++) {
	 fp[k] += fptemp[k];
       }
     }
     else if(size_beta==2) {
       for(size_t k=0;k<n;k++) {
	 fp[k] += fptemp[2*k]+w[k]*fptemp[2*k+1];
       }
     }
   }
   delete[] fptemp;
}
//--------------------------------------------------
void bart::predictb(size_t p, size_t n, double *xp, double *b0p, double *b1p) 
{
  double *fptemp = new double[n*size_beta];
  for(size_t j=0;j<n;j++) {
    b0p[j]=0.0;
    b1p[j]=0.0;
  }
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,xp,size_beta,fptemp);
    if(size_beta==1) {
      for(size_t k=0;k<n;k++) {	
	b0p[k] += fptemp[k];
       }
     }
     else if(size_beta==2) {
       for(size_t k=0;k<n;k++) {
	 b0p[k] += fptemp[2*k];
	 b1p[k] += fptemp[2*k+1];
       }
     }
  }
   delete[] fptemp;
}

//--------------------------------------------------
void bart::draw(double sigma, rn& gen)
{
  for(size_t l=0;l<size_beta;l++)
    std::fill(beta[l].begin(),beta[l].end(),0.0);
  lpp=0.0;
   for(size_t j=0;j<m;j++) {
     fit(t[j],xi,p,n,x,size_beta,ftemp);
     if(size_beta==1){
      for(size_t k=0;k<n;k++) {
         allfit[k] = allfit[k]-ftemp[k];
         r[k] = y[k]-allfit[k];
      }
     }
     else if(size_beta==2){
       for(size_t k=0;k<n;k++) {
         allfit[k] = allfit[k]-ftemp[2*k]-w[k]*ftemp[2*k+1];
         r[k] = y[k]-allfit[k];
       }
     }
     bd(t[j],xi,di,pi,sigma,nv,pv,gen);
     drmu(t[j],xi,di,pi,sigma,gen);
     fit(t[j],xi,p,n,x,size_beta,ftemp);
     if(size_beta==1){
       for(size_t k=0;k<n;k++) {
         allfit[k] += ftemp[k];
         beta[0][k] += ftemp[k];
       }
     }
     else if(size_beta==2){
       for(size_t k=0;k<n;k++) {
         allfit[k] += ftemp[2*k]+w[k]*ftemp[2*k+1];
         beta[0][k] += ftemp[2*k];
         beta[1][k] += ftemp[2*k+1];
       }
     }
     // Compute prior probability of current ensemble
   //   //   tree::npv tpv=t[j].getnodes();
   //   for(size_t b=0;b<tpv.size();b++){
   //     if(tpv[b].ntype()!='b'){
   // 	 lpp+=::log(pi.alpha/pow(1.+tpv[b].depth(),pi.mybeta));
   // 	 lpp+=-::log(p);
   // 	 lpp+=-::log(xi[tpv[b]->v].size());
   //     }
   //     else{
   // 	 lpp+=::log(1.-pi.alpha/pow(1.+tpv[b].depth(),pi.mybeta));
   // 	 double bsm2=0.0;
   // 	 for(size_t i=0;i<size_beta;i++){
   // 	   bsm2+=pow(tpv[b].beta[i],2.0);
   // 	 }
   // 	 lpp+=-::log(tau)-(0.5*bsm2/(pow(pi.tau,2.0));
   //     }
   //   }
   }
}
//--------------------------------------------------
//public functions
void bart::pr() //print to screen
{
   cout << "*****bart object:\n";
   cout << "m: " << m << std::endl;
   cout << "t[0]:\n " << t[0] << std::endl;
   cout << "t[m-1]:\n " << t[m-1] << std::endl;
   cout << "prior and mcmc info:\n";
   pi.pr();
   if(p) cout << "data set: n,p: " << n << ", " << p << std::endl;
   else cout << "data not set\n";
}

#endif
