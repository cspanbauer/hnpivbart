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

template <class T>
class bart {
public:
  //------------------------------
  //friends
  friend bool bd(tree<double>& x, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen);
  friend bool bd(tree<std::vector<double>>& x, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen);
  //------------------------------
  //constructor/destructor
  bart();
  bart(size_t m);
  bart(const bart&);
  ~bart();
  //------------------------------
  //operators
  bart& operator=(const bart&);
  //------------------------------
  //get,set
  size_t getm() {return m;}
  void setm(size_t m);
  void setdata(size_t p,size_t n, double *x, double *y, size_t nc=100);
  void setdata(size_t p,size_t n, double *x, double *y, double *w, size_t nc=100);
  void setpi(pinfo& pi) {this->pi = pi;}
  void setprior(double alpha, double beta, double tau, double *pv, double tau2=1.0) 
  {pi.alpha=alpha; pi.mybeta = beta; pi.tau=tau; pi.tau2=tau2; pi.pv=pv;}
  void settau(double tau) {pi.tau=tau;}
  tree<T>& gettree(size_t i ) { return t[i];}
  xinfo& getxinfo() {return xi;}
  //------------------------------
  //public methods
  void birth(size_t i, size_t nid,size_t v, size_t c, double ml, double mr)
  {t[i].birth(nid,v,c,ml,mr);}
  void death(size_t i,size_t nid, double mu)
  {t[i].death(nid,mu);}
  void pr();
  void tonull() {for(size_t i=0;i!=t.size();i++) t[i].tonull();}
  void count_v_rules(size_t *nv);
  void predict(size_t p, size_t n, double *x, double *fp);
  void predict(size_t p, size_t n, double *x, double *b0p, double *b1p);
  void draw(double sigma, rn& gen);
  double f(size_t i);
  double getbeta0(size_t i);
  double getbeta1(size_t i);
protected:
  size_t m;  //number of trees
  std::vector<tree<T>> t; //the trees 
  pinfo pi; //prior and mcmc info
  //data
  size_t p,n; //x has dim p, n obserations
  double *x,*y,*w;  //x is column stack, pxn
  xinfo xi; //cutpoint info
  //working
  double *allfit; //if the data is set, should be f(x)
  double *r;
  double *ftemp;
  double *bfit;
  dinfo di;
};

#include <iostream>
//#include "bart.h"

//--------------------------------------------------
//constructor
template <>
bart<double>::bart():m(200),t(m,tree<double>(0.)),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),bfit(0),di() {}
template <>
bart<double>::bart(size_t im):m(im),t(m,tree<double>(0.)),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),bfit(0),di() {}
template <>
bart<double>::bart(const bart& ib):m(ib.m),t(m,tree<double>(0.)),pi(ib.pi),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),bfit(0),di() {this->t = ib.t;}
template <>
bart<std::vector<double>>::bart():m(200),t(m,tree<std::vector<double>>(std::vector<double>(2,0.))),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),bfit(0),di() {}
template <>
bart<std::vector<double>>::bart(size_t im):m(im),t(m,tree<std::vector<double>>(std::vector<double>(2,0.))),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),bfit(0),di() {}
template <>
bart<std::vector<double>>::bart(const bart& ib):m(ib.m),t(ib.t),pi(ib.pi),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),bfit(0),di() {this->t = ib.t;}
template <typename T>
bart<T>::~bart() 
{
  if(allfit) delete[] allfit;
  if(r) delete[] r;
  if(ftemp) delete[] ftemp;
  if(bfit) delete[] bfit;
}

//--------------------------------------------------
//operators
template <typename T>
bart<T>& bart<T>::operator=(const bart<T>& rhs)
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
    if(bfit) {delete[] bfit; bfit=0;}

  }
  return *this;
}
//--------------------------------------------------
//get,set
template <typename T>
void bart<T>::setm(size_t m)
{
  t.resize(m);
  this->m = t.size();

  if(allfit && (xi.size()==p)) {
    for(size_t j=0;j<n;j++) allfit[j]=0.0;
  }
}
//--------------------------------------------------
template <>
void bart<double>::setdata(size_t p,size_t n, double *x, double *y, size_t nc)
{

  this->p=p; this->n=n; this->x=x; this->y=y;
  if(xi.size()==0) makexinfo(p,n,&x[0],xi,nc);

  if(allfit) delete[] allfit;
  allfit = new double[n];
  for(size_t j=0;j<n;j++) allfit[j]=0.0;
   
  if(r) delete[] r;
  r = new double[n];

  if(ftemp) delete[] ftemp;
  ftemp = new double[n];
   
  di.n=n; di.p=p; di.x = &x[0]; di.y=r; di.w=0;

}
//--------------------------------------------------
template <>
void bart<std::vector<double>>::setdata(size_t p,size_t n, double *x, double *y, double *w, size_t nc)
{
  
  this->p=p; this->n=n; this->x=x; this->y=y; this->w=w;
  
  if(xi.size()==0) makexinfo(p,n,&x[0],xi,nc);

  if(allfit) delete[] allfit;
  allfit = new double[n];   
  for(size_t j=0;j<n;j++) allfit[j]=0.0;

  if(r) delete[] r;
  r = new double[n];
   
  if(ftemp) delete[] ftemp;
  ftemp = new double[2*n];

  if(bfit) delete[] bfit;
  bfit = new double[2*n];

  di.n=n; di.p=p; di.x = &x[0]; di.y=r; di.w=w;

}
 
template <typename T>
void bart<T>::count_v_rules(size_t *nv){
  for(size_t j=0;j<p;j++) nv[j]=0;
  for(size_t j=0;j<t.size();j++){
    std::vector<tree<T>*> thistree;
    t[j].getnodes(thistree);
    for(size_t node=0;node<thistree.size();node++){
      char ndtyp = thistree[node]->ntype();
      if(ndtyp=='t' && thistree[node]->getl()){ // nodes that are marked top node 't' AND have child should count
          nv[thistree[node]->getv()]++;
        }
      else if(ndtyp=='t' && !(thistree[node]->getl())){ }// nodes that are marked 't' AND have no child should NOT count because they are a bot node in truth
      else if(ndtyp!='b'){ // all other nodes not marked 'b' (bot) should count
        nv[thistree[node]->getv()]++;
      }
    }
  }
}

//--------------------------------------------------
template <>
void bart<double>::predict(size_t p, size_t n, double *x, double *fp) 
//uses: m,t,xi
{
  double *fptemp = new double[n];

  for(size_t j=0;j<n;j++) fp[j]=0.0;
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,fptemp);
    for(size_t k=0;k<n;k++) fp[k] += fptemp[k];
  }

  delete[] fptemp;
}
//--------------------------------------------------
template <>
void bart<std::vector<double>>::predict(size_t p, size_t n, double *x, double *b0p, double *b1p) 
//uses: m,t,xi
{
  double *fptemp = new double[2*n];

  for(size_t j=0;j<n;j++) {
    b0p[j]=0.0;
    b1p[j]=0.0;
  }
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,fptemp);
    for(size_t k=0;k<n;k++) {
      b0p[k] += fptemp[2*k];
      b1p[k] += fptemp[2*k+1];
    }
  }
  delete[] fptemp;
}
//--------------------------------------------------
template <>
void bart<double>::draw(double sigma, rn& gen)
{
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] = allfit[k] - ftemp[k];
      r[k] = y[k] - allfit[k];
    }
    bd(t[j],xi,di,pi,sigma,gen);
    drtheta(t[j],xi,di,pi,sigma,gen);
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
  }
}
//--------------------------------------------------
template <>
void bart<std::vector<double>>::draw(double sigma, rn& gen)
{
  for(size_t k=0;k<(2*n);k++) bfit[k]=0.0;
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] = allfit[k] - ftemp[2*k] - ftemp[2*k+1]*w[k];
      r[k] = y[k] - allfit[k];
    }
    bd(t[j],xi,di,pi,sigma,gen);
    drtheta(t[j],xi,di,pi,sigma,gen);
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] += (ftemp[2*k] + ftemp[2*k+1]*w[k]);
      bfit[2*k] += ftemp[2*k];
      bfit[2*k+1] += ftemp[2*k+1];
    }
  }
}
//--------------------------------------------------
// return value of sum of trees for constant terminal node
template <>
double bart<double>::f(size_t i) {return allfit[i];}
//--------------------------------------------------
// return value of sum of trees for linear terminal node
template <>
double bart<std::vector<double>>::f(size_t i) {return allfit[i];}
//--------------------------------------------------
// return value of sum of trees for intercept of linear terminal node
template <>
double bart<std::vector<double>>::getbeta0(size_t i) {return bfit[2*i];}
//--------------------------------------------------
// return value of sum of trees for slope of linear terminal node
template <>
double bart<std::vector<double>>::getbeta1(size_t i) {return bfit[2*i+1];}

//--------------------------------------------------
//public functions
template <typename T>
void bart<T>::pr() //print to screen
{
  cout << "*****bart object:\n";
  cout << "m: " << m << std::endl;
  //cout << "t[0]:\n " << t[0] << std::endl;
  //cout << "t[m-1]:\n " << t[m-1] << std::endl;
  cout << "prior and mcmc info:\n";
  pi.pr();
  if(p) cout << "data set: n,p: " << n << ", " << p << std::endl;
  else cout << "data not set\n";
}

#endif
