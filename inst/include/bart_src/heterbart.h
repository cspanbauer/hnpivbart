#ifndef GUARD_heterbart_h
#define GUARD_heterbart_h
#include "bart.h"

class heterbart : public bart
{
public:
   heterbart():bart() {}
  heterbart(size_t m, size_t size_beta):bart(m,size_beta) {}
   void pr();
   void draw(double *sigma, rn& gen);
};

//#include "heterbart.h"
#include "heterbartfuns.h"
#include "heterbd.h"

//--------------------------------------------------
void heterbart::pr()
{
   cout << "+++++heterbart object:\n";
   bart::pr();
}
//--------------------------------------------------
void heterbart::draw(double *sigma, rn& gen)
{
  for(size_t l=0;l<size_beta;l++)
    std::fill(beta[l].begin(),beta[l].end(),0.0);
  for(size_t j=0;j<m;j++) {
    //    if(size_beta==2) cout << "Tree: " << j << std::endl;
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
    heterbd(t[j],xi,di,pi,sigma,nv,pv,gen);
    heterdrmu(t[j],xi,di,pi,sigma,gen);
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
  }
}
#endif
