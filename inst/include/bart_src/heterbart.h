#ifndef GUARD_heterbart_h
#define GUARD_heterbart_h
#include "bart.h"

template <class T>
class heterbart : public bart<T>
{
public:
  heterbart():bart<T>() {}
  heterbart(size_t m):bart<T>(m) {}
  void pr();
  void draw(double *sigma, rn& gen);
// private:
//   using bart<T>::m;
//   using bart<T>::t;
//   using bart<T>::p;
//   using bart<T>::n;
//   using bart<T>::x;
//   using bart<T>::ftemp;
//   using bart<T>::allfit;
//   using bart<T>::bfit;
//   using bart<T>::r;
//   using bart<T>::y;
//   using bart<T>::w;
//   using bart<T>::di;
//   using bart<T>::pi;
//   using bart<T>::xi;    
};

//#include "heterbart.h"
#include "heterbartfuns.h"
#include "heterbd.h"

//--------------------------------------------------
template <typename T>
void heterbart<T>::pr()
{
  cout << "+++++heterbart object:\n";
  bart<T>::pr();
}
//--------------------------------------------------
template <>
void heterbart<double>::draw(double *sigma, rn& gen)
{
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] = allfit[k] - ftemp[k];
      r[k] = y[k] - allfit[k];
    }
    heterbd(t[j],xi,di,pi,sigma,gen);
    heterdrtheta(t[j],xi,di,pi,sigma,gen);
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<(n);k++) allfit[k] += ftemp[k];
  }
}
//--------------------------------------------------
template <>
void heterbart<std::vector<double>>::draw(double *sigma, rn& gen)
{
  for(size_t k=0;k<(2*n);k++) bfit[k]=0.0;
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] = allfit[k] - ftemp[2*k] - ftemp[2*k+1]*w[k];
      r[k] = y[k] - allfit[k];
    }    
    heterbd(t[j],xi,di,pi,sigma,gen);
    heterdrtheta(t[j],xi,di,pi,sigma,gen);
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] += (ftemp[2*k] + ftemp[2*k+1]*w[k]);
      bfit[2*k] += ftemp[2*k];
      bfit[2*k+1] += ftemp[2*k+1];
    }
  }
}
#endif
