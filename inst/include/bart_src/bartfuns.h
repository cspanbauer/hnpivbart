#ifndef GUARD_bartfuns_h
#define GUARD_bartfuns_h
#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "rn.h"


//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
template <typename T>
double getpb(tree<T>& t, xinfo& xi, pinfo& pi, typename tree<T>::npv& goodbots);
//--------------------------------------------------
//compute n and \sum y_i with v,c
void getsuff(tree<double>& x,  tree<double>::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& syl, size_t& nr, double& syr);
//--------------------------------------------------
//compute n and \sum y_i for left and right bots
void getsuff(tree<double>& x, tree<double>::tree_p l, tree<double>::tree_p r, xinfo& xi, dinfo& di, size_t& nl, double& syl, size_t& nr, double& syr);
//--------------------------------------------------
//compute n and \sum y_i \sum w_i \sum w_i^2 \sum w_i*y_i and v,c
void getsuff(tree<std::vector<double>>& x,  tree<std::vector<double>>::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& syl, double& swl, double& swyl, double& swwl, size_t& nr, double& syr, double& swr, double& swyr, double& swwr);
//--------------------------------------------------
//compute n and \sum y_i \sum w_i \sum w_i^2 \sum w_i*y_i for left and right bots
void getsuff(tree<std::vector<double>>& x,  tree<std::vector<double>>::tree_p left,  tree<std::vector<double>>::tree_p right, xinfo& xi, dinfo& di, size_t& nl, double& syl, double& swl, double& swyl, double& swwl, size_t& nr, double& syr, double& swr, double& swyr, double& swwr);
//--------------------------------------------------
//lh, replacement for lil that only depends on sum y.
double lh(size_t n, double sy, double sigma, double tau);
//--------------------------------------------------
//draw one mu from post 
double drawnodetheta(size_t n, double sy, double tau, double sigma, rn& gen);
//--------------------------------------------------
//draw two mu (linear) from post
std::vector<double> drawnodetheta(size_t n, double sy, double sw, double swy, double sww, double tau, double tau2, double sigma, rn& gen);
//--------------------------------------------------
// draw all the bottom node mu's
void drtheta(tree<double>& t, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen);
//--------------------------------------------------
// draw all the bottom node theta's
void drtheta(tree<std::vector<double>>& t, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen);
//--------------------------------------------------
//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void allsuff(tree<double>& x, xinfo& xi, dinfo& di, tree<double>::npv& bnv, std::vector<size_t>& nv, std::vector<double>& syv);
//--------------------------------------------------
//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void allsuff(tree<std::vector<double>>& x, xinfo& xi, dinfo& di, tree<std::vector<double>>::npv& bnv, std::vector<size_t>& nv, std::vector<double>& syv, std::vector<double>& swv, std::vector<double>& swyv, std::vector<double>& swwv);
//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc);
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else a/(1+d)^b
template <typename T>
double pgrow(typename tree<T>::tree_p n, xinfo& xi, pinfo& pi);
//--------------------------------------------------
//compute n and \sum y_i for left and right bots
void getsuff(tree<double>& x,  tree<double>::tree_p l,  tree<double>::tree_p r, xinfo& xi, dinfo& di, size_t& nl, double& syl, size_t& nr, double& syr);
//--------------------------------------------------
//birth proposal
template <typename T>
void bprop(tree<T>& x, xinfo& xi, pinfo& pi, typename tree<T>::npv& goodbots, double& PBx, typename  tree<T>::tree_p& nx, size_t& v, size_t& c, double& pr, rn& gen);
//--------------------------------------------------
// death proposal
template <typename T>
void dprop(tree<T>& x, xinfo& xi, pinfo& pi, typename tree<T>::npv& goodbots, double& PBx, typename tree<T>::tree_p& nx, double& pr, rn& gen);



//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
template <typename T>
double getpb(tree<T>& t, xinfo& xi, pinfo& pi, typename tree<T>::npv& goodbots)
{
   double pb;  //prob of birth to be returned
   typename tree<T>::npv bnv; //all the bottom nodes
   t.getbots(bnv);
   for(size_t i=0;i!=bnv.size();i++)
     if(cansplit<T>(bnv[i],xi)) goodbots.push_back(bnv[i]);
   if(goodbots.size()==0) { //are there any bottom nodes you can split on?
      pb=0.0;
   } else {
      if(t.treesize()==1) pb=1.0; //is there just one node?
      else pb=pi.pb;
   }
   return pb;
}
//--------------------------------------------------
//compute n and \sum y_i for left and right give bot and v,c
void getsuff(tree<double>& x,  tree<double>::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& syl, size_t& nr, double& syr)
{
   //cout << "***into getsuff\n";

   double *xx;//current x
   nl=0; syl=0.0;
   nr=0; syr=0.0;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node
         if(xx[v] < xi[v][c]) {
               nl++;
               syl += di.y[i];
          } else {
               nr++;
               syr += di.y[i];
          }
      }
   }
}
//--------------------------------------------------
//compute n and \sum y_i for left and right bots
void getsuff(tree<double>& x, tree<double>::tree_p l, tree<double>::tree_p r, xinfo& xi, dinfo& di, size_t& nl, double& syl, size_t& nr, double& syr)
{
  double *xx;//current x
  nl=0; syl=0.0;
  nr=0; syr=0.0;
  
  for(size_t i=0;i<di.n;i++) {
    xx = di.x + i*di.p;
    tree<double>::tree_cp bn = x.bn(xx,xi);
    if(bn==l) {
      nl++;
      syl += di.y[i];
    }
    if(bn==r) {
      nr++;
      syr += di.y[i];
    }
  }
}
//--------------------------------------------------
//compute n and \sum y_i \sum w_i \sum w_i^2 \sum w_i*y_i for left and right bots
void getsuff(tree<std::vector<double>>& x,  tree<std::vector<double>>::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& syl, double& swl, double& swyl, double& swwl, size_t& nr, double& syr, double& swr, double& swyr, double& swwr)
{
   double *xx;//current x
   nl=0; syl=0.0; swl=0.0; swyl=0.0; swwl=0.0;
   nr=0; syr=0.0; swr=0.0; swyr=0.0; swwr=0.0;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx=x.bn(xx,xi)){ //does the bottom node = xx's bottom node
        if(xx[v] < xi[v][c]) {          
          nl++;
          syl += di.y[i];
          swl += di.w[i];
          swyl += di.w[i]*di.y[i];
          swwl += di.w[i]*di.w[i];
        } else {
          nr++;
          syr += di.y[i];
          swr += di.w[i];
          swyr += di.w[i]*di.y[i];
          swwr += di.w[i]*di.w[i];
        }
      }
   }
}
//--------------------------------------------------
//compute n and \sum y_i \sum w_i \sum w_i^2 \sum w_i*y_i for left and right bots
void getsuff(tree<std::vector<double>>& x, tree<std::vector<double>>::tree_p left,  tree<std::vector<double>>::tree_p right, xinfo& xi, dinfo& di, size_t& nl, double& syl, double& swl, double& swyl, double& swwl, size_t& nr, double& syr, double& swr, double& swyr, double& swwr)
{
   double *xx;//current x
   nl=0; syl=0.0; swl=0.0; swyl=0.0; swwl=0.0;
   nr=0; syr=0.0; swr=0.0; swyr=0.0; swwr=0.0;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tree<std::vector<double>>::tree_cp bn = x.bn(xx,xi);
      if(bn==left) {
         nl++;
         syl += di.y[i];
         swl += di.w[i];
         swyl += di.w[i]*di.y[i];
         swwl += di.w[i]*di.w[i];
      }
      if(bn==right) {
         nr++;
         syr += di.y[i];
         swr += di.w[i];
         swyr += di.w[i]*di.y[i];
         swwr += di.w[i]*di.w[i];
      }
   }
}
//lh, replacement for lil that only depends on sum y.
double lh(size_t n, double sy, double sigma, double tau)
{
   double s2 = sigma*sigma;
   double t2 = tau*tau;
   double k = n*t2+s2;
   return -.5*log(k) + ((t2*sy*sy)/(2.0*s2*k));
}
//--------------------------------------------------
//draw one mu from post 
double drawnodetheta(size_t n, double sy, double tau, double sigma, rn& gen)
{
   double s2 = sigma*sigma;
   double b = n/s2;
   double a = 1.0/(tau*tau);
   return (sy/s2)/(a+b) + gen.normal()/sqrt(a+b);
}
//--------------------------------------------------
//draw two mu (linear) from post
std::vector<double> drawnodetheta(size_t n, double sy, double sw, double swy, double sww, double tau, double tau2, double sigma, rn& gen)
{
  double t2 = tau*tau;
  double s2 = sigma*sigma;
  double inv_V00 = (double)n/s2 + 1./t2;
  double inv_V01 = sw/s2;
  double inv_V11 = sww/s2 + 1./(tau2*tau2);
  double det_invV = inv_V00*inv_V11-inv_V01*inv_V01;
  double V00 = (1./det_invV)*inv_V11;
  double V11 = (1./det_invV)*inv_V00;
  double V01 = -(1./det_invV)*inv_V01;
  double wty0 = sy/s2;
  double wty1 = swy/s2;
  double m0 = wty0*inv_V00+wty1*inv_V01;
  double m1 = wty0*inv_V01+wty1*inv_V11;
  double chol_V00 = sqrt(V00);
  double chol_V01 = V01/chol_V00;
  double chol_V11 = sqrt(V11-chol_V01*chol_V01);
  std::vector<double> out(2,0.);
  double Z0 = gen.normal();
  double Z1 = gen.normal();
  out[0] = chol_V00*Z0 + m0;
  out[1] = chol_V01*Z0 + chol_V11*Z1 + m1;
  return out;
}
//--------------------------------------------------
// draw all the bottom node theta's
void drtheta(tree<double>& t, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen)
{
  typename tree<double>::npv bnv;
   std::vector<size_t> nv;
   std::vector<double> syv;
   allsuff(t,xi,di,bnv,nv,syv);

   for(tree<double>::npv::size_type i=0;i!=bnv.size();i++) 
      bnv[i]->settheta(drawnodetheta(nv[i],syv[i],pi.tau,sigma,gen));
}
//--------------------------------------------------
// draw all the bottom node theta's
void drtheta(tree<std::vector<double>>& t, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen)
{
  tree<std::vector<double>>::npv bnv;
  std::vector<size_t> nv;
  std::vector<double> syv;
  std::vector<double> swv;
  std::vector<double> swyv;
  std::vector<double> swwv;     
  allsuff(t,xi,di,bnv,nv,syv,swv,swyv,swwv);

   for(tree<std::vector<double>>::npv::size_type i=0;i!=bnv.size();i++) 
     bnv[i]->settheta(drawnodetheta(nv[i],syv[i],swv[i],swyv[i],swwv[i],pi.tau,pi.tau2,sigma,gen));
}

//--------------------------------------------------
//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void allsuff(tree<double>& x, xinfo& xi, dinfo& di, tree<double>::npv& bnv, std::vector<size_t>& nv, std::vector<double>& syv)
{
    tree<double>::tree_cp tbn; //the pointer to the bottom node for the current observations
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x

   bnv.clear();
   x.getbots(bnv);

   typedef tree<double>::npv::size_type bvsz;
   bvsz nb = bnv.size();
   nv.resize(nb);
   syv.resize(nb);

   std::map< tree<double>::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) {bnmap[bnv[i]]=i;nv[i]=0;syv[i]=0.0;}

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tbn = x.bn(xx,xi);
      ni = bnmap[tbn];

      ++(nv[ni]);
      syv[ni] += di.y[i];
   }
}
//--------------------------------------------------
//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void allsuff(tree<std::vector<double>>& x, xinfo& xi, dinfo& di, tree<std::vector<double>>::npv& bnv, std::vector<size_t>& nv, std::vector<double>& syv, std::vector<double>& swv, std::vector<double>& swyv, std::vector<double>& swwv)
{
  tree<std::vector<double>>::tree_cp tbn; //the pointer to the bottom node for the current observations
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x

   bnv.clear();
   x.getbots(bnv);

   typedef tree<std::vector<double>>::npv::size_type bvsz;
   bvsz nb = bnv.size();
   nv.resize(nb);
   syv.resize(nb);
   swv.resize(nb);
   swwv.resize(nb);
   swwv.resize(nb);

   std::map< tree<std::vector<double>>::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) {
     bnmap[bnv[i]]=i;
     nv[i]=0;
     syv[i]=0.0;
     swv[i]=0.0;
     swyv[i]=0.0;
     swwv[i]=0.0;}

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tbn = x.bn(xx,xi);
      ni = bnmap[tbn];

      ++(nv[ni]);
      syv[ni] += di.y[i];
      syv[ni] += di.w[i];
      swyv[ni] += di.y[i]*di.w[i];
      swwv[ni] += di.w[i]*di.w[i];
   }
}
//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc)
{
   double xinc;

   //compute min and max for each x
   std::vector<double> minx(p,INFINITY);
   std::vector<double> maxx(p,-INFINITY);
   double xx;
   for(size_t i=0;i<p;i++) {
      for(size_t j=0;j<n;j++) {
         xx = *(x+p*j+i);
         if(xx < minx[i]) minx[i]=xx;
         if(xx > maxx[i]) maxx[i]=xx;
      }
   }
   //make grid of nc cutpoints between min and max for each x.
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
      xinc = (maxx[i]-minx[i])/(nc+1.0);
      xi[i].resize(nc);
      for(size_t j=0;j<nc;j++) xi[i][j] = minx[i] + (j+1)*xinc;
   }
}
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
template <class T>
double pgrow(typename tree<T>::tree_p n, xinfo& xi, pinfo& pi)
{
  if(cansplit<T>(n,xi)) {
      return pi.alpha/pow(1.0+n->depth(),pi.mybeta);
   } else {
      return 0.0;
   }
}
//--------------------------------------------------
//bprop: function to generate birth proposal
template <class T>
void bprop(tree<T>& x, xinfo& xi, pinfo& pi, typename tree<T>::npv& goodbots, double& PBx, typename tree<T>::tree_p& nx, size_t& v, size_t& c, double& pr, rn& gen)
{
      //draw bottom node, choose node index ni from list in goodbots
      size_t ni = floor(gen.uniform()*goodbots.size());
      nx = goodbots[ni]; //the bottom node we might birth at

      //draw v,  the variable
      std::vector<size_t> goodvars; //variables nx can split on
      getgoodvars<T>(nx,xi,goodvars);
      double uu = gen.uniform();
      size_t vi=0; //index of chosen split variable
      double csm = pi.pv[0]; // cumulative sum of pv
      while(csm <= uu){
        vi++;
        csm+=pi.pv[vi];
      }
      //     cout << "Hello " << pi.pv[0] << '\n';
      //size_t vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable
      v = vi;//goodvars[vi];
      //draw c, the cutpoint
      int L,U;
      L=0; U = xi[v].size()-1;
      nx->rg(v,&L,&U);
      c = L + floor(gen.uniform()*(U-L+1)); //U-L+1 is number of available split points

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
      size_t dnx = nx->depth();
      double PGnx = pi.alpha/pow(1.0 + dnx,pi.mybeta); //prior prob of growing at nx

      double PGly, PGry; //prior probs of growing at new children (l and r) of proposal
      if(goodvars.size()>1) { //know there are variables we could split l and r on
         PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.mybeta); //depth of new nodes would be one more
         PGry = PGly;
      } else { //only had v to work with, if it is exhausted at either child need PG=0
         if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
            PGly = 0.0;
         } else {
            PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.mybeta);
         }
         if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
            PGry = 0.0;
         } else {
            PGry = pi.alpha/pow(1.0 + dnx+1.0,pi.mybeta);
         }
      }
      double PDy; //prob of proposing death at y
      if(goodbots.size()>1) { //can birth at y because splittable nodes left
         PDy = 1.0 - pi.pb;
      } else { //nx was the only node you could split on
         if((PGry==0) && (PGly==0)) { //cannot birth at y
            PDy=1.0;
         } else { //y can birth at either l or r
            PDy = 1.0 - pi.pb;
         }
      }

      double Pnogy; //death prob of choosing the nog node at y
      size_t nnogs = x.nnogs();
      typename tree<T>::tree_p nxp = nx->getp();
      if(nxp==0) { //no parent, nx is the top and only node
         Pnogy=1.0;
      } else {
         if(nxp->ntype() == 'n') { //if parent is a nog, number of nogs same at x and y
            Pnogy = 1.0/nnogs;
         } else { //if parent is not a nog, y has one more nog.
           Pnogy = 1.0/(nnogs+1.0);
         }
      }

      pr = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*Pbotx*PBx);
}
//--------------------------------------------------
// death proposal
template <class T>
void dprop(tree<T>& x, xinfo& xi, pinfo& pi, typename tree<T>::npv& goodbots, double& PBx, typename tree<T>::tree_p& nx, double& pr, rn& gen)
{
      //draw nog node, any nog node is a possibility
      typename tree<T>::npv nognds; //nog nodes
      x.getnogs(nognds);
      size_t ni = floor(gen.uniform()*nognds.size());
      nx = nognds[ni]; //the nog node we might kill children at

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double PGny; //prob the nog node grows
      size_t dny = nx->depth();
      PGny = pi.alpha/pow(1.0+dny,pi.mybeta);

      //better way to code these two?
      double PGlx = pgrow<T>(nx->getl(),xi,pi);
      double PGrx = pgrow<T>(nx->getr(),xi,pi);

      double PBy;  //prob of birth move at y
      if(nx->ntype()=='t') { //is the nog node nx the top node
         PBy = 1.0;
      } else {
         PBy = pi.pb;
      }

      double Pboty;  //prob of choosing the nog as bot to split on when y
      int ngood = goodbots.size();
      if(cansplit<T>(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
      if(cansplit<T>(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
      ++ngood;  //know you can split at nx
      Pboty=1.0/ngood;

      double PDx = 1.0-PBx; //prob of a death step at x
      double Pnogx = 1.0/nognds.size();

      pr =  ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
}

#endif
