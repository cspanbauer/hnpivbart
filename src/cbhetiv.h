/*
 *  hnpivBART: Heterogenous Non-Parametric IV Bayesian Additive Regression Trees
 *  Copyright (C) 2021 Charles Spanbauer
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

RcppExport SEXP cbhetiv(
   SEXP _z,
   SEXP _zp,
   SEXP _x,
   SEXP _xp,
   SEXP _T,
   SEXP _Y,
   SEXP _burn,
   SEXP _nd,
   SEXP _burnf1,
   SEXP _burnf2,
   SEXP _m1,
   SEXP _m2,
   SEXP _nc,
   SEXP _power,
   SEXP _base,
   SEXP _tauf1,
   SEXP _tauf2,
   SEXP _doDP,
   SEXP _v,
   SEXP _nu,
   SEXP _a,
   SEXP _ag,
   SEXP _priag,
   SEXP _centermeans,
   SEXP _include_output,
   SEXP _f1s,
   SEXP _f2s,
   SEXP _mTs,
   SEXP _mYs,
   SEXP _sTs,
   SEXP _gammas,
   SEXP _sYs,
   SEXP _nullf1,
   SEXP _nullf2,
   SEXP _printevery,
   SEXP _quiet
)
{
   //-----------------------------------------------------------
   //random number generation
   GetRNGstate();
   rn gen;

   //--------------------------------------------------
   //process arguments

   //z
   Rcpp::NumericMatrix zm(_z);
   size_t nz = zm.ncol();
   size_t pz = zm.nrow();
   double *z = &zm[0];

   //zp
   Rcpp::NumericMatrix zpm(_zp);
   size_t nzp = zpm.ncol();
   size_t pzp = zpm.nrow();
   double *zp = &zm[0];
   
   //x
   Rcpp::NumericMatrix xm(_x);
   size_t nx = xm.ncol();
   size_t px = xm.nrow();
   double *x = &xm[0];
   
   //xp
   Rcpp::NumericMatrix xpm(_xp);
   size_t nxp = xpm.ncol();
   size_t pxp = xpm.nrow();
   double *xp = &xpm[0];

   // T,Y ------------
   Rcpp::NumericVector Tv(_T);
   double *T = &Tv[0];
   size_t nT = Tv.size();

   Rcpp::NumericVector Yv(_Y);
   double *Y = &Yv[0];
   size_t nY = Yv.size();

   // burn, nd
   size_t burn = Rcpp::as<size_t>(_burn);
   size_t nd = Rcpp::as<size_t>(_nd);
   size_t burnf1 = Rcpp::as<size_t>(_burnf1);
   size_t burnf2 = Rcpp::as<size_t>(_burnf2);

   //bart prior
   size_t m1 = Rcpp::as<size_t>(_m1);
   size_t m2 = Rcpp::as<size_t>(_m2);
   size_t nc = Rcpp::as<size_t>(_nc);
   double mybeta = Rcpp::as<double>(_power);
   double alpha = Rcpp::as<double>(_base);
   double tauf1 = Rcpp::as<double>(_tauf1);
   double tauf2 = Rcpp::as<double>(_tauf2);

   //base prior -------------
   double v = Rcpp::as<double>(_v);
   double nu = Rcpp::as<double>(_nu);
   double a = Rcpp::as<double>(_a);

   //alpha prior--------
   Rcpp::NumericVector agv(_ag);
   Rcpp::NumericVector priagv(_priag);
   size_t nag=agv.size();
   Dp::dv ag(nag);
   Dp::dv priag(nag);
   for(size_t i=0;i<nag;i++) { // should use STL
      ag[i]=agv[i];
      priag[i]=priagv[i];
   }

   //should the means be centered after each DpMuSigma draw
   bool centermeans = Rcpp::as<bool>(_centermeans);
   //   size_t include_output = Rcpp::as<size_t>(_include_output);
   size_t include_output = 0;
   bool doDP = Rcpp::as<bool>(_doDP);

   //starting values
   //fs,hs
   Rcpp::NumericVector f1sv(_f1s);
   double *f1s = &f1sv[0];
   Rcpp::NumericVector f2sv(_f2s);
   double *f2s = &f2sv[0];

   double mTs = Rcpp::as<double>(_mTs);
   double mYs = Rcpp::as<double>(_mYs);
   double sTs = Rcpp::as<double>(_sTs);
   double gammas = Rcpp::as<double>(_gammas);
   double sYs = Rcpp::as<double>(_sYs);

   bool nullf1 = Rcpp::as<bool>(_nullf1);
   bool nullf2 = Rcpp::as<bool>(_nullf2);
   
   //other
   size_t printevery = Rcpp::as<size_t>(_printevery);
   bool quiet = Rcpp::as<bool>(_quiet);
   
   size_t n = nT;

   if(!quiet){
     Rprintf("*****************************************************************\n");
     Rprintf("*****Into main of cbhetiv\n");
   }


   
   //--------------------------------------------------
   // print args
   if(!quiet){
     Rprintf("***burn, nd, burnf1, burnf2: %ld, %ld, %ld\n",burn,nd,burnf1,burnf2);

     Rprintf("*** Prior:\n");
     Rprintf("m1 (num trees stage 1), m2 (num trees stage 2), nc (num cut points), %ld, %ld, %ld\n",m1,m2,nc);
     Rprintf("****power, base, tauf1, tauf2: %lf, %lf, %lf, %lf\n",mybeta,alpha,tauf1,tauf2);
     Rprintf("v: %lf\n",v);
     Rprintf("nu: %lf\n",nu);
     Rprintf("a: %lf\n",a);
     Rprintf("alpha prior: grid size: %ld\n",ag.size());
     Rprintf("\tag, first, last: %lf, %lf\n",ag[0],ag[nag-1]);
     Rprintf("\tpriag, first, last: %lf, %lf\n",priag[0],priag[nag-1]);
     
     Rprintf("*** Data:\n");
     Rprintf("nz,pz: %ld, %ld\n",nz,pz);
     Rprintf("nx,px: %ld, %ld\n",nx,px);
     Rprintf("nT: %ld\n",nT);
     Rprintf("first and last T: %lf, %lf\n",T[0],T[nT-1]);
     Rprintf("nY: %ld\n",nY);
     Rprintf("first and last Y: %lf, %lf\n",Y[0],Y[nY-1]);
     
     Rprintf("*** starting values, n is: %ld\n",n);
     Rprintf("\t first and last fs: %lf, %lf\n",f1s[0],f1s[n-1]);
     Rprintf("\t first and last hs: %lf, %lf\n",f2s[0],f2s[n-1]);
     Rprintf(" starting for mT, mY, sT, gamma, sY: %lf, %lf, %lf, %lf, %lf\n",mTs,mYs,sTs,gammas,sYs);
     
     Rprintf("***other\n");
     Rprintf("printevery: %ld\n",printevery);
   }
     

   //-------------------------------------------------
   //-------------------------------------------------
   // bart f setup
   //--------------------------------------------------
   double *z2 = new double[n*pz*2];
   size_t start=0;
   for(size_t i=0;i<n;i++) {
      start = i*2*pz;
      for(size_t j=0;j<pz;j++) {
         z2[start+j] = z[i*pz+j];
      }
      start += pz;
      for(size_t j=0;j<pz;j++) {
         z2[start+j] = z[i*pz+j];
      }
   }
   heterbart<double> bmf1(m1);
   bmf1.setprior(alpha,mybeta,tauf1);
   double *ytempf1 = new double[n*2];  //y for h bart
   double *svecf1 = new double[n*2];   // sigma_i for h bart
   bmf1.setdata(pz,n*2,z2,ytempf1,nc);
   Rcpp::NumericMatrix df1burn(1,1);   //h draws on train
   if(include_output==1) df1burn(burnf1,n);
   Dp::dv fhatb(n,0.0);
   if(!nullf1){
     if(burnf1) {
       for(size_t i=0; i<n; i++) {
         ytempf1[2*i] = T[i] - mTs;
         svecf1[2*i] = sTs;
         ytempf1[2*i+1] = T[i] - mTs;
         svecf1[2*i+1] = sTs;
       }
       for(size_t i=0;i<burnf1;i++) {
         if(i%printevery==0&!quiet) Rprintf("burnf1: done %d (out of %d)\n",i,burnf1);
         bmf1.draw(svecf1,gen);
         if(include_output==1) {for(size_t j=0;j<n;j++) df1burn(i,j) = bmf1.f(2*j);}
         for(size_t j=0;j<n;j++) fhatb[j] = fhatb[j] + bmf1.f(2*j);
       }
       for(size_t j=0;j<n;j++) fhatb[j] = fhatb[j]/burnf1;
     }
   }
   // bmf1 ouput storage
   Rcpp::NumericMatrix df1(1,1); //f draws on train
   if(include_output==1) df1(nd,n);

   //-------------------------------------------------
   //-------------------------------------------------
   // bart f2 setup
   //--------------------------------------------------
   heterbart<std::vector<double>> bmf2(m2);
   bmf2.setprior(alpha,mybeta,tauf2);
   double *ytempf2 = new double[n];  //y for h bart
   double *svecf2 = new double[n];   // sigma_i for h bart
   bmf2.setdata(px,n,x,ytempf2,T,nc);

   //h burn-in   
   Rcpp::NumericMatrix df2burn(1,1); //h draws on train
   if(include_output==1) df2burn(burnf2,n);
   double ZZ1=0.0;
   if(!nullf2){
     for(size_t i=0;i<burnf2;i++) {
       if(i%printevery==0&!quiet) Rprintf("burnf2: done %d (out of %d)\n",i,burnf2);

       // f2 conditional -----------------------------------
       // update ----- 
       for(size_t j=0;j<n;j++) {
         ytempf2[j] = Y[j] - mYs - f2s[j] - gammas * ZZ1;
         svecf2[j] = sYs; 
       }
      //draw -----
       bmf2.draw(svecf2,gen);
       if(include_output==1){
         for(size_t j=0;j<n;j++) {
           df2burn(i,j) = bmf2.f(j);
         }
       }
     }
   }

   
   Rcpp::NumericMatrix df22(1,1);
   Rcpp::NumericMatrix df21(1,1);
   if(include_output==1){
     df22(nd,n);
     df21(nd,n);
   }
   Rcpp::NumericMatrix df21p(nd,nxp);
   Rcpp::NumericMatrix df22p(nd,nxp);
   Rcpp::NumericMatrix df1p(nd,nzp);


   //--------------------------------------------------
   //--------------------------------------------------
   // DpMuSigma setup
   if(!quiet) cout << "\n*** making a DpMuSigma object:\n";
   
   Dp::dv itheta(5,0);
   itheta[0] = mTs; itheta[1] = mYs; itheta[2] = sTs; itheta[3] = gammas; itheta[4] = sYs;
   if(!quiet){
     cout << "itheta:\n";
     printdv(itheta);
   }

   //data for errors
   double *yS = new double[2*n];
   //intialize using starting values
   for(size_t i=0;i<n;i++) {
      yS[2*i] = T[i] - f1s[i];
      yS[2*i+1] = Y[i] - f2s[i];
   }
   
   DpMuSigma dpmS(n,itheta,doDP);
   dpmS.setData(yS);
   dpmS.setPrior(v,nu,a);
   dpmS.setAlphaPrior(ag,priag);
   dpmS.setalpha(1.0);

   // to screen
   if(!quiet)
   dpmS.toscreen();

   //output storage for DpMuSigma
   Rcpp::IntegerVector dnpart(nd);
   Rcpp::NumericVector dalpha(nd);
   Rcpp::NumericMatrix dsigma1(1,1);
   Rcpp::NumericMatrix dsigma2(1,1);
   Rcpp::NumericMatrix dcov(1,1);
   Rcpp::NumericMatrix dcor(1,1);
   Rcpp::NumericMatrix dgamma(1,1);
   Rcpp::NumericMatrix dsY(1,1);
   Rcpp::NumericMatrix dLL(1,1);
   if(include_output==1){
     dcor(nd,n);
     dgamma(nd,n); 
     dsY(nd,n);    
     dLL(nd,n);
     dsigma1(nd,n);
     dsigma2(nd,n);     
   }

   bool doprdz = pzp>0;
   bool doprdx = pxp>0;

   double *f1p = 0;
   double *f22p = 0;
   double *f21p = 0;
   if(doprdx) {
     f21p = new double[nxp];
     f22p = new double[nxp];
   }
   if(doprdz) {
     f1p = new double[nzp];
   }

   //--------------------------------------------------
   //MCMC
   time_t tp;
   int time1 = time(&tp);
   Dp::dvv tmat; 
   tmat = dpmS.thetaMatrix();

   double mT=mTs,mY=mYs,sT=sTs,gamma=gammas,sY=sYs; //temporary values for mean and Lchol of Sigma
   double B0,B1;
   for(size_t i=0;i<(nd+burn);i++) {
      if(i%printevery==0&!quiet) Rprintf("done %d (out of %d)\n",i,nd+burn);
      // f2 conditional -----------------------------------
      // update ----- 
      if(!nullf2){
        for(size_t j=0;j<n;j++) {
          mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];
          if(!nullf1) ZZ1 = (T[j] - mT - bmf1.f(2*j))/sT;
          else ZZ1 = 0.;
          ytempf2[j] = Y[j] - mY - gamma * ZZ1;
          svecf2[j] = sY; 
        }
      //draw -----
      bmf2.draw(svecf2,gen);
      }
      
      // f1 conditional --------------------------------      
      // update -----
      if(!nullf1){
        double R;
        for(size_t j=0; j<n; j++) {
          mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];
          B0 = bmf2.getbeta0(j); B1 = bmf2.getbeta1(j);
          ytempf1[2*j] = T[j] - mT;
          svecf1[2*j] = sT;
          R = (B1*sT+gamma)*(T[j]-mT)-sT*(Y[j]-mY-B0-B1*mT);
          ytempf1[2*j+1] = R/gamma;
          svecf1[2*j+1] = (sT*sY)/gamma;
          //         cout << "YT: " << ytempf[2*j] << std::endl;
        }
      }
      // draws -----
      bmf1.draw(svecf1,gen);
      // Sigma conditional -----------------------------
      //cout << "\n&&&&&&about to do Sigma\n";
      //update -----
      for(size_t j=0;j<n;j++) {
        yS[2*j] = T[j] - bmf1.f(2*j);
        yS[2*j+1] = Y[j] - bmf2.f(j);
      }
      
      //draw -----
      dpmS.draw(gen);
      if(centermeans) dpmS.center();
      tmat = dpmS.thetaMatrix();
      
      if(doprdx){
        bmf2.predict(pxp,nxp,xp,f22p,f21p);
      }
      if(doprdz){
        bmf1.predict(pzp,nzp,zp,f1p);
      }
      
      if(i >= burn) {
         size_t j=i-burn;

         dnpart[j] = dpmS.npart();
         dalpha[j] = dpmS.getalpha();
           for(size_t k=0;k<n;k++) {
             if(include_output==1){               
               dsigma1(j,k) = tmat[k][2];
               dsigma2(j,k) = sqrt(tmat[k][3]*tmat[k][3] + tmat[k][4]*tmat[k][4]);
               dcov(j,k) = tmat[k][3]*tmat[k][2];
               dcor(j,k) = dcov(j,k)/(dsigma1(j,k)*dsigma2(j,k));
               dgamma(j,k) = tmat[k][3];
               dsY(j,k) = tmat[k][4];
               dLL(j,k) = bvnorm(T[k],Y[k],bmf1.f(2*k),df22(j,k)+df21(j,k)*T[k],dsigma1(j,k),dsigma2(j,k),dcov(j,k));
             }
           }

         if(include_output==1){
           for(size_t k=0;k<n;k++) {
             df21(j,k) = bmf2.getbeta1(k);
             df22(j,k) = bmf2.getbeta0(k);
           }
         }
         
         if(doprdx) {
           for(size_t k=0;k<nxp;k++) {
             df22p(j,k) = f22p[k];
             df21p(j,k) = f21p[k];
           }
         }
         
         if(doprdz) {
           for(size_t k=0;k<nzp;k++) {
             df1p(j,k) = f1p[k];
           }
         }

         
         if(include_output==1){
           for(size_t k=0;k<n;k++) {
             df1(j,k) = bmf1.f(2*k);
           }
         }
      }
   }
   int time2 = time(&tp);
   if(!quiet) Rprintf("time: %d\n",time2-time1);
   
      
   //--------------------------------------------------
   // clear RNG at end
   PutRNGstate();
      
   Rcpp::List ret;
   ret["check"] = "biv";
   ret["dnpart"]=dnpart;
   ret["dalpha"]=dalpha;
   if(include_output==1){
     ret["df22"] = df22;
     ret["df21"] = df21;
     ret["dcov"]=dcov;   
     ret["dsigma1"]=dsigma1;
     ret["dsigma2"]=dsigma2;
     ret["dcor"]=dcor;
     ret["dgamma"]=dgamma;
     ret["dsY"]=dsY;
     ret["dLL"]=dLL;
     ret["df1"] = df1;
     ret["df1burn"] = df1burn;
     ret["df2burn"] = df2burn;
   }
   if(doprdz) ret["df1.test"] = df1p;
   if(doprdx){
     ret["df22.test"] = df22p;
     ret["df21.test"] = df21p;
   }
      
   //-------------------------------------------------
   // free
   if(yS) delete [] yS;
   if(ytempf1) delete [] ytempf1;
   if(svecf1) delete [] svecf1;
   if(z2) delete [] z2;
   if(ytempf2) delete [] ytempf2;
   if(svecf2) delete [] svecf2;
   if(f1p) delete [] f1p;
   if(f21p) delete [] f21p;
   if(f22p) delete [] f22p;
   return ret;
}

