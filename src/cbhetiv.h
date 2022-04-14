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
   SEXP _burnf,
   SEXP _burnh,
   SEXP _m1,
   SEXP _m2,
   SEXP _nc,
   SEXP _power,
   SEXP _base,
   SEXP _tauf,
   SEXP _tauh,
   SEXP _doDP,
   SEXP _v,
   SEXP _nu,
   SEXP _a,
   SEXP _ag,
   SEXP _priag,
   SEXP _centermeans,
   SEXP _include_output,
   SEXP _fs,
   SEXP _hs,
   SEXP _mTs,
   SEXP _mYs,
   SEXP _sTs,
   SEXP _gammas,
   SEXP _sYs,
   SEXP _nullf,
   SEXP _nullh,
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
   size_t burnf = Rcpp::as<size_t>(_burnf);
   size_t burnh = Rcpp::as<size_t>(_burnh);

   //bart prior
   size_t m1 = Rcpp::as<size_t>(_m1);
   size_t m2 = Rcpp::as<size_t>(_m2);
   size_t nc = Rcpp::as<size_t>(_nc);
   double mybeta = Rcpp::as<double>(_power);
   double alpha = Rcpp::as<double>(_base);
   double tauf = Rcpp::as<double>(_tauf);
   double tauh = Rcpp::as<double>(_tauh);

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
   Rcpp::NumericVector fsv(_fs);
   double *fs = &fsv[0];
   Rcpp::NumericVector hsv(_hs);
   double *hs = &hsv[0];

   double mTs = Rcpp::as<double>(_mTs);
   double mYs = Rcpp::as<double>(_mYs);
   double sTs = Rcpp::as<double>(_sTs);
   double gammas = Rcpp::as<double>(_gammas);
   double sYs = Rcpp::as<double>(_sYs);

   bool nullf = Rcpp::as<bool>(_nullf);
   bool nullh = Rcpp::as<bool>(_nullh);
   
   //other
   size_t printevery = Rcpp::as<size_t>(_printevery);
   bool quiet = Rcpp::as<bool>(_quiet);
   
   size_t n = nT;

   //Zero vector for first BART model
   Rcpp::NumericVector zeroV(2*n);
   double *zeroVp = &zeroV[0];

   if(!quiet){
     Rprintf("*****************************************************************\n");
     Rprintf("*****Into main of cbhetiv\n");
   }


   
   //--------------------------------------------------
   // print args
   if(!quiet){
     Rprintf("***burn, nd, burnf, burnh: %ld, %ld, %ld\n",burn,nd,burnf,burnh);

     Rprintf("*** Prior:\n");
     Rprintf("m1 (num trees stage 1), m2 (num trees stage 2), nc (num cut points), %ld, %ld, %ld\n",m1,m2,nc);
     Rprintf("****power, base, tauf, tauh: %lf, %lf, %lf, %lf\n",mybeta,alpha,tauf,tauh);
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
     Rprintf("\t first and last fs: %lf, %lf\n",fs[0],fs[n-1]);
     Rprintf("\t first and last hs: %lf, %lf\n",hs[0],hs[n-1]);
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
   heterbart bmf(m1,1);
   bmf.setprior(alpha,mybeta,tauf);
   double *ytempf = new double[n*2];  //y for h bart
   double *svecf = new double[n*2];   // sigma_i for h bart
   double *txf = new double[n*2];
   for(size_t i=0;i<(2*n);i++) txf[i]=0.0;
   bmf.setdata(pz,n*2,z2,ytempf,zeroVp,nc);
   Rcpp::NumericMatrix dfburn(1,1);   //h draws on train
   if(include_output==1) dfburn(burnf,n);
   Dp::dv fhatb(n,0.0);
   if(!nullf){
     if(burnf) {
       for(size_t i=0; i<n; i++) {
         ytempf[2*i] = T[i] - mTs;
         svecf[2*i] = sTs;
         ytempf[2*i+1] = T[i] - mTs;
         svecf[2*i+1] = sTs;
       }
       for(size_t i=0;i<burnf;i++) {
         if(i%printevery==0&!quiet) Rprintf("burnf: done %d (out of %d)\n",i,burnf);
         bmf.draw(svecf,gen);
         if(include_output==1) {for(size_t j=0;j<n;j++) dfburn(i,j) = bmf.f(2*j);}
         for(size_t j=0;j<n;j++) fhatb[j] = fhatb[j] + bmf.f(2*j);
       }
       for(size_t j=0;j<n;j++) fhatb[j] = fhatb[j]/burnf;
     }
   }
   // bmf ouput storage
   Rcpp::NumericMatrix df(1,1); //f draws on train
   if(include_output==1) df(nd,n);

   //-------------------------------------------------
   //-------------------------------------------------
   // bart h setup
   //--------------------------------------------------
   heterbart bmh(m2,2);
   bmh.setprior(alpha,mybeta,tauh);
   double *ytemp = new double[n];  //y for h bart
   double *svec = new double[n];   // sigma_i for h bart
   bmh.setdata(px,n,x,ytemp,T,nc);


   //h burn-in   
   Rcpp::NumericMatrix dhburn(1,1); //h draws on train
   if(include_output==1) dhburn(burnh,n);
   double ZZ1=0.0;
   if(!nullh){
     for(size_t i=0;i<burnh;i++) {
       if(i%printevery==0&!quiet) Rprintf("burnh: done %d (out of %d)\n",i,burnh);

       // h conditional -----------------------------------
       // update ----- 
       for(size_t j=0;j<n;j++) {
         ZZ1 = (T[j] - mTs - bmf.f(2*j))/sTs;
         ytemp[j] = Y[j] - mYs - hs[j] - gammas * ZZ1;
         svec[j] = sYs; 
       }
      //draw -----
       bmh.draw(svec,gen);

       if(include_output==1){
         for(size_t j=0;j<n;j++) {
           dhburn(i,j) = bmh.f(j);
         }
       }
     }
   }

   
   Rcpp::NumericMatrix dbeta0(1,1);
   Rcpp::NumericMatrix dbeta1(1,1);
   if(include_output==1){
     dbeta0(nd,n);
     dbeta1(nd,n);
   }
   Rcpp::NumericMatrix dbeta0p(nd,nxp);
   Rcpp::NumericMatrix dbeta1p(nd,nxp);
   Rcpp::NumericMatrix dfp(nd,nzp);


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
      yS[2*i] = T[i] - fs[i];
      yS[2*i+1] = Y[i] - hs[i];
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

   double *fp = 0;
   double *beta0p = 0;
   double *beta1p = 0;
   if(doprdx) {
     beta0p = new double[nxp];
     beta1p = new double[nxp];
   }
   if(doprdz) {
     fp = new double[nzp];
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
      // h conditional -----------------------------------
      // update ----- 
      if(!nullh){
        for(size_t j=0;j<n;j++) {
          mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];
          ZZ1 = (T[j] - mT - bmf.f(2*j))/sT;
          ytemp[j] = Y[j] - mY - gamma * ZZ1;
          svec[j] = sY; 
        }
      //draw -----
      bmh.draw(svec,gen);
      }
      
      // f conditional --------------------------------      
      // update -----
      if(!nullf){
        double R;
        for(size_t j=0; j<n; j++) {
          mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];
          B0 = bmh.getbeta0(j); B1 = bmh.getbeta1(j);
          ytempf[2*j] = T[j] - mT;
          svecf[2*j] = sT;
          R = (B1*sT+gamma)*(T[j]-mT)-sT*(Y[j]-mY-B0-B1*mT);
          ytempf[2*j+1] = R/gamma;
          svecf[2*j+1] = (sT*sY)/gamma;
          //         cout << "YT: " << ytempf[2*j] << std::endl;
        }
      }
      // draws -----
      bmf.draw(svecf,gen);
      // Sigma conditional -----------------------------
      //cout << "\n&&&&&&about to do Sigma\n";
      //update -----
      for(size_t j=0;j<n;j++) {
        yS[2*j] = T[j] - bmf.f(2*j);
        yS[2*j+1] = Y[j] - bmh.f(j);
      }
      //draw -----
      dpmS.draw(gen);
      if(centermeans) dpmS.center();
      tmat = dpmS.thetaMatrix();
      
      if(doprdx){
        //        bmh.predict(px,nx,x,T,hp);
        bmh.predictb(pxp,nxp,xp,beta0p,beta1p);
      }
      if(doprdz){
        bmf.predict(pzp,nzp,zp,zeroVp,fp);
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
               dLL(j,k) = bvnorm(T[k],Y[k],bmf.f(2*k),dbeta0(j,k)+dbeta1(j,k)*T[k],dsigma1(j,k),dsigma2(j,k),dcov(j,k));
             }
           }

         if(include_output==1){
           for(size_t k=0;k<n;k++) {
             dbeta0(j,k) = bmh.getbeta0(k);
             dbeta1(j,k) = bmh.getbeta1(k);
           }
         }
         
         if(doprdx) {
           for(size_t k=0;k<nxp;k++) {
             dbeta0p(j,k) = beta0p[k];
             dbeta1p(j,k) = beta1p[k];
             //      dhp(j,k) = hp[k];
           }
         }
         
         if(doprdz) {
           for(size_t k=0;k<nzp;k++) {
             dfp(j,k) = fp[k];
             //      dhp(j,k) = hp[k];
           }
         }

         
         if(include_output==1){
           for(size_t k=0;k<n;k++) {
             df(j,k) = bmf.f(2*k);
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
     ret["dbeta0"] = dbeta0;
     ret["dbeta1"] = dbeta1;
     ret["dcov"]=dcov;   
     ret["dsigma1"]=dsigma1;
     ret["dsigma2"]=dsigma2;
     ret["dcor"]=dcor;
     ret["dgamma"]=dgamma;
     ret["dsY"]=dsY;
     ret["dLL"]=dLL;
     ret["df1"] = df;
     ret["dfburn"] = dfburn;
     ret["dhburn"] = dhburn;
   }
   if(doprdz) ret["df1.test"] = dfp;
   if(doprdx){
     ret["dbeta0.test"] = dbeta0p;
     ret["dbeta1.test"] = dbeta1p;
   }
      
   //-------------------------------------------------
   // free
   if(yS) delete [] yS;
   if(ytemp) delete [] ytemp;
   if(svec) delete [] svec;
   if(z2) delete [] z2;
   if(ytempf) delete [] ytempf;
   if(svecf) delete [] svecf;
   if(txf) delete [] txf;
   if(beta0p) delete [] beta0p;
   if(beta1p) delete [] beta1p;
   if(fp) delete [] fp;
   return ret;
}

