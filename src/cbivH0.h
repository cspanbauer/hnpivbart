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

RcppExport SEXP cbivH0(
   SEXP _zx,
   SEXP _zxp,
   SEXP _x,
   SEXP _xp,
   SEXP _T,
   SEXP _Tp,
   SEXP _Y,
   SEXP _burn,
   SEXP _nd,
   SEXP _burnf,
   SEXP _burnh0,
   SEXP _burnh1,
   SEXP _m,
   SEXP _nc,
   SEXP _power,
   SEXP _base,
   SEXP _tauf,
   SEXP _tauh0,
   SEXP _tauh1,
   SEXP _v,
   SEXP _nu,
   SEXP _a,
   SEXP _ag,
   SEXP _priag,
   SEXP _centermeans,
   SEXP _doDP,
   SEXP _fs,
   SEXP _h0s,
   SEXP _h1s,
   SEXP _mTs,
   SEXP _mYs,
   SEXP _sTs,
   SEXP _gammas,
   SEXP _sYs,
   SEXP _printevery,
   SEXP _quiet
)
{
   Rprintf("*****************************************************************\n");
   Rprintf("*****Into main of cbivH0\n");

   //-----------------------------------------------------------
   //random number generation
   GetRNGstate();
   rrn gen;

   //--------------------------------------------------
   //process arguments

   //zx
   Rcpp::NumericMatrix zxm(_zx);
   size_t nzx = zxm.ncol();
   size_t pzx = zxm.nrow();
   double *zx = &zxm[0];

   //zp
   Rcpp::NumericMatrix zpm(_zxp);
   size_t nzp = zpm.ncol();
   size_t pzp = zpm.nrow();
   double *zp = &zpm[0];
   
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

   // Tp   ------------
   Rcpp::NumericVector Tpv(_Tp);
   double *Tp = &Tpv[0];
   size_t nTp = Tpv.size();
   
   // burn, nd
   size_t burn = Rcpp::as<size_t>(_burn);
   size_t nd = Rcpp::as<size_t>(_nd);
   size_t burnf = Rcpp::as<size_t>(_burnf);
   size_t burnh0 = 0.;//Rcpp::as<size_t>(_burnh0);
   size_t burnh1 = Rcpp::as<size_t>(_burnh1);

   //bart prior
   size_t m = Rcpp::as<size_t>(_m);
   size_t nc = Rcpp::as<size_t>(_nc);
   double mybeta = Rcpp::as<double>(_power);
   double alpha = Rcpp::as<double>(_base);
   double tauf = Rcpp::as<double>(_tauf);
   double tauh0 = Rcpp::as<double>(_tauh0);
   double tauh1 = Rcpp::as<double>(_tauh1);


   
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
   bool doDP = Rcpp::as<bool>(_doDP);

   
   //starting values
   //fs,hs
   Rcpp::NumericVector fsv(_fs);
   double *fs = &fsv[0];
   Rcpp::NumericVector h0sv(_h0s);
   double *h0s = &h0sv[0];
   Rcpp::NumericVector h1sv(_h1s);
   double *h1s = &h1sv[0];



   double mTs = Rcpp::as<double>(_mTs);
   double mYs = Rcpp::as<double>(_mYs);
   double sTs = Rcpp::as<double>(_sTs);
   double gammas = Rcpp::as<double>(_gammas);
   double sYs = Rcpp::as<double>(_sYs);
   //other
   size_t printevery = Rcpp::as<size_t>(_printevery);
   bool quiet = Rcpp::as<bool>(_quiet);
   
   size_t n = nT;

   //Zero vector for first BART model
   Rcpp::NumericVector zeroV(2*n);
   double *zeroVp = &zeroV[0];

   
   
   //--------------------------------------------------
   // print args

   // doh1
   bool doh1 = (nx!=0);
   if(!doh1) burnh1=1;
  
   if(!quiet){
     Rprintf("***burn, nd, burnf, burnh0, burnh1: %ld, %ld, %ld %ld\n",burn,nd,burnf,burnh0,burnh1);

     Rprintf("*** Prior:\n");
     Rprintf("m (num trees), nc (num cut points), %ld, %ld\n",m,nc);
     Rprintf("****power, base, tauf, tauh: %lf, %lf, %lf, %lf, %lf\n",mybeta,alpha,tauf,tauh0,tauh1);
     Rprintf("v: %lf\n",v);
     Rprintf("nu: %lf\n",nu);
     Rprintf("a: %lf\n",a);
     Rprintf("alpha prior: grid size: %ld\n",ag.size());
     Rprintf("\tag, first, last: %lf, %lf\n",ag[0],ag[nag-1]);
     Rprintf("\tpriag, first, last: %lf, %lf\n",priag[0],priag[nag-1]);

     Rprintf("*** Data:\n");
     Rprintf("nzx,pzx: %ld, %ld\n",nzx,pzx);
     Rprintf("nx,px: %ld, %ld\n",nx,px);
     Rprintf("nT: %ld\n",nT);
     Rprintf("first and last T: %lf, %lf\n",T[0],T[nT-1]);
     Rprintf("nY: %ld\n",nY);
     Rprintf("first and last Y: %lf, %lf\n",Y[0],Y[nY-1]);

     Rprintf("*** starting values, n is: %ld\n",n);
     Rprintf("\t first and last fs: %lf, %lf\n",fs[0],fs[n-1]);
     Rprintf("\t first and last h0s: %lf, %lf\n",h0s[0],h0s[n-1]);
     Rprintf("\t first and last h1s: %lf, %lf\n",h1s[0],h1s[n-1]);
     Rprintf(" starting for mT, mY, sT, gamma, sY: %lf, %lf, %lf, %lf, %lf\n",mTs,mYs,sTs,gammas,sYs);

     Rprintf("***other\n");
     Rprintf("printevery: %ld\n",printevery);
   }
  

   //-------------------------------------------------
   //-------------------------------------------------
   // bart f setup
   //--------------------------------------------------
   // first make zx2 which is zx with the columns repeated
   double *zx2 = new double[n*pzx*2];
   size_t start=0;
   for(size_t i=0;i<n;i++) {
      start = i*2*pzx;
      for(size_t j=0;j<pzx;j++) {
         zx2[start+j] = zx[i*pzx+j];
      }
      start += pzx;
      for(size_t j=0;j<pzx;j++) {
         zx2[start+j] = zx[i*pzx+j];
      }
   }
   heterbart bmf(m,1);
   bmf.setprior(alpha,mybeta,tauf);
   double *ytempf = new double[2*n];  //y for h bart
   double *svecf = new double[2*n];   // sigma_i for h bart
   bmf.setdata(pzx,2*n,zx2,ytempf,zeroVp,nc);
   
   Rcpp::NumericMatrix dfburn(burnf,n); //h draws on train
   Dp::dv fhatb(n,0.0);
   if(burnf) {
   for(size_t i=0; i<n; i++) {
      ytempf[2*i] = T[i] - mTs;
      svecf[2*i] = sTs;
      ytempf[2*i+1] = T[i] - mTs;
      svecf[2*i+1] = sTs;
   }
   for(size_t i=0;i<burnf;i++) {
      if(i%printevery==0) Rprintf("burnf: done %d (out of %d)\n",i,burnf);
      bmf.draw(svecf,gen);
      for(size_t j=0;j<n;j++) dfburn(i,j) = bmf.f(2*j);
      for(size_t j=0;j<n;j++) fhatb[j] = fhatb[j] + bmf.f(2*j);
   }
   for(size_t j=0;j<n;j++) fhatb[j] = fhatb[j]/burnf;
   }
   //bmf data from starting
   double R;
   for(size_t i=0; i<n; i++) {
      ytempf[2*i] = T[i] - mTs;
      svecf[2*i] = sTs;
      if(doh1) R = (sTs/gammas)*(Y[i] - mYs - h0s[i] - h1s[i]) - T[i] + mTs;
      else R = (sYs/gammas)*(Y[i] - mYs - h0s[i]) - T[i] + mTs;
      ytempf[2*i+1] = -R;
      svecf[2*i+1] = (sTs*sYs)/gammas;
   }
   // bmf ouput storage
   Rcpp::NumericMatrix df(nd,n); //f draws on train
   Rcpp::NumericMatrix dfp(nd,nzp);


   //-------------------------------------------------
   //-------------------------------------------------
   // bart h0 setup
   //--------------------------------------------------
   heterbart bmh0(m,1);
   bmh0.setprior(alpha,mybeta,tauh0);
   double *ytemp0 = new double[n];  //y for h bart
   double *svec0 = new double[n];   // sigma_i for h bart
   bmh0.setdata(1,n,T,ytemp0,zeroVp,nc);


   //h0 burn-in
   Rcpp::NumericMatrix dh0burn(burnh0,n); //h draws on train
   double ZZ10=0.0;
   for(size_t i=0;i<burnh0;i++) {
      if(i%printevery==0) Rprintf("burnh0: done %d (out of %d)\n",i,burnh0);

      // h0 conditional -----------------------------------
      // update ----- 
      for(size_t j=0;j<n;j++) {
         ZZ10 = (T[j] - mTs - bmf.f(2*j))/sTs;
         if(doh1) ytemp0[j] = Y[j] - mYs - h1s[j] - gammas * ZZ10;
	 else ytemp0[j] = Y[j] - mYs - gammas * ZZ10;
         svec0[j] = sYs; 
      }
      //draw -----
      //bmh0.draw(svec0,gen);

      for(size_t j=0;j<n;j++) {
         dh0burn(i,j) = bmh0.f(j);
      }
   }

   //bmh0 data from starting
   for(size_t i=0;i<n;i++) {
      ZZ10 = (T[i]-mTs-fs[i])/sTs;
      if(doh1) ytemp0[i] = Y[i] - mYs - h1s[i] - gammas * ZZ10;
      else ytemp0[i] = Y[i] - mYs - gammas * ZZ10;
      svec0[i] = sYs; 
   }

   //Rcpp::NumericMatrix dh0(nd,n); //h draws on train
   //Rcpp::NumericMatrix dh0p(nd,nTp);

   //-------------------------------------------------
   //-------------------------------------------------
   // bart h1 setup
   //--------------------------------------------------
 
   heterbart bmh1(m,1);
   bmh1.setprior(alpha,mybeta,tauh1);
   double *ytemp1 = new double[n];  //y for h bart
   double *svec1 = new double[n];   // sigma_i for h bart
   bmh1.setdata(px,n,x,ytemp1,zeroVp,nc);

   //h1 burn-in
   Rcpp::NumericMatrix dh1burn(burnh1,n); //h draws on train
   double ZZ11 = 0.0;
   for(size_t i=0;i<burnh1;i++) {
      if(i%printevery==0) Rprintf("burnh1: done %d (out of %d)\n",i,burnh1);

      // h1 conditional -----------------------------------
      // update -----
      for(size_t j=0;j<n;j++) {
         ZZ11 = (T[j] - mTs - bmf.f(2*j))/sTs;
         ytemp1[j] = Y[j] - mYs -  gammas * ZZ11;
         svec1[j] = sYs; 
      }
      //draw -----
      if(doh1) bmh1.draw(svec1,gen);

 if(doh1) {
      for(size_t j=0;j<n;j++) {
        dh1burn(i,j) = bmh1.f(j);
      }
 }
   }
   

  if(doh1){
   //bmh1 data from starting
   for(size_t i=0;i<n;i++) {
      ZZ11 = (T[i]-mTs-fs[i])/sTs;
      ytemp1[i] = Y[i] - mYs - gammas * ZZ11;
      svec1[i] = sYs; 
   }
  }
  
   Rcpp::NumericMatrix dh1(nd,n); //h draws on train
   Rcpp::NumericMatrix dh1p(nd,nxp);

   //--------------------------------------------------
   //--------------------------------------------------
   // DpMuSigma setup
   
   Dp::dv itheta(5,0);
   itheta[0] = mTs; itheta[1] = mYs; itheta[2] = sTs; itheta[3] = gammas; itheta[4] = sYs; 
   if(!quiet) {
     cout << "itheta:\n";
     printdv(itheta);
   }

   //data for errors
   double *yS = new double[2*n];
   //intialize using starting values
   for(size_t i=0;i<n;i++) {
      yS[2*i] = T[i] - fs[i];
      if(doh1) yS[2*i+1] = Y[i] - h1s[i];
      else yS[2*i+1] = Y[i];
   }
   if(!quiet) cout << "check yS: " << yS[0] << ", " << yS[2*n-1] << std::endl;
   
   DpMuSigma dpmS(n,itheta,doDP);
   dpmS.setData(yS);
   dpmS.setPrior(v,nu,a);
   dpmS.setAlphaPrior(ag,priag);
   dpmS.setalpha(1.0);

   // to screen
   if(!quiet){
     cout << "dpmS.toscreen()\n";
     dpmS.toscreen();
   }

   //output storage for DpMuSigma
   Rcpp::IntegerVector dnpart(nd);
   Rcpp::NumericVector dalpha(nd);
   Rcpp::NumericMatrix dmu1(nd,n);
   Rcpp::NumericMatrix dsigma1(nd,n);
   Rcpp::NumericMatrix dmu2(nd,n);
   Rcpp::NumericMatrix dsigma2(nd,n);
   Rcpp::NumericMatrix dcov(nd,n);
   Rcpp::NumericMatrix dcor(nd,n);
   Rcpp::NumericMatrix dgamma(nd,n);
   Rcpp::NumericMatrix dsY(nd,n);
   Rcpp::NumericMatrix dLL(nd,n);
   


   bool dof = true;
   bool doh0 = false;
   bool dodp = true;

   bool doup = true;

   bool doprdx = nxp>0;
   bool doprdT = false;
   bool doprdz = nzp>0;

   double *h0p = 0;
   double *h1p = 0;
   double *fp = 0;
   if(doprdT){
     h0p = new double[nTp];
   }
   if(doprdx){
     h1p = new double[nxp];
   }
   if(doprdz){
     fp = new double[nzp];
   }
   //--------------------------------------------------
   //MCMC
   time_t tp;
   int time1 = time(&tp);
   Dp::dvv tmat;
   Dp::dv logLikelihood;
   tmat = dpmS.thetaMatrix();


   double mT=mTs,mY=mYs,sT=sTs,gamma=gammas,sY=sYs; //temporary values for mean and Lchol of Sigma
   for(size_t i=0;i<(nd+burn);i++) {
      if(i%printevery==0) Rprintf("done %d (out of %d)\n",i,nd+burn);

      // h0 conditional -----------------------------------
      if(doh0) {
      // update ----- 
      if(doup) {
      for(size_t j=0;j<n;j++) {
         mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];
         ZZ10 = (T[j] - mT - bmf.f(2*j))/sT;
         ytemp0[j] = Y[j] - mY - bmh1.f(j) - gamma * ZZ10;
         svec0[j] = sY; 
      }
      } else {
      for(size_t j=0;j<n;j++) {
         ZZ10 = (T[j] - mTs - fs[j])/sTs;
         ytemp0[j] = Y[j] - mYs - h1s[j] - gammas * ZZ10;
         svec0[j] = sYs; 
      }
      }
      //draw -----
      bmh0.draw(svec0,gen);
      }
      // h1 conditional -----------------------------------
      if(doh1) {
      // update ----- 
      if(doup) {
      for(size_t j=0;j<n;j++) {
         mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];
         ZZ11 = (T[j] - mT - bmf.f(2*j))/sT;
         if(doh0) ytemp1[j] = Y[j] - mY - bmh0.f(j) - gamma * ZZ11;
	 else ytemp1[j] = Y[j] - mY - gamma * ZZ11;
         svec1[j] = sY; 
      }
      } else {
      for(size_t j=0;j<n;j++) {
         ZZ11 = (T[j] - mTs - fs[j])/sTs;
         if(doh0) ytemp1[j] = Y[j] - mYs - h0s[j] - gammas * ZZ11;
	 else ytemp1[j] = Y[j] - mYs - gammas * ZZ11;
         svec1[j] = sYs; 
      }
      }
      //draw -----
      bmh1.draw(svec1,gen);
      }

      // f conditional --------------------------------
      if(dof) {
      // update -----
      if(doup) {
      for(size_t j=0; j<n; j++) {
         mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];
         ytempf[2*j] = T[j] - mT;
         svecf[2*j] = sT;
         if(doh1) R = (sT/gamma)*(Y[j] - mY - bmh1.f(j)) - T[j] + mT;
	 else R = (sT/gamma)*(Y[j] - mY) - T[j] + mT;
         ytempf[2*j+1] = -R;
         svecf[2*j+1] = (sT*sY)/gamma;
      }
      } else {
      for(size_t j=0; j<n; j++) {
         ytempf[2*j] = T[j] - mTs;
         svecf[2*j] = sTs;
         if(doh1) R = (sTs/gammas)*(Y[j]-mYs-h1s[j]) - T[j] + mTs;
	 else R = (sTs/gammas)*(Y[j]-mYs) - T[j] + mTs;
         ytempf[2*j+1] = -R;
         svecf[2*j+1] = (sTs*sYs)/gammas;
      }
      }
      // draws -----
      bmf.draw(svecf,gen);
      }
      
      // Sigma conditional -----------------------------
      if(dodp) {
      //update -----
      for(size_t j=0;j<n;j++) {
        yS[2*j] = T[j] - bmf.f(2*j);
        if(doh1) yS[2*j+1] = Y[j] - bmh1.f(j);
	else yS[2*j+1] = Y[j];
      }
      //draw -----
      dpmS.draw(gen);
      if(centermeans) dpmS.center();
      tmat = dpmS.thetaMatrix();
      logLikelihood = dpmS.logLikeVector();
      
      }
      if(doprdT){
	bmh0.predict(1,nTp,Tp,zeroVp,h0p);
      }
      if(doprdx){
	bmh1.predict(pxp,nxp,xp,zeroVp,h1p);
      }
      if(doprdz){
	bmf.predict(pzp,nzp,zp,zeroVp,fp);
      }
      
      if(i >= burn) {
         size_t j=i-burn;

         if(dodp) {
         dnpart[j] = dpmS.npart();
         dalpha[j] = dpmS.getalpha();
         for(size_t k=0;k<n;k++) {
            dmu1(j,k) = tmat[k][0];
            dsigma1(j,k) = tmat[k][2];
            dmu2(j,k) = tmat[k][1];
            dsigma2(j,k) = sqrt(tmat[k][3]*tmat[k][3] + tmat[k][4]*tmat[k][4]);
	    dcov(j,k) = tmat[k][3]*tmat[k][2];
	    dcor(j,k) = dcov(j,k)/(dsigma1(j,k)*dsigma2(j,k));
	    dgamma(j,k) = tmat[k][3];
	    dsY(j,k) = tmat[k][4];
	    dLL(j,k) = logLikelihood[k];
         }
         }

         /* if(doh0) { */
         /* for(size_t k=0;k<n;k++) { */
         /*   dh0(j,k) = bmh0.f(k); */
         /* } */
         /* } */

         if(doh1) {
         for(size_t k=0;k<n;k++) {
           dh1(j,k) = bmh1.f(k);
         }
         }
         
         if(dof) {
         for(size_t k=0;k<n;k++) {
            df(j,k) = bmf.f(2*k);
         }
         }

	 if(doprdz) {
	   for(size_t k=0;k<nzp;k++) {
	     dfp(j,k) = fp[k];
	   }
	 }
	 
	 if(doprdx) {
	   for(size_t k=0;k<nxp;k++) {
	     dh1p(j,k) = h1p[k];
	   }
	 }
      }
   }
   int time2 = time(&tp);
   Rprintf("time: %d\n",time2-time1);


   //--------------------------------------------------
   // clear RNG at end
   PutRNGstate();

   Rcpp::List ret;
   ret["check"] = "bnpiv";
   ret["dnpart"]=dnpart;
   ret["dalpha"]=dalpha;
   ret["dmu1"]=dmu1;
   ret["dsigma1"]=dsigma1;
   ret["dmu2"]=dmu2;
   ret["dsigma2"]=dsigma2;
   ret["dcov"]=dcov;
   ret["dcor"]=dcor;
   ret["dgamma"]=dgamma;
   ret["dsY"]=dsY;
   ret["df"] = df;
   ret["df.test"] = dfp;
   //   ret["dh0"] = dh0;
   ret["dh1"] = dh1;
   //   ret["dh0.test"] = dh0p;
   ret["dh1.test"] = dh1p;
   ret["dLL"] = dLL;
   ret["dfburn"] = dfburn;
   ret["dh0burn"] = dh0burn;
   ret["dh1burn"] = dh1burn;
   ret["dnu"] = nu;
   ret["da"] = a;
   ret["dv"] = v;


   //-------------------------------------------------
   // free
   if(yS) delete [] yS;
   //   if(yb) delete [] yb;
   //   if(xb) delete [] xb;
   if(zx2) delete [] zx2;
   if(ytempf) delete [] ytempf;
   if(svecf) delete [] svecf;
   if(ytemp0) delete [] ytemp0;
   if(svec0) delete [] svec0;
   if(ytemp1) delete [] ytemp1;
   if(svec1) delete [] svec1;
   if(h0p) delete [] h0p;
   if(h1p) delete [] h1p;
   if(fp) delete [] fp;
   
   return ret;
}
