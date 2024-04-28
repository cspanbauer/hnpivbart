/*
 *  hnpivBART: Bayesian Additive Regression Trees
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

RcppExport SEXP cbiv2s(
   SEXP _zx1,
   SEXP _zx2,
   SEXP _zx12,
   SEXP _zxp,
   SEXP _x1,
   SEXP _x2,
   SEXP _xp,
   SEXP _T1,
   SEXP _Y,
   SEXP _burn,
   SEXP _nd,
   SEXP _burnf,
   SEXP _burnh,
   SEXP _m,
   SEXP _nc,
   SEXP _power,
   SEXP _base,
   SEXP _tauf,
   SEXP _tauh,
   SEXP _betabar,
   SEXP _Abeta,
   SEXP _v,
   SEXP _nu,
   SEXP _a,
   SEXP _ag,
   SEXP _priag,
   SEXP _agB,
   SEXP _priagB,
   SEXP _centermeans,
   SEXP _fs,
   SEXP _hs,
   SEXP _betas,
   SEXP _mTs,
   SEXP _mYs,
   SEXP _sTs,
   SEXP _gammas,
   SEXP _sYs,
   SEXP _T2s,
   SEXP _printevery
)
{
   Rprintf("*****************************************************************\n");
   Rprintf("*****Into main of cbiv\n");

   //-----------------------------------------------------------
   //random number generation
   GetRNGstate();
   rrn gen;

   //--------------------------------------------------
   //process arguments

   //zx1
   Rcpp::NumericMatrix zxm1(_zx1);
   size_t nzx1 = zxm1.ncol();
   size_t pzx1 = zxm1.nrow();
   double *zx1 = &zxm1[0];

   //zx2
   Rcpp::NumericMatrix zxm2(_zx2);
   size_t nzx2 = zxm2.ncol();
   size_t pzx2 = zxm2.nrow();
   double *zx2 = &zxm2[0];

   //zx12
   Rcpp::NumericMatrix zxm12(_zx12);
   size_t nzx12 = zxm12.ncol();
   size_t pzx12 = zxm12.nrow();
   double *zx12 = &zxm12[0];
   
   //zxp
   Rcpp::NumericMatrix zpm(_zxp);
   size_t nzp = zpm.ncol();
   size_t pzp = zpm.nrow();
   double *zp = &zpm[0];

   //x1
   Rcpp::NumericMatrix xm1(_x1);
   size_t nx1 = xm1.ncol();
   size_t px1 = xm1.nrow();
   double *x1 = &xm1[0];

   //x2
   Rcpp::NumericMatrix xm2(_x2);
   size_t nx2 = xm2.ncol();
   size_t px2 = xm2.nrow();
   double *x2 = &xm2[0];
   
   //xp
   Rcpp::NumericMatrix xpm(_xp);
   size_t nxp = xpm.ncol();
   size_t pxp = xpm.nrow();
   double *xp = &xpm[0];
   
   // T1,Y ------------
   Rcpp::NumericVector Tv1(_T1);
   double *T1 = &Tv1[0];
   size_t nT1 = Tv1.size();

   Rcpp::NumericVector Yv(_Y);
   double *Y = &Yv[0];
   size_t nY = Yv.size();

   // burn, nd
   size_t burn = Rcpp::as<size_t>(_burn);
   size_t nd = Rcpp::as<size_t>(_nd);
   size_t burnf = Rcpp::as<size_t>(_burnf);
   size_t burnh = Rcpp::as<size_t>(_burnh);

   //bart prior
   size_t m = Rcpp::as<size_t>(_m);
   size_t nc = Rcpp::as<size_t>(_nc);
   double mybeta = Rcpp::as<double>(_power);
   double alpha = Rcpp::as<double>(_base);
   double tauf = Rcpp::as<double>(_tauf);
   double tauh = Rcpp::as<double>(_tauh);

   //beta prior
   double betabar = Rcpp::as<double>(_betabar);
   double Abeta = Rcpp::as<double>(_Abeta);

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
   for(size_t i=0;i<nag;i++) { // should use STL */
       ag[i]=agv[i]; 
       priag[i]=priagv[i]; 
   }
   Rcpp::NumericVector agBv(_agB);
   Rcpp::NumericVector priagBv(_priagB);
   size_t nagB=agBv.size();
   Dp::dv agB(nagB);
   Dp::dv priagB(nagB);
   for(size_t i=0;i<nagB;i++) { // should use STL */
       agB[i]=agBv[i]; 
       priagB[i]=priagBv[i]; 
   }

   //should the means be centered after each DpMuSigma draw
   bool centermeans = Rcpp::as<bool>(_centermeans);

   //starting values
   //fs,hs
   Rcpp::NumericVector fsv(_fs);
   double *fs = &fsv[0];
   Rcpp::NumericVector hsv(_hs);
   double *hs = &hsv[0];
   Rcpp::NumericVector T2sv(_T2s);
   double *T2s = &T2sv[0];
   Rcpp::NumericVector T2v(_T2s);
   double *T2 = &T2v[0];
   Rcpp::NumericVector T1sv(_T2s);
   double *T1s = &T1sv[0];
   
   
   
   double betas = Rcpp::as<double>(_betas);
   double mTs = Rcpp::as<double>(_mTs);
   double mYs = Rcpp::as<double>(_mYs);
   double sTs = Rcpp::as<double>(_sTs);
   double gammas = Rcpp::as<double>(_gammas);
   double sYs = Rcpp::as<double>(_sYs);


   
   //other
   size_t printevery = Rcpp::as<size_t>(_printevery);

   size_t n1 = nT1;
   size_t n2 = nY;

   //Zero vector for first BART model
   Rcpp::NumericVector zeroVp2(n2);
   Rcpp::NumericVector zeroVp12(n1+n2);
   double *zeroV2 = &zeroVp2[0];
   double *zeroV12 = &zeroVp12[0];
   
   //--------------------------------------------------
   // print args

   Rprintf("***burn, nd, burnf, burnh: %ld, %ld, %ld\n",burn,nd,burnf,burnh);

   Rprintf("*** Prior:\n");
   Rprintf("m (num trees), nc (num cut points), %ld, %ld\n",m,nc);
   Rprintf("****power, base, tauf, tauh: %lf, %lf, %lf, %lf\n",mybeta,alpha,tauf,tauh);
   Rprintf("betabar: %lf\n",betabar);
   Rprintf("Abeta: %lf\n",Abeta);
   Rprintf("v: %lf\n",v);
   Rprintf("nu: %lf\n",nu);
   Rprintf("a: %lf\n",a);
   Rprintf("alpha prior: grid size: %ld\n",ag.size()); 
   Rprintf("\tag, first, last: %lf, %lf\n",ag[0],ag[nag-1]); 
   Rprintf("\tpriag, first, last: %lf, %lf\n",priag[0],priag[nag-1]); 

   Rprintf("*** Data:\n");
   Rprintf("nzx1,nzx2,pzx1,pzx2: %ld, %ld, %ld, %ld\n",nzx1,nzx2,pzx1,pzx2);
   Rprintf("nx1,nx2,px1,px2: %ld, %ld, %ld, %ld\n",nx1,nx2,px1,px2);
   Rprintf("nT1: %ld\n",nT1);
   Rprintf("first and last T1: %lf, %lf\n",T1[0],T1[nT1-1]);
   Rprintf("nY: %ld\n",nY);
   Rprintf("first and last Y: %lf, %lf\n",Y[0],Y[nY-1]);

   Rprintf("*** starting values, n1 and n2 are: %ld and %ld\n",n1,n2);
   Rprintf("\t first and last fs: %lf, %lf\n",fs[0],fs[n2-1]);
   Rprintf("\t first and last hs: %lf, %lf\n",hs[0],hs[n2-1]);
   Rprintf(" starting for beta, mT, mY, sT, gamma, sY: %lf, %lf, %lf, %lf, %lf, %lf\n",betas,mTs,mYs,sTs,gammas,sYs);
   //   Rprintf(" starting for T2: %lf\n",T2s[0]);
   
   Rprintf("***other\n");
   Rprintf("printevery: %ld\n",printevery);


   //-------------------------------------------------
   //-------------------------------------------------
   // bart f setup
   //--------------------------------------------------
   // first make zx2 which is zx with the columns repeated

   heterbart bmf(m,1);
   bmf.setprior(alpha,mybeta,tauf);
   double *ytempf = new double[n1+n2];  //y for h bart
   double *svecf = new double[n1+n2];   // sigma_i for h bart
   
   bmf.setdata(pzx12,n1+n2,zx12,ytempf,zeroV12,nc);
   Rcpp::NumericMatrix dfburn(burnf,n1+n2); //h draws on train
   Dp::dv fhatb(n2,0.0);
   if(burnf) {
   for(size_t i=0; i<n1; i++) {
      ytempf[i] = T1s[i] - mTs;
      svecf[i] = sTs;
   }
   for(size_t i=0; i<n2; i++) {
     ytempf[n1+i] = T2s[i] - mTs;
     svecf[n1+i] = sTs;
   }
   for(size_t i=0;i<burnf;i++) {
      if(i%printevery==0) Rprintf("burnf: done %d (out of %d)\n",i,burnf);
      bmf.draw(svecf,gen);
      for(size_t j=0;j<n1+n2;j++) dfburn(i,j) = bmf.f(j);
      for(size_t j=n1+1;j<n1+n2;j++) fhatb[j-n1] = fhatb[j-n1] + bmf.f(j);
   }
   for(size_t j=0;j<n2;j++) fhatb[j] = fhatb[j]/burnf;
   }
   //bmf data from starting
   double R;
   for(size_t i=0; i<n1; i++) {
      ytempf[i] = T1s[i] - mTs;
      svecf[i] = sTs;
   }
   for(size_t i=0; i<n2; i++) {
      R = (betas*sTs + gammas)*(T2s[i-n1] - mTs) - sTs*(Y[i-n1] - mYs - hs[i-n1] - betas*mTs);
      ytempf[n1+i] = R/gammas;
      svecf[n1+i] = (sTs*sYs)/gammas;
   }
   // bmf ouput storage
   Rcpp::NumericMatrix df(nd,n2); //f draws on train
   Rcpp::NumericMatrix dfp(nd,nzp);

   //-------------------------------------------------
   //-------------------------------------------------
   // bart h setup
   //--------------------------------------------------
   heterbart bmh(m,1);
   bmh.setprior(alpha,mybeta,tauh);
   double *ytemp = new double[n2];  //y for h bart
   double *svec = new double[n2];   // sigma_i for h bart
   bmh.setdata(px2,n2,x2,ytemp,zeroV2,nc);


   //h burn-in
   Rcpp::NumericMatrix dhburn(burnh,n2); //h draws on train
   double ZZr=0.0;
   for(size_t i=0;i<burnh;i++) {
      if(i%printevery==0) Rprintf("burnh: done %d (out of %d)\n",i,burnh);

      // h conditional -----------------------------------
      // update ----- 
      for(size_t j=0;j<n2;j++) {
         ZZr = (T2s[j] - mTs - bmf.f(n1+j))/sTs;
         ytemp[j] = Y[j] - mYs - betas*T2s[j] - gammas * ZZr;
         svec[j] = sYs; 
      }
      //draw -----
      bmh.draw(svec,gen);

      for(size_t j=0;j<n2;j++) {
         dhburn(i,j) = bmh.f(j);
      }
   }

   //bmh data from starting
   for(size_t i=0;i<n2;i++) {
      ZZr = (T2s[i]-mTs-fs[i])/sTs;
      ytemp[i] = Y[i] - mYs - betas*T2s[i] - gammas * ZZr;
      svec[i] = sYs; 
   }


   Rcpp::NumericMatrix dh(nd,n2); //h draws on train
   Rcpp::NumericMatrix dhp(nd,nxp);

   //--------------------------------------------------
   //--------------------------------------------------
   // beta setup
   double *yb = new double[n2];
   double *xb = new double[n2];
   double Zr=0.0;
   for(size_t i=0;i<n2;i++) {
      Zr = (T2s[i]-mTs-fs[i])/sTs;
      yb[i] = (Y[i] - mYs - hs[i] - gammas * Zr)/sYs;
      xb[i] = T2s[i]/sYs;
   }
   double betad = betas;
   Rcpp::NumericVector dbeta(nd); //storage for output beta draws

   //--------------------------------------------------
   //--------------------------------------------------
   // DpMuSigma setup for (imputed) stage 1 and 2
   cout << "\n*** making a DpMuSigma object:\n";
   
   Dp::dv itheta(5,0); 
   itheta[0] = mTs; itheta[1] = mYs; itheta[2] = sTs; itheta[3] = gammas; itheta[4] = sYs; 
   cout << "itheta:\n"; 
   printdv(itheta); 

   //data for errors
   double *yS = new double[2*n2];
   //intialize using starting values
   for(size_t i=0;i<n2;i++) {
      yS[2*i] = T2s[i] - fs[i];
      yS[2*i+1] = Y[i] - betas*T2s[i] - hs[i];
   }
   cout << "check yS: " << yS[0] << ", " << yS[2*n2-1] << std::endl;
   cout << "Hi0\n";
   DpMuSigma dpmS(n2,itheta,true);
   dpmS.setData(yS);
   dpmS.setPrior(v,nu,a);
   dpmS.setAlphaPrior(ag,priag);
   dpmS.setalpha(1.0);

   //-------------------------------------------------
   //-------------------------------------------------
   // DpMuSigma setup for stage 1

   Dp::dv ithetaB(1,0);
   ithetaB[0]=sTs;
   printdv(ithetaB);

   //data for errors
   double *ySB = new double[n1];
   for(size_t i=0;i<n1;i++) {
     ySB[i]=T1s[i]-fs[i];
   }

   DpMuTau dpmSB(n1,ithetaB,true);
   dpmSB.setData(ySB);
   //   dpmSB.setPrior(0.0,a,1.0,vnu,a);
   dpmSB.setAlphaPrior(agB,priagB);
   dpmSB.setalpha(1.0);
   
   for(size_t j=0;j<n2;j++)
     T2[j]=sTs*gen.normal()+mTs+bmf.f(n1+j);

   cout << "Before output setup\n";
   /* // to screen */
   /* cout << "dpmS.toscreen()\n"; */
   /* dpmS.toscreen(); */
   //output storage for DpMuSigma
   Rcpp::IntegerVector dnpart(nd);
   cout << "A\n";
   Rcpp::NumericVector dalpha(nd);
   cout << "A\n";
   Rcpp::NumericMatrix dmu1(nd,n2);
   cout << "A\n";
   Rcpp::NumericMatrix dsigma1(nd,n2);
   cout << "A\n";
   Rcpp::NumericMatrix dmu2(nd,n2);
   cout << "A\n";
   Rcpp::NumericMatrix dsigma2(nd,n2);
   cout << "A\n";
   Rcpp::NumericMatrix dgamma(nd,n2);
   cout << "A\n";
   Rcpp::NumericMatrix dcov(nd,n2);
   cout << "A\n";
   Rcpp::NumericMatrix dcor(nd,n2);
   cout << "A\n";
   Rcpp::NumericMatrix dsY(nd,n2);

   cout << "After output setup\n";
   
   bool dobeta = true;
   bool dof = true;
   bool doh = true;
   bool dodp = true;

   bool doup = true;
   cout << nxp << '\n';
   bool doprdx = nxp>0;
   bool doprdz = nzp>0;

   double *hp = 0;
   double *fp = 0;
   if(doprdx)
     {
       hp = new double[nxp];
     }
   if(doprdz)
     {
       fp = new double[nzp];
     }
   cout << "After prediction setup\n";
   //--------------------------------------------------
   //MCMC
   time_t tp;
   int time1 = time(&tp);
   Dp::dvv tmat; 
   tmat = dpmS.thetaMatrix();
   Dp::dvv tmat2;
   tmat2 = dpmSB.thetaMatrix();
   cout << tmat2[1][0] << '\n';

   cout << "Before MCMC\n";
   double mT=mTs,mY=mYs,sT=sTs,gamma=gammas,sY=sYs; //temporary values for mean and Lchol of Sigma
   for(size_t i=0;i<(nd+burn);i++) {
      if(i%printevery==0) Rprintf("done %d (out of %d)\n",i,nd+burn);
      for(size_t j=0;j<n2;j++)
	T2[j]=tmat2[j][0]*gen.normal()+0.0+bmf.f(n1+j);

      // beta conditional -----------------------------------
      if(dobeta) {
      // update -----
      if(doup) {
      for(size_t j=0;j<n2;j++) {
         mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];
         Zr = (T2[j] - mT - bmf.f(n1+j))/sT;
         yb[j] = (Y[j] - mY - bmh.f(j) - gamma * Zr)/sY;
         xb[j] = T2[j]/sY;
      }
      } else {
      for(size_t j=0;j<n2;j++) {
         Zr = (T2s[j] - mTs - fs[j])/sTs;
         yb[j] = (Y[j] - mYs - hs[j] - gammas * Zr)/sYs;
         xb[j] = T2s[j]/sYs;
      }
      }
      //draw -----
      betad = lreg(n2,xb,yb,1.0,betabar,Abeta,gen);
      }

      // h conditional -----------------------------------
      if(doh) {
      // update ----- 
      if(doup) {
      for(size_t j=0;j<n2;j++) {
         mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];
         ZZr = (T2[j] - mT - bmf.f(n1+j))/sT;
         ytemp[j] = Y[j] - mY - betad*T2[j] - gamma * ZZr;
         svec[j] = sY;
      }
      } else {
      for(size_t j=0;j<n2;j++) {
         ZZr = (T2s[j] - mTs - fs[j])/sTs;
         ytemp[j] = Y[j] - mYs - betas*T2[j] - gammas * ZZr;
         svec[j] = sYs; 
      }
      }
      //draw -----
      bmh.draw(svec,gen);
      }

      // f conditional --------------------------------
      if(dof) {
      // update -----
      if(doup) {
	for(size_t j=0; j<n1; j++) {
	  ytempf[j] = T1[j] - 0.0;
	  svecf[j] = tmat2[j][0];
	}
	for(size_t j=0; j<n2; j++) {
	  mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];
	  R = (betad*sT+gamma)*(T2[j]-mT) - sT*(Y[j] - mY - bmh.f(j) - betad*mT);
	  ytempf[n1+j] = R/gamma;
	  svecf[n1+j] = (sT*sY)/gamma;
	}
      } else {
      for(size_t j=0; j<n1; j++) {
         ytempf[j] = T1[j] - 0.0;
         svecf[j] = sTs;
      }
      for(size_t j=0; j<n2; j++) {
         R = (betas*sTs+gammas)*(T2[j]-mTs) - sTs*(Y[j] - mYs - hs[j] - betas*mTs);
         ytempf[n1+j] = R/gammas;
         svecf[n1+j] = (sTs*sYs)/gammas;
      }
      }
      // draws -----
      bmf.draw(svecf,gen);
      //      cout << "1: " << bmf.f(0) << '\n' << "2: " << bmf.f(2) << '\n' << "3: " << bmf.f(4) << '\n';
      }

      // Sigma conditional -----------------------------
      //cout << "\n&&&&&&about to do Sigma\n";
      if(dodp) {
      //update -----
      for(size_t j=0;j<n2;j++) {
         yS[2*j] = T2[j] - bmf.f(n1+j);
         yS[2*j+1] = Y[j] - betad*T2[j] - bmh.f(j);
      }
      //draw -----
      dpmS.draw(gen); 
      if(centermeans) dpmS.center(); 
      tmat = dpmS.thetaMatrix();
      for(size_t j=0;j<n2;j++){
	tmat[j][0]=0.; tmat[j][1]=0.;
      }

      if(dodp){
	for(size_t j=0;j<n1;j++) {
	  ySB[j] = T1[j]-bmf.f(j);
	}
      }
      dpmSB.draw(gen);
      //      if(centermeans) dpmSB.center();
      tmat2 = dpmS.thetaMatrix();

	//cout << "h(0): " << bmh.f(0) << '\n'; */
      }

      if(doprdx)
	bmh.predict(pxp,nxp,xp,T2,hp);
      if(doprdz)
	bmf.predict(pzp,nzp,zp,zeroV12,fp);
      
      if(i >= burn) {
         size_t j=i-burn;
	 if(dodp) {
	   dnpart[j] = dpmS.npart();
	   dalpha[j] = dpmS.getalpha();
         
           for(size_t k=0;k<n2;k++) {
             dmu1(j,k) = tmat[k][0];
             dsigma1(j,k) = tmat[k][2];
             dmu2(j,k) = tmat[k][1];
             dsigma2(j,k) = sqrt(tmat[k][3]*tmat[k][3] + tmat[k][4]*tmat[k][4]);
	     dgamma(j,k) = tmat[k][3];
	     dcov(j,k) = tmat[k][3]*tmat[k][2];
	     dcor(j,k) = dcov(j,k)/(dsigma1(j,k)*dsigma2(j,k));
	     dsY(j,k) = tmat[k][4];
	     df(j,k) = bmf.f(n1+k);
           }
	 }

         if(dobeta) {dbeta[j] = betad;}

         if(doh) {
	   for(size_t k=0;k<n2;k++) {
            dh(j,k) = bmh.f(k);
	   }
         }
	 if(doprdx) {
	   for(size_t k=0;k<nxp;k++) {
	     dhp(j,k) = hp[k];
	   }
	 }
	 if(doprdz) {
	   for(size_t k=0;k<nzp;k++) {
	     dfp(j,k) = fp[k];
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
   ret["check"] = "biv";
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
   ret["dbeta"] = dbeta;
   ret["dh"] = dh;
   ret["dh.test"] = dhp;
   ret["df"] = df;
   ret["df.test"] = dfp;
   ret["dfburn"] = dfburn;
   ret["dhburn"] = dhburn;
   ret["dnu"] = nu;
   ret["da"] = a;
   ret["dv"] = v;

   
   //-------------------------------------------------
   // free
   if(yS) delete [] yS;
if(ySB) delete [] ySB;
   if(yb) delete [] yb;
   if(xb) delete [] xb;
   if(ytemp) delete [] ytemp;
   if(svec) delete [] svec;
   if(zx2) delete [] zx2;
   if(ytempf) delete [] ytempf;
   if(svecf) delete [] svecf;
   if(hp) delete [] hp;
   if(fp) delete [] fp;

   return ret;
}
