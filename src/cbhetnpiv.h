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

RcppExport SEXP cbhetnpiv(
   SEXP _z,
   SEXP _zp,
   SEXP _T,
   SEXP _TX,
   SEXP _TXp,
   SEXP _Y,
   SEXP _type1,
   SEXP _type2,
   SEXP _burn,
   SEXP _nd,
   SEXP _keepevery,
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
   SEXP _offset,
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

   //z
   Rcpp::NumericMatrix zpm(_zp);
   size_t nzp = zpm.ncol();
   size_t pzp = zpm.nrow();
   double *zp = &zpm[0];

   // TX ------------
   Rcpp::NumericMatrix TXm(_TX);
   size_t nTx = TXm.ncol();
   size_t pTx = TXm.nrow();
   double *Tx = &TXm[0];

   //TXp
   Rcpp::NumericMatrix TXpm(_TXp);
   size_t nTXp = TXpm.ncol();
   size_t pTXp = TXpm.nrow();
   double *TXp = &TXpm[0];

   // T,Y ------------
   Rcpp::NumericVector Tv(_T);
   double *T = &Tv[0];
   size_t nT = Tv.size();
   Rcpp::NumericVector Yv(_Y);
   double *Y = &Yv[0];
   size_t nY = Yv.size();

   // type2 ----------
   size_t type1 = Rcpp::as<size_t>(_type1);
   size_t type2 = Rcpp::as<size_t>(_type2);
   // burn, nd
   size_t burn = Rcpp::as<size_t>(_burn);
   size_t nd = Rcpp::as<size_t>(_nd);
   size_t ke = Rcpp::as<size_t>(_keepevery);
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
   double offset = Rcpp::as<double>(_offset);
   size_t include_output = 0;//Rcpp::as<size_t>(_include_output);
   bool doDP = Rcpp::as<bool>(_doDP);

   //starting values
   //fs,hs
   Rcpp::NumericVector fsv(_fs);
   double *fs = &fsv[0];
   Rcpp::NumericVector h0sv(_hs);
   double *hs = &h0sv[0];

   double mTs = Rcpp::as<double>(_mTs);
   double mYs = Rcpp::as<double>(_mYs);
   double sTs = Rcpp::as<double>(_sTs);
   double gammas = Rcpp::as<double>(_gammas);
   double sYs = Rcpp::as<double>(_sYs);

   //null model bools
   //nullf, nullh
   bool nullf = Rcpp::as<bool>(_nullf);
   bool nullh = Rcpp::as<bool>(_nullh);

   //other
   size_t printevery = Rcpp::as<size_t>(_printevery);
   bool quiet = Rcpp::as<bool>(_quiet);

   size_t n = nTx;
   
   if(!quiet){
     Rprintf("*****************************************************************\n");
     Rprintf("*****Into main of cbhetnpiv\n");
   }

   
   //--------------------------------------------------
   // print args
   if(!quiet){
     Rprintf("***burn, nd, burnf, burnh0, burnh1: %ld, %ld, %ld\n",burn,nd,burnf,burnh);
     Rprintf("***offset: %lf\n",offset);

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
     Rprintf("nTx: %ld\n",nTx);
     //Rprintf("first and last T: %lf, %lf\n",Tx[0],T[nT-1]);
     Rprintf("nY: %ld\n",nY);
     Rprintf("first and last Y: %lf, %lf\n",Y[0],Y[nY-1]);
     
     Rprintf("*** starting values, n is: %ld\n",n);
     Rprintf("\t first and last fs: %lf, %lf\n",fs[0],fs[n-1]);
     Rprintf("\t first and last hs: %lf, %lf\n",hs[0],hs[n-1]);
     Rprintf(" starting for mT, mY, sT, gamma, sY: %lf, %lf, %lf, %lf, %lf\n",mTs,mYs,sTs,gammas,sYs);
     
     Rprintf("***other\n");
     Rprintf("printevery: %ld\n",printevery);
   }
   
   double *YY = new double[n];
   for(size_t j=0;j<n;j++){
     if(type2==1) YY[j] = Y[j]-offset; // continuous Y
     else if(type2==2) { // binary Y with probit link
       double sign;
       if(Y[j]==1) sign=1.;
       else sign=-1.;
       YY[j] = sign*rtnorm(sign*(hs[j]+(T[j]-fs[j])*gammas/sTs), -sign*offset, sqrt(1.-gammas*gammas), gen);
     }
   }

   //-------------------------------------------------
   //-------------------------------------------------
   // bart f1 setup
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
   heterbart<double> bmf(m1);
   bmf.setprior(alpha,mybeta,tauf);
   double *ytempf = new double[2*n];  //y for h bart
   double *svecf = new double[2*n];   // sigma_i for h bart
   bmf.setdata(pz,2*n,z2,ytempf,nc);
   for(size_t i=0; i<n; i++) {
     ytempf[2*i] = T[i] - mTs;
     svecf[2*i] = sTs;
     ytempf[2*i+1] = T[i] - mTs;
     svecf[2*i+1] = sTs;
   }
   bmf.draw(svecf,gen);
   // bmf ouput storage
   Rcpp::NumericMatrix df(nd,n); //f draws on train
   Rcpp::NumericMatrix dfp(nd,nzp);

   //-------------------------------------------------
   //-------------------------------------------------
   // bart h setup
   //--------------------------------------------------
   heterbart<double> bmh(m2);
   bmh.setprior(alpha,mybeta,tauh);
   double *ytemp = new double[n];  //y for h bart
   double *svec = new double[n];   // sigma_i for h bart
   bmh.setdata(pTx,nTx,Tx,ytemp,nc);
   //h burn-in
   Rcpp::NumericMatrix dhburn(1,1); //h draws on train
   if(include_output==1) dhburn(nd,n);
   double ZZ10=0.0;
   //  if(!nullh){
      for(size_t j=0;j<n;j++) {
        //else { } // binary Y with logistic link, not implemented yet
        ytemp[j] = YY[j] - mYs - gammas * ZZ10;
        svec[j] = sYs; 
      }
      //  }
   bmh.draw(svec,gen);
   Rcpp::NumericMatrix dh(nd,n); //h draws on train
   //   cout << nTXp << '\n';
   //cout << nTXn << '\n';
   Rcpp::NumericMatrix dhp(nd,nTXp);

   //--------------------------------------------------
   //--------------------------------------------------
   // DpMuSigma setup
   
   Dp::dv itheta(5,0);
   itheta[0] = mTs; itheta[1] = mYs; itheta[2] = sTs; itheta[3] = gammas; itheta[4] = sYs; 

   //data for errors
   double *yS = new double[2*n];
   //intialize using starting values
   // cout << "AAA\n";
   //cout << betas << '\n';
   for(size_t i=0;i<n;i++) {
     yS[2*i] = T[i] - fs[i];
     yS[2*i+1] = YY[i] - hs[i];
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
   Rcpp::NumericVector dcovMean(n);
   Rcpp::NumericVector dcorMean(n);
   Rcpp::NumericVector dcovSampMean(nd);
   Rcpp::NumericVector dcorSampMean(nd);   
   //   Rcpp::NumericMatrix dclusters(nd,n);
   Rcpp::NumericVector dsT(nd);
   Rcpp::NumericMatrix dsigma1(1,1);
   Rcpp::NumericMatrix dsigma2(1,1);
   Rcpp::NumericMatrix dcov(1,1);   
   Rcpp::NumericMatrix dcor(1,1);   
   Rcpp::NumericVector dgamma(nd); 
   Rcpp::NumericMatrix dsY(1,1);    
   Rcpp::NumericMatrix dLL(1,1);       
   if(include_output==1) {
      dsigma1(nd,n);
      dsigma2(nd,n);
      dcov(nd,n);
      dcor(nd,n);
      dgamma(nd,n);
      dsY(nd,n);
      // dLL(nd,n);
   }
   bool doprdx = nTXp>0;
   bool doprdz = nzp>0;
   
   double *hp = 0;
   if(doprdx){
     hp = new double[nTXp];
   }
   double *fp = 0;
   if(doprdz){
     fp = new double[nzp];
   }

   //--------------------------------------------------
   //MCMC
   time_t tp;
   int time1 = time(&tp);
   Dp::dvv tmat;
   dpmS.draw(gen);
   if(centermeans) dpmS.center();
   tmat = dpmS.thetaMatrix();
   double mT=mTs,mY=mYs,sT=sTs,gamma=gammas,sY=sYs; //temporary values for mean and Lchol of Sigma
   for(size_t i=0;i<(nd*ke+burn);i++) {
     //     cout << "i: " << i << '\n';
      if(i%printevery==0&!quiet) Rprintf("done %d (out of %d)\n",i,nd*ke+burn);
      if(type2==2) { // probit link draw auxiliary
        double sign;
        for(size_t j=0;j<n;j++){
          if(i>0) {mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];}        
          if(Y[j]==1) sign=1.;
          else sign=-1.;
          YY[j] = sign*rtnorm(sign*(bmh.f(j)+(T[j]-bmf.f(2*j))*gamma/sT), -sign*offset, sqrt(1.-gamma*gamma), gen);
          //   if(j==0) cout << "YY[0]: " << YY[0] << '\n';
        }
      }
      // h conditional -----------------------------------
      //   if(!nullh) {
      // update ----- 
      for(size_t j=0;j<n;j++) {
        if(i>0) {mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];}
         if(!nullf) ZZ10 = (T[j] - mT - bmf.f(2*j))/sT;
         else ZZ10 = 0.;
         ytemp[j] = YY[j] - mY - gamma * ZZ10;
         svec[j] = sY; 
      }
      //draw -----
      bmh.draw(svec,gen);
      //   }
      //     cout << "f2: " << bmh.f(0) << '\n';
      
      // f conditional --------------------------------
      if(!nullf){
        double R;
        for(size_t j=0; j<n; j++) {
          if(i>0){mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];}
          ytempf[2*j] = T[j] - mT;
          svecf[2*j] = sT;
          R = (sT/gamma)*(YY[j] - mY - bmh.f(j)) - T[j] + mT;
          ytempf[2*j+1] = -R;
          svecf[2*j+1] = (sT*sY)/gamma;
        }
        // draws -----
        bmf.draw(svecf,gen);
	  }
      //    cout << "f1: " << bmf.f(0) << '\n';

      if(type2==1){
        // Sigma conditional -----------------------------
        //update -----
        for(size_t j=0;j<n;j++) {
          yS[2*j] = T[j] - bmf.f(2*j);
          yS[2*j+1] = YY[j] - bmh.f(j);
        }
        //draw -----
        dpmS.draw(gen); 
        if(centermeans) dpmS.center(); 
        tmat = dpmS.thetaMatrix();
      }
      else{
        //update -----
        for(size_t j=0;j<n;j++) {
          yS[2*j] = T[j] - bmf.f(2*j);
          // cout << yS[2*j] << '\n';
          yS[2*j+1] = YY[j] - bmh.f(j);
        }
        // Sufficient statistics
        double sumY1sq=0.0,sumY2sq=0.0,sumY1Y2=0.0;
        for(size_t j=0;j<n;j++) {
          sumY1sq += yS[2*j]*yS[2*j];
          sumY2sq += yS[2*j+1]*yS[2*j+1];
          sumY1Y2 += yS[2*j]*yS[2*j+1];
        }   
        // Draw rho
        double sig;
        if(i>0) sig=tmat[0][2];
        else sig=sTs;
        size_t rho_grid_size=2000;
        Rcpp::NumericVector rho_grid(rho_grid_size-1);
        Rcpp::NumericVector log_post_rho(rho_grid_size-1);
        Rcpp::NumericVector prob_rho(rho_grid_size-1);
        for(size_t g=0;g<rho_grid_size-1;g++) {
          rho_grid[g]=0.001*((double)(g-999.));
        }
        size_t rho_index;
        double max_log_post_rho=log_rho_posterior(rho_grid[0],sumY1sq,sumY2sq,sumY1Y2,sig,1,1,n);
        for(size_t g=0;g<rho_grid_size-1;g++) {
          log_post_rho[g]=log_rho_posterior(rho_grid[g],sumY1sq,sumY2sq,sumY1Y2,sig,1.,1.,n);
          if(log_post_rho[g]>max_log_post_rho) max_log_post_rho=log_post_rho[g];
        }
        double sum_exp_max=0.;
        for(size_t g=0;g<rho_grid_size-1;g++){
          sum_exp_max += exp(log_post_rho[g]-max_log_post_rho);
        }
        double lse=max_log_post_rho+log(sum_exp_max);
        for(size_t g=0;g<rho_grid_size-1;g++){
          prob_rho[g] = exp(log_post_rho[g]-lse);
        }
        rho_index = 0;
        double uu = gen.uniform();
        double sum_prob = prob_rho[0];
        //      cout << "UU: " << uu << '\n';
        //        cout << "SP: " << sum_prob << '\n';
        while(uu>sum_prob){
          rho_index++;
          sum_prob += prob_rho[rho_index];
        }
        //      cout << "Rho index: " << rho_index << '\n';
        // Draw sigma_t
        size_t sig_grid_size=1000;
        Rcpp::NumericVector sig_grid(sig_grid_size-1);
        Rcpp::NumericVector log_post_sig(sig_grid_size-1);
        Rcpp::NumericVector prob_sig(sig_grid_size-1);
        for(size_t g=0;g<sig_grid_size-1;g++) {
          sig_grid[g]=0.01*((double)g+1.1);
          //    cout << sig_grid[g] << '\n';
        }
        double max_log_post_sig=log_sig_posterior(sig_grid[0],sumY1sq,sumY2sq,sumY1Y2,rho_grid[rho_index],0.,sqrt(10.),n);
        //       cout << log_sig_posterior(sig_grid[50],sumY1sq,sumY2sq,sumY1Y2,rho_grid[rho_index],0.,sqrt(10.),n) << '\n';
        //       cout << sumY1sq << '\n' << sumY2sq << '\n' << sumY1Y2 << '\n';
        for(size_t g=0;g<sig_grid_size-1;g++) {
          log_post_sig[g]=log_sig_posterior(sig_grid[g],sumY1sq,sumY2sq,sumY1Y2,rho_grid[rho_index],0.,sqrt(10.),n);
          if(log_post_sig[g]>max_log_post_sig) max_log_post_sig=log_post_sig[g];
        }
        //     cout << max_log_post_sig << '\n';
        sum_exp_max=0.;
        for(size_t g=0;g<sig_grid_size-1;g++){
          sum_exp_max += exp(log_post_sig[g]-max_log_post_sig);
        }
        lse=max_log_post_sig+log(sum_exp_max);
        for(size_t g=0;g<sig_grid_size-1;g++){
          prob_sig[g] = exp(log_post_sig[g]-lse);
          //  cout << prob_sig[g] << '\n';
        }
        size_t sig_index=0;
        uu = gen.uniform();
        sum_prob=prob_sig[0];
        while(uu>sum_prob){
          sig_index++;
          sum_prob += prob_sig[sig_index];
          //       cout << sig_index << '\n';
        }
        for(size_t j=0;j<n;j++){
          tmat[j][2]=sig_grid[sig_index];
          //     cout << sig_index << '\n';
          tmat[j][3]=rho_grid[rho_index];
          if(tmat[j][3]==0.) tmat[j][3]=.001;
          tmat[j][4]=sqrt(1.-tmat[j][3]*tmat[j][3]);
        }
        //             cout << "Chol(sigma):\n" << tmat[0][2] << '\n' << tmat[0][3] << '\n' << tmat[0][4] << '\n';
      }

      if(doprdx){
        bmh.predict(pTXp,nTXp,TXp,hp);
        //      cout << hp[0] << '\n';
      }

      if(doprdz){
        bmf.predict(pzp,nzp,zp,fp);
      }

      
      if(i >= burn) {
        size_t j=i-burn;
        if(j % ke == 0){
          size_t draw=j/ke;

          dsT[draw] = tmat[0][0];
          dgamma[draw] = tmat[0][3];
          dnpart[draw] = dpmS.npart();
          dalpha[draw] = dpmS.getalpha();
          for(size_t k=0;k<n;k++) {
            dcovMean(k) += (1./(double)nd)*tmat[k][3]*tmat[k][2];
            dcorMean(k) += (1./(double)nd)*(tmat[k][3]*tmat[k][2])/(tmat[k][2]*sqrt(tmat[k][3]*tmat[k][3]+tmat[k][4]*tmat[k][4]));
              dcovSampMean(draw) += (1./(double)n)*tmat[k][3]*tmat[k][2];
              dcorSampMean(draw) += (1./(double)n)*(tmat[k][3]*tmat[k][2])/(tmat[k][2]*sqrt(tmat[k][3]*tmat[k][3]+tmat[k][4]*tmat[k][4]));
              if(include_output==1){
                dsigma1(draw,k) = tmat[k][2];
                dsigma2(draw,k) = sqrt(tmat[k][3]*tmat[k][3] + tmat[k][4]*tmat[k][4]);
                dcov(draw,k) = tmat[k][3]*tmat[k][2];
                dcor(draw,k) = dcov(draw,k)/(dsigma1(draw,k)*dsigma2(draw,k));
                //   dgamma(draw,k) = tmat[k][3];
                dsY(draw,k) = tmat[k][4];
                //  dLL(draw,k) = bvnorm(T[k],Y[k],bmf.f(2*k),bmh.f(k),dsigma1(draw,k),dsigma2(draw,k),dcov(draw,k));
              }
            }
         
              for(size_t k=0;k<n;k++) {
                dh(draw,k) = bmh.f(k)+offset;
              }
           
              for(size_t k=0;k<n;k++) {
                df(draw,k) = bmf.f(2*k);
              }
            
            if(doprdx){
              for(size_t k=0;k<nTXp;k++) {          
                dhp(draw,k) = hp[k]+offset;
                //                if(k==0) cout << hp[0] << '\n';
              }
            }
         
            if(doprdz){
              for(size_t k=0;k<nzp;k++) {
                dfp(draw,k) = fp[k];
              }
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
   ret["check"] = "bhetnpiv";
   ret["dnpart"]=dnpart;
   ret["dalpha"]=dalpha;
   ret["dgamma"]=dgamma;
   ret["dsT"]=dsT;
   if(include_output==1){
     ret["dsigma1"]=dsigma1;
     ret["dsigma2"]=dsigma2;
     ret["dcov"]=dcov;
     ret["dcor"]=dcor;
     ret["dsY"]=dsY;
     // ret["dLL"]=dLL;
   }
   ret["df1"] = df;
   ret["df2"] = dh;
   ret["dcov.mean"] = dcovMean;
   ret["dcor.mean"] = dcorMean;
   ret["dcov.samp.mean"] = dcovSampMean;
   ret["dcor.samp.mean"] = dcorSampMean;
   ret["df1.test"] = dfp;
   ret["df2.test"] = dhp;
   ret["dnu"] = nu;
   ret["da"] = a;
   ret["dv"] = v;
   ret["time"] = time2-time1;

   //-------------------------------------------------
   // free
   if(yS) delete [] yS;
   //   if(yb) delete [] yb;
   //   if(xb) delete [] xb;
   if(z2) delete [] z2;
   if(ytempf) delete [] ytempf;
   if(svecf) delete [] svecf;
   if(ytemp) delete [] ytemp;
   if(svec) delete [] svec;
   if(YY) delete [] YY;
   if(doprdx) {if(hp) delete [] hp;}
   if(doprdz) {if(fp) delete [] fp;}
   
   return ret;
}
