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

RcppExport SEXP cbiv(
   SEXP _z,
   SEXP _zp,
   SEXP _x,
   SEXP _xp,
   SEXP _T,
   SEXP _Y,
   SEXP _type1,
   SEXP _type2,
   SEXP _burn,
   SEXP _nd,
   SEXP _keepevery,
   SEXP _burnf1,
   SEXP _burnf22,
   SEXP _m1,
   SEXP _m2,
   SEXP _nc,
   SEXP _power,
   SEXP _base,
   SEXP _tauf1,
   SEXP _tauf22,
   SEXP _betabar,
   SEXP _Abeta,
   SEXP _doDP,
   SEXP _v,
   SEXP _nu,
   SEXP _a,
   SEXP _sparse1,
   SEXP _sparse2,
   SEXP _ag,
   SEXP _priag,
   SEXP _centermeans,
   SEXP _offset,
   SEXP _in_sample,
   SEXP _aggregate_sigma,
   SEXP _f1s,
   SEXP _f22s,
   SEXP _betas,
   SEXP _mTs,
   SEXP _mYs,
   SEXP _sTs,
   SEXP _gammas,
   SEXP _sYs,
   SEXP _nullf1,
   SEXP _nullf22,
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

   //zx
   Rcpp::NumericMatrix zm(_z);
   size_t nz = zm.ncol();
   size_t pz = zm.nrow();
   double *z = &zm[0];

   //zxp
   Rcpp::NumericMatrix zpm(_zp);
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

   // type2 ----------
   size_t type1 = Rcpp::as<size_t>(_type1);
   size_t type2 = Rcpp::as<size_t>(_type2);
   
   // burn, nd
   size_t burn = Rcpp::as<size_t>(_burn);
   size_t nd = Rcpp::as<size_t>(_nd);
   size_t ke = Rcpp::as<size_t>(_keepevery);
   size_t burnf1 = Rcpp::as<size_t>(_burnf1);
   size_t burnf22 = Rcpp::as<size_t>(_burnf22);

   //bart prior
   size_t m1 = Rcpp::as<size_t>(_m1);
   size_t m2 = Rcpp::as<size_t>(_m2);
   size_t nc = Rcpp::as<size_t>(_nc);
   double mybeta = Rcpp::as<double>(_power);
   double alpha = Rcpp::as<double>(_base);
   double tauf1 = Rcpp::as<double>(_tauf1);
   double tauf22 = Rcpp::as<double>(_tauf22);
   size_t sparse1 = Rcpp::as<size_t>(_sparse1);
   size_t sparse2 = Rcpp::as<size_t>(_sparse2);
     
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

   //should the means be centered after each DpMuSigma draw
   bool centermeans = Rcpp::as<bool>(_centermeans);
   bool in_sample = Rcpp::as<bool>(_in_sample);
   bool aggregate_sigma = Rcpp::as<bool>(_aggregate_sigma);
   double offset = Rcpp::as<double>(_offset);
   bool doDP = Rcpp::as<bool>(_doDP);
   if(type1==2||type2==2) doDP=false;

   //starting values
   //fs,hs
   Rcpp::NumericVector f1sv(_f1s);
   double *f1s = &f1sv[0];
   Rcpp::NumericVector f22sv(_f22s);
   double *f22s = &f22sv[0];

   double betas = Rcpp::as<double>(_betas);
   double mTs = Rcpp::as<double>(_mTs);
   double mYs = Rcpp::as<double>(_mYs);
   double sTs = Rcpp::as<double>(_sTs);
   double gammas = Rcpp::as<double>(_gammas);
   double sYs = Rcpp::as<double>(_sYs);

   //null model bools
   //nullf, nullf22
   bool nullf1 = Rcpp::as<bool>(_nullf1);
   bool nullf22 = Rcpp::as<bool>(_nullf22);
   
   //other
   size_t printevery = Rcpp::as<size_t>(_printevery);
   bool quiet = Rcpp::as<bool>(_quiet);

   size_t n = nT;

   // doh1
   bool dof22 = (nx!=0);
   if(dof22) burnf22=1;

   if(!quiet){
     Rprintf("*****************************************************************\n");
     Rprintf("*****Into main of cbiv\n");
   }
   
   //--------------------------------------------------
   // print args
   if(!quiet){
     Rprintf("***burn, nd, burnf1, burnf22: %ld, %ld, %ld\n",burn,nd,burnf1,burnf22);
     Rprintf("*** Prior:\n");
     Rprintf("m1 (num trees stage 1), m2 (num trees stage 2), nc (num cut points), %ld, %ld, %ld\n",m1,m2,nc);
     Rprintf("****power, base, tauf1, tauf22: %lf, %lf, %lf, %lf\n",mybeta,alpha,tauf1,tauf22);
     Rprintf("betabar: %lf\n",betabar);
     Rprintf("Abeta: %lf\n",Abeta);
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
     Rprintf("\t first and last f1s: %lf, %lf\n",f1s[0],f1s[n-1]);
     Rprintf("\t first and last f22s: %lf, %lf\n",f22s[0],f22s[n-1]);
     Rprintf(" starting for beta, mT, mY, sT, gamma, sY: %lf, %lf, %lf, %lf, %lf, %lf\n",betas,mTs,mYs,sTs,gammas,sYs);

     Rprintf("***other\n");
     Rprintf("printevery: %ld\n",printevery);
   }

   //-------------------------------------------------
   //-------------------------------------------------
   // bart probit auxiliary setup (if stage2==2)
   //--------------------------------------------------
   double *YY = new double[n];
   if(type2==1)
     for(size_t j=0;j<n;j++) YY[j]=Y[j]-offset;
   else if(type2==2){
     double sign;
     for(size_t j=0;j<n;j++){
       if(Y[j]==1) sign=1.;
       else sign=-1.;
       YY[j] = sign*rtnorm(sign*(f22s[j]+(T[j]-f1s[j])*gammas/sTs), -sign*offset, sqrt(1-gammas*gammas), gen);
     }
   }
   
   //-------------------------------------------------
   //-------------------------------------------------
   // bart f1 setup
   //--------------------------------------------------   
   // first make z2 which is z with the columns repeated
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
   double *ytempf1 = new double[n*2];  //y for h bart
   double *svecf1 = new double[n*2];   // sigma_i for h bart
   size_t *vcf1 = new size_t[pz];     // variable count for f1 bart
   double *lpvf1 = new double[pz];    // log probability variable split
   double *pvf1 = new double[pz];
   for(size_t i=0; i<n; i++) {
     ytempf1[2*i] = T[i] - mTs;
     svecf1[2*i] = sTs;
     ytempf1[2*i+1] = T[i] - mTs;
     svecf1[2*i+1] = sTs;
   }
   for(size_t j=0; j<pz; j++) {
     vcf1[j]=0;
     lpvf1[j]=-1*::log(pz);
     pvf1[j]=::exp(lpvf1[j]);
   }
   bmf1.setdata(pz,n*2,z2,ytempf1,nc);
   bmf1.setprior(alpha,mybeta,tauf1,pvf1);
   bmf1.draw(svecf1,gen);
   // bmf1 ouput storage
   Rcpp::NumericMatrix df1(1,1);
   if(in_sample) df1=Rcpp::NumericMatrix(nd,n); //f draws on train
   Rcpp::NumericMatrix df1p(nd,nzp);
   Rcpp::IntegerMatrix dvcf1(nd,pz);
     
   //-------------------------------------------------
   //-------------------------------------------------
   // bart f22 setup
   //--------------------------------------------------
   heterbart<double> bmf22(m2);
  
   double *ytempf22 = new double[n];  //y for h bart
   double *svecf22 = new double[n];   // sigma_i for h bart
   size_t *vcf22 = new size_t[px]; // variable count for f2 bart
   double *lpvf22 = new double[px];
   double *pvf22 = new double[px];
      
   double ZZ1=0.0;
   for(size_t j=0;j<n;j++) {
     ZZ1 = (T[j] - mTs - bmf1.f(2*j))/sTs;
     ytempf22[j] = YY[j] - mYs - betas*T[j] - gammas * ZZ1;
     svecf22[j] = sYs; 
   }
   for(size_t j=0; j<px; j++) {
     vcf22[j]=0;
     lpvf22[j]=-1.*::log(px);
     pvf22[j]=::exp(lpvf22[j]);
   }
   bmf22.setdata(px,n,x,ytempf22,nc);
   bmf22.setprior(alpha,mybeta,tauf22,pvf22);
   bmf22.draw(svecf22,gen);
   Rcpp::NumericMatrix df22(1,1);
   if(in_sample) df22=Rcpp::NumericMatrix(nd,n); //h draws on train
   Rcpp::NumericMatrix df22p(nd,nxp);
   Rcpp::IntegerMatrix dvcf22(nd,px);

   //--------------------------------------------------
   //--------------------------------------------------
   // beta setup
   double *yb = new double[n];
   double *xb = new double[n];
   double Z1=0.0;
   double betad = betas;
   if(!nullf22){
     for(size_t i=0;i<n;i++) {
       Z1 = (T[i]-mTs-f1s[i])/sTs;
       yb[i] = (YY[i] - mYs - f22s[i] - gammas * Z1)/sYs;       
       xb[i] = T[i]/sYs;
     }
     betad = lreg(n,xb,yb,1.,betabar,Abeta,gen);
   }
   Rcpp::NumericVector dbeta(nd); //storage for output beta draws

   //--------------------------------------------------
   //--------------------------------------------------
   // DpMuSigma setup 
   if(!quiet) std::cout << "\n*** making a DpMuSigma object:\n";
   Dp::dv itheta(5,0); 
   itheta[0] = mTs; itheta[1] = mYs; itheta[2] = sTs; itheta[3] = gammas; itheta[4] = sYs; 
   if(!quiet){
     std::cout << "itheta:\n"; 
     printdv(itheta);
   }

   //data for errors
   double *yS = new double[2*n];
   //intialize using starting values
   for(size_t i=0;i<n;i++) {
      yS[2*i] = T[i] - f1s[i];
      yS[2*i+1] = YY[i] - betas*T[i] - f22s[i];
   }
   if(!quiet)
     std::cout << "check yS: " << yS[0] << ", " << yS[2*n-1] << std::endl;
   
   DpMuSigma dpmS(n,itheta,doDP); 
   dpmS.setData(yS); 
   dpmS.setPrior(v,nu,a); 
   dpmS.setAlphaPrior(ag,priag); 
   dpmS.setalpha(1.0); 

   /* // to screen */
   if(!quiet) {
     dpmS.toscreen();
   }
   
   //output storage
   Rcpp::IntegerVector dnpart(nd);
   Rcpp::NumericVector dalpha(nd);
   Rcpp::NumericMatrix dsigma1(nd,n);
   Rcpp::NumericMatrix dsigma2(nd,n);
   Rcpp::NumericMatrix dcov(nd,n);
   
   // test prediction setup
   bool doprdx = nxp>0;
   bool doprdz = nzp>0;

   double *f22p = 0;
   double *fp = 0;
   if(doprdx)
     {
       f22p = new double[nxp];
     }
   if(doprdz)
     {
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
      if(i%printevery==0&!quiet) Rprintf("done %d (out of %d)\n",i,nd+burn);
      // auxiliary y_tilde (if probit)
         if(type2==2){
           double sign;
           for(size_t j=0;j<n;j++){
             if(i>0) {mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];}
             if(Y[j]==1) sign=1.;
             else sign=-1.;
             YY[j] = sign*rtnorm(sign*(bmf22.f(j)+betad*T[j]+(T[j]-bmf1.f(2*j))*gamma/sT), -sign*offset, sqrt(1-gamma*gamma), gen);
           }
         }
      // beta conditional -----------------------------------
      // update -----
      if(!nullf22){
        for(size_t j=0;j<n;j++) {
          if(i>0) {mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];}
          Z1 = (T[j] - mT - bmf1.f(2*j))/sT;
          yb[j] = (YY[j] - mY - bmf22.f(j) - gamma * Z1)/sY;	  
          xb[j] = T[j]/sY;
        }
        //draw -----
        betad = lreg(n,xb,yb,1.0,betabar,Abeta,gen);
      }

      // f22 conditional -----------------------------------
      // update ----- 
          for(size_t j=0;j<n;j++) {
            if(i>0) {mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];}
            ZZ1 = (T[j] - mT - bmf1.f(2*j))/sT;
            ytempf22[j] = YY[j] - mY - betad*T[j] - gamma * ZZ1;
            svecf22[j] = sY;
          }
        //draw -----
        bmf22.draw(svecf22,gen);      

      // f conditional --------------------------------
      if(!nullf1) {
        // update -----
        double R;
        for(size_t j=0; j<n; j++) {
          if(i>0) {mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];}
          ytempf1[2*j] = T[j] - mT;
          svecf1[2*j] = sT;
          R = (betad*sT+gamma)*(T[j]-mT) - sT*(YY[j] - mY - bmf22.f(j) - betad*mT);	  
          ytempf1[2*j+1] = R/gamma;
          svecf1[2*j+1] = (sT*sY)/gamma;      
        }
        // draws -----
        bmf1.draw(svecf1,gen);      
      }

      if(type2==1){
      // Sigma conditional -----------------------------
        for(size_t j=0;j<n;j++) {
          yS[2*j] = T[j] - bmf1.f(2*j);
          yS[2*j+1] = YY[j] - betad*T[j] - bmf22.f(j);    
        }
        //draw -----
        dpmS.draw(gen); 
        if(centermeans) dpmS.center(); 
        tmat = dpmS.thetaMatrix();
      }
      else{
        //update -----
        for(size_t j=0;j<n;j++) {
          yS[2*j] = T[j] - bmf1.f(2*j);
          yS[2*j+1] = YY[j] - betad*T[j] - bmf22.f(j);
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
        size_t rho_index = 0;
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
        }
        double max_log_post_sig=log_sig_posterior(sig_grid[0],sumY1sq,sumY2sq,sumY1Y2,rho_grid[rho_index],0.,sqrt(10.),n);
        for(size_t g=0;g<sig_grid_size-1;g++) {
          log_post_sig[g]=log_sig_posterior(sig_grid[g],sumY1sq,sumY2sq,sumY1Y2,rho_grid[rho_index],0.,sqrt(10.),n);
          if(log_post_sig[g]>max_log_post_sig) max_log_post_sig=log_post_sig[g];
        }
        sum_exp_max=0.;
        for(size_t g=0;g<sig_grid_size-1;g++){
          sum_exp_max += exp(log_post_sig[g]-max_log_post_sig);
        }
        lse=max_log_post_sig+log(sum_exp_max);
        for(size_t g=0;g<sig_grid_size-1;g++){
          prob_sig[g] = exp(log_post_sig[g]-lse);
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
          tmat[j][0]=0.;
          tmat[j][1]=0.;
          tmat[j][2]=sig_grid[sig_index];
          tmat[j][3]=rho_grid[rho_index];
          if(tmat[j][3]==0.) tmat[j][3]=0.001;
          tmat[j][4]=sqrt(1.-tmat[j][3]*tmat[j][3]);
        }      
      }
      
      // Predictions     
      if(doprdx)
        bmf22.predict(pxp,nxp,xp,f22p);      
      if(doprdz)
        bmf1.predict(pzp,nzp,zp,fp);

      // Get variable counts within trees
      bmf1.count_v_rules(vcf1);
      bmf22.count_v_rules(vcf22);
      double conc1, conc2;
      if(sparse1==1){
        // Draw concentration parameter and then log Dirichlet
        conc1=draw_concentration(lpvf1,pz);
        draw_log_dirichlet(vcf1,lpvf1,conc1,pz,gen);
        // Convert to probability
        for(size_t j=0;j<pz;j++) pvf1[j]=exp(lpvf1[j]);
        // Reset splitting rules
        bmf1.setprior(alpha,mybeta,tauf1,pvf1);
      }
      if(sparse2==1){
        // Draw concentration parameter and then log Dirichlet
        conc2=draw_concentration(lpvf22,px);
        draw_log_dirichlet(vcf22,lpvf22,conc2,px,gen);
        // Convert to probability
        for(size_t j=0;j<px;j++) pvf22[j]=exp(lpvf22[j]);
        //Reset splitting rules
        bmf22.setprior(alpha,mybeta,tauf22,pvf22);
      }

      
      // Record posterior samples
      if(i >= burn) {
        size_t j=i-burn;
        if(j % ke == 0){
          size_t draw=j/ke;
          
          dnpart[draw] = dpmS.npart();
          dalpha[draw] = dpmS.getalpha();
	  for(size_t k=0;k<n;k++) {
	    dsigma1(draw,k) = tmat[k][2];
	    dsigma2(draw,k) = sqrt(tmat[k][3]*tmat[k][3] + tmat[k][4]*tmat[k][4]);
	    dcov(draw,k) = tmat[k][3]*tmat[k][2];
	  }  

	  if(in_sample){
	    for(size_t k=0;k<n;k++) {
	      df22(draw,k) = bmf22.f(k)+offset;
	    }	    
	    for(size_t k=0;k<n;k++) {
	      df1(draw,k) = bmf1.f(2*k);
	    }
	  }	           
        	 
          dbeta[draw] = betad;
        
          if(doprdx) {
            for(size_t k=0;k<nxp;k++) {
              df22p(draw,k) = f22p[k]+offset;
            }
          }
        
          if(doprdz) {
            for(size_t k=0;k<nzp;k++) {
              df1p(draw,k) = fp[k];
            }
          }
          for(size_t k=0;k<pz;k++){
            dvcf1(draw,k)=vcf1[k];
          }            
          for(size_t k=0;k<px;k++){
            dvcf22(draw,k)=vcf22[k];
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
   if(in_sample){
     ret["df1"]=df1;
     ret["df22"]=df22;
   }
   ret["dbeta"] = dbeta;
   ret["dsigma1"]=dsigma1;
   ret["dsigma2"]=dsigma2;
   ret["dcov"]=dcov;
   ret["df1.test"] = df1p;
   ret["df22.test"] = df22p;
   ret["df1.varcount"] = dvcf1;
   ret["df22.varcount"] = dvcf22;

   
   //-------------------------------------------------
   // free
   if(yS) delete [] yS;
   if(yb) delete [] yb;
   if(xb) delete [] xb;
   if(ytempf22) delete [] ytempf22;
   if(svecf22) delete [] svecf22;
   if(vcf1) delete [] vcf1;
   if(lpvf1) delete [] lpvf1;
   if(pvf1) delete [] pvf1;
   if(z2) delete [] z2;
   if(ytempf1) delete [] ytempf1;
   if(svecf1) delete [] svecf1;
   if(vcf22) delete [] vcf22;
   if(lpvf22) delete [] lpvf22;
   if(pvf22) delete [] pvf22;
   if(YY) delete [] YY;
   if(doprdx) {if(f22p) delete [] f22p;}
   if(doprdz) {if(fp) delete [] fp;}

   return ret;
   }
