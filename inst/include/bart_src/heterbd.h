#ifndef GUARD_heterbd_h
#define GUARD_heterbd_h

#include "rn.h"
#include "info.h"
#include "tree.h"
#include "treefuns.h"
#include "bartfuns.h"
#include "heterbartfuns.h"

bool heterbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, std::vector<size_t>& nv, std::vector<double>& pv, rn& gen);

#include <iostream>
//#include "heterbd.h"
#include "heterbartfuns.h"

//using std::cout;
using std::endl;

bool heterbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, std::vector<size_t>& nv, std::vector<double>& pv, rn& gen)
{
   tree::npv goodbots;  //nodes we could birth at (split on)
   double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x

   if(gen.uniform() < PBx) { //do birth or death

      //--------------------------------------------------
      //draw proposal
      tree::tree_p nx; //bottom node
      size_t v,c; //variable and cutpoint
      double pr; //part of metropolis ratio from proposal and prior
      bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,pv,gen);
      //--------------------------------------------------
      //compute sufficient statistics
      size_t nr,nl; //counts in proposed bots
      double srl, sryl, srwl, srwyl, srwwl;
      double srr, sryr, srwr, srwyr, srwwr;
      hetergetsuff(x,nx,v,c,xi,di,nl,srl,sryl,srwl,srwyl,srwwl,nr,srr,sryr,srwr,srwyr,srwwr,sigma);
      //--------------------------------------------------
      //compute alpha
      double alpha=0.0, lalpha=0.0;
      double lhl, lhr, lht;
      if((nl>=5) && (nr>=5)) { //cludge?
        lhl = heterlh(srl,sryl,pi.tau);
        lhr = heterlh(srr,sryr,pi.tau);
        lht = heterlh(srl+srr,sryl+sryr,pi.tau);
        
        alpha=1.0;
        lalpha = log(pr) + (lhl+lhr-lht); 
        lalpha = std::min(0.0,lalpha);
      }
      //--------------------------------------------------
      //try metrop
      std::vector<double> betal(di.size_beta,0.0); std::vector<double> betar(di.size_beta,0.0); //means for new bottom nodes, left and right
      double uu = gen.uniform();
      bool dostep = (alpha > 0) && (log(uu) < lalpha);
        if(dostep) {
          if(di.size_beta==1){
            betal.push_back(heterdrawnodemu0(srl,sryl,pi.tau,gen));
            betar.push_back(heterdrawnodemu0(srr,sryr,pi.tau,gen));
          }
          else if(di.size_beta==2){
            std::vector<double> betalpr,betarpr;
            betal=heterdrawnodemu1(srl,sryl,srwl,srwyl,srwwl,pi.tau,gen);
            betar=heterdrawnodemu1(srr,sryr,srwr,srwyr,srwwr,pi.tau,gen);
          }
          x.birthp(nx,v,c,betal,betar);
	  nv[v]++;
          return true;
        } else {
          return false;
        }
   } else {
     //--------------------------------------------------
     //draw proposal
     double pr;  //part of metropolis ratio from proposal and prior
     tree::tree_p nx; //nog node to death at
     dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);
     
     //--------------------------------------------------
     //compute sufficient statistics
     double srl, sryl, srwl, srwyl, srwwl;
     double srr, sryr, srwr, srwyr, srwwr;
     hetergetsuff(x, nx->getl(), nx->getr(), xi, di, srl, sryl, srwl, srwyl, srwwl, srr, sryr, srwr, srwyr, srwwr, sigma);
     
     //--------------------------------------------------
     //compute alpha
     double lhl, lhr, lht;
     lhl = heterlh(srl,sryl,pi.tau);
     lhr = heterlh(srr,sryr,pi.tau);
     lht = heterlh(srl+srr,sryl+sryr,pi.tau);

     double lalpha = log(pr) + (lht - lhl - lhr);
     lalpha = std::min(0.0,lalpha);

     //--------------------------------------------------
     //try metrop
     //double a,b,s2,yb;
     std::vector<double> beta;
     if(log(gen.uniform()) < lalpha) {
       if(di.size_beta==1){
         beta.push_back(heterdrawnodemu0(srl+srr,sryl+sryr,pi.tau,gen));
       }
       else if(di.size_beta==2){
         beta.resize(2);
         beta=heterdrawnodemu1(srl+srr,sryl+sryr,srwl+srwr,srwyl+srwyr,srwwl+srwwr,pi.tau,gen);
       }
       nv[nx->getv()]--;
       x.deathp(nx,beta);
       return true;
     } else {
       return false;
     }
   }
}

#endif
