#ifndef GUARD_heterbd_h
#define GUARD_heterbd_h

#include "rn.h"
#include "info.h"
#include "tree.h"
#include "treefuns.h"
#include "bartfuns.h"
#include "heterbartfuns.h"

bool heterbd(tree<double>& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, rn& gen);

#include <iostream>
//#include "heterbd.h"
#include "heterbartfuns.h"

//using std::cout;
using std::endl;

bool heterbd(tree<double>& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, rn& gen)
{
   tree<double>::npv goodbots;  //nodes we could birth at (split on)
   double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x
   if(gen.uniform() < PBx) { //do birth or death

      //--------------------------------------------------
      //draw proposal
      tree<double>::tree_p nx; //bottom node
      size_t v,c; //variable and cutpoint
      double pr; //part of metropolis ratio from proposal and prior
      bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,gen);
      //--------------------------------------------------
      //compute sufficient statistics
      size_t nr,nl; //counts in proposed bots
      double bl,br; //sums of weights
      double Ml, Mr; //weighted sum of y in proposed bots
      hetergetsuff(x,nx,v,c,xi,di,nl,bl,Ml,nr,br,Mr,sigma);

      //--------------------------------------------------
      //compute alpha
      double alpha=0.0, lalpha=0.0;
      double lhl, lhr, lht;
      if((nl>=5) && (nr>=5)) { //cludge?
         lhl = heterlh(bl,Ml,pi.tau);
         lhr = heterlh(br,Mr,pi.tau);
         lht = heterlh(bl+br,Ml+Mr,pi.tau);
   
         alpha=1.0;
         lalpha = log(pr) + (lhl+lhr-lht); 
         lalpha = std::min(0.0,lalpha);
      }

      //--------------------------------------------------
      //try metrop
      double thetal,thetar; //means for new bottom nodes, left and right
      double uu = gen.uniform();
      bool dostep = (alpha > 0) && (log(uu) < lalpha);
      if(dostep) {
         thetal = heterdrawnodetheta(bl,Ml,pi.tau,gen);
         thetar = heterdrawnodetheta(br,Mr,pi.tau,gen);
         x.birthp(nx,v,c,thetal,thetar);
         return true;
      } else {
         return false;
      }
   } else {
      //--------------------------------------------------
      //draw proposal
      double pr;  //part of metropolis ratio from proposal and prior
      tree<double>::tree_p nx; //nog node to death at
      dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      double srr,srl; //sums of weights
      double sryr, sryl; //weighted sums of y
      hetergetsuff(x, nx->getl(), nx->getr(), xi, di, srl, sryl, srr, sryr, sigma);

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
      double theta;
      if(log(gen.uniform()) < lalpha) {
         theta = heterdrawnodetheta(srl+srr,sryl+sryr,pi.tau,gen);
         x.deathp(nx,theta);
         return true;
      } else {
         return false;
      }
   }
}

bool heterbd(tree<std::vector<double>>& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, rn& gen);

bool heterbd(tree<std::vector<double>>& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, rn& gen)
{
   tree<std::vector<double>>::npv goodbots;  //nodes we could birth at (split on)
   double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x

   if(gen.uniform() < PBx) { //do birth or death

      //--------------------------------------------------
      //draw proposal
      tree<std::vector<double>>::tree_p nx; //bottom node
      size_t v,c; //variable and cutpoint
      double pr; //part of metropolis ratio from proposal and prior
      bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,gen);

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
      std::vector<double>thetal,thetar; //means for new bottom nodes, left and right
      double uu = gen.uniform();
      bool dostep = (alpha > 0) && (log(uu) < lalpha);
      if(dostep) {
        thetal = heterdrawnodetheta(srl,sryl,srwl,srwyl,srwwl,pi.tau,pi.tau2,gen);
        thetar = heterdrawnodetheta(srr,sryr,srwr,srwyr,srwwr,pi.tau,pi.tau2,gen);
         x.birthp(nx,v,c,thetal,thetar);
         return true;
      } else {
         return false;
      }
   } else {
      //--------------------------------------------------
      //draw proposal
      double pr;  //part of metropolis ratio from proposal and prior
      tree<std::vector<double>>::tree_p nx; //nog node to death at
      dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      double srl,sryl,srwl,srwyl,srwwl;
      double srr,sryr,srwr,srwyr,srwwr;
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
      std::vector<double> theta;
      if(log(gen.uniform()) < lalpha) {
        theta = heterdrawnodetheta(srl+srr,sryl+sryr,srwl+srwr,srwyl+srwyr,srwwl+srwwr,pi.tau,pi.tau2,gen);
         x.deathp(nx,theta);
         return true;
      } else {
         return false;
      }
   }
}

#endif
