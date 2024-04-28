#ifndef GUARD_bd_h
#define GUARD_bd_h

#include "rn.h"
#include "info.h"
#include "tree.h"
#include "treefuns.h"
#include "bartfuns.h"

bool bd(tree<double>& x, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen);

#include <iostream>
//#include "bd.h"
//#include "bartfuns.h"

//using std::cout;
using std::endl;

bool bd(tree<double>& x, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen)
{
   typename tree<double>::npv goodbots;  //nodes we could birth at (split on)
   double PBx = getpb<double>(x,xi,pi,goodbots); //prob of a birth at x

   if(gen.uniform() < PBx) { //do birth or death

      //--------------------------------------------------
      //draw proposal
      tree<double>::tree_p nx; //bottom node
      size_t v,c; //variable and cutpoint
      double pr; //part of metropolis ratio from proposal and prior
      bprop<double>(x,xi,pi,goodbots,PBx,nx,v,c,pr,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      size_t nr,nl; //counts in proposed bots
      double syl, syr; //sum of y in proposed bots
      getsuff(x,nx,v,c,xi,di,nl,syl,nr,syr);

      //--------------------------------------------------
      //compute alpha
      double alpha=0.0, lalpha=0.0;
      double lhl, lhr, lht;
      if((nl>=5) && (nr>=5)) { //cludge?
         lhl = lh(nl,syl,sigma,pi.tau);
         lhr = lh(nr,syr,sigma,pi.tau);
         lht = lh(nl+nr,syl+syr,sigma,pi.tau);
   
         alpha=1.0;
         lalpha = log(pr) + (lhl+lhr-lht) + log(sigma);
         lalpha = std::min(0.0,lalpha);
      }

      //--------------------------------------------------
      //try metrop
      double thetal,thetar; //means for new bottom nodes, left and right
      double uu = gen.uniform();
      bool dostep = (alpha > 0) && (log(uu) < lalpha);
      if(dostep) {
        thetal = drawnodetheta(nl,syl,pi.tau,sigma,gen);
        thetar = drawnodetheta(nr,syr,pi.tau,sigma,gen);
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
      size_t nr,nl; //counts at bots of nx
      double syl, syr; //sum at bots of nx
      getsuff(x, nx->getl(), nx->getr(), xi, di, nl, syl, nr, syr);

      //--------------------------------------------------
      //compute alpha
      double lhl, lhr, lht;
      lhl = lh(nl,syl,sigma,pi.tau);
      lhr = lh(nr,syr,sigma,pi.tau);
      lht = lh(nl+nr,syl+syr,sigma,pi.tau);

      double lalpha = log(pr) + (lht - lhl - lhr) - log(sigma);
      lalpha = std::min(0.0,lalpha);

      //--------------------------------------------------
      //try metrop
      double mu;
      if(log(gen.uniform()) < lalpha) {
        mu = drawnodetheta(nl+nr,syl+syr,pi.tau,sigma,gen);
         x.deathp(nx,mu);
         return true;
      } else {
         return false;
      }
   }
}

bool bd(tree<std::vector<double>>& x, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen);

bool bd(tree<std::vector<double>>& x, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen)
{
  tree<std::vector<double>>::npv goodbots;  //nodes we could birth at (split on)
  double PBx = getpb<std::vector<double>>(x,xi,pi,goodbots); //prob of a birth at x

  if(gen.uniform() < PBx) { //do birth or death

    //--------------------------------------------------
    //draw proposal
    tree<std::vector<double>>::tree_p nx; //bottom node
    size_t v,c; //variable and cutpoint
    double pr; //part of metropolis ratio from proposal and prior
    bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,gen);

    //--------------------------------------------------
    //compute sufficient statistics (Constant)
    size_t nr,nl; //counts in proposed bots
    double syl, swl, swyl, swwl;
    double syr, swr, swyr, swwr; //sum of y in proposed bots
    getsuff(x,nx,v,c,xi,di,nl,syl,swl,swyl,swwl,nr,syr,swr,swyr,swwr);
      
    //--------------------------------------------------
    //compute alpha
    double alpha=0.0, lalpha=0.0;
    double lhl, lhr, lht;
    if((nl>=5) && (nr>=5)) { //cludge?
      lhl = lh(nl,syl,sigma,pi.tau);
      lhr = lh(nr,syr,sigma,pi.tau);
      lht = lh(nl+nr,syl+syr,sigma,pi.tau);
      alpha=1.0;
      lalpha = log(pr) + (lhl+lhr-lht) + log(sigma);
      lalpha = std::min(0.0,lalpha);
    }       
    //--------------------------------------------------
    //try metrop
    std::vector<double> thetal, thetar; //means for new bottom nodes, left and right
    double uu = gen.uniform();
    bool dostep = (alpha > 0) && (log(uu) < lalpha);
    if(dostep) {
      thetal = drawnodetheta(nl,syl,swl,swyl,swwl,pi.tau,pi.tau2,sigma,gen);
      thetar = drawnodetheta(nr,syr,swr,swyr,swwr,pi.tau,pi.tau2,sigma,gen);
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
    size_t nr,nl;
    double syl, swl, swyl, swwl;
    double syr, swr, swyr, swwr;
    getsuff(x, nx->getl(), nx->getr(), xi, di, nl, syl, swl, swyl, swwl, nr, syr, swr, swyr, swwr);
     
    //--------------------------------------------------
    //compute alpha
    double lhl, lhr, lht;
    lhl = lh(nl,syl,sigma,pi.tau);
    lhr = lh(nr,syr,sigma,pi.tau);
    lht = lh(nl+nr,syl+syr,sigma,pi.tau);
     
    double lalpha = log(pr) + (lht - lhl - lhr) - log(sigma);
    lalpha = std::min(0.0,lalpha);
     
    //--------------------------------------------------
    //try metrop
    std::vector<double> theta;
    if(log(gen.uniform()) < lalpha) {
      theta = drawnodetheta(nl+nr,syl+syr,swl+swr,swyl+swyr,swwl+swwr,pi.tau,pi.tau2,sigma,gen);
      x.deathp(nx,theta);
      return true;
    } else {
      return false;
    }
  }
}

#endif
