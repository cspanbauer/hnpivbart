#ifndef GUARD_tree_h
#define GUARD_tree_h

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstddef>
#include <vector>
#include "rn.h"

//--------------------------------------------------
//xinfo xi, then xi[v][c] is the c^{th} cutpoint for variable v.
//left if x[v] < xi[v][c]
typedef std::vector<double> vec_d; //double vector
typedef std::vector<vec_d> xinfo; //vector of vectors, will be split rules

//--------------------------------------------------
template <typename T>
class tree {
public:
   //typedefs--------------------
   typedef tree<T>* tree_p;
   typedef const tree<T>* tree_cp;
  typedef std::vector<tree_p> npv; 
  typedef std::vector<tree_cp> cnpv;
   //contructors,destructors--------------------
   tree(): v(0),c(0),p(0),l(0),r(0) {}
   tree(const tree& n): v(0),c(0),p(0),l(0),r(0) {cp(this,&n);}
   tree(T itheta): theta(itheta),v(0),c(0),p(0),l(0),r(0) {}
   void tonull(); //like a "clear", null tree has just one node
   ~tree() {tonull();}
   //operators----------
   tree& operator=(const tree&);
   //interface--------------------
   //set
   void settheta(T _theta) {this->theta=_theta;}
  //  void settheta(std::vector<double> _theta) {this->theta=_theta;}//for(size_t k=0;k<_theta.size();k++) this->theta[k]=_theta[k];}
   void setv(size_t v) {this->v = v;}
   void setc(size_t c) {this->c = c;}
   //get
   T gettheta() const {return theta;}
   size_t getv() const {return v;}
   size_t getc() const {return c;}
   tree_p getp() {return p;}  
   tree_p getl() {return l;}
   tree_p getr() {return r;}
   //tree functions--------------------
   tree_p getptr(size_t nid); //get node pointer from node id, 0 if not there
   void pr(bool pc=true); //to screen, pc is "print children"
   size_t treesize(); //number of nodes in tree
   size_t nnogs();    //number of nog nodes (no grandchildren nodes)
   size_t nbots();    //number of bottom nodes
   bool birth(size_t nid, size_t v, size_t c, T thetal, T thetar);
   bool death(size_t nid, T theta);
   void birthp(tree_p np,size_t v, size_t c, T thetal, T thetar);
   void deathp(tree_p nb, T theta);
  void getbots(npv& bv);         //get bottom nodes
   void getnogs(npv& nv);         //get nog nodes (no granchildren)
   void getnodes(npv& v);         //get vector of all nodes
   void getnodes(cnpv& v) const;  //get vector of all nodes (const)
   tree_p bn(double *x,xinfo& xi); //find Bottom Node
   void rg(size_t v, int* L, int* U); //recursively find region [L,U] for var v
   //node functions--------------------
   size_t nid() const; //nid of a node
   size_t depth();  //depth of a node
   char ntype(); //node type t:top, b:bot, n:no grandchildren i:interior (t can be b)
   bool isnog();
private:
   T theta; //template to allow multiple types
   //rule: left if x[v] < xinfo[v][c]
   size_t v;
   size_t c;
   //tree structure
   tree_p p; //parent
   tree_p l; //left child
   tree_p r; //right child
   //utiity functions
   void cp(tree_p n,  tree_cp o); //copy tree
};

#include <string>
#include <map>
//#include "tree.h"
//--------------------
// node id
template <typename T>
size_t tree<T>::nid() const 
{
   if(!p) return 1; //if you don't have a parent, you are the top
   if(this==p->l) return 2*(p->nid()); //if you are a left child
   else return 2*(p->nid())+1; //else you are a right child
}
//--------------------
template <typename T>
typename tree<T>::tree_p tree<T>::getptr(size_t nid)
{
   if(this->nid() == nid) return this; //found it
   if(l==0) return 0; //no children, did not find it
   tree_p lp = l->getptr(nid);
   if(lp) return lp; //found on left
   tree_p rp = r->getptr(nid);
   if(rp) return rp; //found on right
   return 0; //never found it
}
//--------------------
//add children to  bot node nid
template <typename T>
bool tree<T>::birth(size_t nid,size_t v, size_t c, T thetal, T thetar)
{
   tree_p np = getptr(nid);
   if(np==0) {
      cout << "error in birth: bottom node not found\n";
      return false; //did not find note with that nid
   }
   if(np->l!=0) {
      cout << "error in birth: found node has children\n";
      return false; //node is not a bottom node
   }

   //add children to bottom node np
   tree_p l = new tree<T>;
   l->theta=thetal;
   tree_p r = new tree<T>;
   r->theta=thetar;
   np->l=l;
   np->r=r;
   np->v = v; np->c=c;
   l->p = np;
   r->p = np;

   return true;
}
//--------------------
//depth of node
template <typename T>
size_t tree<T>::depth()
{
   if(!p) return 0; //no parents
   else return (1+p->depth());
}
//--------------------
//tree size
template <typename T>
size_t tree<T>::treesize()
{
   if(l==0) return 1;  //if bottom node, tree size is 1
   else return (1+l->treesize()+r->treesize());
}
//--------------------
//node type
template <typename T>
char tree<T>::ntype()
{
   //t:top, b:bottom, n:no grandchildren, i:internal
   if(!p) return 't';
   if(!l) return 'b';
   if(!(l->l) && !(r->l)) return 'n';
   return 'i';
}
//--------------------
//print out tree(pc=true) or node(pc=false) information
template <typename T>
void tree<T>::pr(bool pc) 
{
   size_t d = depth();
   size_t id = nid();

   size_t pid;
   if(!p) pid=0; //parent of top node
   else pid = p->nid();

   std::string pad(2*d,' ');
   std::string sp(", ");
   if(pc && (ntype()=='t'))
      cout << "tree size: " << treesize() << std::endl;
   cout << pad << "(id,parent): " << id << sp << pid;
   cout << sp << "(v,c): " << v << sp << c;
   // cout << sp << "theta: " << theta;
   cout << sp << "type: " << ntype();
   cout << sp << "depth: " << depth();
   cout << sp << "pointer: " << this << std::endl;

   if(pc) {
      if(l) {
         l->pr(pc);
         r->pr(pc);
      }
   }
}
//--------------------
//kill children of  nog node nid
template <typename T>
bool tree<T>::death(size_t nid, T theta)
{
   tree_p nb = getptr(nid);
   if(nb==0) {
      cout << "error in death, nid invalid\n";
      return false;
   }
   if(nb->isnog()) {
      delete nb->l;
      delete nb->r;
      nb->l=0;
      nb->r=0;
      nb->v=0;
      nb->c=0;
      nb->theta=theta;
      return true;
   } else {
      cout << "error in death, node is not a nog node\n";
      return false;
   }
}
//--------------------
//is the node a nog node
template <typename T>
bool tree<T>::isnog() 
{
   bool isnog=true;
   if(l) {
      if(l->l || r->l) isnog=false; //one of the children has children.
   } else {
      isnog=false; //no children
   }
   return isnog;
}
//--------------------
template <typename T>
size_t tree<T>::nnogs() 
{
   if(!l) return 0; //bottom node
   if(l->l || r->l) { //not a nog
      return (l->nnogs() + r->nnogs());
   } else { //is a nog
      return 1;
   }
}
//--------------------
template <typename T>
size_t tree<T>::nbots() 
{
   if(l==0) { //if a bottom node
      return 1;
   } else {
      return l->nbots() + r->nbots();
   }
}
//--------------------
//get bottom nodes
template <typename T>
void tree<T>::getbots(npv& bv)
{
   if(l) { //have children
      l->getbots(bv);
      r->getbots(bv);
   } else {
      bv.push_back(this);
   }
}
//--------------------
//get nog nodes
template <typename T>
void tree<T>::getnogs(npv& nv)
{
   if(l) { //have children
      if((l->l) || (r->l)) {  //have grandchildren
         if(l->l) l->getnogs(nv);
         if(r->l) r->getnogs(nv);
      } else {
         nv.push_back(this);
      }
   }
}
//--------------------
//get all nodes
template <typename T>
void tree<T>::getnodes(npv& v)
{
   v.push_back(this);
   if(l) {
      l->getnodes(v);
      r->getnodes(v);
   }
}
//--------------------
template <typename T>
void tree<T>::getnodes(cnpv& v)  const
{
   v.push_back(this);
   if(l) {
      l->getnodes(v);
      r->getnodes(v);
   }
}
//--------------------
template <typename T>
typename tree<T>::tree_p tree<T>::bn(double *x, xinfo& xi)
{
   if(l==0) return this; //no children
   if(x[v] < xi[v][c]) {
      return l->bn(x,xi);
   } else {
      return r->bn(x,xi);
   }
}
//--------------------
//find region for a given variable
template <typename T>
void tree<T>::rg(size_t v, int* L, int* U)
{
   if(this->p==0)  {
      return;
   }
   if((this->p)->v == v) { //does my parent use v?
      if(this == p->l) { //am I left or right child
         if((int)(p->c) <= (*U)) *U = (p->c)-1;
         p->rg(v,L,U);
      } else {
         if((int)(p->c) >= *L) *L = (p->c)+1;
         p->rg(v,L,U);
      }
   } else {
      p->rg(v,L,U);
   }
}
//--------------------
//cut back to one node
template <>
void tree<double>::tonull()
{
   size_t ts = treesize();
   //loop invariant: ts>=1
   while(ts>1) { //if false ts=1
      npv nv;
      getnogs(nv);
      for(size_t i=0;i<nv.size();i++) {
         delete nv[i]->l;
         delete nv[i]->r;
         nv[i]->l=0;
         nv[i]->r=0;
      }
      ts = treesize(); //make invariant true
   }
   theta=0.0;
   v=0;c=0;
   p=0;l=0;r=0;
}
//--------------------
//cut back to one node
template <>
void tree<std::vector<double>>::tonull()
{
   size_t ts = treesize();
   //loop invariant: ts>=1
   while(ts>1) { //if false ts=1
      npv nv;
      getnogs(nv);
      for(size_t i=0;i<nv.size();i++) {
         delete nv[i]->l;
         delete nv[i]->r;
         nv[i]->l=0;
         nv[i]->r=0;
      }
      ts = treesize(); //make invariant true
   }
   std::fill(theta.begin(),theta.end(),0.0);
   v=0;c=0;
   p=0;l=0;r=0;
}

//--------------------
//copy tree tree o to tree n
template <typename T>
void tree<T>::cp(tree_p n, tree_cp o)
//assume n has no children (so we don't have to kill them)
//recursion down
{
   if(n->l) {
      cout << "cp:error node has children\n";
      return;
   }

   n->theta = o->theta;
   n->v = o->v;
   n->c = o->c;

   if(o->l) { //if o has children
      n->l = new tree<T>;
      (n->l)->p = n;
      cp(n->l,o->l);
      n->r = new tree<T>;
      (n->r)->p = n;
      cp(n->r,o->r);
   }
}
//--------------------------------------------------
//operators
template <typename T>
tree<T>& tree<T>::operator=(const tree<T>& rhs)
{
   if(&rhs != this) {
      tonull(); //kill left hand side (this)
      cp(this,&rhs); //copy right hand side to left hand side
   }
   return *this;
}
//--------------------
//add children to bot node *np
template <typename T>
void tree<T>::birthp(tree_p np,size_t v, size_t c, T thetal, T thetar)
{
   tree_p l = new tree<T>;
   l->theta=thetal;
   tree_p r = new tree<T>;
   r->theta=thetar;
   np->l=l;
   np->r=r;
   np->v = v; np->c=c;
   l->p = np;
   r->p = np;
}
//--------------------
//kill children of  nog node *nb
template <typename T>
void tree<T>::deathp(tree_p nb, T theta)
{
   delete nb->l;
   delete nb->r;
   nb->l=0;
   nb->r=0;
   nb->v=0;
   nb->c=0;
   nb->theta=theta;
}

#endif
