#ifndef TREE
#define TREE

#include <vector>

using namespace std;

template <class V>
class BSTree {
 public:
  
  class Node {
  public:
    Node (V v, Node *p) : v(v), p(p), l(0), r(0) {}

    void list (vector<V> &res) const {
      if (l) l->list(res);
      res.push_back(v);
      if (r) r->list(res);
    }

    int depth () const {
      int dl = l ? l->depth() : 0, dr = r ? r->depth() : 0, d = dl < dr ? dr : dl;
      return 1 + d;
    }
  
    int size () const {
      int sl = l ? l->size() : 0, sr = r ? r->size() : 0;
      return 1 + sl + sr;
    }
  
    V v;
    Node *p, *l, *r;
  };

  BSTree () : r(0) {}
  
  ~BSTree () { erase(r); }
  
  void erase (Node *x) {
    if (x) {
      erase(x->l);
      erase(x->r);
      delete x;
    }
  }	

  Node * predecessor (Node *x) {
    if (x->l) {
      Node *y = x->l;
      while (y->r)
	y = y->r;
      return y;
    }
    Node *y = x->p;
    while (y && x == y->l) {
      x = y;
      y = y->p;
    }
    return y;
  }
  
  Node * successor (Node *x) {
    if (x->r) {
      Node *y = x->r;
      while (y->l)
	y = y->l;
      return y;
    }
    Node *y = x->p;
    while (y && x == y->r) {
      x = y;
      y = y->p;
    }
    return y;
  }

  Node * insert (V v) {
    Node *x = r, *y = 0;
    int s;
    while (x) {
      s = v->order(x->v);
      if (s == 0) {
	//splay(x);
	return x;
      }
      y = x;
      x = s == 1 ? x->l : x->r;
    }
    Node *z = new Node(v, y);
    if (!y)
      r = z;
    else if (s == 1)
      y->l = z;
    else
      y->r = z;
    //splay(z);
    return z;
  }

  void remove (Node *z) {
    if (!z->l)
      transplant(z, z->r);
    else if (!z->r)
      transplant(z, z->l);
    else {
      Node *y = successor(z);
      if (y->p != z) {
	transplant(y, y->r);
	y->r = z->r;
	y->r->p = y;
      }
      transplant(z, y);
      y->l = z->l;
      y->l->p = y;
    }
    delete z;
  }

  void transplant (Node *u, Node *v) {
    if (!u->p)
      r = v;
    else if (u == u->p->l)
      u->p->l = v;
    else
      u->p->r = v;
    if (v)
      v->p = u->p;
  }
    
  void splay (Node *x) {
    while (x->p)
      splayAux(x);
  }

  void splayAux (Node *x) {
    Node *y = x->p, *z = y->p, *xl = x->l, *xr = x->r, *yl = y->l, *yr = y->r;
    if (!z && x == y->r) {
      transplant(y, x);
      x->l = y;
      y->p = x;
      y->r = xl;
      if (xl) xl->p = y;
    }
    else if (!z) {
      transplant(y, x);
      x->r = y;
      y->p = x;
      y->l = xr;
      if (xr) xr->p = y;
    }
    else if (x == y->r && y == z->r) {
      transplant(z, x);
      x->l = y;
      y->p = x;
      y->l = z;
      y->r = xl;
      if (xl) xl->p = y;
      z->p = y;
      z->r = yl;
      if (yl) yl->p = z;
    }
    else if (x == y->l && y == z->l) {
      transplant(z, x);
      x->r = y;
      y->p = x;
      y->r = z;
      y->l = xr;
      if (xr) xr->p = y;
      z->p = y;
      z->l = yr;
      if (yr) yr->p = z;
    }
    else if (x == y->l) {
      transplant(z, x);
      x->l = z;
      x->r = y;
      y->p = x;
      y->l = xr;
      if (xr) xr->p = y;
      z->p = x;
      z->r = xl;
      if (xl) xl->p = z;
    }
    else {
      transplant(z, x);
      x->r = z;
      x->l = y;
      y->p = x;
      y->r = xl;
      if (xl) xl->p = y;
      z->p = x;
      z->l = xr;
      if (xr) xr->p = z;
    }
  }
    
  void list (vector<V> &res) { if (r) r->list(res); }
  int depth () const { return r ? r->depth() : 0; }
  int size () const { return r ? r->size() : 0; }

  Node *r;
};

#endif
