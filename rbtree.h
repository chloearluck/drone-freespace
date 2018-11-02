#ifndef RBTREE
#define RBTREE

#include <vector>

using namespace std;

template <class V>
class RBTree {
 public:
  
  class Node {
  public:
    Node () : v(0), p(0), l(0), r(0), red(false) {}
    
    Node (V v, Node *p, Node *l, Node *r, bool red)
      : v(v), p(p), l(l), r(r), red(red) {}

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
    bool red;
  };

  Node *r, *nil;

  RBTree () : nil(new Node) { r = nil; }
  
  ~RBTree () { erase(r); delete nil; }
    
  void erase (Node *x) {
    if (x != nil) {
      erase(x->l);
      erase(x->r);
      delete x;
    }
  }

  Node * min (Node *x) {
    while (x->l != nil)
      x = x->l;
    return x;
  }
  
  bool find (V v) const {
    Node *x = r, *y = nil;
    int s;
    while (x != nil) {
      s = v->order(x->v);
      if (s == 0) {
	return true;
      }
      y = x;
      x = s == 1 ? x->l : x->r;
    }
    return false;
  }

  Node * insert (V v) {
    Node *x = r, *y = nil;
    int s;
    while (x != nil) {
      s = v->order(x->v);
      if (s == 0) {
	return x;
      }
      y = x;
      x = s == 1 ? x->l : x->r;
    }
    Node *z = new Node(v, y, nil, nil, true);
    if (y == nil)
      r = z;
    else if (s == 1)
      y->l = z;
    else
      y->r = z;
    insertFixup(z);
    return z;
  }

  void insertFixup (Node *z) {
    while (z->p->red)
      if (z->p == z->p->p->l) {
	Node *y = z->p->p->r;
	if (y->red) {
	  z->p->red = false;
	  y->red = false;
	  z->p->p->red = true;
	  z = z->p->p;
	}
	else {
	  if (z == z->p->r) {
	    z = z->p;
	    leftRotate(z);
	  }
	  z->p->red = false;
	  z->p->p->red = true;
	  rightRotate(z->p->p);
	}
      }
      else {
	Node *y = z->p->p->l;
	if (y->red) {
	  z->p->red = false;
	  y->red = false;
	  z->p->p->red = true;
	  z = z->p->p;
	}
	else {
	  if (z == z->p->l) {
	    z = z->p;
	    rightRotate(z);
	  }
	  z->p->red = false;
	  z->p->p->red = true;
	  leftRotate(z->p->p);
	}
      }
    r->red = false;
  }
  
  void leftRotate (Node *x) {
    Node *y = x->r;
    x->r = y->l;
    if (y->l != nil)
      y->l->p = x;
    y->p = x->p;
    if (x->p == nil)
      r = y;
    else if (x == x->p->l)
      x->p->l = y;
    else
      x->p->r = y;
    y->l = x;
    x->p = y;
  }

  void rightRotate (Node *x) {
    Node *y = x->l;
    x->l = y->r;
    if (y->r != nil)
      y->r->p = x;
    y->p = x->p;
    if (x->p == nil)
      r = y;
    else if (x == x->p->l)
      x->p->l = y;
    else
      x->p->r = y;
    y->r = x;
    x->p = y;
  }

  void remove (Node *z) {
    Node *y = z, *x = nil;
    bool yred = y->red;
    if (z->l == nil) {
      x = z->r;
      transplant(z, z->r);
    }
    else if (z->r == nil) {
      x = z->l;
      transplant(z, z->l);
    }
    else {
      y = min(z->r);
      yred = y->red;
      x = y->r;
      if (y->p == z)
	x->p = y;
      else {
	transplant(y, y->r);
	y->r = z->r;
	y->r->p = y;
      }
      transplant(z, y);
      y->l = z->l;
      y->l->p = y;
      y->red = z->red;
    }
    delete z;
    if (!yred)
      deleteFixup(x);
  }

  void transplant (Node *u, Node *v) {
    if (u->p == nil)
      r = v;
    else if (u == u->p->l)
      u->p->l = v;
    else
      u->p->r = v;
    v->p = u->p;
  }

  void deleteFixup (Node *z) {} // to do

  void list (vector<V> &res) { if (r) r->list(res); }
  
  int depth () const { return r ? r->depth() : 0; }
  
  int size () const { return r ? r->size() : 0; }
};

#endif
