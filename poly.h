/*
  ACP (Adaptive Controlled Precision) Library
  for robust computational geometry

  Copyright (c) 2013-07-15
  Victor Joseph Milenkovic
  University of Miami
  vjm@miami.edu
  Elisha Sacks
  Purdue University
  eps@cs.purdue.edu

  This file contains code described in

Robust Complete Path Planning in the Plane
Victor Milenkovic, Elisha Sacks, and Steven Trac
Proceedings of the Workshop on the Algorithmic Foundations of Robotics (WAFR)
pages 37-52, 2012

   This code is under development.
   It is free for use for academic and research purposes.
   Permission from the authors is required for commercial use.
*/

#ifndef POLY_H
#define POLY_H

#include "object.h"
#include "ppoly.h"

namespace acp {

class SubKnot : public Object<Parameter> {
  PTR<Object<Parameter>> l, r;
  unsigned int i, n;
  
  Parameter calculate () {
    return (l->get() * (n - i) + r->get() * i) / n;
  }
public:
  SubKnot (Object<Parameter>* l, Object<Parameter>* r, unsigned int i, unsigned int n)
    : l(l), r(r), i(i), n(n) {}
};

class DInt {
 public:
  PTR<Object<Parameter>> l, u;
  unsigned int n;
   
  DInt (Object<Parameter>* l, Object<Parameter>* u, unsigned int n) : l(l), u(u), n(n) {}
};

class Root : public Object<Parameter> {
protected:
  PTR<Object<PPoly>> p;
public:
  Root (Object<PPoly> *p) : p(p) {}
  PTR<Object<PPoly>> getPoly () { return p; }
};

typedef vector<PTR<Root>> Roots;

class PolySolver {
  Roots linearRoots (Object<Parameter>* l, Object<Parameter>* u);
  Roots quadraticRoots ();
  Roots quadraticRoots (Object<Parameter>* l, Object<Parameter>* u);
  Roots descartesRoots (Object<Parameter>* l, Object<Parameter>* u);
  bool descartes1 (vector<DInt> &st, const DInt &i, int v, bool lflag);
  bool descartes2 (vector<DInt> &st, const DInt &i, int v);
  int descartes2int (const DInt &i, int v);
  PTR<Object<PPoly>> poly;
  PPoly get () { return poly->get(); }
  PPoly getApprox (double e) { return poly->getApprox(e); }
  Primitive1(Degree, Object<PPoly>*, poly);
 public:
  PolySolver (Object<PPoly> *poly) : poly(poly) {}
  int deg () { return Degree(poly); }
  Roots getRoots ();
  Roots getRoots (Object<Parameter>* l, Object<Parameter>* u);
  Roots getRoots (double l, double u) {
    PTR<Object<Parameter>> pl = new Object<Parameter>(Parameter::constant(l));
    PTR<Object<Parameter>> pr = new Object<Parameter>(Parameter::constant(u));
    return getRoots(pl, pr);              
  }
};

class LinearRoot : public Root {
  Parameter calculate () {
    PPoly q = p->get();
    return - q[0]/q[1];
  }
 public:
  LinearRoot (Object<PPoly> *p) : Root(p) {}
};

class QuadraticRoot : public Root {
  bool flag;
  
  Parameter calculate () {
    PPoly q = p->get();
    const Parameter &a = q[2], &b = q[1], &c = q[0];
    Parameter d = (b*b - 4.0*a*c).sqrt();
    if (b.sign() == -1)
      return flag ? 2.0*c/(d-b) : 0.5*(d-b)/a;
    return flag ? -0.5*(d+b)/a : -2.0*c/(d+b);
  }
  
 public:
  QuadraticRoot (Object<PPoly> *p, bool flag) : Root(p), flag(flag) {}
};

class PolyRoot : public Root {
  PTR<Object<Parameter>> l, u;
  Parameter calculate () {
    const PPoly &q = p->get();
    int slb = q.value(l->get()).sign(); 
    int sub = q.value(u->get()).sign();
    Parameter lu = l->get().ubP();
    Parameter ul = u->get().lbP();
    if (uninitialized())
      return q.shrinkBracket(lu.interval(ul), sub);
    Parameter x = getCurrentP();
    Parameter xl(x.lbP());
    Parameter xu(x.ubP());

    return q.shrinkBracket(Parameter((xl-lu).sign(false) > 0 ? xl : ul,
                             (xu-ul).sign(false) < 0 ? xu : lu), sub);
  }
 public:
  PolyRoot (Object<PPoly> *p, Object<Parameter>* l, Object<Parameter>* u) : Root(p), l(l), u(u) {}
};   

class CauchyBound : public Object<Parameter> {
  PTR<Object<PPoly>> p;
  
  Parameter calculate () {
    PPoly q = p->get();
    q.checkDegree();
    Parameter b = Parameter::constant(1.0), qd = q.lc();
    for (int i = 0; i < q.deg(); ++i) {
      Parameter bi = (q.a[i]/qd).abs();
      if (b < bi)
	b = bi;
    }
    return Parameter(1.0 + b);
  }
  
 public:
  CauchyBound (Object<PPoly> *p) : p(p) {}
};

Primitive1(QuadraticRoots, Object<PPoly> *, p);

Primitive2(Order, Object<Parameter> *, a, Object<Parameter> *, b);

Primitive3(Descartes, Object<PPoly> *, p, Object<Parameter> *, l, Object<Parameter> *, u);

typedef PTR<Root> RootPTR;

class Root2;

class Root2PTR : public PTR<Root2> {
public:
  Root2PTR () {}
  Root2PTR (Root2 *r) : PTR<Root2>(r) {}
  operator PTR<Object<PV2>> () const;
};

typedef Object<PPoly2> Poly2;

class Root2 : public Object<PV2> {
  PV2 r;
  PTR<Poly2> f, g;
public:
  Root2(Poly2 *f, Poly2 *g, PV2 &r) : r(r), f(f), g(g) {}
  PTR<Poly2> getF () { return f; }
  PTR<Poly2> getG () { return g; }
private:
  PV2 calculate () { 
    if (uninitialized())
      return polish(r);
    else
      return polish(getCurrentP());
  }

  //order of roots in v: bot, right, top, left
  // 0-3 for f, 4-7 for g
  class RootBoundary {
  public:
    PV2 I;
    vector<vector<RootPTR>> v;
  };

  /*
  //order of roots in v: bot, right, top, left
  typedef struct {
    PV2 I;
    vector< vector<RootPTR> > v;
  } RootBoundary;
  */

  static bool polishing;

  PV2 polish (PV2 p);
  void newton(RootBoundary &I);
  bool subdivide(RootBoundary &I);
  std::vector<int> intersections(vector<vector<RootPTR>> &rbv);
  int parity(const std::vector<int> &alt);
  void setRoots(RootBoundary &rb);
  RootBoundary splitHoriz(RootBoundary &rb, Parameter c);
  RootBoundary splitVert(RootBoundary &rb, Parameter c);
  bool solve(PV2, PV2, PV2, PV2 &);
};

inline Root2PTR::operator PTR<Object<PV2>> () const {
  return PTR<Object<PV2>>(operator Root2 *());
}

}
#endif
