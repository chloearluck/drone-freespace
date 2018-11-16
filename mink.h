#ifndef MINK
#define MINK

#include "polyhedron.h"

Polyhedron * minkowskiSum (Polyhedron *a, Polyhedron *b);

Shells minkowskiShells (Polyhedron *a, Polyhedron *b, const Shells &sh);

Shells minkowskiShellsAux (Polyhedron *a, Polyhedron *b, const Shells &sh);

void minkowskiShellsT (void *ptr);

class MSData {
 public:
  unsigned int i, is, ie;
  Polyhedron *a, *b;
  const Shells *sh;
  Shells res;
};

bool minkowskiShell (Polyhedron *a, Polyhedron *b, Shell *s);

class UnitVector : public Point {
  Point *a;

  PV3 calculate () { return a->get().unit(); }
 public:
  UnitVector (Point *a) : a (a) {}
};

class GreatCircleN : public Point {
  Point *t, *h;

  PV3 calculate () { return t->get().cross(h->get()); }
 public:
  GreatCircleN (Point *t, Point *h) : t(t), h(h) {}
};

class GreatCircleMinX : public Point {
  Point *n;
  
  PV3 calculate () {
    PV3 p = n->get();
    return PV3(- p.y*p.y - p.z*p.z, p.y*p.x, p.z*p.x);
  }
 public:
  GreatCircleMinX (Point *n) : n(n) {}
};

class GreatCircleMinY : public Point {
  Point *n;
  
  PV3 calculate () {
    PV3 p = n->get();
    return PV3(p.x*p.y, - p.x*p.x - p.z*p.z, p.z*p.y);
  }
 public:
  GreatCircleMinY (Point *n) : n(n) {}
};

class GreatCircleMinZ : public Point {
  Point *n;
  
  PV3 calculate () {
    PV3 p = n->get();
    return PV3(p.x*p.z, p.y*p.z, - p.x*p.x - p.y*p.y);
  }
 public:
  GreatCircleMinZ (Point *n) : n(n) {}
};

class EENormal : public Point {
  HEdge *e, *f;

  PV3 calculate () { return e->getU().cross(f->getU()); }
 public:
  EENormal (HEdge *e, HEdge *f) : e(e), f(f) {}
};

class FNormal : public Point {
  Face *f;
  
  PV3 calculate () { return f->getP()->get().n; }
 public:
  FNormal (Face *f) : f(f) {}
};

class TripleProduct :public Primitive {
  Point *a, *b, *c;

  Parameter calculate () {
    return a->get().tripleProduct(b->get(), c->get());
  }
 public:
  TripleProduct (Point *a, Point *b, Point *c) : a(a), b(b), c(c) {}
};

void sphereBBox (Edge *e, double *bbox);

void sphereBBox (Point *t, Point *h, double *bbox);

void sphereBBox (const HEdges &ed, double *bbox);

int coordinate (Point *a, int c);

class BSPElt {
 public:
  BSPElt (Vertex *v, const HEdges &ed, bool convex = true);
  BSPElt (Edge *e, bool convex = true);
  BSPElt (Face *f);
  BSPElt (const BSPElt &x);
  bool compatible (const BSPElt &e) const { return (l & e.l) == 0u; }
  int side (Point *r) const;

  union {
    Vertex *v;
    Edge *e;
    Face *f;
  } d;
  HEdges ed;
  ID l;
  double bbox[6];
  bool convex;
};

typedef vector<BSPElt> BSPElts;

void BSPTree (BSPElts &aelts, BSPElts &belts, BSPElts &ea, BSPElts &eb,
	      int nmax = 40, int dmax = 20, ID c = 0u);

void BSPPartition (BSPElts &elts, Point *r, ID c, BSPElts &elts1, BSPElts &elts2);

void BSPLeaf (const BSPElts &aelts, const BSPElts &belts, BSPElts &ea, 
	      BSPElts &eb);

class MinkHullFace {
 public:
  MinkHullFace (HEdge *e, MinkHullFace *prev, MinkHullFace *next)
    : e(e), prev(prev), next(next) {}
  bool inCset (HEdge *f) const
  { return find(cset.begin(), cset.end(), f) != cset.end(); }
  void updateCset (MinkHullFace *h, HEdge *f);
  bool conflict (HEdge *f) const;
  void cone (HEdges &hedges) const;

  HEdges cset;
  HEdge *e;
  MinkHullFace *prev, *next;
};

class DegenerateConflict : public Primitive {
  Point *a, *b, *c, *d;

  Parameter calculate () {
    PV3 u = b->get() - a->get(), v = c->get() - a->get(),
      w = d->get() - a->get(), x = u.cross(v);
    Parameter k = x.tripleProduct(v, w);
    return k.sign() == 1 ? k : x.tripleProduct(w, u);
  }
 public:
  DegenerateConflict (Point *a, Point *b, Point *c, Point *d)
    : a(a), b(b), c(c), d(d) {}
};

bool convexCone (Vertex *v, HEdges &hedges);

MinkHullFace * initHull (HEdges &hedges);

int circulationEEE (HEdge *e, HEdge *f, HEdge *g);

MinkHullFace * updateHull (MinkHullFace *hull, HEdge *e, bool &flag);

MinkHullFace * updateHullAux (MinkHullFace *fs, MinkHullFace *fe, HEdge *e);

void deleteHull (MinkHullFace *hull);

bool convexOrder (const HEdges &hedges);

typedef map<pair<Vertex *, Vertex *>, Vertex *> VVVMap;

Polyhedron * convolution (Polyhedron *a, Polyhedron *b);

void sumVF (Polyhedron *a, Polyhedron *b, bool avflag, VVVMap &vmap,
	    Polyhedron *con);

void convexVertices (Polyhedron *a, BSPElts &elts);

bool compatibleVF (HEdges &ed, Face *f);

void sumVF (Vertex *v, Face *f, bool avflag, VVVMap &vmap, Polyhedron *con);

class InnerProductEF : public Primitive {
  HEdge *e;
  Face *f;

  Parameter calculate () {
    return e->getU().dot(f->getP()->get().n);
  }
 public:
  InnerProductEF (HEdge *e, Face *f) : e(e), f(f) {}
};

Vertex * sumVV (Vertex *a, Vertex *b, bool aflag, VVVMap &vmap, Polyhedron *con);

void sumEE (Polyhedron *a, Polyhedron *b, VVVMap &vmap, Polyhedron *con);

void convexEdges (Polyhedron *a, BSPElts &elts);

int convexEdge (Edge *e);

class ConvexEdge : public Primitive {
  HEdge *e1, *e2;

  Parameter calculate () {
    return e1->getU().tripleProduct(e1->getF()->getP()->get().n,
				    e2->getF()->getP()->get().n);
  }
 public:
  ConvexEdge (HEdge *e1, HEdge *e2) : e1(e1), e2(e2) {}
};

bool compatibleEE (Edge *e, Edge *f, bool &aflag);

void sumEE (Edge *e, Edge *f, bool aflag, VVVMap &vmap, Polyhedron *con);

#endif
