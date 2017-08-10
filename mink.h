#ifndef MINK
#define MINK

#include "polyhedron.h"

class Convolution : public Polyhedron {
  EFVMap efvmap;
  FFFVMap fffvmap;
 public:
  Convolution () { oneWayInt = true; }
  Face * minkowskiInit (FaceSet &fdone, Octree<Face *> *octree,
			Polyhedron *a, VVMap &vvmap);
  Face * rmaxFace (Vertex *&v);
  void subdivide (Face *f, FaceSet &fdone, Octree<Face *> *octree,
		  Polyhedron *a, VVMap &vvmap);
  void intersectFF (Face *f, const FaceSet &fdone, Octree<Face *> *octree);
  void expand (FaceSet &fdone, Octree<Face *> *octree,
	       Polyhedron *a, VVMap &vvmap, Face *f, Faces &st);
  void neighborFaces (HEdge *e, Faces &fa) const;
};

Primitive2(NormalOrderR, Plane *, p, Plane *, q);

Polyhedron * minkowskiSum (Polyhedron *a, Polyhedron *b);

void sortHEdges (Polyhedron *a, int nf);

class UnitVector : public Point {
  Point *a;
  PV3 calculate () { return a->getP().unit(); }
 public:
  UnitVector (Point *a) : a (a) {}
};

class GreatCircleN : public Point {
  Point *t, *h;
  PV3 calculate () { return t->getP().cross(h->getP()); }
 public:
  GreatCircleN (Point *t, Point *h) : t(t), h(h) {}
};

class GreatCircleMinX : public Point {
  Point *n;
  PV3 calculate () {
    PV3 p = n->getP();
    return PV3(- p.y*p.y - p.z*p.z, p.y*p.x, p.z*p.x);
  }
 public:
  GreatCircleMinX (Point *n) : n(n) {}
};

class GreatCircleMinY : public Point {
  Point *n;
  PV3 calculate () {
    PV3 p = n->getP();
    return PV3(p.x*p.y, - p.x*p.x - p.z*p.z, p.z*p.y);
  }
 public:
  GreatCircleMinY (Point *n) : n(n) {}
};

class GreatCircleMinZ : public Point {
  Point *n;
  PV3 calculate () {
    PV3 p = n->getP();
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
  PV3 calculate () { return f->getP()->getN(); }
 public:
  FNormal (Face *f) : f(f) {}
};

void sphereBBox (Edge *e, double *bbox);

void sphereBBox (Point *t, Point *h, double *bbox);

void sphereBBox (const HEdges &ed, double *bbox);

int coordinate (Point *a, int c);

void merge (double *bbox1, double *bbox2);

class BSPElt {
 public:
  BSPElt (Vertex *v, const HEdges &ed);
  BSPElt (Edge *e);
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

Primitive4(DegenerateConflict, Point *, a, Point *, b, Point *, c, Point *, d);

bool convexCone (Vertex *v, HEdges &hedges);

MinkHullFace * initHull (HEdges &hedges);

int circulationEEE (HEdge *e, HEdge *f, HEdge *g);

MinkHullFace * updateHull (MinkHullFace *hull, HEdge *e, bool &flag);

MinkHullFace * updateHullAux (MinkHullFace *fs, MinkHullFace *fe, HEdge *e);

void deleteHull (MinkHullFace *hull);

bool convexOrder (const HEdges &hedges);

Polyhedron * triangulate (const FaceSet &fs);

Polyhedron * minkowskiSumFull (Polyhedron *a, Polyhedron *b);

bool minkowskiShell (Polyhedron *a, Polyhedron *b, Shell *s);

typedef map<VVPair, Vertex *> VVPairVMap;

Convolution * convolution (Polyhedron *a, Polyhedron *b);

void sumVF (Polyhedron *a, Polyhedron *b, bool avflag, FaceDescSet &fds,
	    VVPairVMap &vmap, Polyhedron *con);

void convexVertices (Polyhedron *a, BSPElts &elts);

bool compatibleVF (HEdges &ed, Face *f);

void sumVF (Vertex *v, Face *f, bool avflag, FaceDescSet &fds,
	    VVPairVMap &vmap, Polyhedron *con);

Primitive2(InnerProductEF, HEdge *, e, Face *, f);

Vertex * sumVV (Vertex *a, Vertex *b, bool aflag, VVPairVMap &vmap,
		Polyhedron *con);

void sumEE (Polyhedron *a, Polyhedron *b, FaceDescSet &fds,
	    VVPairVMap &vmap, Polyhedron *con);

void convexEdges (Polyhedron *a, BSPElts &elts);

int convexEdge (Edge *e);

Primitive2(ConvexEdge, HEdge *, e1, HEdge *, e2);

bool compatibleEE (Edge *e, Edge *f, bool &aflag);

Primitive2(CompatibleEdge, Edge *, e, Face *, f);

void sumEE (Edge *e, Edge *f, bool aflag, FaceDescSet &fds,
	    VVPairVMap &vmap, Polyhedron *con);

#endif
