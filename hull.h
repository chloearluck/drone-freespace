#ifndef HULL
#define HULL

#include "io.h"

class HullEdge;

class HullFace;

class HullVertex {
 public:
  HullVertex (Point *p) : p(p), edge(0) {}
  PTR<Point> p;
  HullEdge *edge;
  vector<HullFace *> conflicts;
};

typedef vector<HullVertex *> HullVertices;

class HullEdge {
 public:
  HullEdge (HullVertex *tail, HullEdge *twin, HullEdge *next) 
    : tail(tail), twin(twin), next(next), face(0), flag(false) {}
  HullVertex * head () const { return twin->tail; }
  bool horizonEdge (HullVertex *v) const;

  HullVertex *tail;
  HullEdge *twin, *next;
  HullFace *face;
  bool flag;
};

typedef vector<HullEdge *> HullEdges;

class HullFace {
 public:
  HullFace (HullEdge *boundary) : boundary(boundary) {}
  bool visible (HullVertex *v) const;
  bool contains (HullVertex *v) const;
  bool conflict (HullVertex *v) const;

  HullEdge *boundary;
  HullVertices conflicts;
};

typedef vector<HullFace *> HullFaces;

class Hull {
 public:
  ~Hull ();
  HullVertex * addVertex (Point *p);
  HullEdge * addEdge (HullVertex *t, HullVertex *h);
  HullFace * addFace (HullEdge *e);
  void formHull ();
  void initHull ();
  void expandHull (int i);
  void findHorizon (HullVertex *v, HullEdges &h) const;
  HullEdge * findHorizonEdge (HullVertex *v) const;
  void addCone (HullVertex *v, HullEdges &h);
  void addConeFace (HullVertex *v, HullEdge *e);
  void removeVisible (int i) const;
  Polyhedron * polyhedron () const;

  HullVertices vertices;
  HullEdges edges;
  HullFaces faces;
};

Polyhedron * convexHull (Polyhedron *a);

Polyhedron * convexHull (const Points &pts);

void randomPermutation (int n, int *p);

int randomInteger (int lb, int ub);

#endif
