#ifndef HULL
#define HULL

#include "polyhedron.h"

Polyhedron * convexHull (Polyhedron *a);

Polyhedron * convexHull (Points &pts, bool perturbed);

int * randomPermutation (int n);

int randomInteger (int lb, int ub);

Polyhedron * hullInit (Points &pts, bool perturbed);

class ConflictGraph {
  map<Vertex *, FaceSet *> vcon;
  map<Face *, VertexSet *> fcon;
 public:
  ~ConflictGraph ();
  void insert (Vertex *v, Face *f);
  void erase (Face *f);
  FaceSet * conflicts (Vertex *v);
  VertexSet * conflicts (Face *f);
  void update (Vertex *v, Face *f, VertexSet *vs);
};

ConflictGraph conflictGraphInit (Polyhedron *a);

bool visible (Vertex *v, Face *f);

void expandHull (Polyhedron *a, Vertex *v, ConflictGraph &cg);

HEdges horizon (const FaceSet &fs);

void expandHull (Polyhedron *a, Vertex *v, ConflictGraph &cg, HEdge *h);

#endif
