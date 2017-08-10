#ifndef IO
#define IO

#include "polyhedron.h"
#include <sstream>

class FaceRecord {
 public:
  FaceRecord () : m1(0), m2(0), p(0), o1(false), o2(false) {}
  IVectors b;
  ID m1, m2;
  Plane *p;
  bool o1, o2;
};

typedef vector<FaceRecord> FaceRecords;

Polyhedron * readPolyhedron (istream &istr, bool tflag = false, bool check = false);

void readVertices (istream &istr, Polyhedron *a);

void readAttributes (istream &istr, Polyhedron *a);

void skipComments (istream &istr);

void readFaceRecords (istream &istr, FaceRecords &frs, Polyhedron *a);

FaceRecord readFaceRecord (istream &istr);

void formFaces (const FaceRecords &frs, Polyhedron *a, bool tflag, bool check);

void faceVertices (const IVectors &iv, Polyhedron *a, VVertices &reg);

bool facePlane (const VVertices &reg, Polyhedron *a);

int projectionCoordinate (const Vertices &ve);

Polyhedron * readPolyhedronVTK (istream &istr, bool perturb = true);

typedef pair<Vertex *, ID> VIPair;

typedef map<Vertex *, ID> VIMap;

void writePolyhedron (Polyhedron *a, ostream &ostr);

void writeHFaces (Polyhedron *a, ostream &ostr);

void writePolyhedronVTK (const Faces &fa, ostream &ostr);

int getPoint (VIMap &vimap, PV3s &pts, Vertex *v);

void outputVTK (const PV3s &pts, const IVector &data, bool pflag, 
		ostream &ostr);

void writePolyhedronVTK (const HFaces &fa, ostream &ostr);

Polyhedron * readPolyhedronSTL (istream &istr);

void ptriangles (const Triangles &tr, ostream &ostr);

string outi (int i);

void plines (const vector<PV3s> &lines, int i);

#endif
