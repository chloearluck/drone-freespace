#ifndef IO
#define IO

#include "polyhedron.h"
#include <fstream>
#include <sstream>

Polyhedron * readPolyhedron (istream &istr, bool perturbed = true);

void readCells (istream &istr, Polyhedron *a);

Cell * readCell (istream &istr, Polyhedron *a);

Shell * readShell (istream &istr, Polyhedron *a);

typedef pair<Vertex *, ID> VIPair;

typedef map<Vertex *, ID> VIMap;

typedef pair<Face *, ID> FIPair;

typedef map<Face *, ID> FIMap;

void writePolyhedron (Polyhedron *a, ostream &ostr);

void writeCells (Polyhedron *a, ostream &ostr);

void writeShell (Shell *s, FIMap &fimap, ostream &ostr);

Polyhedron * readPolyhedronVTK (istream &istr, bool perturbed = true);

void skipComments (istream &istr);

void writePolyhedronVTK (Polyhedron *a, ostream &ostr);

void writePolyhedronVTK (const Faces &fa, ostream &ostr);

int getPoint (VIMap &vimap, vector<PV3> &pts, Vertex *v);

void outputVTK (const vector<PV3> &pts, const vector<ID> &data,
		int ptype, ostream &ostr);

Polyhedron * readPolyhedronOBJ (istream &istr, bool perturbed = true);

bool readPointOBJ (istream &istr, double &x, double &y, double &z);

bool readTriangleOBJ (istream &istr, int n, int &i, int &j, int &k);

void readIndexOBJ (istream &istr, int n, int &i);

void addTriangleOBJ (Polyhedron *a, Vertex *u, Vertex *v, Vertex *w);

#endif
