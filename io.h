#ifndef IO
#define IO

#include "polyhedron.h"
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

void writePolyhedronOBJ (const Faces &fa, ostream &ostr);

int getPoint (VIMap &vimap, PV3s &pts, Vertex *v);

void outputVTK (const PV3s &pts, const IVector &data, bool pflag, 
		ostream &ostr);

void outputOBJ (const PV3s &pts, const IVector &data, bool pflag, 
    ostream &ostr);

#endif
