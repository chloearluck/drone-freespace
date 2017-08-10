#ifndef PACK
#define PACK

#include "io.h"
#include "mink.h"

double pack3 (Polyhedron *a, Polyhedron *b, Polyhedron *c, double minsep, PTR<Point> t[3]);

double maxBBoxDimension (Polyhedron *a);

bool pack3 (Polyhedron *a, Polyhedron *b, Polyhedron *c, PTR<Point> t[3], double bs);

void freeSpace (Polyhedron *a, double *box, double *res);

void freeSpace (double *a, double *b, double *res);

void translateBox (double *a, PTR<Point> t, double *res);

void intersectBoxes (double *a, double *b, double *res);

PTR<Point> midpointBox (double *a);

PTR<Point> interiorPoint (Polyhedron *a);

Face * largestFace (Polyhedron *a);

double area (Face *f);

void pack3output (Polyhedron *a, Polyhedron *b, Polyhedron *c,
		  PTR<Point> t[3], ostream &ostr);

#endif
