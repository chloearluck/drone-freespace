#include "polyhedron.h"
#include "io.h"
#include <cstring>
#include "simplify.h"
#include "mink.h"

Polyhedron * loadPoly(const char * filename) {
  Polyhedron * poly;
  ifstream infile (filename);
  if (infile.is_open()) {
    poly = readPolyhedronVTK (infile);
    infile.close();
  } else {
    cout<<"could not read from file"<<endl;
    return NULL;
  }

  return poly;
}

void savePoly(Polyhedron * p, const char * filename) {
  ofstream out;
  out.open(filename);
  if (out.is_open()) {
    writePolyhedronVTK (p->faces, out);
    out.close();
  } else {
    cout<<"could not write to file"<<endl;
  }
}

int main (int argc, char *argv[]) {
  Parameter::enable();
  
  unsigned i = atoi(argv[1]);
  srandom(i);

  Polyhedron * p1 = loadPoly("two_room_output/sum18-out.vtk");
  simplify(p1, 1e-6);
  Polyhedron * p2 = loadPoly("two_room_output/sum19-out.vtk");
  simplify(p2, 1e-6);
  Polyhedron * p3 = p1->boolean(p2, Union);
  simplify(p3, 1e-6);
}
