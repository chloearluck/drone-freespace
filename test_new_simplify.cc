#include "polyhedron.h"
#include "io.h"
#include <cstring>
#include "simplify.h"

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
  Polyhedron * test = loadPoly(argv[1]);
  Polyhedron * test2 = loadPoly(argv[2]);
  Polyhedron * test3 = test->boolean(test2, Union);
  simplify(test3, 1e-6);
}