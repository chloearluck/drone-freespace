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
  double d1[6] = {-5, 5, -5, 5, -5, 5};
  Polyhedron * big_box = box(d1);

  double d2[6] = {1, 4, -4, 4, -4, 4};
  Polyhedron * room1 = box(d2);
  
  double d3[6] = {-4, -1, -4, 4, -4, 4};
  Polyhedron * room2 = box(d3);

  double d4[6] = {-2, 2, -.45, .45, -.9, 0};
  Polyhedron * hall1 = box(d4);

  double d5[6] = {-2.95, -2.05, 3, 6, 1, 2.5};
  Polyhedron * hall2 = box(d5);

  Polyhedron * out = big_box;
  out = out->boolean(room1, Complement);
  out = out->boolean(room2, Complement);
  out = out->boolean(hall1, Complement);
  out = out->boolean(hall2, Complement);
  simplify(out, 1e-6);
  savePoly(out, "tworooms.vtk");


}