#include "hull.h"
#include "geometry3d.h"
#include "mink.h"
#include <cstring>
#include "io.h"
#include "simplify.h"

Polyhedron * loadPoly(const char * str) {
  Polyhedron * poly;
  ifstream infile (str);
  if (infile.is_open()) {
    poly = readPolyhedronVTK (infile);
    infile.close();
  } else {
    printf("could not read from file\n");
    return NULL;
  }

  return poly;
}

void savePoly(Polyhedron * p, const char * filename) {
  int n = strlen(filename);
  char str[n+9];
  strncpy(str, filename, n);
  strncpy(str+n, "-out.vtk", 9);

  ofstream out;
  out.open(str);
  if (out.is_open()) {
    writePolyhedronVTK (p->faces, out);
    out.close();
  } else {
    cout<<"could not write to file"<<endl;
  }
}

int main (int argc, char *argv[]) {
  Parameter::enable();

  if (argc == 2) { 
    unsigned seed = atoi(argv[1]);
    srandom(seed);
  }

  Polyhedron * space = loadPoly("sum37-out.vtk");
  simplify(space, 1e-6);
  Polyhedron * ball = loadPoly("sphere.vtk");
  simplify(ball, 1e-6);
  Polyhedron * unit_ball = ball->scale(0.05);
  simplify(unit_ball, 1e-6);
  Polyhedron * msum = minkowskiSumFull(space, unit_ball);
  msum->formCells();
  cout<<msum->cells.size()<<" cells"<<endl;
  savePoly(msum, "test");
}