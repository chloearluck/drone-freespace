#include "hull.h"
#include "geometry3d.h"
#include "mink.h"
#include <cstring>
#include "io.h"

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

//now try point in cell...

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

  Polyhedron * robot = loadPoly("frustum.vtk");
  Polyhedron * room = loadPoly("tworooms.vtk");
  Polyhedron * msum = minkowskiSumFull(robot, room);
  savePoly(msum, "msum");

  return 0;
}