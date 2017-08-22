//minkowski sum of rect-drone and room is incorrect

#include "hull.h"
#include "geometry3d.h"
#include "mink.h"
#include <cstring>

Polyhedron * loadPoly(char * str) {
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


void savePoly(Polyhedron * p, char * filename) {
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

void saveShell(HFaces hf,  char * filename) {
  int n = strlen(filename);
  char str[n+9];
  strncpy(str, filename, n);
  strncpy(str+n, "-out.vtk", 9);

  ofstream out;
  out.open(str);
  if (out.is_open()) {
    writePolyhedronVTK (hf, out);
    out.close();
  } else {
    cout<<"could not write to file"<<endl;
  }
}


void saveWithShells(Polyhedron * poly, char * filename) {
  savePoly(poly, filename);
  poly->formCells();
  for (int i=1; i<poly->cells.size(); i++) {
    Cell * c = poly->cells[i];
    for (int j=0; j<c->nShells(); j++) {
      Shell * s = c->getShell(j);
      char str[50];
      sprintf(str, "%s-%d-%d", filename, i, j);
      saveShell(s->getHFaces(), str);
    }
  }
}

int main (int argc, char *argv[]) {
  if (argc == 2) { 
    unsigned seed = atoi(argv[1]);
    srandom(seed);
  }

  Polyhedron * robot = loadPoly("rect-drone.vtk");
  Polyhedron * room = loadPoly("room.vtk");
  Polyhedron * mSum = minkowskiSumFull(robot, room);
  savePoly(mSum, "mSum");
  
  return 0;
}