#include "hull.h"
#include "geometry3d.h"
#include "mink.h"
#include <cstring>

Polyhedron * loadPoly(char * filename) {
  int n = strlen(filename);
  char str[n+5];
  strncpy(str, filename, n);
  strncpy(str+n, ".vtk", 5);

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


int main (int argc, char *argv[]) {
if (argc < 3) {
    printf("not enough arguments\n");
    return 1;
  }

  char * filename1 =  argv[1];
  char * filename2 =  argv[2];

  if (argc == 4) { 
    unsigned seed = atoi(argv[3]);
    srandom(seed);
  }

  Polyhedron * poly1 = loadPoly(filename1);
  Polyhedron * poly2 = loadPoly(filename2);

  Polyhedron * poly3 = minkowskiSum (poly1, poly2);

  return 0;
}
