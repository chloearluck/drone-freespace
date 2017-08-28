#include "path3d.h"
#include "io.h"

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
    cout<<"could not read from file"<<endl;
    return NULL;
  }

  return poly;
}

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout<<"not enough arguments"<<endl;
    return 1;
  }

  char * filename =  argv[1];

  if (argc == 3) {
    unsigned seed = atoi(argv[2]);
    srandom(seed);
  }

  Polyhedron * blockspace = loadPoly(filename);

  PTR<Point> start = new InputPoint(5, 5, 2);
  PTR<Point> end = new InputPoint(15, 15, 2);
  
  Points path;

  findPath(blockspace, start, end, path);
}
