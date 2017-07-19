#include "freespace.h"

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

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout<<"not enough arguments"<<endl;
    return 1;
  }

  char * filename =  argv[1];

  // if (argc == 3) { 
    // unsigned seed = atoi(argv[2]);
    srandom(10);
  // }

  Polyhedron * poly = loadPoly(filename);
  if (poly == NULL)
    return 1;

  Polyhedron * obstacle;
  if (argc == 3) 
    obstacle = loadPoly(argv[2]); 
  else {
    double bounds[6] = { 5, 7, 5, 7, 5, 7};
    obstacle = box(bounds);
  }
  
  double theta = M_PI / 10;
  double bb_bounds[6] = {-20, 20, -20, 20, -20, 20};

  FreeSpace * fs  = new FreeSpace(poly, obstacle, theta, bb_bounds);
}