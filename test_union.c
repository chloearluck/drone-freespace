#include "hull.h"
#include "geometry3d.h"



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

void savePoly(Polyhedron * p, char * str) {
  ofstream out;
  out.open(str);
  if (out.is_open()) {
    writePolyhedronVTK (p->faces, out);
    out.close();
  } else {
    printf("could not write to file\n");
  }
}

// //take the union of 2 polyhedrons given at vtk files
// //first 2 args are name of vtk files, 3rd is random seed
// int main (int argc, char *argv[]) {
//   if (argc < 3) {
//     printf("not enough arguments\n");
//     return 1;
//   }

//   char * filename1 =  argv[1];
//   char * filename2 =  argv[2];

//   if (argc == 3) { 
//     unsigned seed = atoi(argv[3]);
//     srandom(seed);
//   }

//   Polyhedron * poly1 = loadPoly(filename1);
//   if (poly1 == NULL)
//     return 1;
//   Polyhedron * poly2 = loadPoly(filename2);
//   if (poly2 == NULL)
//     return 1;

//   Polyhedron * out = poly1->boolean(poly2, Union);

//   savePoly(out, "out.vtk");
// }

int main (int argc, char *argv[]) {
  Polyhedron * outerApprox = loadPoly("output/0-out.vtk");
  if (outerApprox == NULL) 
    return 1;

  for (int i=1; i<168; i++) {
    printf("%d of 168\n", i);
    char s[30];
    sprintf(s, "output/%d-out.vtk", i);
    Polyhedron * poly = loadPoly(s);
    outerApprox = outerApprox->boolean(poly, Union);
    sprintf(s, "output/union-%d.vtk", i);
    savePoly(outerApprox, s);
  }

  savePoly(outerApprox, "out.vtk");
}