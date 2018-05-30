#include "freespace.h"
#include "io.h"
#include <queue>

double tanHalfAngle (int n) {
  return tan((1.0 + 1.0e-8) * M_PI / n);
}

int main (int argc, char *argv[]) {
  int numRotations = 40;
  PTR<Object<Parameter> > tan_half_angle = new InputParameter(tanHalfAngle(numRotations));
  
  Parameter::enable();

  if (argc < 2) {
    cout<<"not enough arguments"<<endl;
    return 1;
  }

  char * filename =  argv[1];

  if (argc == 4) { 
    unsigned seed = atoi(argv[3]);
    srandom(seed);
  }

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
  
  FreeSpace * fs  = new FreeSpace(poly, obstacle, tan_half_angle, numRotations);
}
