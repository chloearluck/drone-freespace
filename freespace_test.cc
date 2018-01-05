#include "freespace.h"
#include "io.h"
#include <queue>

int main (int argc, char *argv[]) {
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
  
  double theta = M_PI / 20;

  FreeSpace * fs  = new FreeSpace(poly, obstacle, theta);
}
