#include "freespace-graph.h"

int main (int argc, char *argv[]) {

  std::vector<Polyhedron *> blockspaces;
  char s[50];
  for (int i=0; i<=40; i++) {
    sprintf(s, "output/sum%02d-out.vtk", i);
    Polyhedron * p = loadPoly(s);
    blockspaces.push_back(p);
  }
  double theta = M_PI / 20;
  double clearance_unit = 0.25;
  int num_levels = 3;
  FreeSpaceGraph graph(blockspaces, theta, clearance_unit, num_levels);

  return 0;
}