#include "freespace-graph.h"

void testCreateGraph(const char * blockspace_dir, const char * graph_dir, double theta, double clearance_unit, int num_levels) {
  int numRotations = floor(2*M_PI/theta);
  
  // std::vector<Polyhedron *> close_blockspaces;
  char s[50];
  // for (int i=0; i<=40; i++) {
  //   sprintf(s, "%s/close%02d-out.vtk", blockspace_dir, i);
  //   Polyhedron * p = loadPolyVTK(s);
  //   close_blockspaces.push_back(p);
  // }

  std::vector<Polyhedron *> rough_blockspaces;
  for (int i=0; i<=40; i++) {
    sprintf(s, "%s/rough%02d-out.vtk", blockspace_dir, i);
    Polyhedron * p = loadPolyVTK(s);
    rough_blockspaces.push_back(p);
  }

  FreeSpaceGraph graph(rough_blockspaces, rough_blockspaces, theta, clearance_unit, num_levels, graph_dir);
}

void testSearchGraph(const char * graph_dir, PTR<Point> start, PTR<Point> end, int startRotationIndex, int endRotationIndex) {
  FreeSpaceGraph test(graph_dir);
  std::vector<std::pair<PTR<Point>, double > > path;

  Polyhedron * p = loadPoly((std::string(graph_dir) + "/0-" + std::to_string(startRotationIndex) + ".tri").c_str());
  p->computeWindingNumbers();
  assert(!p->contains(start));
  int cell_index = p->containingCell(start);
  FreeSpaceGraph::Node * startNode = test.graph[0][startRotationIndex]->get(cell_index);
  delete p;

  p = loadPoly((std::string(graph_dir) + "/0-" + std::to_string(endRotationIndex) + ".tri").c_str());
  p->computeWindingNumbers();
  assert(!p->contains(end));
  cell_index = p->containingCell(end);
  FreeSpaceGraph::Node * endNode = test.graph[0][endRotationIndex]->get(cell_index);
  delete p;

  test.getPath(startNode, endNode, start, end, path);
  cout<<"done search graph"<<endl;

  cout<<endl<<endl;
  for (int i=0; i<path.size(); i++)
    cout<<path[i].first->getApprox().getX().mid()<<" "<<path[i].first->getApprox().getY().mid()<<" "<<path[i].first->getApprox().getZ().mid()<<" "<<path[i].second<<endl;
}

int main (int argc, char *argv[]) {
  Parameter::enable();
  time_t start_t,end_t;

  if (argc == 2) { 
    unsigned seed = atoi(argv[1]);
    srandom(seed);
  }

  const char * graph_dir = "skinny-lattice-graph";

  //--------------------------------
  // Test creating a freespace-graph
  const char * blockspace_dir = "skinny-lattice-output";
  double theta = M_PI / 20;
  double clearance_unit = 0.05;
  int num_levels = 1;
  
  time(&start_t);
  testCreateGraph(blockspace_dir, graph_dir, theta, clearance_unit, num_levels);
  time(&end_t);
  cout << "Elapsed time: "<< difftime (end_t,start_t)/60.0 << " minutes" << endl;
  //--------------------------------

  //--------------------------------
  // Test searching a freespace-graph
  PTR<Point> start = new Point(-2.0,0.0, -0.5);
  PTR<Point> end = new Point(5.0,0.0,0.0);
  int startRotationIndex = 0;
  int endRotationIndex = 0;
  
  time(&start_t);
  testSearchGraph(graph_dir, start, end, startRotationIndex, endRotationIndex);
  time(&end_t);
  cout << "Elapsed time: "<< difftime (end_t,start_t) << " seconds" << endl;
  //--------------------------------
}