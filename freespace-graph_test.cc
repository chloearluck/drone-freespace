#include "freespace-graph.h"
#include "freespace.h"

void testPathExists() {
  Polyhedron * robot = loadPolyVTK("spiny.vtk");
  Polyhedron * obstacle = loadPolyVTK("wirecage.vtk");
  int numRotations = 40;
  
  Parameter::disable();
  PTR<Object<Parameter> > tan_half_angle = new InputParameter(tan((1.0 + 1.0e-8) * M_PI / numRotations));
  Parameter::enable();

  FreeSpace * fs  = new FreeSpace(robot, obstacle, tan_half_angle, numRotations);
  FreeSpaceGraph graph(fs->blockspaces_close, fs->blockspaces_rough, M_PI * (2/numRotations), 0.01, 1, "wireGraph");

  //search for a path from a to b
  PTR<Point> a = new Point(-4.5, 0.0, 0.0);
  PTR<Point> b = new Point( 4.5, 0.0, 0.0);

  Polyhedron * p = fs->blockspaces_close[0];
  p->computeWindingNumbers();
  assert(!p->contains(a) && !p->contains(b));
  int a_index = p->containingCell(a);
  int b_index = p->containingCell(b);
  cout << "finding a path from cell "<<a_index<<" to cell "<<b_index<<" in freespace rotation 0"<<endl;
  FreeSpaceGraph::Node * startNode = graph.graph[0][0]->get(a_index);
  FreeSpaceGraph::Node * endNode = graph.graph[0][0]->get(b_index);
  
  //call bfspath and print out the nodes
  std::vector<FreeSpaceGraph::Node*> nodes;
  graph.bfsPath(startNode, endNode, nodes);
  for (int i=0; i<nodes.size(); i++) 
    cout << nodes[i]->level << " " << nodes[i]->blockspace_index << " " << nodes[i]->cell_index << endl;
}

void testCreateGraph(const char * blockspace_dir, const char * graph_dir, double theta, double clearance_unit, int num_levels) {
  int numRotations = floor(2*M_PI/theta);
  
  std::vector<Polyhedron *> close_blockspaces;
  char s[50];
  for (int i=0; i<=40; i++) {
    sprintf(s, "%s/close%02d-out.vtk", blockspace_dir, i);
    Polyhedron * p = loadPolyVTK(s);
    close_blockspaces.push_back(p);
  }

  std::vector<Polyhedron *> rough_blockspaces;
  for (int i=0; i<=40; i++) {
    sprintf(s, "%s/rough%02d-out.vtk", blockspace_dir, i);
    Polyhedron * p = loadPolyVTK(s);
    rough_blockspaces.push_back(p);
  }

  FreeSpaceGraph graph(close_blockspaces, rough_blockspaces, theta, clearance_unit, num_levels, graph_dir);
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

  testPathExists();
  return 0;

  const char * graph_dir = "droneRoomGraph";

  //--------------------------------
  // Test creating a freespace-graph
  const char * blockspace_dir = "droneRoomOutput";
  double theta = M_PI / 20;
  double clearance_unit = 0.05;
  int num_levels = 3;
  
  time(&start_t);
  testCreateGraph(blockspace_dir, graph_dir, theta, clearance_unit, num_levels);
  time(&end_t);
  cout << "Elapsed time: "<< difftime (end_t,start_t)/60.0 << " minutes" << endl;
  //--------------------------------

  // //--------------------------------
  // // Test searching a freespace-graph
  // PTR<Point> start = new Point(-1.0,-1.0,-1.0);
  // PTR<Point> end = new Point(5.0,5.0,5.0);
  // int startRotationIndex = 0;
  // int endRotationIndex = 0;
  
  // time(&start_t);
  // testSearchGraph(graph_dir, start, end, startRotationIndex, endRotationIndex);
  // time(&end_t);
  // cout << "Elapsed time: "<< difftime (end_t,start_t) << " seconds" << endl;
  // //--------------------------------
}