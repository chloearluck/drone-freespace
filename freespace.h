#include "hull.h"
#include "geometry3d.h"
#include "mink.h"
#include <cstring>



class FreeSpace {
public:
  Polyhedron * robot, * obstacle, * bb;
  std::vector<Polyhedron*> cspaces;
  
  class Node {
  public:
    Polyhedron * poly;
    int cell_index;
    Node(Polyhedron * poly, int cell_index) : poly(poly), cell_index(cell_index) {}
  };

  class Edge {
  public:
    Node *a, *b;
    PTR<Point> p;
    Edge(Node * a, Node * b, PTR<Point> p) : a(a), b(b), p(p) {}
  };

  std::vector<std::vector<Node*> > nodes;
  std::vector<Edge*> edges;
  int numRotations;

  FreeSpace(Polyhedron * robot, Polyhedron * obstacle, double theta, double * bounding_box);
};