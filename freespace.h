#include "hull.h"
#include "geometry3d.h"
#include "mink.h"
#include <cstring>



class FreeSpace {
public:
  Polyhedron * robot, * obstacle, * bb;
  std::vector<Polyhedron*> cspaces;
  
  class Edge;
  class Node {
  public:
    Polyhedron * poly;
    int cell_index;
    std::vector<Edge *> edges;
    Node(Polyhedron * poly, int cell_index) : poly(poly), cell_index(cell_index) {}
  };

  class Edge {
  public:
    Node *a, *b;
    PTR<Point> p;
    Edge(Node * a, Node * b, PTR<Point> p) {
      this->a = a;
      this->b = b;
      a->edges.push_back(this);
      b->edges.push_back(this);
    }
  };

  std::vector<std::vector<Node*> > nodes;
  std::vector<Edge*> edges;
  int numRotations;

  FreeSpace(Polyhedron * robot, Polyhedron * obstacle, double theta, double * bounding_box);

private:
  Node * findOrAddNode(int cspace_index, int cell_index);
};

