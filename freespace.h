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
    int cspace_index;
    int cell_index;
    std::vector<Edge *> edges;
    bool discovered;
    Node * parent;
    Node(int cspace_index, int cell_index) : cspace_index(cspace_index), cell_index(cell_index), discovered(false), parent(NULL) {}
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
  Node * findNode(int cspace_index, int cell_index); 

private:
  Node * findOrAddNode(int cspace_index, int cell_index);
};

