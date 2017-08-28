#include "hull.h"
#include "geometry3d.h"
#include "mink.h"
#include <cstring>


class FreeSpace {
public:
  Polyhedron * robot, * obstacle;
  std::vector<Polyhedron*> blockspaces;
  
  class Edge;
  class Node {
  public:
    int space_index;
    int cell_index;
    std::vector<Edge *> edges;
    bool discovered;
    Node * parent;
    Node(int space_index, int cell_index) : space_index(space_index), cell_index(cell_index), discovered(false), parent(NULL) {}
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

  FreeSpace(Polyhedron * robot, Polyhedron * obstacle, double theta);
  Node * findNode(int space_index, int cell_index); 

private:
  Node * findOrAddNode(int space_index, int cell_index);
};

