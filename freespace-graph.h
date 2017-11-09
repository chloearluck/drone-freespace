#include "path3d.h"
#include "io.h"

Polyhedron * loadPoly(const char * filename);

class FreeSpaceGraph {
 public: 
  std::vector<Polyhedron*> blockspaces;
  Polyhedron * unit_ball;
  double theta;

  class Node {
   public:
    int level;
    int blockspace_index;
    int cell_index;
    Node * parent; //component of the previous level which contains this component
    std::vector<Node*> children; //components of the next level which are contained in this component
    std::vector<Node*> neighbors; //components of neighboring angles which intersect this component
    Node(int level, int blockspace_index, int cell_index, Node * parent) 
     : level(level), blockspace_index(blockspace_index), cell_index(cell_index), parent(parent) {}
  };

  class BlockSpaceNode {
   public:
    Polyhedron * blockspace;
    std::vector<Node*> nodes;
    BlockSpaceNode(Polyhedron * blockspace) : blockspace(blockspace) {}
    Node * get(int index) {
      for (int i=0; i<nodes.size(); i++)
        if (nodes[i]->cell_index == index)
          return nodes[i];
      return NULL;
    }
  };

  std::vector<std::vector<BlockSpaceNode *> > graph;
  // graph[level][blockspace_index]
  // graph[level][blockspace_index]->get(cell_index)

  FreeSpaceGraph(std::vector<Polyhedron*> & blockspaces, double theta, double clearance_unit, int num_levels);
};

