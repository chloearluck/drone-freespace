#include "path3d.h"
#include "mink.h"
#include "io.h"

Polyhedron * loadPoly(const char * filename);

class FreeSpaceGraph {
 public: 
  Polyhedron * unit_ball;
  double theta;

  class Node {
   public:
    int level;
    int blockspace_index;
    int cell_index;
    Node * parent;  
    std::vector<Node*> children; 
    std::vector<Node*> neighbors; 
    std::vector<Node*> siblings;
    std::vector< pair< PTR<Point>, PTR<Point> > > siblingPoints;
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

  FreeSpaceGraph(std::vector<Polyhedron*> & original_blockspaces, double theta, double clearance_unit, int num_levels);
};

