#include "path3d.h"
#include "mink.h"
#include "io.h"
#include <sys/stat.h>
#include <errno.h>
#include <string>
#include <iostream>

Polyhedron * loadPoly(const char * filename);
void savePoly(Polyhedron * p, char * filename);

class FreeSpaceGraph {
 public: 
  double theta;

  class Node {
   public:
    int level;
    int blockspace_index;
    int cell_index;
    bool enabled;
    Node * parent;  
    std::vector<Node*> children; 
    std::vector<Node*> neighbors;
    std::vector<int> neighborIntersectionIndex; 
    std::vector<Node*> siblings;
    std::vector< std::pair< PTR<Point>, PTR<Point> > > siblingPoints;
    Node(int level, int blockspace_index, int cell_index, Node * parent)
     : level(level), blockspace_index(blockspace_index), cell_index(cell_index), parent(parent), enabled(true) {}
  };

  class BlockSpaceNode {
   public:
    int level;
    int blockspace_index;
    std::vector<Node*> nodes;
    BlockSpaceNode(int level, int blockspace_index) : level(level), blockspace_index(blockspace_index) {}
    Node * get(int index) {
      for (int i=0; i<nodes.size(); i++)
        if (nodes[i]->cell_index == index)
          return nodes[i];
      return NULL;
    }
    Node * getOrCreate(int index) {
      Node * n = get(index);
      if (n == NULL) {
        n = new Node(level, blockspace_index, index, NULL);
        nodes.push_back(n);
      }
      return n;
    }
  };

  int num_levels;
  int blockspaces_per_level;
  std::string dir;
  std::vector<std::vector<BlockSpaceNode *> > graph;

  void deepestPath(FreeSpaceGraph::Node * n, PTR<Point> p, PTR<Point> q, std::vector<std::pair<FreeSpaceGraph::Node * , PTR<Point> > > & path);
  void nodePointPath(std::vector<FreeSpaceGraph::Node*> & nodes, PTR<Point> a, PTR<Point> b, std::vector<std::pair<FreeSpaceGraph::Node * , PTR<Point> > > & path);
  void bfsPath(FreeSpaceGraph::Node * start, FreeSpaceGraph::Node * end, std::vector<FreeSpaceGraph::Node*> & nodes);
  bool isConnected(FreeSpaceGraph::Node * start, FreeSpaceGraph::Node * end);
  void getPath(Node * start, Node * end, PTR<Point> a, PTR<Point> b, std::vector<std::pair<PTR<Point>, double > > & path);
  double angle(int blockspace_index);

  FreeSpaceGraph(std::vector<Polyhedron*> & original_blockspaces, double theta, double clearance_unit, int num_levels, const char * dir);
  FreeSpaceGraph(const char * dir);
};

