#include "freespace-graph.h"

Polyhedron * loadPoly(const char * filename) {
  Polyhedron * poly;
  ifstream infile (filename);
  if (infile.is_open()) {
    poly = readPolyhedronVTK (infile);
    infile.close();
  } else {
    cout<<"could not read from file"<<endl;
    return NULL;
  }

  return poly;
}

class InputParameter : public Object<Parameter> {
public:
  InputParameter (double x) { set(Parameter::input(x)); }
};

void scale(Polyhedron * p, PTR<Object<Parameter> > unit) {
  for (int i=0; i> p->vertices.size(); i++) {
    Vertex * v = p->vertices[i];
    p->moveVertex(v, new ScalePoint(v->getP(), unit));
  }
}

PTR<Point> pointInCell(Polyhedron * poly, int i) {
  Cell * cell =  poly->cells[i];
  Face * face = cell->getShell(0)->getHFaces()[0]->getF();
  PTR<Point> fp;
  double unit = 1;
  do {
    fp = new FacePoint(cell, unit);
    unit /= 2;
  } while (!face->contains(fp));
  PTR<Point> p;
  unit = 1;
  do {
    p = new CellInternalPoint(cell, fp, unit);
    unit /= 2;
  } while (!cell->contains(p));
  return p;
}

FreeSpaceGraph::FreeSpaceGraph(std::vector<Polyhedron*> & blockspaces, double theta, double clearance_unit, int num_levels) {
  //read in the unit sphere and scale it by unit
  unit_ball = loadPoly("sphere.vtk");
  PTR<Object<Parameter> > unit =  new InputParameter(clearance_unit);
  scale(unit_ball, unit);

  this->blockspaces.insert(this->blockspaces.begin(), blockspaces.begin(), blockspaces.end());
  this->theta = theta;

  //initialize the graph
  for(int i=0; i<num_levels; i++) {
    std::vector<BlockSpaceNode*> v;
    graph.push_back(v);
  }
  for (int i=0; i<num_levels; i++) 
    for (int j=0; j< blockspaces.size(); j++)
      graph[i].push_back(new FreeSpaceGraph::BlockSpaceNode(blockspaces[j]));

  //create the first level nodes
  for (int i=0; i<blockspaces.size(); i++) {
    blockspaces[i]->computeWindingNumbers();
    for (int j=0; j< blockspaces[i]->cells.size(); j++) {
      if (blockspaces[i]->cells[j]->getWN() == 0) {
        graph[0][i]->nodes.push_back(new FreeSpaceGraph::Node(0, i, j, NULL));
      }
    }
  }

  //create the first level edges
  for (int i=0; i<blockspaces.size(); i++) {
    int j = (i+1)%blockspaces.size();
    Polyhedron * block_union = blockspaces[i]->boolean(blockspaces[j], Union);
    block_union->computeWindingNumbers();
    for (int k=0; k< block_union->cells.size(); k++)
      if (block_union->cells[k]->getWN() == 0) {
        cout<<i<<" "<<j<<" "<<k<<endl;
        //generate a point in the cell
        PTR<Point> p = pointInCell(block_union, k);
        //find with components of blockspaces i and j it lies in
        int ci = blockspaces[i]->containingCell(p);
        int cj = blockspaces[j]->containingCell(p);
        //make them neighbors
        graph[0][i]->get(ci)->neighbors.push_back(graph[0][j]->get(cj));
        graph[0][j]->get(cj)->neighbors.push_back(graph[0][i]->get(ci));
      }
    delete block_union;
  }

  //create subsequent layers...

}