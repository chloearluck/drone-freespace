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

pair <PTR<Point>, PTR<Point> > nearestPointPair(Polyhedron * poly, int i, int j) {
  PTR<Point> a = NULL;
  PTR<Point> b = NULL;

  //populate HFace lists
  HFaces hfs1, hfs2;
  for (int k=0; k<poly->cells[i]->nShells(); k++) {
    Shell * s = poly->cells[i]->getShell(k);
    HFaces hfs = s->getHFaces();
    hfs1.insert(hfs1.end(), hfs.begin(), hfs.end());
  }
  for (int k=0; k<poly->cells[j]->nShells(); k++) {
    Shell * s = poly->cells[j]->getShell(k);
    HFaces hfs = s->getHFaces();
    hfs2.insert(hfs2.end(), hfs.begin(), hfs.end());
  }

  //populate vertex lists
  Vertices v1, v2;
  for (int k=0; k<hfs1.size(); k++) {
    HEdge * heFirst = hfs1[k]->getF()->getBoundary(0);
    HEdge * he = heFirst;
    v1.push_back(he->tail());
    do {
      if (std::find(v1.begin(), v1.end(), he->tail()) == v1.end())
        v1.push_back(he->tail());
      he = he->getNext();
    } while (he != heFirst);
  }
  for (int k=0; k<hfs2.size(); k++) {
    HEdge * heFirst = hfs2[k]->getF()->getBoundary(0);
    HEdge * he = heFirst;
    v2.push_back(he->tail());
    do {
      if (std::find(v2.begin(), v2.end(), he->tail()) == v2.end())
        v2.push_back(he->tail());
      he = he->getNext();
    } while (he != heFirst);
  }

  //find closest pair
  for (int k=0; k<v1.size(); k++)
    for (int l=0; l<v2.size(); l++)
      if (a == NULL || CloserPair(v1[k]->getP(), v2[l]->getP(), a, b) > 0) {
        a = v1[k]->getP();
        b = v2[l]->getP();
      }

  return std::make_pair(a,b);
}

FreeSpaceGraph::FreeSpaceGraph(std::vector<Polyhedron*> & original_blockspaces, double theta, double clearance_unit, int num_levels, const char * dir) {
  //read in the unit sphere and scale it by unit
  Polyhedron * unit_ball = loadPoly("sphere.vtk");
  PTR<Object<Parameter> > unit =  new InputParameter(clearance_unit);
  scale(unit_ball, unit);

  this->theta = theta;
  this->num_levels = num_levels;
  this->blockspaces_per_level = original_blockspaces.size();
  std::vector<Polyhedron*> blockspaces;
  blockspaces.insert(blockspaces.begin(), original_blockspaces.begin(), original_blockspaces.end());

  //create directory
  if (mkdir(dir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) == -1 && errno != EEXIST)
    cout<<"could not create directry"<<endl;

  //initialize the graph
  for(int i=0; i<num_levels; i++) 
    graph.push_back(std::vector<BlockSpaceNode*>());
  for (int i=0; i<num_levels; i++) //TO DO: save blockspaces and block unions to file
    for (int j=0; j< blockspaces.size(); j++)
      graph[i].push_back(new FreeSpaceGraph::BlockSpaceNode(i, j));

  std::vector<Polyhedron*> prev_blockspaces;
  for (int level = 0; level<num_levels; level++) {
    cout<<"level "<<level<<endl;

    if (level > 0) {
      for (int i=0; i<prev_blockspaces.size(); i++)
        delete prev_blockspaces[i];
      prev_blockspaces.clear();

      prev_blockspaces.insert(prev_blockspaces.begin(), blockspaces.begin(), blockspaces.end());
      blockspaces.clear();

      for (int i=0; i<prev_blockspaces.size(); i++) 
        blockspaces.push_back(minkowskiSumFull(prev_blockspaces[i], unit_ball));
    }

    cout<<"creating nodes"<<endl;
    //create nodes
    for (int i=0; i<blockspaces.size(); i++) {
      blockspaces[i]->computeWindingNumbers();
      for (int j=0; j<blockspaces[i]->cells.size(); j++) 
        if (blockspaces[i]->cells[j]->getWN() == 0) {
          if (level > 0) {
            //find the parent of this cell
            PTR<Point> p = pointInCell(blockspaces[i], j);
            int comp = prev_blockspaces[i]->containingCell(p);
            FreeSpaceGraph::Node * parent = graph[level-1][i]->get(comp);
            FreeSpaceGraph::Node * newNode = new FreeSpaceGraph::Node(level, i, j, parent);
            graph[level][i]->nodes.push_back(newNode);
            parent->children.push_back(newNode);
          } else 
            graph[level][i]->nodes.push_back(new FreeSpaceGraph::Node(level, i, j, NULL));
        }
    }

    cout<<"creating edges"<<endl;
    //create neighbor edges
    for (int i=0; i<blockspaces.size(); i++) {
      int j = (i+1)%blockspaces.size();
      Polyhedron * block_union = blockspaces[i]->boolean(blockspaces[j], Union);
      block_union->computeWindingNumbers();
      for (int k=0; k< block_union->cells.size(); k++)
        if (block_union->cells[k]->getWN() == 0) {
          cout<<i<<" "<<j<<" "<<k<<endl;
          PTR<Point> p = pointInCell(block_union, k);
          int ci = blockspaces[i]->containingCell(p);
          int cj = blockspaces[j]->containingCell(p);
          graph[level][i]->get(ci)->neighbors.push_back(graph[level][j]->get(cj));
          graph[level][j]->get(cj)->neighbors.push_back(graph[level][i]->get(ci));
        }
      delete block_union;  
    }

    cout<<"finding siblings"<<endl;
    if (level > 0) {
      for (int i=0; i<blockspaces.size(); i++) 
        for (int j=0; j<graph[level-1][i]->nodes.size(); j++) {
          FreeSpaceGraph::Node * parent = graph[level-1][i]->nodes[j];
          if (parent->children.size() > 1)
            for (int c1 = 0; c1 < parent->children.size()-1; c1++)
              for (int c2 = c1+1; c2 < parent->children.size(); c2++) {
                FreeSpaceGraph::Node * child1 = parent->children[c1];
                FreeSpaceGraph::Node * child2 = parent->children[c2];
                child1->siblings.push_back(child2);
                child2->siblings.push_back(child1);
                std::pair< PTR<Point>, PTR<Point> > ab = nearestPointPair(blockspaces[i], child1->cell_index, child2->cell_index);
                child1->siblingPoints.push_back(ab);
                child2->siblingPoints.push_back(std::make_pair(ab.second,ab.first));
              }
        }
    }
  }

  delete unit_ball;

  cout<<"saving"<<endl;

  std::string s(dir);
  s.append("/graph.txt");
  ofstream out;
  out.open(s.c_str());
  out << num_levels << " " << blockspaces_per_level << " " << theta << " " << clearance_unit <<endl;
  if (out.is_open()) {
    for (int i=0; i<num_levels; i++)
      for (int j=0; j<blockspaces_per_level; j++)
        for (int k=0; k<graph[i][j]->nodes.size(); k++) {
          FreeSpaceGraph::Node * node = graph[i][j]->nodes[k];
          out << node->level << " " << node->blockspace_index << " " << node->cell_index << ",";
          if (node->parent != NULL) out << " " << node->parent->level << " " << node->parent->blockspace_index << " " << node->parent->cell_index << ",";
          else out << " NULL,";
          for (int l=0; l<node->children.size(); l++)
            out<< " " << node->children[l]->level << " " << node->children[l]->blockspace_index << " " << node->children[l]->cell_index;
          out<<",";
          for (int l=0; l<node->neighbors.size(); l++)
            out<< " " << node->neighbors[l]->level << " " << node->neighbors[l]->blockspace_index << " " << node->neighbors[l]->cell_index;
          if (node->siblings.size()>0) {
          out<<",";
            for (int l=0; l<node->siblings.size(); l++)
              out<< " " << node->siblings[l]->level << " " << node->siblings[l]->blockspace_index << " " << node->siblings[l]->cell_index;
            out<<",";
            for (int l=0; l<node->siblingPoints.size(); l++)
              out << " " << node->siblingPoints[l].first->getApprox().getX().mid()
                  << " " << node->siblingPoints[l].first->getApprox().getY().mid()
                  << " " << node->siblingPoints[l].first->getApprox().getZ().mid()
                  << " " << node->siblingPoints[l].second->getApprox().getX().mid()
                  << " " << node->siblingPoints[l].second->getApprox().getY().mid()
                  << " " << node->siblingPoints[l].second->getApprox().getZ().mid();
          }
          out<<endl;
        }
    out.close();
  }
}
FreeSpaceGraph::FreeSpaceGraph(const char * dir) {
  std::string s(dir);
  s.append("/graph.txt");
  ifstream infile (s.c_str());
  double clearance_unit;
  if (infile.is_open()) {
    infile >> num_levels >> blockspaces_per_level >> theta >> clearance_unit;

    //initialize graph
    for(int i=0; i<num_levels; i++) 
      graph.push_back(std::vector<BlockSpaceNode*>());
    for (int i=0; i<num_levels; i++)
      for (int j=0; j< blockspaces_per_level; j++)
        graph[i].push_back(new FreeSpaceGraph::BlockSpaceNode(i, j));

    std::string line;
    std::getline(infile, line); //eat new line
    while (std::getline(infile, line)) {
      istringstream ss(line);
      string id, parent, children, neighbors, siblings, siblingPoints;
      FreeSpaceGraph::Node * n;
      if (getline(ss, id, ',')) {
        istringstream ss2(id);
        int i,j,k;
        if (ss2 >> i >> j >> k)
          n = graph[i][j]->getOrCreate(k);
      } else { cout << "cannot read line ("<<line<<")"<<endl; return; }
      if (getline(ss, parent, ',') && parent.find("NULL") == string::npos) {
        istringstream ss2(parent);
        int i,j,k;
        if (ss2 >> i >> j >> k)
          n->parent = graph[i][j]->getOrCreate(k);
      }
      if (getline(ss, children, ',')) {
        istringstream ss2(children);
        int i,j,k;
        while (ss2 >> i >> j >> k)
          n->children.push_back(graph[i][j]->getOrCreate(k));
      }
      if (getline(ss, neighbors, ',')) {
        istringstream ss2(neighbors);
        int i,j,k;
        while (ss2 >> i >> j >> k)
          n->neighbors.push_back(graph[i][j]->getOrCreate(k));
      }
      if (getline(ss, siblings, ',')) {
        istringstream ss2(siblings);
        int i,j,k;
        while (ss2 >> i >> j >> k)
          n->siblings.push_back(graph[i][j]->getOrCreate(k));
      }
      if (getline(ss, siblingPoints, ',')) {
        istringstream ss2(siblingPoints);
        double x1, y1, z1, x2, y2, z2;
        while (ss2 >> x1 >> y1 >> z1 >> x2 >> y2 >> z2) {
          PTR<Point> a = new InputPoint(x1, y1, z1);
          PTR<Point> b = new InputPoint(x2, y2, z2);
          n->siblingPoints.push_back(std::make_pair(a,b));
        }
      }
    }
    infile.close();
  } else {
    cout<<"could not read from file"<<endl;
  }
}