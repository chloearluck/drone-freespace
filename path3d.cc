#include "path3d.h"

struct ComparePointOrder {
  PTR<Point> r;
  ComparePointOrder(PTR<Point> r) : r(r) {}
  bool operator()(PTR<FaceIntersectionPoint> i, PTR<FaceIntersectionPoint> j) {
    return (PointOrder(i, j, r) > 0); 
  }
};

HEdge * commonEdge(HFace * hf1, HFace * hf2) {
  Face * f1 = hf1->getF();
  Face * f2 = hf2->getF();
  HEdge * start_edge = f1->getBoundary(0);
  HEdge * curr_edge = start_edge;
  HEdge * common_edge = NULL;
  do {
    if (curr_edge->getE()->otherFace(f1) == f2) 
      return curr_edge;
    curr_edge = curr_edge->getNext();
  } while (curr_edge != start_edge);
  return NULL;
}

void bfs(PTR<FaceIntersectionPoint> a, PTR<FaceIntersectionPoint> b, Points & path) {
  cout<<"finding a sequence of faces connecting"<<endl;
  pp(a); cout<<"and"<<endl; pp(b);

  HFace * fa  = a->getHFace();
  HFace * fb  = b->getHFace();
  HFaces pathfaces;

  std::map<HFace *, HFace *> parents;
  std::queue<HFace *> q;
  q.push(fa);
  parents[fa] = NULL;

  while (!q.empty()) {
    HFace * current = q.front(); q.pop();
    if (current == fb) {
      cout<<"found path"<<endl;
      while(parents[current] != NULL) {
        pathfaces.push_back(current);
        current = parents[current];
      }
      break;
    }

    HFaces hfaces;
    current->neighbors(hfaces);
    for (int i=0; i<hfaces.size(); i++) {
      HFace * hf = hfaces[i];
      if (parents.find(hf) == parents.end()) {
        parents[hf] = current;
        q.push(hf);
      }
    }
  }
  if (pathfaces.size() == 0) {
    cout<<"no path found"<<endl;
    return;
  }

  path.push_back((PTR<Point>) a);
  for (int i=pathfaces.size()-1; i>0; i--) {
    HEdge * e = commonEdge(pathfaces[i], pathfaces[i-1]);
    assert(e != NULL);
    path.push_back(new MidPoint(e->tail()->getP(), e->head()->getP()));
  }
}

void findPath(Polyhedron * blockspace, PTR<Point> start, PTR<Point> end, Points &path) {
  //make sure start and end are not in blackspace
  if (blockspace->contains(start)) {
    cout<<"invalid start"<<endl;
    return;
  }
  if (blockspace->contains(end)) {
    cout<<"invalid end"<<endl;
    return;
  }
  
  //make sure start and end are in the same component of blockspace
  blockspace->computeWindingNumbers();
  cout<<"blockspace has "<<blockspace->cells.size()<<" cells"<<endl;
  int cell_index = blockspace->containingCell(start);
  if (blockspace->containingCell(end) != cell_index) {
    cout<<"start and end are not in the same component"<<endl;
    return;
  }
  cout<<"start and end are in cell "<<cell_index<<endl;   

  //do raycasting
  Cell * cell = blockspace->cells[cell_index];
  std::vector<PTR<FaceIntersectionPoint> > points;
  PTR<Point> r =  new DiffPoint(end, start);
  for (int i=0; i<cell->nShells(); i++) {
    Shell * shell = cell->getShell(i);
    for (int j=0; j<shell->getHFaces().size(); j++) {
      HFace * hface = shell->getHFaces()[j];
      PTR<FaceIntersectionPoint> q = new FaceIntersectionPoint(start, end, hface);
      bool onFace = hface->getF()->contains(q);
      bool afterStart = (PointOrder(start, q, r) == 1);
      bool beforeEnd = (PointOrder(q, end, r) == 1);
      if (onFace && afterStart && beforeEnd)
        points.push_back(q);
    }
  }
  
  cout<<"intersections: "<<points.size()<<endl;
  assert(points.size() % 2 == 0);

  std::sort(points.begin(), points.end(), ComparePointOrder(r));
  pp(start);
  for (int i=0; i<points.size(); i++)
    pp(points[i]);
  pp(end);

  path.push_back(start);
  for (int i=0; i<points.size(); i+=2) {
    bfs(points[i], points[i+1], path);
  }
  if (points.size()>0)
    path.push_back((PTR<Point>) points[points.size()-1]);
  path.push_back(end);

  for (int i=0; i<path.size(); i++)
    pp(path[i]);
}
