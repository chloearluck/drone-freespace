#include "path3d.h"

class TransformationData {
 public:
  PV3 u, v, w, t;
  TransformationData () {}
  TransformationData (PV3 u, PV3 v, PV3 w, PV3 t) : u(u), v(v), w(w), t(t) {}
  Parameter & operator[](int i) { 
    if (i < 3) return u[i];
    if (i < 6) return v[i%3];
    if (i < 9) return w[i%3];
    return t[i%3];
  }
  int size () const { return 12; }
};

class Transformation : public Object<TransformationData> {
 public:
  PV3 getU () { return get().u; }
  PV3 getV () { return get().v; }
  PV3 getW () { return get().w; }
  PV3 getT () { return get().t; }
};

//transformation which rotates d about bc to lie in the same plane as abc
class UnfoldTriangleTransformation : public Transformation {
  PTR<Point> pa, pb, pc, pd;
  TransformationData calculate () {
    PV3 a = pa->getP();
    PV3 b = pb->getP();
    PV3 c = pc->getP();
    PV3 d = pd->getP();

    PV3 u = (c-b).unit();
    PV3 uprime = u;
    PV3 w = (d-b).cross(u).unit();
    PV3 wprime = u.cross(a-b).unit();
    PV3 v = w.cross(u);
    PV3 vprime = wprime.cross(u);

    PV3 t = b - b.dot(u)*(uprime) - b.dot(v)*(vprime) - b.dot(w)*(wprime);
    PV3 x = uprime.getX()*u + vprime.getX()*v + wprime.getX()*w;
    PV3 y = uprime.getY()*u + vprime.getY()*v + wprime.getY()*w;
    PV3 z = uprime.getZ()*u + vprime.getZ()*v + wprime.getZ()*w;
    return TransformationData(x,y,z,t);
  }
 public:
  UnfoldTriangleTransformation(PTR<Point> pa, PTR<Point> pb, PTR<Point> pc, PTR<Point> pd) : pa(pa), pb(pb), pc(pc), pd(pd) {}
};

class CompositeTransformation : public Transformation {
  PTR<Transformation> f, g;
  TransformationData calculate () {
    PV3 u1 = g->getU();
    PV3 v1 = g->getV();
    PV3 w1 = g->getW();
    PV3 t1 = g->getT();
    PV3 u2 = f->getU();
    PV3 v2 = f->getV();
    PV3 w2 = f->getW();
    PV3 t2 = f->getT();

    PV3 x1 = PV3(u1.getX(), v1.getX(), w1.getX());
    PV3 y1 = PV3(u1.getY(), v1.getY(), w1.getY());
    PV3 z1 = PV3(u1.getZ(), v1.getZ(), w1.getZ());

    PV3 u = PV3(x1.dot(u2), y1.dot(u2), z1.dot(u2)); 
    PV3 v = PV3(x1.dot(v2), y1.dot(v2), z1.dot(v2)); 
    PV3 w = PV3(x1.dot(w2), y1.dot(w2), z1.dot(w2)); 
    PV3 t = PV3(u2.dot(t1), v2.dot(t1), w2.dot(t1)) + t2;

    return TransformationData(u,v,w,t);
  }
 public:
  CompositeTransformation(PTR<Transformation> f, PTR<Transformation> g) : f(f), g(g) {}
};

//apply transformation t to point p
class TransformedPoint : public Point {
  PTR<Point> point;
  PTR<Transformation> t;
  PV3 calculate () {
    PV3 p = point->getP();
    PV3 u = t->getU();
    PV3 v = t->getV();
    PV3 w = t->getW();
    PV3 tt = t->getT();
    return PV3(u.dot(p), v.dot(p), w.dot(p)) + tt;
  }
 public:
  TransformedPoint(PTR<Point> point, PTR<Transformation> t) : point(point), t(t) {}
};

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

void localPath(PTR<FaceIntersectionPoint> a, PTR<FaceIntersectionPoint> b, HFaces & pathfaces, Points & path) {
  path.push_back((PTR<Point>) a);
  for (int i=1; i<pathfaces.size(); i++) {
    HEdge * e = commonEdge(pathfaces[i], pathfaces[i-1]);
    assert(e != NULL);
    path.push_back(new MidPoint(e->tail()->getP(), e->head()->getP()));
  }
}

void bfs(PTR<FaceIntersectionPoint> a, PTR<FaceIntersectionPoint> b, HFaces & pathfaces) {
  cout<<"finding a sequence of faces connecting"<<endl;
  pp(a); cout<<"and"<<endl; pp(b);

  HFace * fa  = a->getHFace();
  HFace * fb  = b->getHFace();
  HFaces pathfaces_rev;

  std::map<HFace *, HFace *> parents;
  std::queue<HFace *> q;
  q.push(fa);
  parents[fa] = NULL;

  while (!q.empty()) {
    HFace * current = q.front(); q.pop();
    if (current == fb) {
      cout<<"found path"<<endl;
      while(current != NULL) {
        pathfaces_rev.push_back(current);
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
  if (pathfaces_rev.size() == 0) {
    cout<<"no path found"<<endl;
    return;
  }

  for (int i=pathfaces_rev.size()-1; i>=0; i--)
    pathfaces.push_back(pathfaces_rev[i]);
}

void findPath(Polyhedron * blockspace, PTR<Point> start, PTR<Point> end, Points &path) {
  //make sure start and end are not in blockspace
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
    HFaces subHfaces;
    Points subPath;
    bfs(points[i], points[i+1], subHfaces);
    localPath(points[i], points[i+1], subHfaces, subPath);
    path.insert(path.end(), subPath.begin(), subPath.end());
  }
  if (points.size()>0)
    path.push_back((PTR<Point>) points[points.size()-1]);
  path.push_back(end);

  for (int i=0; i<path.size(); i++)
    pp(path[i]);
}
