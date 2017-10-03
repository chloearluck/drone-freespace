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

//transformation which rotates triangle abc into the xy-plane centered around b
class XYPlaneTriangleTransfromation : public Transformation {
  PTR<Point> pa, pb, pc, pd;
  TransformationData calculate () {
    PV3 a = pa->getP();
    PV3 b = pb->getP();
    PV3 c = pc->getP();

    PV3 u = (c-b).unit();
    PV3 w = (a-b).cross(u).unit();
    PV3 v = u.cross(w);
    PV3 t = b - PV3(b.dot(u), b.dot(v), b.dot(w));

    return TransformationData(u,v,w,t);
  }
 public:
  XYPlaneTriangleTransfromation(PTR<Point> pa, PTR<Point> pb, PTR<Point> pc) : pa(pa), pb(pb), pc(pc) {}
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

Edge * commonEdge(HFace * hf1, HFace * hf2) {
  Face * f1 = hf1->getF();
  Face * f2 = hf2->getF();
  HEdge * start_edge = f1->getBoundary(0);
  HEdge * curr_edge = start_edge;
  HEdge * common_edge = NULL;
  do {
    if (curr_edge->getE()->otherFace(f1) == f2) 
      return curr_edge->getE();
    curr_edge = curr_edge->getNext();
  } while (curr_edge != start_edge);
  return NULL;
}

class PathTriangle {
 public:
  PTR<Point> p[3];
  PathTriangle();
  PathTriangle(PTR<Point> p0, PTR<Point> p1, PTR<Point> p2) {
    p[0] = p0;  p[1] = p1;  p[2] = p2;
  }
  PathTriangle(const PathTriangle& t) {
    p[0] = t.p[0];  p[1] = t.p[1];  p[2] = t.p[2];
  }
};

void savePathTriangles(std::vector<PathTriangle> ts, const char * filename) {
  ofstream ostr;
  ostr.open(filename);
  if (ostr.is_open()) {
    ostr << setprecision(20) << "# vtk DataFile Version 3.0" << endl
         << "vtk output" << endl << "ASCII" << endl
         << "DATASET POLYDATA" << endl 
         << "POINTS " << ts.size()*3 << " double" << endl;
    for (int i=0; i<ts.size(); i++)
      for (int j=0; j<3; j++)
        ostr << ts[i].p[j]->getApprox().getX().mid() << " " << ts[i].p[j]->getApprox().getY().mid() << " " << ts[i].p[j]->getApprox().getZ().mid() << endl;
  
    ostr << endl << "POLYGONS " << ts.size() << " " << 4*ts.size() << endl;
    for (int i=0; i<ts.size(); i++)
      ostr << "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << endl;
    ostr.close();
  } else { cout<<"could not open file"<<endl; return; }
}

void flattenTriangles(std::vector<PathTriangle> & triangles, std::vector<PTR<Transformation> > & transformations, std::vector<PathTriangle> & flattened, PTR<Transformation> xyplane) {
  PTR<Transformation> cumulative = xyplane;
  flattened.push_back(PathTriangle(new TransformedPoint(triangles[0].p[0], cumulative), 
                                   new TransformedPoint(triangles[0].p[1], cumulative),
                                   new TransformedPoint(triangles[0].p[2], cumulative)));

  for (int i=1; i<triangles.size(); i++) {
    PathTriangle t(triangles[i]);
    t.p[0] = new TransformedPoint(t.p[0], cumulative);
    t.p[1] = new TransformedPoint(t.p[1], cumulative);
    
    cumulative = new CompositeTransformation(cumulative, transformations[i-1]);
    t.p[2] = new TransformedPoint(t.p[2], cumulative);
    flattened.push_back(t);
  }
}


void localPath(PTR<FaceIntersectionPoint> a, PTR<FaceIntersectionPoint> b, HFaces & pathfaces, Points & path) {
  /* ---- unoptimized BFS solution ---- */
  path.push_back((PTR<Point>) a);
  for (int i=1; i<pathfaces.size(); i++) {
    Edge * e = commonEdge(pathfaces[i], pathfaces[i-1]);
    assert(e != NULL);
    path.push_back(new MidPoint(e->getT()->getP(), e->getH()->getP()));
  }
  /* ---------------------------------- */

  std::vector<PathTriangle> triangles;
  std::vector<PTR<Transformation> > transformations;
  for (int i=1; i<pathfaces.size(); i++) {
    HFace * hf1 = pathfaces[i-1];    
    HFace * hf2 = pathfaces[i];
    Edge * ce = commonEdge(hf1, hf2);
    HEdge * he1 = (ce->getHEdge(0)->getF() == hf1->getF() ? ce->getHEdge(0) : ce->getHEdge(1) );
    HEdge * he2 = (he1 == ce->getHEdge(0) ? ce->getHEdge(1) : ce->getHEdge(0));
    PTR<Point> p1 = he1->getNext()->head()->getP();
    PTR<Point> p2 = he2->getNext()->head()->getP();
    PTR<Transformation> t = new UnfoldTriangleTransformation(p1, ce->getT()->getP(), ce->getH()->getP(), p2);
    transformations.push_back(t);
    if (i == 1) triangles.push_back(PathTriangle(p1, ce->getT()->getP(), ce->getH()->getP()));
    triangles.push_back(PathTriangle(ce->getT()->getP(), ce->getH()->getP(), p2));
  }
  std::vector<PathTriangle> flattened;
  PTR<Transformation> xyplane = new XYPlaneTriangleTransfromation(triangles[0].p[0], triangles[0].p[1], triangles[0].p[2]);
  flattenTriangles(triangles, transformations, flattened, xyplane);
  savePathTriangles(triangles, "triangles.vtk");
  savePathTriangles(flattened, "flattened.vtk");
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
}
