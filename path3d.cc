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

class XYComponents : public Object<PV2> {
 protected:
  PTR<Point> p;
  PV2 calculate () { return PV2(p->getP().getX(), p->getP().getY()); }
 public:
  XYComponents (PTR<Point> p) : p(p) {} 
};

class ABintersectCD : public Object<PV2> {
 protected:
  PTR<Object<PV2> > pa,pb,pc,pd;
 public:
  ABintersectCD(PTR<Object<PV2> > pa, PTR<Object<PV2> > pb, PTR<Object<PV2> > pc, PTR<Object<PV2> >pd)
   : pa(pa), pb(pb), pc(pc), pd(pd) {}
  PV2 calculate () {
    PV2 a = pa->get();
    PV2 b = pb->get();
    PV2 c = pc->get();
    PV2 d = pd->get();

    Parameter t = -((a-c).cross(d-c)) / ((b-a).cross(d-c));
    return a + t *(b-a);
  }
};

Primitive3(AreaABC, PTR<Object<PV2> >, pa, PTR<Object<PV2> >, pb, PTR<Object<PV2> >, pc);
int AreaABC::sign() {
  PV2 a = pa->get();
  PV2 b = pb->get();
  PV2 c = pc->get();
  return (b-a).cross(c-a).sign();
}

Primitive2(AEqualB, PTR<Object<PV2> >, pa, PTR<Object<PV2> >, pb);
int AEqualB::sign() {
  PV2 a = pa->get();
  PV2 b = pb->get();
  if ((b-a).length() == 0)
    return 1;
  return 0;
}

class PathVertex {
 public:
  PTR<Point> original;
  PTR<Object<PV2> > transformed2d;
  PathVertex(PTR<Point> p) {
    original = p; transformed2d = 0;
  }
  PathVertex(PTR<Point> original, PTR<Object<PV2> > transformed2d) : original(original), transformed2d(transformed2d) {}
};

//find the 3d equivalent of the 2d intersection between triangle edge p->q and line segment c->d
class ABintersectCDto3D : public Point {
 protected:
  PathVertex *vp, *vq, *vc, *vd;
 public:
  ABintersectCDto3D(PathVertex * vp, PathVertex * vq, PathVertex * vc, PathVertex * vd)
   : vp(vp), vq(vq), vc(vc), vd(vd) {}
  PV3 calculate () {
    PV2 p = vp->transformed2d->get();
    PV2 q = vq->transformed2d->get();
    PV2 c = vc->transformed2d->get();
    PV2 d = vd->transformed2d->get();
    PV3 p3d = vp->original->getP();
    PV3 q3d = vq->original->getP();

    Parameter t = -((p-c).cross(d-c)) / ((q-p).cross(d-c));
    PV3 result = p3d + t *(q3d-p3d);
    return p3d + t *(q3d-p3d);
  }
};

class PathTriangle {
 public:
  PathVertex * p[3];
  HFace * hface;
  PathTriangle(PathVertex * p0, PathVertex * p1, PathVertex * p2, HFace * hface) {
    p[0] = p0; p[1] = p1; p[2] = p2; this->hface = hface;
  }
};

void savePathTriangles(std::vector<PathTriangle> ts, const char * filename, bool saveTranformed) {
  ofstream ostr;
  ostr.open(filename);
  if (ostr.is_open()) {
    ostr << setprecision(20) << "# vtk DataFile Version 3.0" << endl
         << "vtk output" << endl << "ASCII" << endl
         << "DATASET POLYDATA" << endl 
         << "POINTS " << ts.size()*3 << " double" << endl;
    for (int i=0; i<ts.size(); i++)
      for (int j=0; j<3; j++)
        if (saveTranformed)
          ostr << ts[i].p[j]->transformed2d->getApprox().getX().mid() << " " << ts[i].p[j]->transformed2d->getApprox().getY().mid() << " 0.0" << endl;
        else
          ostr << ts[i].p[j]->original->getApprox().getX().mid() << " " << ts[i].p[j]->original->getApprox().getY().mid() << " " << ts[i].p[j]->original->getApprox().getZ().mid() << endl;
  
    ostr << endl << "POLYGONS " << ts.size() << " " << 4*ts.size() << endl;
    for (int i=0; i<ts.size(); i++)
      ostr << "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << endl;
    ostr.close();
  } else { cout<<"could not open file"<<endl; return; }
}

PathVertex * getVertex(PTR<Point> p, std::vector<PathVertex * > & vertices) {
  for (int i=0; i<vertices.size(); i++)
    if (p == vertices[i]->original)
      return vertices[i];
  PathVertex * v = new PathVertex(p);
  vertices.push_back(v);
  return v;
}

void flattenTriangles(std::vector<PathTriangle> & triangles, std::vector<PTR<Transformation> > & transformations, PTR<Transformation> xyplane, PathVertex * start, PathVertex * end) {
  PTR<Transformation> cumulative = xyplane;
  triangles[0].p[0]->transformed2d = new XYComponents(new TransformedPoint(triangles[0].p[0]->original, cumulative));
  triangles[0].p[1]->transformed2d = new XYComponents(new TransformedPoint(triangles[0].p[1]->original, cumulative));
  triangles[0].p[2]->transformed2d = new XYComponents(new TransformedPoint(triangles[0].p[2]->original, cumulative));
  start->transformed2d = new XYComponents(new TransformedPoint(start->original, cumulative));

  for (int i=1; i<triangles.size(); i++) {
    if (triangles[i].p[0]->transformed2d == 0)
      triangles[i].p[0]->transformed2d = new XYComponents(new TransformedPoint(triangles[i].p[0]->original, cumulative));
    if (triangles[i].p[1]->transformed2d == 0)
      triangles[i].p[1]->transformed2d = new XYComponents(new TransformedPoint(triangles[i].p[1]->original, cumulative));
    
    cumulative = new CompositeTransformation(cumulative, transformations[i-1]);
    triangles[i].p[2]->transformed2d = new XYComponents(new TransformedPoint(triangles[i].p[2]->original, cumulative));
  }
  end->transformed2d = new XYComponents(new TransformedPoint(end->original, cumulative));
}

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

void commonEdge(PathTriangle t1, PathTriangle t2, PathVertex *& p, PathVertex *& q) {
  p = 0; q = 0;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      if (t1.p[i] == t2.p[j])
        if (!p)
          p = t1.p[i];
        else {
          q = t1.p[i];
          return;
        }
}

class PathEdge {
 public:
  PathVertex * left;
  PathVertex * right;
  PathEdge(PathVertex * left, PathVertex * right) : left(left), right(right) {}
};

void shortestPath(std::vector<PathTriangle> triangles, PathVertex * start, PathVertex * end, Points & path) {
  std::vector<PathVertex *> left;
  std::vector<int> leftIndices;
  std::vector<PathEdge> edges;
  for (int i=0; i<triangles.size()-1; i++) {
    PathVertex * p, * q;
    commonEdge(triangles[i], triangles[i+1], p, q);
    PathVertex * from;
    for (int j=0; j<3; j++)
      if (triangles[i].p[j] != p && triangles[i].p[j] != q)
        from = triangles[i].p[j];
    PathVertex * l = ((AreaABC(from->transformed2d, p->transformed2d, q->transformed2d) > 0)? p : q);
    PathVertex * r = ((l != p)? p : q);
    edges.push_back(PathEdge(l, r)); 
  }
  edges.push_back(PathEdge(end, end));

  left.push_back(start);
  leftIndices.push_back(-1);
  for (int i=0; i<edges.size(); i++) {
    
    if (left[left.size()-1] == edges[i].left)
      continue;

    //find the shortest path through edges[0]...edges[i-1] to edges[i].left;
    if (left.size() > 1) {
      int bound = left.size()-1;
      while (bound > 0 && AreaABC(left[bound-1]->transformed2d, left[bound]->transformed2d, edges[i].left->transformed2d) > 0)
        bound--;

      if (bound < left.size()-1) {
        //calculate intersections for edges leftIndices[bound+1]
        //  (or right point is intersection is invalid)
        //add them all to newPoints;
        std::vector<PathVertex *> newPoints;
        std::vector<int> newPointsIndices;
        for (int j=leftIndices[bound+1]; j<i; j++) {
          if (left[bound] == edges[j].left)
            continue;
          if (AreaABC(left[bound]->transformed2d, edges[j].right->transformed2d, edges[i].left->transformed2d) > 0) 
            newPoints.push_back(edges[j].right);
          else {
            PTR<Object<PV2> > int2d = new ABintersectCD(edges[j].left->transformed2d, edges[j].right->transformed2d, left[bound]->transformed2d, edges[i].left->transformed2d);
            PTR<Point> int3d = new ABintersectCDto3D(edges[j].left, edges[j].right, left[bound], edges[i].left);
            newPoints.push_back(new PathVertex(int3d, int2d));
          }
          newPointsIndices.push_back(j);
        }
        //remove all elements of left after bound;
        left.erase(left.begin()+bound+1, left.end());
        leftIndices.erase(leftIndices.begin()+bound+1, leftIndices.end());
        //replace them with newPoints
        left.insert(left.end(), newPoints.begin(), newPoints.end());
        leftIndices.insert(leftIndices.end(), newPointsIndices.begin(), newPointsIndices.end());
      }
    }
    left.push_back(edges[i].left);
    leftIndices.push_back(i);
  }

  //add 3d points of left to path
  path.clear();
  for (int i=0; i<left.size(); i++) {
    path.push_back(left[i]->original);
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
  std::vector<PathVertex * > vertices;
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
    //create the 3 vertices
    PathVertex * v0 = getVertex(ce->getT()->getP(), vertices);
    PathVertex * v1 = getVertex(ce->getH()->getP(), vertices);
    PathVertex * v2 = getVertex(p2, vertices);
    if (i ==1) 
      triangles.push_back(PathTriangle(getVertex(p1, vertices), v0, v1, hf1));
    triangles.push_back(PathTriangle(v0, v1, v2, hf2));
  }
  PTR<Transformation> xyplane = new XYPlaneTriangleTransfromation(triangles[0].p[0]->original, triangles[0].p[1]->original, triangles[0].p[2]->original);
  PathVertex * start = new PathVertex((PTR<Point>) a);
  PathVertex * end = new PathVertex((PTR<Point>) b);
  flattenTriangles(triangles, transformations, xyplane, start, end);
  savePathTriangles(triangles, "triangles.vtk", false);
  savePathTriangles(triangles, "flattened.vtk", true);

  shortestPath(triangles, start, end, path);
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
