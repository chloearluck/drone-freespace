#include "path3d.h"

void save(Points path, const char * filename) {
  ofstream ostr;
  ostr.open(filename);
  if (ostr.is_open()) { 
    ostr << setprecision(20) << "# vtk DataFile Version 3.0" << endl
         << "vtk output" << endl << "ASCII" << endl
         << "DATASET POLYDATA" << endl 
         << "POINTS " << path.size() << " double" << endl;
    for (int i=0; i<path.size(); i++)
      ostr << path[i]->getApprox(1e-16).getX().mid() << " " << path[i]->getApprox(1e-16).getY().mid() << " " << path[i]->getApprox(1e-16).getZ().mid() << endl;
    ostr<<endl<<"LINES 1 "<<path.size()+1<<endl<<path.size()<<" ";
    for (int i=0; i<path.size(); i++)
      ostr<<i<<" ";
    ostr.close();
  } else {
    cout<<"could not write to file"<<endl;
  }
}

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
  UnfoldTriangleTransformation(HFace * hf1, HFace * hf2) {
    PTR<Point> p1[3];
    p1[0] = hf1->getF()->getBoundary(0)->tail()->getP();
    p1[1] = hf1->getF()->getBoundary(0)->getNext()->tail()->getP();
    p1[2] = hf1->getF()->getBoundary(0)->getNext()->getNext()->tail()->getP();
  
    PTR<Point> p2[3];
    p2[0] = hf2->getF()->getBoundary(0)->tail()->getP();
    p2[1] = hf2->getF()->getBoundary(0)->getNext()->tail()->getP();
    p2[2] = hf2->getF()->getBoundary(0)->getNext()->getNext()->tail()->getP();

    int j1, j2;
    for (int i=0; i<3; i++)
      if (p1[i] != p2[0] && p1[i] != p2[1] && p1[i] != p2[2])
        j1 = i;
    
    for (int i=0; i<3; i++)
      if (p1[0] != p2[i] && p1[1] != p2[i] && p1[2] != p2[i])
        j2 = i;
   
    this->pa = p1[j1];
    this->pb = p1[(j1+1)%3];      
    this->pc = p1[(j1+2)%3];
    this->pd = p2[j2];      
  }
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

Primitive3(AreaABC, PTR<Object<PV2> >, pa, PTR<Object<PV2> >, pb, PTR<Object<PV2> >, pc);
int AreaABC::sign() {
  PV2 a = pa->get();
  PV2 b = pb->get();
  PV2 c = pc->get();
  return (b-a).cross(c-a).sign();
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
    return p3d + t *(q3d-p3d);
  }
};

class PathTriangle {
 public:
  PathVertex * p[3];
  HFace * hface;
  PTR<Transformation> cumulative;
  PathVertex * getVertex(PTR<Point> q) {
    for (int i=0; i<3; i++)
      if (p[i]->original == q)
        return p[i];
    return NULL;
  }
  PathTriangle(PathVertex * p0, PathVertex * p1, PathVertex * p2, HFace * hface) {
    p[0] = p0; p[1] = p1; p[2] = p2; this->hface = hface;
  }
  PathTriangle(HFace * hf,  PathTriangle prev) {
    PTR<Point> ps[3];
    ps[0] = hf->getF()->getBoundary(0)->tail()->getP();
    ps[1] = hf->getF()->getBoundary(0)->getNext()->tail()->getP();
    ps[2] = hf->getF()->getBoundary(0)->getNext()->getNext()->tail()->getP();
    
    int uncommon = -1;
    for (int i=0; i<3; i++)
      if (ps[i] != prev.p[0]->original  &&  ps[i] != prev.p[1]->original  &&  ps[i] != prev.p[2]->original)
        uncommon = i;
    assert(uncommon > -1);

    PTR<Point> p = ps[uncommon];
    ps[uncommon] = ps[2];
    ps[2] = p;

    this->p[0] = prev.getVertex(ps[0]); 
    this->p[1] = prev.getVertex(ps[1]); 
    this->p[2] = new PathVertex(ps[2]);
    this->hface = hf;
  }
};

void savePathTriangles(std::vector<PathTriangle> & ts, const char * filename, bool saveTranformed) {
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

void flattenTriangles(std::vector<PathTriangle> & triangles, std::vector<PTR<Transformation> > & transformations, PTR<Transformation> xyplane, PathVertex * start, PathVertex * end) {
  PTR<Transformation> cumulative = xyplane;
  triangles[0].p[0]->transformed2d = new XYComponents(new TransformedPoint(triangles[0].p[0]->original, cumulative));
  triangles[0].p[1]->transformed2d = new XYComponents(new TransformedPoint(triangles[0].p[1]->original, cumulative));
  triangles[0].p[2]->transformed2d = new XYComponents(new TransformedPoint(triangles[0].p[2]->original, cumulative));
  start->transformed2d = new XYComponents(new TransformedPoint(start->original, cumulative));
  triangles[0].cumulative = cumulative;

  for (int i=1; i<triangles.size(); i++) {
    if (triangles[i].p[0]->transformed2d == 0)
      triangles[i].p[0]->transformed2d = new XYComponents(new TransformedPoint(triangles[i].p[0]->original, cumulative));
    if (triangles[i].p[1]->transformed2d == 0)
      triangles[i].p[1]->transformed2d = new XYComponents(new TransformedPoint(triangles[i].p[1]->original, cumulative));
    
    cumulative = new CompositeTransformation(cumulative, transformations[i-1]);
    triangles[i].p[2]->transformed2d = new XYComponents(new TransformedPoint(triangles[i].p[2]->original, cumulative));
    triangles[i].cumulative = cumulative;
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
  int uncommon = -1;
  int commonCount = 0;
  for (int i=0; i<3; i++)
    if (t1.p[i] != t2.p[0] && t1.p[i] != t2.p[1] && t1.p[i] != t2.p[2])
      uncommon = i;
    else
      commonCount++;
  assert(uncommon > -1); assert(commonCount == 2);
  p = t1.p[(uncommon+1)%3];
  q = t1.p[(uncommon+2)%3];
}

class PathEdge {
 public:
  PathVertex * left;
  PathVertex * right;
  PathEdge(PathVertex * left, PathVertex * right) : left(left), right(right) {}
};

bool hasVertexPoint(HFace * hf, PTR<Point> p) {
  PTR<Point> p0 = hf->getF()->getBoundary(0)->tail()->getP();
  PTR<Point> p1 = hf->getF()->getBoundary(0)->getNext()->tail()->getP();
  PTR<Point> p2 = hf->getF()->getBoundary(0)->getNext()->getNext()->tail()->getP();
  return ((p0 == p) || (p1 == p) || (p2 == p));
}

//DEBUG
class State {
 public:
  std::vector<PathTriangle> triangles;
  std::vector<PathVertex*> path;
  State(std::vector<PathTriangle> & ts, std::vector<PathVertex*> & p, int n) {
    triangles.insert(triangles.begin(), ts.begin(), ts.begin() + n);
    path.insert(path.begin(), p.begin(), p.end());
  }
  void save(const char * filename) {
    ofstream ostr;
    ostr.open(filename);
    if (ostr.is_open()) {
      ostr << setprecision(20) << triangles.size()<<" 0"<<endl;
      for (int i=0; i<triangles.size(); i++) 
        for (int j=0; j<3; j++) 
          ostr<<triangles[i].p[j]->transformed2d->getApprox(1e-16).getX().mid()<<" "<<triangles[i].p[j]->transformed2d->getApprox(1e-16).getY().mid()<<endl;
      // ostr<<endl;
      for (int i=0; i<path.size(); i++)
        ostr<<path[i]->transformed2d->getApprox(1e-16).getX().mid()<<" "<<path[i]->transformed2d->getApprox(1e-16).getY().mid()<<endl;
      ostr.close();
    }
  }
};
class State2 {
 public:
  std::vector<PathEdge> edges;
  std::vector<PathVertex*> path;
  State2(std::vector<PathEdge> & es, std::vector<PathVertex*> & p, int n) {
    edges.insert(edges.begin(), es.begin(), es.begin() + n);
    path.insert(path.begin(), p.begin(), p.end());
  }
  void save(const char * filename) {
    ofstream ostr;
    ostr.open(filename);
    if (ostr.is_open()) {
      ostr << setprecision(20)<< edges.size()<<" 0"<<endl;
      for (int i=0; i<edges.size(); i++) {
        ostr<<edges[i].left->transformed2d->getApprox(1e-16).getX().mid()<<" "<<edges[i].left->transformed2d->getApprox(1e-16).getY().mid()<<endl;
        ostr<<edges[i].right->transformed2d->getApprox(1e-16).getX().mid()<<" "<<edges[i].right->transformed2d->getApprox(1e-16).getY().mid()<<endl;
      }  
      for (int i=0; i<path.size(); i++)
        ostr<<path[i]->transformed2d->getApprox(1e-16).getX().mid()<<" "<<path[i]->transformed2d->getApprox(1e-16).getY().mid()<<endl;
      ostr.close();
    }
  }
};
//DEBUG

void shortestPathHelper(std::vector<PathEdge> & edges, std::vector<PathVertex*> & newPoints, std::vector<int> & newPointsIndices, int startIndex, int endIndex, PathVertex * a, PathVertex * b) {
  //alway use the lowest index for a point that is the path point for multiple edges
  if (edges[endIndex].right == b || edges[endIndex].left == b) {
    shortestPathHelper(edges, newPoints, newPointsIndices, startIndex, endIndex-1, a, b);
    if (newPoints.size()==0 || newPoints[newPoints.size()-1] !=  b) {
      newPoints.push_back(b);
      newPointsIndices.push_back(endIndex);
    }
    return; 
  }

  if (edges[startIndex].right == a || edges[startIndex].left == a) {
    shortestPathHelper(edges, newPoints, newPointsIndices, startIndex+1, endIndex, a, b);
    return;
  }

  for (int i=endIndex; i>=startIndex; i--) {
    if (AreaABC(a->transformed2d, edges[i].right->transformed2d, b->transformed2d) > 0) {
      shortestPathHelper(edges, newPoints, newPointsIndices, startIndex, i-1, a, edges[i].right);
      if (newPoints.size()==0 || newPoints[newPoints.size()-1] != edges[i].right) {//avoid duplicate points on the path
        newPoints.push_back(edges[i].right);
        newPointsIndices.push_back(i);
      }
      shortestPathHelper(edges, newPoints, newPointsIndices, i+1, endIndex, edges[i].right, b);
      return;
    }

    if (AreaABC(a->transformed2d, edges[i].left->transformed2d, b->transformed2d) < 0) {
      shortestPathHelper(edges, newPoints, newPointsIndices, startIndex, i-1, a, edges[i].left);
      if (newPoints.size()==0 || newPoints[newPoints.size()-1] != edges[i].left) {//avoid duplicate points on the path
        newPoints.push_back(edges[i].left);
        newPointsIndices.push_back(i);
      }
      shortestPathHelper(edges, newPoints, newPointsIndices, i+1, endIndex, edges[i].left, b);
      return;
    }
  }
}

void shortestPath(std::vector<PathTriangle> & triangles, PathVertex * start, PathVertex * end, Points & path, std::vector<int> & candidatesToSwap, std::vector<PathVertex*> & verts) {
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

  int state_num = 0;
  char s[50];

  left.push_back(start);
  leftIndices.push_back(-1);
  for (int i=0; i<edges.size(); i++) {
    //DEBUG
    int n = i+2; if (i == triangles.size()-1) n = i+1;
    sprintf(s, "state/%d.txt", state_num);
    State(triangles, left, i+1).save(s);
    sprintf(s, "state/edge%d.txt", state_num);
    State2(edges, left, i+1).save(s); 
    state_num++;
    //-----

    if (left[left.size()-1] == edges[i].left)
      continue;

    int bound;
    if (left.size() > 1)  {
      bound = left.size()-1;
      while ((bound > 0) && (left[bound] != edges[leftIndices[bound]].right) && (AreaABC(left[bound-1]->transformed2d, left[bound]->transformed2d, edges[i].left->transformed2d) > 0) )
        bound--;
    }

    if (left.size() > 1 && bound < left.size()-1) {
      std::vector<PathVertex *> newPoints;
      std::vector<int> newPointsIndices;
      shortestPathHelper(edges, newPoints, newPointsIndices, leftIndices[bound]+1, i, left[bound], edges[i].left);
      //remove all elements of left after bound;
      left.erase(left.begin()+bound+1, left.end());
      leftIndices.erase(leftIndices.begin()+bound+1, leftIndices.end());
      //replace them with newPoints
      left.insert(left.end(), newPoints.begin(), newPoints.end());
      leftIndices.insert(leftIndices.end(), newPointsIndices.begin(), newPointsIndices.end());
    } else {
      left.push_back(edges[i].left);
      leftIndices.push_back(i);
    }
  }

  //DEBUG
  sprintf(s, "state/%d.txt", state_num);
  State(triangles, left, triangles.size()).save(s);
  sprintf(s, "state/edge%d.txt", state_num);
  State2(edges, left, edges.size()).save(s); 
  state_num++;
  //-----

  //save candidatesToSwap
  for (int i=1; i<left.size()-1; i++) 
      candidatesToSwap.push_back(leftIndices[i]);
      int j=leftIndices[i]+1;
      while(edges[j].left == left[i] || edges[j].right == left[i])
        j++;
      candidatesToSwap.push_back(j);
      verts.push_back(left[i]);
    }

  //populate path with point in left and edge intersections
  path.push_back(start->original);
  int j = 1;
  for (int i=0; i<edges.size(); i++) {
    if (leftIndices[j] == i) {
      path.push_back(left[j]->original);
      j++;
    } else {
      assert(leftIndices[j-1] < i); assert(i < leftIndices[j]);
      if (left[j-1] != edges[i].left && left[j-1] != edges[i].right)
        path.push_back(new ABintersectCDto3D(edges[i].left, edges[i].right, left[j-1], left[j]));
    }
  }
}

//populate newPath with a sequence of PathTriangles from oldPath[startIndex] to oldPathpendIndex] with goes
//the other way around pv
void otherWay(vector<PathTriangle> & oldPath, vector<PathTriangle> & newPath, int startIndex, int endIndex, PathVertex * pv) {
  std::vector<HFace*> neighbors;
  HFace * nextHF = NULL;
  HFace * hfaceStart = oldPath[startIndex].hface;
  hfaceStart->neighbors(neighbors);
  assert(neighbors.size() == 3);
  
  for (int i=0; i< neighbors.size(); i++)
    if (neighbors[i] != oldPath[startIndex-1].hface && neighbors[i] != oldPath[startIndex+1].hface && hasVertexPoint(neighbors[i],pv->original))
      nextHF = neighbors[i];
  assert(nextHF != NULL);
  PathTriangle next(nextHF, oldPath[startIndex]);
  newPath.push_back(next);
  while (true) {
    neighbors.clear();
    HFace * prev = ((newPath.size()==1)? oldPath[startIndex].hface : newPath[newPath.size()-2].hface);
    newPath[newPath.size()-1].hface->neighbors(neighbors);
    assert(neighbors.size() == 3);
    for (int i=0; i<neighbors.size(); i++)
      if (neighbors[i] != prev && hasVertexPoint(neighbors[i],pv->original))
        nextHF = neighbors[i];
    if (nextHF == oldPath[endIndex].hface)
      return;
    PathTriangle next(nextHF, newPath[newPath.size()-1]);
    newPath.push_back(next); 
  }
}

void localPath(PTR<FaceIntersectionPoint> a, PTR<FaceIntersectionPoint> b, HFaces & pathfaces, Points & path) {
  std::vector<PathTriangle> triangles;
  std::vector<PTR<Transformation> > transformations;
  HFace * hf = pathfaces[0];
  PathVertex * v0 = new PathVertex(hf->getF()->getBoundary(0)->tail()->getP());
  PathVertex * v1 = new PathVertex(hf->getF()->getBoundary(0)->getNext()->tail()->getP());
  PathVertex * v2 = new PathVertex(hf->getF()->getBoundary(0)->getNext()->getNext()->tail()->getP());
  triangles.push_back(PathTriangle(v0,v1,v2,hf));
  for (int i=1; i<pathfaces.size(); i++) {
    PathTriangle prev= triangles[i-1];
    hf = pathfaces[i];

    PTR<Point> p[3];
    p[0] = hf->getF()->getBoundary(0)->tail()->getP();
    p[1] = hf->getF()->getBoundary(0)->getNext()->tail()->getP();
    p[2] = hf->getF()->getBoundary(0)->getNext()->getNext()->tail()->getP();

    int uncommon1 = -1;
    for (int j=0; j<3; j++)
      if (prev.p[j]->original != p[0] && prev.p[j]->original != p[1] && prev.p[j]->original != p[2])
        uncommon1 = j;
    assert(uncommon1 > -1);

    int uncommon2 = -1;
    for (int j=0; j<3; j++)
      if (p[j] != prev.p[0]->original && p[j] != prev.p[1]->original && p[j] != prev.p[2]->original)
        uncommon2 = j;
    assert(uncommon2 > -1);

    PathTriangle t(prev.p[(uncommon1+1)%3], prev.p[(uncommon1+2)%3], new PathVertex(p[uncommon2]), hf);
    triangles.push_back(t);
    transformations.push_back(new UnfoldTriangleTransformation(prev.p[uncommon1]->original, t.p[0]->original, t.p[1]->original, t.p[2]->original));  
  } 
  PTR<Transformation> xyplane = new XYPlaneTriangleTransfromation(triangles[0].p[0]->original, triangles[0].p[1]->original, triangles[0].p[2]->original);
  PathVertex * start = new PathVertex((PTR<Point>) a);
  PathVertex * end = new PathVertex((PTR<Point>) b);
  flattenTriangles(triangles, transformations, xyplane, start, end);

  std::vector<int> candidatesToSwap;
  std::vector<PathVertex*> verts;
  shortestPath(triangles, start, end, path, candidatesToSwap, verts);

  save(path, "path0.vtk");
  char s[50];  

  bool trueFlip;
  int iter_num  = 0;
  do {
    iter_num++;
    trueFlip = false;
    cout<<"new iteration------------------"<<endl;
    for (int i=1; i<triangles.size(); i++) {
      int j = find(candidatesToSwap.begin(), candidatesToSwap.end(), i) - candidatesToSwap.begin();
      if (j == candidatesToSwap.size() || j%2 == 1)
        continue;

      int startIndex = candidatesToSwap[j]; 
      int endIndex = candidatesToSwap[j+1];

      std:vector<PathTriangle> newPath;
      cout<<"startIndex: "<<startIndex<<" endIndex: "<<endIndex<<endl;

      PathVertex * replacingVertex = verts[j/2];
      otherWay(triangles, newPath, startIndex, endIndex, replacingVertex);

      PTR<Transformation> cumulative = triangles[startIndex].cumulative;
      newPath[0].p[0]->transformed2d = new XYComponents(new TransformedPoint(newPath[0].p[0]->original, cumulative));
      newPath[0].p[1]->transformed2d = new XYComponents(new TransformedPoint(newPath[0].p[1]->original, cumulative));
      PTR<Transformation> t = new UnfoldTriangleTransformation(triangles[startIndex].hface, newPath[0].hface);
      cumulative = new CompositeTransformation(cumulative, t);
      newPath[0].p[2]->transformed2d = new XYComponents(new TransformedPoint(newPath[0].p[2]->original, cumulative));
      newPath[0].cumulative = cumulative;
      for (int j=1; j<newPath.size(); j++) {
        newPath[j].p[0]->transformed2d = new XYComponents(new TransformedPoint(newPath[j].p[0]->original, cumulative));
        newPath[j].p[1]->transformed2d = new XYComponents(new TransformedPoint(newPath[j].p[1]->original, cumulative));
        t = new UnfoldTriangleTransformation(newPath[j-1].hface, newPath[j].hface);
        cumulative = new CompositeTransformation(cumulative, t);
        newPath[j].p[2]->transformed2d = new XYComponents(new TransformedPoint(newPath[j].p[2]->original, cumulative));
        newPath[j].cumulative = cumulative;
      }
      //give triangles[endIndex] new PathVertices to match the new intermediate path
      //also replace the one that doesn't match with a new vertex, it may match something earlier in the path
      int uncommon = -1;
      int commonCount = 0;
      for (int j=0; j<3; j++) {
        int matchingIndex = -1;
        for (int k=0; k<3; k++)
          if (triangles[endIndex].p[j]->original == newPath[newPath.size()-1].p[k]->original)
            matchingIndex = k;
        if (matchingIndex == -1) {
          uncommon = j;
          triangles[endIndex].p[j] = new PathVertex(triangles[endIndex].p[j]->original);
        } else {
          triangles[endIndex].p[j] = newPath[newPath.size()-1].p[matchingIndex];
          commonCount++;
        }
      }
      assert(uncommon > -1); assert(commonCount == 2);
      //swap the the triangle[endIndex]'s vertices so that the uncommon 1 is at position 2
      PathVertex * v_tmp = triangles[endIndex].p[uncommon];
      triangles[endIndex].p[uncommon] = triangles[endIndex].p[2];
      triangles[endIndex].p[2] = v_tmp;
      //replace these vertices in an subsequent triangles 
      for (int j=0; j<3; j++) {
        PathVertex * v = triangles[endIndex].p[j];
        for (int k=endIndex+1; k<triangles.size(); k++) {
          if (triangles[k].p[0]->original == v->original)
            triangles[k].p[0] = v;
          else if (triangles[k].p[1]->original == v->original)
            triangles[k].p[1] = v;
          else if (triangles[k].p[2]->original == v->original)
            triangles[k].p[2] = v;
          else
            break;
        }
      }

      //re-transform the rest of the triangles or else triangle[endIndex] wont necessarily share a face with the last triangle in newPath
      triangles[endIndex].p[0]->transformed2d = new XYComponents(new TransformedPoint(triangles[endIndex].p[0]->original, cumulative));
      triangles[endIndex].p[1]->transformed2d = new XYComponents(new TransformedPoint(triangles[endIndex].p[1]->original, cumulative));
      t = new UnfoldTriangleTransformation(newPath[newPath.size()-1].hface, triangles[endIndex].hface);
      cumulative = new CompositeTransformation(cumulative, t);
      triangles[endIndex].p[2]->transformed2d = new XYComponents(new TransformedPoint(triangles[endIndex].p[2]->original, cumulative));
      triangles[endIndex].cumulative = cumulative;
      for (int j=endIndex+1; j<triangles.size(); j++) {
        triangles[j].p[0]->transformed2d = new XYComponents(new TransformedPoint(triangles[j].p[0]->original, cumulative));
        triangles[j].p[1]->transformed2d = new XYComponents(new TransformedPoint(triangles[j].p[1]->original, cumulative));
        t = new UnfoldTriangleTransformation(triangles[j-1].hface, triangles[j].hface);
        cumulative = new CompositeTransformation(cumulative, t);
        triangles[j].p[2]->transformed2d = new XYComponents(new TransformedPoint(triangles[j].p[2]->original, cumulative));
        triangles[j].cumulative = cumulative;
      }
      end->transformed2d = new XYComponents(new TransformedPoint(end->original, cumulative));

      triangles.erase(triangles.begin()+startIndex+1, triangles.begin()+endIndex);
      triangles.insert(triangles.begin()+startIndex+1, newPath.begin(), newPath.end());

      candidatesToSwap.clear();
      verts.clear();
      path.clear();

      shortestPath(triangles, start, end, path, candidatesToSwap, verts);

      //does this new path still contain replacingVertex at edge startIndex?
      //if it doesn't, set trueFlip to true
      //how can we tell.  If theres a flush vertex, it will be in candidatesToSwap
      //meaning we should find a j where j%2 ==0 
      // and candidatesToSwap[j]<= startIndex and candidatesToSwap[j+1] > startIndex
      // and verts[j] = replacingVertex
      for (int j=0; j<candidatesToSwap.size(); j++) {
        if (candidatesToSwap[j] > startIndex)
          break;
        if (candidatesToSwap[j+1] > startIndex && verts[j/2] != replacingVertex) {
          cout<<"trueFlip"<<endl;
          trueFlip = true;
        }
      }

      //DEBUG save path
      sprintf(s, "path%d-%02d.vtk", iter_num, i);
      save(path, s);
      //DEBUG
    }
  } while(trueFlip);
}

void bfs(HFace * fa, HFace * fb, HFaces & pathfaces) {
  HFaces pathfaces_rev;

  std::map<HFace *, HFace *> parents;
  std::queue<HFace *> q;
  q.push(fa);
  parents[fa] = NULL;

  while (!q.empty()) {
    HFace * current = q.front(); q.pop();
    if (current == fb) {
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
    cout<<"no path found by bfs"<<endl;
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
  int cell_index = blockspace->containingCell(start);
  if (blockspace->containingCell(end) != cell_index) {
    cout<<"start and end are not in the same component"<<endl;
    return;
  }

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
    bfs(points[i]->getHFace(), points[i+1]->getHFace(), subHfaces);
    localPath(points[i], points[i+1], subHfaces, subPath);
    path.insert(path.end(), subPath.begin(), subPath.end());
  }
  if (points.size()>0)
    path.push_back((PTR<Point>) points[points.size()-1]);
  path.push_back(end);
}
