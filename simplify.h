#ifndef SIMPLIFY
#define SIMPLIFY

#include <queue>
#include "polyhedron.h"
#include "io.h"
#include "expander2.h"

void simplify (Polyhedron *a, double d, bool opt2 = false);

void simplify1 (Polyhedron *a, double d);

FOctree * faceOctree (Polyhedron *a, double s = 0.0);

typedef pair<Vertex *, PTR<Point> > VPPair;

typedef map<Vertex *, PTR<Point> > VPMap;

bool simplify1 (Polyhedron *a, double d, FOctree *octree, VPMap &vpmap,
		int &n1, int &n2);

class DistancePP : public Primitive {
  Point *a, *b;
  double d;
  int sign () {
    PV3 u = a->getP() - b->getP();
    return (d*d - u.dot(u)).sign();
  }
 public:
  DistancePP (Point *a, Point *b, double d) : a(a), b(b), d(d) {}
};

class DistancePL1 : public Primitive {
  Point *a, *t, *h;
  double d;
  int sign () {
    PV3 tp = t->getP(), u = h->getP() - tp,
      p = tp + ((a->getP() - tp).dot(u)/u.dot(u))*u, w = a->getP() - p;
    return (d*d - w.dot(w)).sign();
  }
 public:
  DistancePL1 (Point *a, Point *t, Point *h, double d) : a(a), t(t), h(h), d(d) {}
};

bool closeVE (Vertex *v, Vertex *w, Vertex *x, double d);

class IFeature {
 public:
  Edge *e;
  Vertex *v;
  bool forward;

  IFeature (Edge *e) : e(e), v(0), forward(false) {}
  IFeature (HEdge *h) : e(h->getE()), v(h->getNext()->head()),
        forward(h->getForward()) {}
  
  HEdge * getH () const {
    HEdge *h = e->getHEdge(0);
    return h->getForward() == forward ? h : h->ccw();
  }

  bool valid () const {
    return e->HEdgesN() == 2 && (!v || getH()->getNext()->head() == v);
  }
  
  bool operator< (const IFeature &x) const {
    if (!x.v && v)
      return true;
    if (x.v && !v)
      return false;
    if (v == x.v && e == x.e && forward == x.forward)
      return false;
    int s = CloserPair(x.e->getT()->getP(), x.e->getH()->getP(),
		       e->getT()->getP(), e->getH()->getP()); 
    return s == 1 || s == 0 && e < x.e;
  }
};

typedef priority_queue<IFeature> IFeatureQ;

void addCollapseFlips (Edge *e, double d, IFeatureQ &fq);

void addFlips (Edge *e, double d, IFeatureQ &fq);

bool collapse (Polyhedron *a, double d, Edge *e, IFeatureQ &fq, FOctree *octree,
	       VPMap &vpmap, bool &bad);

bool collapsible (Edge *e);

double collapseDistance (Edge *e);

bool badCollapse (Vertex *v, FOctree *octree, double dc);

void newEdges (Face *f,  Vertex *v, bool *ef);

bool intersects (Face *f, bool *evf, Face *g);

bool intersectsFFP (Face *f, bool *ef, Face *g, bool *eg);

bool intersectsPE (Face *f, bool *evf, Face *g, bool *evg, Edges &iedges,
		   Faces &ifaces, bool &pflag);

bool intersectsFEP (Face *f, bool *ef, Edge *e, bool ee);

bool intersectsEE (Edge *e, Edge *f, int pc);

void addStar (Vertex *v, double d, IFeatureQ &fq);

bool flip (Polyhedron *a, double d, HEdge *e, IFeatureQ &fq, FOctree *octree,
	   VPMap &vpmap, bool &bad);

bool flippable (HEdge *e, double d);

HEdge * flip (Polyhedron *a, HEdge *e);

bool badFlip (HEdge *e, FOctree *octree);

void newEdges (Face *f, Edge *e, bool *ef);

void describe1 (double t, int n1, int n2, const VPMap &vpmap);

double displacement (const VPMap &vpmap, double &dmax);

class LinePoint : public Point {
  Point *a, *t, *h;
  PV3 calculate () {
    PV3 ap = a->getP(), tp = t->getP(), u = h->getP() - tp;
    return tp + (u.dot(ap - tp)/u.dot(u))*u;
  }
 public:
  LinePoint (Point *a, Point *t, Point *h) : a(a), t(t), h(h) {
    IDSet psth = intersection(t->getps(), h->getps());
    for (IDSet::iterator i = psth.begin(); i != psth.end(); ++i)
      ps.insert(*i);
  }
};

class PlanePoint : public Point {
  Plane *p;
  Point *a;
  PV3 calculate () {
    PV3 n = p->getN(), q = a->getP();
    return  q - ((q.dot(n) + p->getK())/n.dot(n))*n;
  }
 public:
  PlanePoint (Plane *p, Point *a) : p(p), a(a) {
    ps.insert(p->getid());
  }
};

class LineLinePoint : public Point {
  Point *a, *b, *c, *d;
  PV3 calculate () {
    PV3 ap = a->getP(), u = b->getP() - ap, cp = c->getP(), v = d->getP() - cp,
      w = cp - ap;
    Parameter k1 = u.dot(u), k2 = u.dot(v), k4 = - v.dot(v), k5 = u.dot(w),
      k6 = v.dot(w), den = k1*k4 + k2*k2, k = (k5*k4 + k2*k6)/den;
      return ap + k*u;
  }
 public:
  LineLinePoint (Point *a, Point *b, Point *c, Point *d)
    : a(a), b(b), c(c), d(d) {
    IDSet psab = intersection(a->getps(), b->getps());
    for (IDSet::iterator i = psab.begin(); i != psab.end(); ++i)
      ps.insert(*i);
  }
};

void simplify2 (Polyhedron *a, double d, bool opt2);

bool closeVF (Vertex *v, Face *f, double d);

bool closeEE (Edge *e, Edge *f, double d);

enum FeatureType {VV, VE, VF, EE};

class Feature {
 public:
  Feature () : v(0), e(0), f(0), fa(0), p(0), q(0) {}
  
  Feature (Vertex *vi, Vertex *wi, bool flag = true) : type(VV),
    v(vi < wi ? vi : wi), w(vi < wi ? wi : vi), e(0), f(0), fa(0),
    flag(flag), p(0), q(0) {}
  
  Feature (Vertex *v, Edge *e, bool flag = true)
    : type(VE), v(v), w(0), e(e), f(0), fa(0), flag(flag), p(0), q(0) {}
  
  Feature (Vertex *v, Face *fa, bool flag = true)
    : type(VF), v(v), w(0), e(0), f(0), fa(fa), flag(flag), p(0), q(0) {}
  
  Feature (Edge *ei, Edge *fi, bool flag = true) : type(EE), v(0), w(0),
    e(ei < fi ? ei : fi), f(ei < fi ? fi : ei), fa(0), flag(flag), p(0), q(0) {}
  
  bool operator< (const Feature &x) const {
    if (type != x.type)
      return type < x.type;
    switch (type) {
    case VV: return v < x.v || v == x.v && w < x.w;
    case VE: return v < x.v || v == x.v && e < x.e;
    case VF: return v < x.v || v == x.v && fa < x.fa;
    case EE: return e < x.e || e == x.e && f < x.f;
    }
  }

  bool small (double d) const {
    switch (type) {
    case VV: return DistancePP(v->getP(), w->getP(), d) == 1;
    case VE: return closeVE(v, e->getT(), e->getH(), d);
    case VF: return closeVF(v, fa, d);
    case EE: return closeEE(e, f, d);
    }
  }
  
  void closestPoints (PV3 &pp, PV3 &qp) {
    if (type == EE) {
      Point *et = e->getT()->getP(), *eh = e->getH()->getP(),
	*ft = f->getT()->getP(), *fh = f->getH()->getP();
      p = new LineLinePoint(et, eh, ft, fh);
      q = new LineLinePoint(ft, fh, et, eh);
    }
    else {
      p = v->getP();
      if (type == VV)
	q = w->getP();
      else if (type == VE)
	q = new LinePoint(v->getP(), e->getT()->getP(), e->getH()->getP());
      else
	q = new PlanePoint(fa->getP(), v->getP());
    }
    pp = p->getP();
    qp = q->getP();
  }

  bool point (Point *a) const {
    VertexSet vs;
    vertices(vs);
    for (VertexSet::iterator i = vs.begin(); i != vs.end(); ++i)
      if (a == (*i)->getP())
	return true;
    return false;
  }
  
  void vertices (VertexSet &vs) const {
    if (type == EE) {
      vs.insert(e->getT());
      vs.insert(e->getH());
      vs.insert(f->getT());
      vs.insert(f->getH());
    }
    else {
      vs.insert(v);
      if (type == VV)
	vs.insert(w);
      else if (type == VE) {
	vs.insert(e->getT());
	vs.insert(e->getH());
      }
      else {
	Vertices ve;
	fa->boundaryVertices(ve);
	vs.insert(ve.begin(), ve.end());
      }
    }
  }

  FeatureType type;
  Vertex *v, *w;
  Edge *e, *f;
  Face *fa;
  bool flag;
  PTR<Point> p, q;
};

typedef set<Feature> FeatureSet;

FeatureSet smallFeatures (Polyhedron *a, double d, VPMap *vpmap);

class TData {
 public:
  FOctree *octree;
  double d;
  VPMap *vpmap;
  FeatureSet fs;
};

void smallCandidatesT (void *ptr);

bool moved (Face *f, const VPMap &vpmap);

void smallCandidates (Face *f, Face *g, double d, FeatureSet &fs);

bool closeVVT (Vertex *v, Vertex *w, double d);

bool closeVET (Vertex *v, Edge *e, double d);

bool closeVFT (Vertex *v, Face *f, double d);

bool closeEET (Edge *e, Edge *f, double d);

FeatureSet smallFeatures (const FeatureSet &fs, double d);

PV3 orthogonal (const PV3 &u);

double separation (const FeatureSet &fs);

void simplify2s (Polyhedron *a, double d, FeatureSet &fs, VPMap &vpmap,
		 bool vobj = true, double vbound = 1.0);

double * solveLP (const FeatureSet &fs, double d, const VPMap &vpmap,
		  bool vobj, double vbound, VIMap &vimap);

void setupLP (const FeatureSet &fs, double d, VIMap &vimap, Expander2 &exp);

Primitive4(InnerProduct4, Point *, a, Point *, b, Point *, c, Point *, d);

void featuresLPVV (Vertex *v, Vertex *w, FeatureSet &fs);

void featuresLPVE (Vertex *v, Edge *e, FeatureSet &fs);

void setupLP (const Feature &f, double d, VIMap &vimap, Expander2 &exp,
	      bool minimal);

Feature minimalFeature (const Feature &f);

Feature minimalFeature (Vertex *v, Face *f);

Feature minimalFeature (Edge *e, Edge *f);

class SeparatorData {
 public:
  SeparatorData () {}
  
  SeparatorData (const PV3 &p, const PV3 &q, const PV3 &u, const PV3 &v,
		 const PV3 &w, const Parameter &pq)
    : p(p), q(q), u(u), v(v), w(w), pq(pq) {}
  
  int size () const { return 16; }
  
  Parameter & operator[](int i) {
    if (i < 6)
      return i < 3 ? p[i] : q[i-3];
    if (i < 12)
      return i < 9 ? u[i-6] : v[i-9];
    return i < 15 ? w[i-12] : pq;
  }
  
  PV3 p, q, u, v, w;
  Parameter pq;
};

class Separator : public Object<SeparatorData> {
  Feature f;
  double d;
  
  SeparatorData calculate () {
    PV3 p, q;
    f.closestPoints(p, q);
    PV3 u = f.flag ? q - p : p - q,
      v = f.type == VE ? f.e->getU() : orthogonal(u),
      w = u.cross(v);
#ifdef UNIT_U
    Parameter pq = u.dot(u).sqrt();
    return SeparatorData(p, q, u/pq, v.unit(), w.unit(), pq);
#else
    Parameter pq = u.dot(u);
    return SeparatorData(p, q, u, v, w, pq);
#endif
  }
 public:
  Separator (const Feature &f, double d) : f(f), d(d) {}
  
  PV3 getU () { return getApprox(1e-7).u; }
  
  double getPQ () { return (getApprox(1e-7).pq/d).mid(); }

  class VWR : public Object<PV3> {
  public:
    Separator *s;
    Point *a;
    bool flag;
    VWR (Separator *s, Point *a, bool flag) : s(s), a(a), flag(flag) {}
    PV3 calculate () {
      if (a == s->f.p || a == s->f.q)
	return PV3::constant(0, 0, 0);
      PV3 p = a->getP() - (flag ? s->get().p : s->get().q);
      double d = s->d;
      Parameter zero = Parameter::constant(0);
      PV3 vwr(s->get().v.dot(p)/d, 
              s->f.type == VE && s->f.point(a) ? zero : s->get().w.dot(p)/d,
	      s->f.point(a) ? zero : s->get().u.dot(p)/d);
      return vwr;
    }
  };      
  
  void getVWR (Point *a, bool flag, double *vwr) {
    VWR vwrObj(this, a, flag);
    PV3 p = vwrObj.getApprox(1e-7);
#ifdef UNIT_U
    vwr[0] = p.x.mid();
    vwr[1] = p.y.mid();
    vwr[2] = p.z.mid();
#else
    SeparatorData sd = getApprox(1e-7);
    vwr[0] = p.x.mid() / sd.v.length().mid();
    vwr[1] = p.y.mid() / sd.w.length().mid();
    vwr[2] = p.z.mid() / sd.u.length().mid();
#endif
  }
};

int getVertex (Vertex *v, VIMap &vimap);

void getDV (Vertex *v, double d, const VPMap &vpmap, double *dv);

bool check (Separator *s, const Vertices &v, const Vertices &w);

Primitive3(SeparatorOrder, Separator *, s, Point *, a, Point *, b);

void moveVertices (Polyhedron *a, VPMap &vpmap, const VIMap &vimap, double *dvw);

FeatureSet smallFeatures (Polyhedron *a, double d, const FeatureSet &fs,
			  VPMap *vpmap);

VPMap savePos (const FeatureSet &fs);

void restorePos (Polyhedron *a, const VPMap &vpmap);

void describe2 (int n, const FeatureSet &fs, double t, const VPMap &vpmap);

#endif
