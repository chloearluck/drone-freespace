#ifndef SIMPLIFY
#define SIMPLIFY

#include <queue>
#include "polyhedron.h"
#include "io.h"
#include "expander2.h"
#include "poly.h"

void simplify (Polyhedron *a, double d, bool opt2 = false);

void simplify1 (Polyhedron *a, double d);

class DistancePP : public Primitive {
  Point *a, *b;
  double d;

  Parameter calculate () {
    PV3 u = a->get() - b->get();
    return (d*d - u.dot(u));
  }
 public:
  DistancePP (Point *a, Point *b, double d) : a(a), b(b), d(d) {}
};

class DistancePL : public Primitive {
  Point *a, *t, *h;
  double d;

  Parameter calculate () {
    PV3 tp = t->get(), u = h->get() - tp,
      p = tp + ((a->get() - tp).dot(u)/u.dot(u))*u, w = a->get() - p;
    return d*d - w.dot(w);
  }
 public:
  DistancePL (Point *a, Point *t, Point *h, double d) : a(a), t(t), h(h), d(d) {}
};

bool closeVE (Vertex *v, Vertex *w, Vertex *x, double d);

Octree<Face *> * faceOctree (Polyhedron *a, double s = 0.0);

typedef pair<Vertex *, PTR<Point> > VPPair;

typedef map<Vertex *, PTR<Point> > VPMap;

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
  
  bool operator< (const IFeature &x) const { return !x.v && v; }
};

typedef priority_queue<IFeature> IFeatureQ;

void addCollapseFlips (Edge *e, double d, IFeatureQ &fq);

void addFlips (Edge *e, double d, IFeatureQ &fq);

void simplify1 (Polyhedron *a, double d, IFeatureQ &fq, IFeatureQ &nfq,
		Octree<Face *> *octree, VPMap &vpmap, int &n1, int &n2);

bool collapse (Polyhedron *a, double d, Edge *e, IFeatureQ &fq, IFeatureQ &nfq,
	       Octree<Face *> *octree, VPMap &vpmap);

bool collapsible (Edge *e);

bool badCollapse (Vertex *v, Octree<Face *> *octree);

void addStar (Vertex *v, double d, IFeatureQ &fq);

bool flip (Polyhedron *a, double d, HEdge *e, IFeatureQ &fq, IFeatureQ &nfq,
	   Octree<Face *> *octree, VPMap &vpmap);

bool flippable (HEdge *e, double d);

HEdge * flip (Polyhedron *a, HEdge *e);

bool badFlip (HEdge *e, Octree<Face *> *octree);

void describe1 (double t, int n1, int n2, const VPMap &vpmap);

double displacement (const VPMap &vpmap, double &dmax);

class PlanePoint : public Point {
  Plane *p;
  Point *a;

  PV3 calculate () {
    PV3 n = p->get().n, q = a->get();
    return  q - ((q.dot(n) + p->get().k)/n.dot(n))*n;
  }
 public:
  PlanePoint (Plane *p, Point *a) : p(p), a(a) {
    ps.insert(p->getid());
  }
};

class LinePoint : public Point {
  Point *a, *t, *h;
  
  PV3 calculate () {
    PV3 ap = a->get(), tp = t->get(), u = h->get() - tp;
    return tp + (u.dot(ap - tp)/u.dot(u))*u;
  }
 public:
  LinePoint (Point *a, Point *t, Point *h) : a(a), t(t), h(h) {
    set<ID> psth = intersection(t->getps(), h->getps());
    for (set<ID>::iterator i = psth.begin(); i != psth.end(); ++i)
      ps.insert(*i);
  }
};

class LineLinePoint : public Point {
  Point *a, *b, *c, *d;
  
  PV3 calculate () {
    PV3 ap = a->get(), u = b->get() - ap, cp = c->get(), v = d->get() - cp,
      w = cp - ap;
    Parameter k1 = u.dot(u), k2 = u.dot(v), k4 = - v.dot(v), k5 = u.dot(w),
      k6 = v.dot(w), den = k1*k4 + k2*k2, k = (k5*k4 + k2*k6)/den;
    return ap + k*u;
  }
 public:
  LineLinePoint (Point *a, Point *b, Point *c, Point *d)
    : a(a), b(b), c(c), d(d) {
    set<ID> psab = intersection(a->getps(), b->getps());
    for (set<ID>::iterator i = psab.begin(); i != psab.end(); ++i)
      ps.insert(*i);
  }
};

void simplify2 (Polyhedron *a, double d, bool opt2);

enum FeatureType {VV, VE, VF, EE, EF};

class Feature {
 public:
  Feature () : v(0), w(0), e(0), f(0), fa(0) {}
  
  Feature (Vertex *vi, Vertex *wi, bool flag = true) : type(VV),
    v(vi < wi ? vi : wi), w(vi < wi ? wi : vi), e(0), f(0), fa(0),
    flag(flag) {}
  
  Feature (Vertex *v, Edge *e, bool flag = true)
    : type(VE), v(v), w(0), e(e), f(0), fa(0), flag(flag) {}
  
  Feature (Vertex *v, Face *fa, bool flag = true)
    : type(VF), v(v), w(0), e(0), f(0), fa(fa), flag(flag) {}
  
  Feature (Edge *ei, Edge *fi, bool flag = true) : type(EE), v(0), w(0),
    e(ei < fi ? ei : fi), f(ei < fi ? fi : ei), fa(0), flag(flag) {}

  Feature (Edge *e, Face *fa) : type(EF), v(0), w(0), e(e), f(0), fa(fa) {}

  bool operator< (const Feature &x) const {
    if (type != x.type)
      return type < x.type;
    switch (type) {
    case VV: return v < x.v || v == x.v && w < x.w;
    case VE: return v < x.v || v == x.v && e < x.e;
    case VF: return v < x.v || v == x.v && fa < x.fa;
    case EE: return e < x.e || e == x.e && f < x.f;
    case EF: return e < x.e || e == x.e && fa < x.fa;
    default: assert(0);
    }
  }
  
  bool point (Point *a) const {
    set<Vertex *> vs;
    vertices(vs);
    for (set<Vertex *>::iterator i = vs.begin(); i != vs.end(); ++i)
      if (a == (*i)->getP())
	return true;
    return false;
  }
  
  void vertices (set<Vertex *> &vs) const {
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
	Vertices ve = fa->getBoundary()->loop();
	vs.insert(ve.begin(), ve.end());
      }
    }
  }

  FeatureType type;
  Vertex *v, *w;
  Edge *e, *f;
  Face *fa;
  bool flag;
};

typedef set<Feature> FeatureSet;

FeatureSet smallFeatures (Polyhedron *a, double d, VPMap *vpmap);

void smallFeaturesT (void *ptr);

class SFData {
 public:
  unsigned int i, is, ie;
  vector<pair<Face *, Face *> > *ff;
  double d;
  VPMap *vpmap;
  FeatureSet fs;
};

bool moved (Face *f, const VPMap &vpmap);

void smallCandidates (Face *f, Face *g, double d, FeatureSet &fs);

bool firstFace (Vertex *v, Face *f);

bool closeVVT (Vertex *v, Vertex *w, double d);

bool closeVET (Vertex *v, Edge *e, double d);

bool closeVFT (Vertex *v, Face *f, double d);

bool closeEET (Edge *e, Edge *f, double d);

FeatureSet smallFeatures (const FeatureSet &fs, double d);

double separation (const FeatureSet &fs);

void simplify2s (Polyhedron *a, double d, FeatureSet &fs, VPMap &vpmap,
		 bool vobj = true, double vbound = 1.0);

Expander2 * solveLP (const FeatureSet &fs, double d, const VPMap &vpmap,
		     bool vobj, double vbound, set<Vertex *> &vs);

void setupLP (const FeatureSet &fs, double d, set<Vertex *> &vs, Expander2 *exp);

Feature minimalFeature (const Feature &f, PTR<Point> &p, PTR<Point> &q);

Feature minimalFeature (Vertex *v, Face *f, PTR<Point> &p, PTR<Point> &q);

Feature minimalFeature (Edge *e, Edge *f, PTR<Point> &p, PTR<Point> &q);

bool small (const Feature &f, double d);

PV3 orthogonal (const PV3 &u);

class SeparatorData {
 public:
  PV3 p, q, u, v, w;
  Parameter pq;
  
  int size () const { return 16; }

  Parameter & operator[](int i) {
    if (i < 6)
      return i < 3 ? p[i] : q[i-3];
    if (i < 12)
      return i < 9 ? u[i-6] : v[i-9];
    return i < 15 ? w[i-12] : pq;
  }

  const Parameter & operator[](int i) const {
    if (i < 6)
      return i < 3 ? p[i] : q[i-3];
    if (i < 12)
      return i < 9 ? u[i-6] : v[i-9];
    return i < 15 ? w[i-12] : pq;
  }

  SeparatorData () {}
  SeparatorData (const PV3 &p, const PV3 &q, const PV3 &u, const PV3 &v,
		 const PV3 &w, const Parameter &pq)
    : p(p), q(q), u(u), v(v), w(w), pq(pq) {}
};

class Separator : public Object<SeparatorData> {
  Feature f;
  double d;
  PTR<Point> ac, bc;

  SeparatorData calculate () {
    PV3 p = ac->get(), q = bc->get(), u = f.flag ? q - p : p - q,
      v = f.type == VE ? f.e->getU() : orthogonal(u), w = u.cross(v);
    Parameter pq = u.dot(u);
    return SeparatorData(p, q, u, v, w, pq);
  }
  
 public:
  Separator (const Feature &fin, double d) : d(d) {
    f = minimalFeature(fin, ac, bc);
  }
  
  PV3 getU () { return getApprox(1e-7).u; }
  
  double getPQ () { return (getApprox(1e-7).pq/d).mid(); }

  class VWR : public Object<PV3> {
  public:
    Separator *s;
    Point *a;
    bool flag;
    VWR (Separator *s, Point *a, bool flag) : s(s), a(a), flag(flag) {}
    PV3 calculate () {
      if (a == s->ac || a == s->bc)
	return PV3::constant(0.0, 0.0, 0.0);
      PV3 p = a->get() - (flag ? s->get().p : s->get().q);
      double d = s->d;
      Parameter zero = Parameter::constant(0.0);
      PV3 vwr(s->get().v.dot(p)/d, 
	      s->f.type == VE && s->f.point(a) ? zero : s->get().w.dot(p)/d,
	      s->f.point(a) ? zero : s->get().u.dot(p)/d);
      return vwr;
    }
  };      
  
  void getVWR (Point *a, bool flag, double *vwr) {
    VWR vwrObj(this, a, flag);
    PV3 p = vwrObj.getApprox(1e-7);
    SeparatorData sd = getApprox(1e-7);
    vwr[0] = p.x.mid() / sd.v.length().mid();
    vwr[1] = p.y.mid() / sd.w.length().mid();
    vwr[2] = p.z.mid() / sd.u.length().mid();
  }
};

void getDV (Vertex *v, double d, const VPMap &vpmap, double *dv);

void moveVertices (Polyhedron *a, double d, VPMap &vpmap, const set<Vertex *> &vs,
		   Expander2 *exp);

FeatureSet smallFeatures (Polyhedron *a, double d, const FeatureSet &fs,
			  VPMap *vpmap);

VPMap savePos (const FeatureSet &fs);

void restorePos (Polyhedron *a, const VPMap &vpmap);

void describe2 (int n, const FeatureSet &fs, double t, const VPMap &vpmap);

bool intersects (Expander2 *e, Expander2::Pair *p);

bool intersects (const Points &pts, const Points &dsp, bool vf);

class FIPoly : public Object<PPoly> {
  Points pts, dsp;

  PPoly calculate () {
    PPoly p;
    addTerms(p, 0, 1, 2);
    addTerms(p, 0, 2, 3);
    addTerms(p, 0, 3, 1);
    addTerms(p, 2, 1, 3);
    return p;
  }
  
  void addTerms (PPoly &p, int i, int j, int k) {
    p.add(pts[i]->get().tripleProduct(pts[j]->get(), pts[k]->get()), 0);
    p.add(dsp[i]->get().tripleProduct(pts[j]->get(), pts[k]->get()), 1);
    p.add(pts[i]->get().tripleProduct(dsp[j]->get(), pts[k]->get()), 1);
    p.add(pts[i]->get().tripleProduct(pts[j]->get(), dsp[k]->get()), 1);
    p.add(pts[i]->get().tripleProduct(dsp[j]->get(), dsp[k]->get()), 2);
    p.add(dsp[i]->get().tripleProduct(pts[j]->get(), dsp[k]->get()), 2);
    p.add(dsp[i]->get().tripleProduct(dsp[j]->get(), pts[k]->get()), 2);
    p.add(dsp[i]->get().tripleProduct(dsp[j]->get(), dsp[k]->get()), 3);
  }

 public:
  FIPoly (const Points &pts, const Points &dsp) : pts(pts), dsp(dsp) {}
};

class DisplacedPoint : public Point {
  PTR<Point> a, u;
  Object<Parameter> *s;

  PV3 calculate () {
    return a->get() + s->get()*u->get();
  }
 public:
  DisplacedPoint (Point *a, Point *u, Object<Parameter> *s) : a(a), u(u), s(s) {}
};

bool intersectsVF (const Points &pts);

bool intersectsEE (const Points &pts);

Polyhedron * round (Polyhedron *a);

#endif
