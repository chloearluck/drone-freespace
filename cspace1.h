#include "mink.h"

#define PI 3.141592653589793

int quadraticRoots (Parameter a, Parameter b, Parameter c, Parameter *x);

int quarticRoots (Parameter a4, Parameter a3, Parameter a2, Parameter a1,
		  Parameter a0, Parameter *x);

VertexSet intersection (const VertexSet &a, const VertexSet &b);

class Angle : public Object<Parameter> {
 public:
  static Angle *mpi, *ppi;
  virtual VertexSet getVS () const { return VertexSet(); }
  Parameter getT () { return get(); }
  
  bool piAngle () const { return this == mpi || this == ppi; }

  Parameter cos () {
    if (this == mpi || this == ppi)
      return Parameter::constant(-1.0);
    Parameter t = get();
    return (1.0 - t*t)/(1.0 + t*t);
  }
  
  Parameter sin () {
    if (this == mpi || this == ppi)
      return Parameter::constant(0.0);
    Parameter t = get();
    return 2.0*t/(1.0 + t*t);
  }

  PV3 rotateZ (const PV3 &p) {
    Parameter c = cos(), s = sin();
    return PV3(c*p.x - s*p.y, s*p.x + c*p.y, p.z);
  }

  PV2 rotate (const PV2 &p) {
    Parameter c = cos(), s = sin();
    return PV2(c*p.x - s*p.y, s*p.x + c*p.y);
  }

  PV2 rotateDt (const PV2 &p) {
    Parameter t = get(), dt = 2.0/(1.0 + t*t);
    PV2 q = rotate(p);
    return PV2(- dt*q.y, dt*q.x);
  }

  int sharedVertices (Angle *a) const {
    return intersection(getVS(), a->getVS()).size();
  }
  
  double theta () {
    if (this == mpi)
      return - PI;
    if (this == ppi)
      return PI;
    return 2.0*atan(get().mid());
  }

  double lb () {
    if (this == mpi)
      return -1e308;
    if (this == ppi)
      return 1e308;
    return get().lb();
  }

  double ub () {
    if (this == mpi)
      return -1e308;
    if (this == ppi)
      return 1e308;
    return get().ub();  
  }
};

typedef vector<PTR<Angle> > Angles;

Primitive2(AngleOrder, Angle *, a, Angle *, b);

bool angleOrder (Angle *a, Angle *b);

bool inInterval (Angle *a, Angle *s, Angle *e);

bool intervalOverlap (Angle *s1, Angle *e1, Angle *s2, Angle *e2);

typedef pair<Angle *, Angle *> AngleInterval;

typedef vector<AngleInterval> AngleIntervals;

AngleIntervals intersect (const AngleIntervals &a, const AngleIntervals &b);

class InputAngle : public Angle {
 public:
  InputAngle (double t, bool flag = true) {
    double tt = tan(0.5*t);
    set(flag ? Parameter::input(tt) : Parameter::constant(tt));
  }
};

int sinCosRoots (Parameter k1, Parameter k2, Parameter k3, Parameter *x);

// k1*s + k2*c + k3 = 0
class AnglesCS : public Object<PV2> {
  friend class AngleCS;
  friend class Cspace;
 protected:
  bool roots, min;
  PTR<Angle> a[2];
  VertexSet vs;
  virtual void coeffs (Parameter &k1, Parameter &k2, Parameter &k3) = 0;
  
 Parameter getP (bool flag) { return flag ? get().x : get().y; }

 PV2 calculate () {
   Parameter k1, k2, k3, x[2];
   coeffs(k1, k2, k3);
   min = (k3 - k2).sign() == 1;
   int n = sinCosRoots(k1, k2, k3, x);
   if (n == 0) {
     roots = false;
     return PV2();
   }
   roots = true;
   return PV2(x[0], x[1]);
 }
 
 public:
 AnglesCS () { a[0] = a[1] = 0; }
 Angle * getA (bool flag);
};

class AngleCS : public Angle {
  AnglesCS *a;
  bool flag;
  Parameter calculate () { return a->getP(flag); }
 public:
  AngleCS (AnglesCS *a, bool flag) : a(a), flag(flag) {}
  VertexSet getVS () const { return a->vs; }
};

class EFKey {
 public:
  Vertex *v[5];
  EFKey (Edge *e, Face *f) {
    v[0] = e->getT();
    v[1] = e->getH();
    sort(v, v + 2);
    HEdge *b = f->getBoundary(0);
    v[2] = b->tail();
    v[3] = b->head();
    v[4] = b->getNext()->head();
    sort(v + 2, v + 5);
  }

  EFKey (Vertex *v1, Vertex *v2, Vertex *v3, Vertex *v4, Vertex *v5) {
    v[0] = v1; v[1] = v2; v[2] = v3; v[3] = v4; v[4] = v5;
    sort(v, v + 2);
    sort(v + 2, v + 5);
  }

  bool operator< (const EFKey &x) const {
    for (int i = 0; i < 5; ++i)
      if (v[i] < x.v[i])
	return true;
      else if (x.v[i] < v[i])
	return false;
    return false;
  }
};  

class AnglesEF : public AnglesCS {
  Edge *e;
  Face *f;
  bool aflag;
  Vertex **v;
  void coeffs (Parameter &k1, Parameter &k2, Parameter &k3) {
    PV3 u, n;
    if (e) {
      u =  e->getU();
      n = f->getP()->getN();
    }
    else {
      PV3 p[5];
      for (int i = 0; i < 5; ++i)
	p[i] = v[i]->getP()->getP();
      u = p[1] - p[0];
      n = (p[4] - p[3]).cross(p[2] - p[3]);
    }
    k1 = aflag ? u.x*n.y - u.y*n.x : n.x*u.y - n.y*u.x;
    k2 = u.x*n.x + u.y*n.y;
    k3 = u.z*n.z;
  }
 public:
  AnglesEF (Edge *e, Face *f, bool aflag) : e(e), f(f), aflag(aflag), v(0) {
    vs.insert(e->getT()); vs.insert(e->getH());
    HEdge *b = f->getBoundary(0);
    vs.insert(b->tail());
    vs.insert(b->head());
    vs.insert(b->getNext()->head());
    get();
  }

  AnglesEF (const EFKey &k, bool aflag)
    : e(0), f(0), aflag(aflag), v(new Vertex * [5]) {
    for (int i = 0; i < 5; ++i)
      v[i] = k.v[i];
    vs.insert(k.v, k.v + 5);
    get();
  }

  ~AnglesEF () { delete [] v; }

};

class CPoint : public Object<PV2> {
  friend class CFace;
  friend class CEdge;
  friend class Cspace;
 protected:
  PTR<Angle> a;
 public:
  CPoint (PTR<Angle> a) : a(a) {}
  PTR<Angle> getA () const { return a; }
  PV2 getP () { return get(); }
  void getBBox (double *bbox) {
    PV2 xy = get();
    bbox[0] = xy.x.lb();
    bbox[1] = xy.x.ub();
    bbox[2] = xy.y.lb();
    bbox[3] = xy.y.ub();
    bbox[4] = a->lb();
    bbox[5] = a->ub();
  }
};

typedef vector<CPoint *> CPoints;

Primitive3(COrder, CPoint *, a, CPoint *, t, CPoint *, h);

bool onEdge (CPoint *a, CPoint *t, CPoint *h);

Primitive2(CPointOrderX, CPoint *, a, CPoint *, b);

class CPointOrderXO {
 public:
  bool operator() (CPoint *a, CPoint *b) const {
    return a != b && CPointOrderX(a, b) == 1;
  }
};

class CEdge;

class Spiral {
 public:
  Vertex *v;
  Edge *e;
  bool aflag;
  VertexSet vs;
  vector<CEdge *> edges;

  Spiral (Vertex *v, Edge *e, bool aflag) : v(v), e(e), aflag(aflag) {
    vs.insert(v);
    vs.insert(e->getT());
    vs.insert(e->getH());
  }

  void vertices (VertexSet &vsa, VertexSet &vsb) const {
    if (aflag) {
      vsa.insert(v);
      vsb.insert(e->getT());
      vsb.insert(e->getH());
    }
    else {
      vsa.insert(e->getT());
      vsa.insert(e->getH());
      vsb.insert(v);
    }
  }

  int sharedVertices (Angle *a) const { return intersection(vs, a->getVS()).size(); }

  void coeffs (PV2 &a, PV2 &b) const;
  void bboxTP (Angle *as, Angle *ae, double *bbox);

  PV2 xy (Angle *a) {
    PV2 u, v;
    coeffs(u, v);
    return a->rotate(u) + v;
  }

  PV2 dt (Angle *a) {
    PV2 u, v;
    coeffs(u, v);
    return a->rotateDt(u);
  }
};

class SpiralOrder {
 public:
  bool operator() (Spiral *a, Spiral *b) const {
    if (a == b) return false;
    return a->v < b->v || a->v == b->v && a->e < b->e;
  }
};

typedef set<Spiral *, SpiralOrder> SpiralSet;

class SpiralPoint : public CPoint {
  friend class CEdge;
  Spiral *s;
  PV2 calculate () { return s->xy(a); }
 public:
  SpiralPoint (PTR<Angle> a, Spiral *s) : CPoint(a), s(s) {}
  Spiral * getS () const { return s; }
};

class CEdge;

class HCEdge;

class CVertex {
  friend class CEdge;
  friend class CHEdge;
  friend class CFace;
  friend class Cspace;
 protected:
  CPoint *p;
  vector<CEdge *> edges;
  double bbox[6];
 public:
  CVertex (CPoint *p) : p(p) { p->getBBox(bbox); }
  ~CVertex () { delete p; }
  CPoint * getP () { return p; }
  int EdgesN () const { return edges.size(); }
  CEdge * getEdge (int i) const { return edges[i]; }
  double * getBBox () { return bbox; }
  void outgoingHEdges (vector<HCEdge *> &ed) const;
};

typedef vector<CVertex *> CVertices;

typedef set<CVertex *> CVertexSet;

CVertexSet intersection (const CVertexSet &a, const CVertexSet &b);

class VertexAngleOrder {
 public:
  bool operator() (CVertex *a, CVertex *b) const {
    return a != b && angleOrder(a->getP()->getA(), b->getP()->getA());
  }
};

class CFace;

class CEdge {
  friend class CPoint;
  friend class CVertex;
  friend class HCEdge;
  friend class CFace;
  friend class Cspace;
 protected:
  CVertex *t, *h;
  bool bflag;
  Spiral *s;
  double bbox[6];
  vector<HCEdge *> hedges;
  CVertices vertices;
 public:
  CEdge (CVertex *t, CVertex *h, bool bflag);
  ~CEdge ();
  CVertex * getT () const { return t; }
  CVertex * getH () const { return h; }
  Spiral * getSpiral () const { return s; }
  int HEdgesN () const { return hedges.size(); }
  HCEdge * getHEdge (int i) const { return hedges[i]; }
  double * getBBox () { return bbox; }
  void addVertex (CVertex *v);
  HCEdge * addHEdge (bool forward);
  void removeHEdge (HCEdge *e);
  void edgeVertices (CVertices &ve) const;
  CFace * otherFace (CFace *f) const;
  bool piEdge () const;
  bool horizontal () const { return t->p->a == h->p->a; }
  bool increasing () const { return angleOrder(t->p->a, h->p->a); }
  bool decreasing () const { return angleOrder(h->p->a, t->p->a); }
  AngleInterval angleInterval () const;
  bool contains (Angle *a) const { return inInterval(a, t->p->a, h->p->a); }
  bool contains (CPoint *p) const { return onEdge(p, t->p, h->p); }
  PV2 xy (Angle *a) const;
  PV3 getU (CPoint *p) const;
  void sortHEdges ();
};

typedef vector<CEdge *> CEdges;

typedef set<CEdge *> CEdgeSet;

Primitive3(CEdgeOrderP, CEdge *, e, HCEdge *, f, HCEdge *, g);

class CEdgeOrder {
 public:
  CEdgeOrder (CEdge *e) : e(e) {}
  bool operator() (HCEdge *f, HCEdge *g) const {
    return f != g && CEdgeOrderP(e, f, g) == 1;
  }

  CEdge *e;
};

class CPointEdge : public CPoint {
  CEdge *e;
  PV2 calculate () {
    PV2 p = e->xy(a);
    if ((p.x.intervalWidth() > .001 || p.y.intervalWidth() > .001) &&
	Parameter::highPrecision == 53u) {
      Parameter::highPrecision = 212u;
      throw signException;
    }
    return p;
  }
  
 public:
  CPointEdge (CEdge *e, Angle *a) : e(e), CPoint(a) {}
};

class HCEdge {
  friend class CEdge;
  friend class CFace;
  friend class Cspace;
 protected:
  CEdge *e;
  bool forward, flag;
  CFace *f;
  HCEdge *next, *hp;
 public:
  HCEdge (CEdge *e, bool forward) 
    : e(e), forward(forward), flag(false), next(0), hp(0), f(0) {}
  CVertex * tail () const { return forward ? e->t : e->h; }
  CVertex * head () const { return forward ? e->h : e->t; }
  CEdge * getE () const { return e; }
  bool getForward () const { return forward; }
  CFace * getF () const { return f; }
  HCEdge * getNext () const { return next; }
  void setNext (HCEdge *h) { next = h; }
  HCEdge * getHP () const { return hp; }
  HCEdge * cw () const;
  HCEdge * ccw () const;
  void loop (vector<HCEdge *> &ed);
  bool onLoop ();
  bool increasing () const { return forward ? e->increasing() : e->decreasing(); }
  bool decreasing () const { return forward ? e->decreasing() : e->increasing(); }
  PV3 getU (CPoint *p) const;
  PV3 getN (CPoint *p) const;
};

typedef vector<HCEdge *> HCEdges;

typedef set<HCEdge *> HCEdgeSet;

enum HHCEdgeType { RightEnd, LeftEnd, LeftStart, RightStart};

class HHCEdge {
 public:
  HCEdge *e;
  bool f;
  HHCEdge (HCEdge *e, bool f) : e(e), f(f) {}
  CVertex * tail () const { return f ? e->tail() : e->head(); }
  CVertex * head () const { return f ? e->head() : e->tail(); }
  
  HHCEdgeType type () const {
    if (f)
      return e->increasing() ? RightStart : LeftEnd;
    return e->increasing() ? RightEnd : LeftStart;
  }
  
  bool operator< (const HHCEdge &x) const {
    Angle *a = tail()->getP()->getA(), *xa = x.tail()->getP()->getA();
    return a == xa ? type() < x.type() : angleOrder(a, xa);
  }
};

typedef vector<HHCEdge> HHCEdges;

class Slab {
 public:
  CEdge *l, *r;
  Angle *s, *e;
  CPoints pb, pt;
  CFace *f;

  Slab (CEdge *l, CEdge *r, Angle *s, Angle *e, CFace *f)
    : l(l), r(r), s(s), e(e), f(f) {}
  void addPoint (CPoint *p);
};

typedef vector<Slab *> Slabs;

Primitive2(InSlab, CPoint *, p, Slab *, s);

class CShell;

enum CFaceType { CFaceVF, CFaceFV, CFaceEEP, CFaceEEM };

class CFace {
  friend class AnglesSF;
  friend class HCEdge;
  friend class CShell;
  friend class Cspace;
 protected:
  CFaceType type;
  int *feature1, *feature2;
  HCEdges boundary;
  double bbox[6];
  CEdges edges;
  VertexSet vs;
  CShell *s;
  Slabs slabs;
 public:
  CFace (Vertex *v, Face *f, bool aflag);
  CFace (Edge *e, Edge *f, bool pflag);
  CFace (CFaceType type, int *feature1, int *feature2)
    : type(type), feature1(feature1), feature2(feature2), s(0) {}
  ~CFace () {
    for (Slabs::iterator s = slabs.begin(); s != slabs.end(); ++s)
      delete *s;
  }
  CFaceType getType () const { return type; }
  int * getFeature1 () const { return feature1; }
  int * getFeature2 () const { return feature2; }
  HCEdge * getBoundary (int i) const { return boundary[i]; }
  const HCEdges &getBoundary () { return boundary; } 
  double * getBBox () { return bbox; }
  const CEdges & getEdges () const { return edges; }
  const Slabs & getSlabs () const { return slabs; }
  Angle * startAngle () const { return boundary[0]->tail()->p->a; }
  Angle * endAngle () const { return boundary[0]->next->head()->p->a; }
  CEdge * leftEdge () const { return boundary[0]->next->next->next->e; }
  CEdge * rightEdge () const { return boundary[0]->next->e; }
  void addBoundary (HCEdge *e);
  void vertices (VertexSet &vsa, VertexSet &vsb) const;
  Spiral * sharedSpiral (CFace *f) const;
  void boundaryVertices (CVertices &ve) const;
  bool boundaryEdge (CEdge *e) const;
  void boundaryHEdges (HCEdges &ed) const;
  void sharedBoundaryVertices (CFace *f, CVertices &vfg) const;
  void edgeVertices (CVertexSet &vs) const;
  int sharedVertices (Angle *a) const { return intersection(vs, a->getVS()).size(); }
  bool sameBoundary (CVertex *a, CVertex *b) const;
  bool sharedBoundaryVertex (CFace *f, CFace *g) const;
  bool bboxOverlap (CFace *f, CFace *g) const;
  bool angleOverlap (CFace *f, CFace *g) const;
  bool inInterval (Angle *a) const;
  bool contains (CPoint *p) const;
  bool contains2(CPoint *p) const;
  PV3 getN (CPoint *p) const;
  void coeffs (PV2 &a1, PV2 &a2, Parameter &k1, Parameter &k2,
	       Parameter &k3) const;
  void coeffsVF (PV2 &a1, PV2 &a2, Parameter &k1, Parameter &k2,
		 Parameter &k3) const;
  void coeffsFV (PV2 &a1, PV2 &a2, Parameter &k1, Parameter &k2,
		 Parameter &k3) const;
  void coeffsEE (PV2 &a1, PV2 &a2, Parameter &k1, Parameter &k2,
		 Parameter &k3) const;
  void line (Angle *a, Parameter &u, Parameter &v, Parameter &w) const;
  PV2 intersectionXY (CFace *f, Angle *a) const;
  CFace * neighbor (HCEdge *e) const;
  void formSlabs();
  void initSlabs (HHCEdges &hhe, CPoints &hpts) const;
  void updateRE (Slabs &sl, const HHCEdge &h);
  void updateLE (Slabs &sl, const HHCEdge &h);
  void updateLS (Slabs &sl, const HHCEdge &h);
  void updateRS (Slabs &sl, const HHCEdge &h);
  void splitSlab (Slab *s, const Angles &an, Slabs &sl);
};

typedef vector<CFace *> CFaces;

Primitive3(IndependentFaces, CFace *, f, CFace *, g, Angle *, a);

typedef pair<CEdge *, CEdge *> CEEPair;

typedef set<CEEPair> CEEPairSet;

typedef pair<EFKey, AnglesCS *> EFAPair;

typedef map<EFKey, AnglesCS *> EFAMap;

typedef pair<Angle *, Spiral *> ASPair;

typedef pair<ASPair, CVertex *> ASVPair;

typedef map<ASPair, CVertex *> ASVMap;

typedef pair<CEdge *, CFace *> CEFPair;

typedef pair<CEFPair, CVertices> CEFVPair;

typedef map<CEFPair, CVertices> CEFVMap;

typedef pair<Spiral *, CFace *> SFPair;

typedef pair<Angle *, CFace *> AFPair;

typedef pair<AFPair, CVertex *> AFVPair;

typedef map<AFPair, CVertex *> AFVMap;

typedef pair<SFPair, AnglesCS *> SFAPair;

typedef map<SFPair, AnglesCS *> SFAMap;

class LinePatchPoint : public CPoint {
 protected:
  CEdge *e;
  CFace *f;
  PV2 calculate () {
    CPoint *et = e->getT()->getP(), *eh = e->getH()->getP();
    Spiral *l = f->leftEdge()->getSpiral(), *r = f->rightEdge()->getSpiral();
    PV2 pl = l->xy(a), v = r->xy(a) - pl, u = eh->getP() - et->getP(),
      w = et->getP() - pl;
    Parameter k = u.cross(w)/u.cross(v);
    return pl + k*v;    
  }
 public:
  LinePatchPoint (CEdge *e, CFace *f)
    : CPoint(e->getT()->getP()->getA()), e(e), f(f) {}
};

class AnglesSF : public AnglesCS {
 protected:
  Spiral *s;
  CFace *f;
  void coeffs (Parameter &k1, Parameter &k2, Parameter &k3) {
    PV2 a1, b1, a2, b2, a3, b3;
    s->coeffs(a1, b1);
    f->leftEdge()->getSpiral()->coeffs(a2, b2);
    f->rightEdge()->getSpiral()->coeffs(a3, b3);
    PV2 a12 = a1 - a2, b12 = b1 - b2, a23 = a2 - a3, b23 = b2 - b3;
    k1 = a23.dot(b12) - a12.dot(b23);
    k2 = a12.cross(b23) - a23.cross(b12);
    k3 = a12.cross(a23) + b12.cross(b23);
  }
 public:
  AnglesSF (Spiral *s, CFace *f) : s(s), f(f) { get(); }
};

typedef pair<CFace *, CFace *> CFFPair;

typedef pair<CFFPair, CEdges> CFFEPair;

typedef map<CFFPair, CEdges> CFFEMap;

class FFFData {
 public:
  Parameter p[4];
  int n;
  FFFData () : n(0) {}
  int size () const { return n; }
  Parameter & operator[](int i) { return p[i]; }
};

class AnglesFFF : public Object<FFFData>
{
  friend class AngleFFF;
  friend class Cspace;
  CFace *f1, *f2, *f3;
  int n;
  PTR<Angle> a[4];
  Parameter getP (int i) { return get().p[i]; }
  FFFData calculate ();
  void coeffs (Parameter &k1, Parameter &k2, Parameter &k3, Parameter &k4,
	       Parameter &k5, Parameter &k6);
  void coeffs (Parameter &k1, Parameter &k2, Parameter &k3, PV2 &a1, PV2 &b1,
	       PV2 &a2, PV2 &b2, Parameter &x1, Parameter &x2, Parameter &x3,
	       Parameter &x4, Parameter &x5, Parameter &x6);
 public:
  AnglesFFF (CFace *f1, CFace *f2, CFace *f3) : f1(f1), f2(f2), f3(f3) {
    a[0] = a[1] = a[2] = a[3] = 0;
    get();
  }
  
  Angle * getA (int i);
};

class AngleFFF : public Angle {
  AnglesFFF *s;
  int i;
  Parameter calculate () { return s->getP(i); }
 public:
  AngleFFF (AnglesFFF *s, int i) : s(s), i(i) {}
};

int sinCosRoots (Parameter k1, Parameter k2, Parameter k3, Parameter k4,
		 Parameter k5, Parameter k6, Parameter *x);

class FFFPoint : public CPoint {
  CFace *f1, *f2, *f3;
  PV2 calculate () {
    PV2 p12 = f1->intersectionXY(f2, a);
    if (f3) {
      PV2 p13 = f1->intersectionXY(f3, a);
      if (p13.dot(p13).intervalWidth() < p12.dot(p12).intervalWidth())
	return p13;
    }
    return p12;
  }
  
 public:
  FFFPoint (PTR<Angle> a, CFace *f1, CFace *f2, CFace *f3)
    : CPoint(a), f1(f1), f2(f2), f3(f3) {}
};    

typedef pair<CVertex *, CVertex *> CVVPair;

typedef map<CVertex *, CVertex *> CVVMap;

class CCell;

class CShell {
 public:
  CFaces faces;
  CVertex *vm;
  bool outer;
  CCell *c;

  CShell () : vm(0), outer(0), c(0) {}
  void init ();
};

typedef vector<CShell *> CShells;
  
class CShellOrder {
 public:
  bool operator() (CShell *a, CShell *b) const {
    return a != b && CPointOrderX(a->vm->getP(), b->vm->getP()) == 1;
  }
};

class CCell {
 public:
  CShells shells;
  ~CCell ();
};

typedef vector<CCell *> CCells;

typedef pair<CPoint *, Vertex *> CPVPair;

typedef map<CPoint *, Vertex *> CPVMap;

class Cspace {
 public:
  CVertices vertices;
  CEdges edges;
  CFaces faces;
  CCells cells;
  double bbox[6];
  EFAMap efamap;
  SpiralSet ss;
  vector<AnglesFFF *> fff;
  ASVMap asvmap;
  CEFVMap efvmap;
  AFVMap afvmap;
  SFAMap sfamap;
  Cspace *par;

  Cspace (Cspace *par = 0) : par(par) {}
  ~Cspace ();
  AnglesCS * angles (Edge *e, Face *f, bool aflag);
  AnglesCS * angles (const EFKey &k, bool aflag);
  void angleIntervals (Edge *e, Face *f, bool aflag, AngleIntervals &pos,
		       AngleIntervals &neg);
  void angleIntervals (HEdge *e, Face *f, bool aflag, AngleIntervals &pos,
		       AngleIntervals &neg);
  Spiral * getSpiral (Vertex *v, Edge *e, bool aflag);
  CVertex * getVertex (CPoint *p);
  CVertex * getVertex (Angle *a, Spiral *s);
  HCEdge * addHEdge (CVertex *a, CVertex *b, bool bflag, HCEdge *hp = 0);
  CEdge * getEdge (CVertex *a, CVertex *b, bool bflag, HCEdge *hp = 0);
  HCEdge * addLoop (const CVertices &ve);
  void patches (Polyhedron *a, Polyhedron *b);
  void patchesVF (Polyhedron *a, Polyhedron *b, bool aflag);
  void patchVF (Vertex *v, Face *f, bool aflag);
  void patchVF (Vertex *v, Face *f, bool aflag, Angle *s, Angle *e);
  void patch (Angle *s, Angle *e, Spiral *l, Spiral *r, CFace *f);
  void patchesEE (Polyhedron *a, Polyhedron *b);
  void patchEE (Edge *e, Edge *f);
  void patchEE (Edge *e1, Edge *e2, Angle *s, Angle *e, bool flag);
  void intersectFF ();
  void intersectEE (CFace *f, CFace *g, CEEPairSet &eps);
  void intersectEE (CEdge *e, CEdge *f);
  void intersectSS (CEdge *e, CEdge *f) const;
  void intersectSL (CEdge *e, CEdge *f);
  void intersectLL (CEdge *e, CEdge *f);
  void intersectFF (CFace *f, CFace *g);
  void intersectFF (CFace *f, CFace *g, CVertices &vfg);
  void intersectEF (HCEdge *e, CFace *f, CVertices &vfg);
  void intersectEF (CEdge *e, CFace *f, CVertices &ve);
  CVertex * intersectEFLine (CEdge *e, CFace *f);
  AnglesCS * intersectSF (Spiral *s, CFace *f);
  AnglesCS * intersectSFLine (Spiral *s, CFace *f);
  void formFF (CFace *f, CFace *g, CVertices &vfg);
  bool intersectsFFLine (CFace *f, CFace *g) const;
  void formFFLine (CFace *f, CFace *g, CVertices &vfg);
  CVertex * spiralVertex (Angle *a, Spiral *s);
  void formFFCurve (CFace *f, CFace *g, CVertices &vfg);
  void formFF (CFace *f, CFace *g, CVertex *v, CVertex *w);
  void intersectFFF ();
  void formFFEMap (CFFEMap &ffemap) const;
  void intersectFFF (CFace *f, const CFFEMap &ffemap);
  void intersectFFF (CFace *f1, CFace *f2, CFace *f3, CEdge *e12, CEdge *e13,
		     const CFFEMap &ffemap);
  void intersectFFFH (CFace *f1, CFace *f2, CFace *f3, CEdge *eh,
		      CEdge *e1, CEdge *e2);
  void intersectFFFG (CFace *f1, CFace *f2, CFace *f3, CEdge *e12, CEdge *e13,
		      const CEdges &ed);
  void sortVertices ();
  Cspace * subfaces ();
  void subfaces (CFace *f, Cspace *a, CVVMap &vvmap) const;
  void subedges (CFace *f, Cspace *a, CVVMap &vvmap, HCEdges &he) const;
  void subedges (HCEdge *e, Cspace *a, CVVMap &vvmap, HCEdges &he) const;
  CVertex * getVertex (CVertex *v, CVVMap &vmap);
  void setNext (CFace *f, const HCEdges &he) const;
  CFace * addFace (HCEdge *e, CFace *f);
  void removeBad (const HCEdges &ed);
  void formCells (Polyhedron *a, Polyhedron *b);
  void formShells (Polyhedron *a, Polyhedron *b, CShells &sh);
  void formShells (CShells &sh) const;
  CEdge * maxAngleEdge (CShell *s) const;
  CShell * enclosingShell (CShell *s, Octree<CFace *> *octree);
  void removeFace (CFace *f) const;
  void describe () const;
};

Primitive2(PointOrderZ, Point *, a, Point *, b);

bool inIntervalZ (Vertex *v, Edge *e);

class OrientationVF : public Primitive {
  bool aflag;
  Face *f;
  Angle *a;
  Spiral *l, *r;
  int sign () {
    PV3 n3 = aflag ? f->getP()->getN() : - a->rotateZ(f->getP()->getN());
    PV2 n(n3.x, n3.y), t = r->xy(a) - l->xy(a);
    return n.cross(t).sign();
  }
 public:
  OrientationVF (bool aflag, Face *f, Angle *a, Spiral *l, Spiral *r)
    : aflag(aflag), f(f), a(a), l(l), r(r) {}
};

class OrientationEE : public Primitive {
  bool flag;
  Edge *e, *f;
  Angle *a;
  Spiral *l, *r;
  int sign () {
    PV3 n3 = a->rotateZ(e->getU()).cross(f->getU());
    PV2 n(n3.x, n3.y), t = r->xy(a) - l->xy(a);
    return flag ? t.cross(n).sign() : n.cross(t).sign();
  }
 public:
  OrientationEE (bool flag, Edge *e, Edge *f, Angle *a, Spiral *l, Spiral *r)
    : flag(flag), e(e), f(f), a(a), l(l), r(r) {}
};

Primitive4(CFFOrder, CFace *, f, CFace *, g, CVertex *, v, CVertex *, w);

Primitive4(EdgeOrderL, CPoint *, t, CPoint *, h, CPoint *, v, CPoint *, w);

class EdgeVertexOrderL {
 public:
  EdgeVertexOrderL (CEdge *e) : t(e->getT()), h(e->getH()) {}
  bool operator() (CVertex *v, CVertex *w) const {
    return v != w && EdgeOrderL(t->getP(), h->getP(), v->getP(), w->getP()) == 1;
  }

  CVertex *t, *h;
};

class EdgeVertexOrderS {
 public:
  EdgeVertexOrderS (CEdge *e) : inc(e->increasing()) {}
  bool operator() (CVertex *v, CVertex *w) const {
    return v != w && inc == angleOrder(v->getP()->getA(), w->getP()->getA());
  }
  bool inc;
};

Primitive4(CVertexHHEdgeOrder, CVertex *, v, HCEdge *, e, HCEdge *, f, CFace *, g);

class HHCEdgeOrder {
 public:
  HHCEdgeOrder (CFace *g) : g(g) {}
  bool operator() (const HHCEdge &e, const HHCEdge &f) const {
    CVertex *u = e.tail(), *v = f.tail();
    return u != v ? u < v : CVertexHHEdgeOrder(u, e.e, f.e, g) == 1;
  }

  CFace *g;
};

class RayPatchPointX : public CPoint {
 protected:
  CPoint *p;
  CFace *f;
  PV2 calculate () {
    HCEdge *e = f->getBoundary(0);
    CFace *g = e->getHP()->getF();

    Spiral *l = g->leftEdge()->getSpiral(), *r = g->rightEdge()->getSpiral();
    PV2 pl = l->xy(a), v = r->xy(a) - pl, w = p->getP() - pl;
    Parameter k = w.y/v.y;
    return pl + k*v;
  }
 public:
  RayPatchPointX (CPoint *p, CFace *f) : CPoint(p->getA()), p(p), f(f) {}
};

Cspace * cspace (Polyhedron *a, Polyhedron *b);

class TransformedPoint : public Point {
  CPoint *c;
  Point *a;
  
  PV3 calculate () {
    PV3 p = a->getP();
    PV2 q = c->getA()->rotate(PV2(p.x, p.y)) + c->getP();
    return PV3(q.x, q.y, p.z);
  }

 public:
  TransformedPoint (CPoint *c, Point *a) : c(c), a(a) {}
};

Polyhedron * transform (Polyhedron *p, CPoint *c);

class EdgePointOrderL {
 public:
  EdgePointOrderL (CPoint *t, CPoint *h) : t(t), h(h) {}
  bool operator() (CPoint *v, CPoint *w) const {
    return v != w && EdgeOrderL(t, h, v, w) == 1;
  }

  CPoint *t, *h;
};

typedef pair<CEdge *, CPoints> EPPair;

typedef map<CEdge *, CPoints> EPMap;

typedef pair<Face *, Slab *> FSPair;

typedef map<Face *, Slab *> FSMap;

typedef pair<Slab *, Faces> SFSPair;

typedef map<Slab *, Faces> SFSMap;

Polyhedron * discretize (Cspace *c, double d);

Polyhedron * discretize (const CFaces &fa, double d, CPVMap &pvmap);

EPMap discretizeEdges (const CFaces &fa, double d);

CPoints discretize (CEdge *e, double d);

void discretize (CEdge *e, Angle *as, Angle *ae, double d, CPoints &pts);

bool close (CPoint *a, CPoint *t, CPoint *h, double d);

class DistancePL : public Primitive {
  Point *a, *t, *h;
  double d;
  int sign () {
    PV3 tp = t->getP(), u = h->getP() - tp,
      p = tp + ((a->getP() - tp).dot(u)/u.dot(u))*u, w = a->getP() - p;
    return (d*d - w.dot(w)).sign();
  }
 public:
  DistancePL (Point *a, Point *t, Point *h, double d) : a(a), t(t), h(h), d(d) {}
};

void delentil (const CFaces &fa, EPMap &epmap);

void delentil (Slab *s, EPMap &epmap);

class CTriangle {
 public:
  CPoint *a, *b, *c;
  CTriangle (CPoint *a, CPoint *b, CPoint *c) : a(a), b(b), c(c) {}
};

typedef vector<CTriangle> CTriangles;

void discretize (Slab *s, CPVMap &pvmap, const EPMap &epmap, Polyhedron *a,
		 FSMap &fsmap, SFSMap &sfmap);

CTriangles discretize (Slab *s, const EPMap &epmap);

CPoints getPoints (CEdge *f, Angle *s, Angle *e, const EPMap &epmap);

CTriangles discretize (const CPoints &l, const CPoints &r, Slab *s);

void ctriangle (CPoint *a, CPoint *b, CPoint *c, CTriangles &tr);

void splitB (const CTriangle &t, const CPoints &pb, bool lflag, CTriangles &tr);

void splitT (const CTriangle &t, const CPoints &pt, bool lflag, CTriangles &tr);

Vertex * getVertex (CPoint *p, CPVMap &cpvmap, Polyhedron *a);

class PointCPoint : public Point {
  CPoint *p;
  bool flag;
  PV3 calculate () {
    PV2 xy = p->getP();
    PTR<Angle> a = p->getA();
    Parameter th;
    if (a == Angle::mpi)
      th = - Parameter::constant(PI);
    else if (a == Angle::ppi)
      th = Parameter::constant(PI);
    else
      th = a->getT().sign()*a->cos().acos();
    return PV3(xy.x, xy.y, th);
  }
 public:
  PointCPoint (CPoint *p, bool flag = false) : p(p), flag(flag) {}
  ~PointCPoint () { if (flag) delete p; }
  CPoint * getCP () const { return p; }
};

void removeIntersections (Polyhedron *a, CPVMap &pvmap, EPMap &epmap,
			  FSMap &fsmap, SFSMap &sfmap);

void splitSlabs (Face *f, Face *g, CPVMap &pvmap, EPMap &epmap, Polyhedron *a,
		 FSMap &fsmap, SFSMap &sfmap, FOctree *octree);

void splitSlab (Slab *s, const Angles &as, CPVMap &pvmap, EPMap &epmap,
		Polyhedron *a, FSMap &fsmap, SFSMap &sfmap, FOctree *octree);

void neighbors (Slab *s, CEdge *e, Slabs &sl);

void remove (Slab *s, Polyhedron *a, FSMap &fsmap, SFSMap &sfmap);

void addPoint (CEdge *e, Angle *a, EPMap &epmap);

void rediscretize (Slab *s, CPVMap &pvmap, const EPMap &epmap, Polyhedron *a,
		   FSMap &fsmap, SFSMap &sfmap, FOctree *octree);

void facesPI (Cspace *c, Polyhedron *a, CPVMap &pvmap, bool flag);

Vertices * loopPI (Polyhedron *a, CPVMap &pvmap, CEdgeSet &es, CEdge *e);

