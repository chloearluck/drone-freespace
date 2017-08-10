#include "cspace.h"

int quadraticRoots (Parameter a, Parameter b, Parameter c, Parameter *x)
{
  Parameter d = b*b - 4.0*a*c;
  if (d.sign() == -1)
    return 0;
  Parameter k = b.sign() == -1 ? d.sqrt() - b : - d.sqrt() - b,
    x1 = 0.5*k/a, x2 = 2.0*c/k;
  if ((x2 - x1).sign() == 1) {
    x[0] = x1;
    x[1] = x2;
  }
  else {
    x[0] = x2;
    x[1] = x1;
  }    
  return 2;
}

int quarticRoots (Parameter a, Parameter b, Parameter c, Parameter d, Parameter e,
		  Parameter *x)
{
  Parameter p = c/a - 3.0*b*b/(8.0*a*a),
    q = (b*b*b - 4.0*a*b*c + 8.0*a*a*d)/(8.0*a*a*a),
    delta0 = c*c - 3.0*b*d + 12.0*a*e,
    delta1 = 2.0*c*c*c - 9.0*b*c*d + 27.0*b*b*e + 27.0*a*d*d - 72.0*a*c*e,
    k1 = delta1*delta1 - 4.0*delta0*delta0*delta0, s;
  if (k1.sign() > -1) {
    if (k1.sign() == 0)
      return 0;
    Parameter kq = 0.5*(delta1 + k1.sqrt()),
      Q = kq.sign()*(kq.abs().root(3)),      
      k2 = -2.0*p/3.0 + (Q + delta0/Q)/(3.0*a);
    s = 0.5*k2.sqrt();
  }
  else {
    if (p.sign() == 1)
      return 0;
    Parameter D = 64.0*a*a*a*e - 16.0*a*a*c*c + 16.0*a*b*b*c - 16.0*a*a*b*d
      - 3.0*b*b*b*b;
    if (D.sign() == 1)
      return 0;
    Parameter ct = 0.5*delta1/(delta0*delta0*delta0).sqrt(), phi3 = ct.acos()/3.0,
      ks = -2.0*p/3.0 + 2.0*delta0.sqrt()*phi3.cos()/(3.0*a);
    s = ks.lb() > 0.0 ? 0.5*ks.sqrt() : Parameter::constant(0.0);
  }
  int n = 0;
  for (int t = -1; t < 2; t += 2) {
    Parameter k3 = -4.0*s*s - 2.0*p - t*q/s;
    if (k3.sign() == 1) {
      Parameter k4 = 0.5*k3.sqrt();
      for (int tt = -1; tt < 2; tt += 2)
	x[n++] = - b/(4.0*a) + t*s + tt*k4;
    }
  }
  return n;
}

VertexSet intersection (const VertexSet &a, const VertexSet &b)
{
  VertexSet c;
  set_intersection(a.begin(), a.end(), b.begin(), b.end(),
		   inserter(c, c.begin()));
  return c;
}

Angle *Angle::mpi = new InputAngle(0.0, false);
Angle *Angle::ppi = new InputAngle(0.0, false);

int AngleOrder::sign ()
{
  return (b->getT() - a->getT()).sign();
}

bool angleOrder (Angle *a, Angle *b)
{
  if (a == b || a == Angle::ppi || b == Angle::mpi)
    return false;
  return a == Angle::mpi || b == Angle::ppi || AngleOrder(a, b) == 1;
};

bool inInterval (Angle *a, Angle *s, Angle *e)
{
  return angleOrder(s, e) ? angleOrder(s, a) && angleOrder(a, e)
    : angleOrder(e, a) && angleOrder(a, s);
}

bool intervalOverlap (Angle *s1, Angle *e1, Angle *s2, Angle *e2)
{
  return angleOrder(s2, e1) && angleOrder(s1, e2);
}

void intersect (const AngleInterval &x, const AngleInterval &y,
		AngleIntervals &res)
{
  Angle *ns = angleOrder(x.first, y.first) ? y.first : x.first,
    *ne = angleOrder(x.second, y.second) ? x.second : y.second;
  if (angleOrder(ns, ne))
    res.push_back(AngleInterval(ns, ne));
}

AngleIntervals intersect (const AngleIntervals &a, const AngleIntervals &b)
{
  AngleIntervals res;
  for (AngleIntervals::const_iterator ai = a.begin(); ai != a.end(); ++ai)
    for (AngleIntervals::const_iterator bi = b.begin(); bi != b.end(); ++bi)
      intersect(*ai, *bi, res);
  return res;
}

int COrder::sign ()
{
  return (a->getP() - t->getP()).dot(h->getP() - t->getP()).sign();
}

bool onEdge (CPoint *a, CPoint *t, CPoint *h)
{
  if ( a == t || a == h)
    return false;
  Angle *ta = t->getA(), *ha = h->getA();
  return ta == ha ? COrder(a, t, h) == 1 && COrder(a, h, t) == 1
    : inInterval(a->getA(), ta, ha);
}

// k1*s + k2*c + k3 = 0
int sinCosRoots (Parameter k1, Parameter k2, Parameter k3, Parameter *x)
{
  return quadraticRoots(k3 - k2, 2.0*k1, k2 + k3, x);
}

Angle * AnglesCS::getA (bool flag)
{
  if (!a[0]) {
    a[0] = new AngleCS(this, true);
    a[1] = new AngleCS(this, false);
  }
  return flag ? a[0] : a[1];
}

// xy(th) = th*a + b
void Spiral::coeffs (PV2 &a, PV2 &b) const
{
  PV3 p = v->getP()->getP(), t = e->getT()->getP()->getP(),
    h = e->getH()->getP()->getP(), u = h - t, w = t.cross(u);
  if (aflag) {
    a.x = - p.x;
    a.y = - p.y;
    b.x = (p.z*u.x - w.y)/u.z;
    b.y = (p.z*u.y + w.x)/u.z;
  }
  else {
    a.x = (w.y - p.z*u.x)/u.z; 
    a.y = - (w.x + p.z*u.y)/u.z;
    b.x = p.x;
    b.y = p.y;
  }
}

void Spiral::bboxTP (Angle *as, Angle *ae, double *bbox)
{
  PV2 u, v;
  coeffs(u, v);
  double t = - atan2(u.y.mid(), u.x.mid());
  for (int i = 0; i < 4; ++i) {
    InputAngle ti(t, false);
    if (inInterval(&ti, as, ae)) {
      PV2 p = xy(&ti);
      bbox[0] = min(bbox[0], p.x.lb());
      bbox[1] = max(bbox[1], p.x.ub());
      bbox[2] = min(bbox[2], p.y.lb());
      bbox[3] = max(bbox[3], p.y.ub());
    }
    t += 0.5*PI;
    if (t > PI)
      t -= 2.0*PI;
  }
}

void CVertex::outgoingHEdges (HCEdges &ed) const
{
  for (CEdges::const_iterator e = edges.begin(); e != edges.end(); ++e)
    for (HCEdges::iterator f = (*e)->hedges.begin(); f != (*e)->hedges.end(); ++f)
      if ((*f)->tail() == this)
	ed.push_back(*f);
}

CVertexSet intersection (const CVertexSet &a, const CVertexSet &b)
{
  CVertexSet c;
  set_intersection(a.begin(), a.end(), b.begin(), b.end(),
		   inserter(c, c.begin()));
  return c;
}

CEdge::CEdge (CVertex *t, CVertex *h, bool bflag) : t(t), h(h), bflag(bflag), s(0)
{ 
  copyBBox(t->bbox, bbox);
  mergeBBox(h->bbox, bbox);
  SpiralPoint *tp = dynamic_cast<SpiralPoint *>(t->p),
    *hp = dynamic_cast<SpiralPoint *>(h->p);
  s = (tp && hp && tp->getS() == hp->getS()) ? tp->getS() : 0;
  if (s)
    if (angleOrder(tp->a, hp->a))
      s->bboxTP(tp->a, hp->a, bbox);
    else
      s->bboxTP(hp->a, tp->a, bbox);
}

CEdge::~CEdge ()
{
  for (HCEdges::iterator e = hedges.begin(); e != hedges.end(); ++e)
    delete *e;
}

void CEdge::addVertex (CVertex *v)
{
  if (find(vertices.begin(), vertices.end(), v) == vertices.end())
    vertices.push_back(v);
}

HCEdge * CEdge::addHEdge (bool forward)
{
  HCEdge *e = new HCEdge(this, forward);
  hedges.push_back(e);
  return e;
}

void CEdge::removeHEdge (HCEdge *e)
{
  HCEdges::iterator j = remove(hedges.begin(), hedges.end(), e);
  hedges.erase(j, hedges.end());
  delete e;
}

void CEdge::edgeVertices (CVertices &ve) const
{
  ve.push_back(t);
  ve.insert(ve.end(), vertices.begin(), vertices.end());
  ve.push_back(h);
}

CFace * CEdge::otherFace (CFace *f) const
{
  return hedges[0]->f == f ? hedges[1]->f : hedges[0]->f;
}

bool CEdge::piEdge () const
{
  return t->p->a == h->p->a && t->p->a->piAngle();
}

AngleInterval CEdge::angleInterval () const
{
  Angle *at = t->p->getA(), *ah = h->p->getA();
  if (angleOrder(at, ah))
    return AngleInterval(at, ah);
  return AngleInterval(ah, at);
}

PV2 CEdge::xy (Angle *a) const
{
  if (s)
    return s->xy(a);
  if (hedges[0]->hp)
    return hedges[0]->hp->e->xy(a);
  return hedges[0]->f->intersectionXY(hedges[1]->f, a);
}

PV3 CEdge::getU (CPoint *p) const
{
  if (horizontal()) {
    PV2 uxy = h->getP()->getP() - t->getP()->getP();
    return PV3(uxy.x, uxy.y, Parameter::constant(0.0));
  }
  if (s) {
    PV2 uxy = s->dt(p->a);
    Parameter uz = Parameter::constant(increasing() ? 1.0 : -1.0);
    return PV3(uz*uxy.x, uz*uxy.y, uz);
  }
  return hedges[1]->f->getN(p).cross(hedges[0]->f->getN(p));
}

void CEdge::sortHEdges ()
{
  if (hedges.size() > 2)
    sort(hedges.begin(), hedges.end(), CEdgeOrder(this));
}

int CEdgeOrderP::sign ()
{
  CPoint *p = e->getT()->getP();
  PV3 nf = f->getN(p), ng = g->getN(p), r = Rdir.getP();
  int yf = nf.dot(r).sign(), yg = ng.dot(r).sign();
  if (yf != yg)
    return yf;
  PV3 u = e->getU(p);
  return u.tripleProduct(ng, nf).sign();
}

HCEdge * HCEdge::cw () const
{
  int n = e->hedges.size();
  for (int i = 0; i < n; ++i)
    if (e->hedges[i] == this)
      return e->hedges[(i+1)%n];
}

HCEdge * HCEdge::ccw () const
{
  int n = e->hedges.size();
  for (int i = 0; i < n; ++i)
    if (e->hedges[i] == this)
      return e->hedges[i == 0 ? n - 1 : i - 1];
}

void HCEdge::loop (HCEdges &ed)
{
  HCEdge *h = this;
  do {
    ed.push_back(h);
    h = h->next;
  }
  while (h != this);
}

bool HCEdge::onLoop ()
{
  HCEdge *h = this;
  do {
    h->flag = true;
    h = h->next;
    if (!h)
      return false;
  }
  while (h != this);
  return true;
}

PV3 HCEdge::getU (CPoint *p) const
{
  if (hp)
    return hp->getU(p);
  return forward ? e->getU(p) : - e->getU(p);
}

PV3 HCEdge::getN (CPoint *p) const
{
  return forward ? f->getN(p) : - f->getN(p);
}

void Slab::addPoint (CPoint *p)
{
  if (p->getA() == s || p->getA() == e)
    if (InSlab(p, this) == 1)
      if (p->getA() == s)
	pb.push_back(p);
      else
	pt.push_back(p);
}

int InSlab::sign ()
{
  Angle *a = p->getA();
  if (p == s->l->getT()->getP() || p == s->l->getH()->getP() ||
      p == s->r->getT()->getP() || p == s->r->getH()->getP())
    return 0;
  PV2 pl = s->l->xy(a), pr = s->r->xy(a), u = pr - pl;
  Parameter k = u.dot(p->getP() - pl);
  int ks = k.sign();
  return ks < 1 ? ks : (u.dot(u) - k).sign();
}

CFace::CFace (Vertex *v, Face *f, bool aflag) : s(0)
{
  type = aflag ? CFaceVF : CFaceFV;
  feature1 = (int *) v;
  feature2 = (int *) f;
  HEdge *e = f->getBoundary(0);
  vs.insert(v);
  vs.insert(e->tail());
  vs.insert(e->head());
  vs.insert(e->getNext()->head());
}

CFace::CFace (Edge *e, Edge *f, bool pflag) : s(0)
{
  type = pflag ? CFaceEEP : CFaceEEM;
  feature1 = (int *) e;
  feature2 = (int *) f;
  vs.insert(e->getT());
  vs.insert(e->getH());
  vs.insert(f->getT());
  vs.insert(f->getH());
}

void CFace::addBoundary (HCEdge *e)
{
  HCEdges ed;
  e->loop(ed);
  if (boundary.empty()) {
    copyBBox(ed[0]->e->bbox, bbox);
    for (int i = 1; i < ed.size(); ++i)
      mergeBBox(ed[i]->e->bbox, bbox);
  }
  boundary.push_back(e);
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    (*e)->f = this;
}

void CFace::vertices (VertexSet &vsa, VertexSet &vsb) const
{
  if (type == CFaceVF) {
    vsa.insert((Vertex *) feature1);
    HEdge *b = ((Face *) feature2)->getBoundary(0);
    vsb.insert(b->tail());
    vsb.insert(b->head());
    vsb.insert(b->getNext()->head());
  }
  else if (type == CFaceFV) {
    vsb.insert((Vertex *) feature1);
    HEdge *b = ((Face *) feature2)->getBoundary(0);
    vsa.insert(b->tail());
    vsa.insert(b->head());
    vsa.insert(b->getNext()->head());
  }
  else {
    Edge *ea = (Edge *) feature1, *eb = (Edge *) feature2;
    vsa.insert(ea->getT());
    vsa.insert(ea->getH());
    vsb.insert(eb->getT());
    vsb.insert(eb->getH());
  }
}

Spiral * CFace::sharedSpiral (CFace *f) const
{
  Spiral *sl = leftEdge()->s, *sr = rightEdge()->s,
    *fsl = f->leftEdge()->s, *fsr = f->rightEdge()->s;
  if (sl == fsl || sl == fsr)
    return sl;
  if (sr == fsl || sr == fsr)
    return sr;
  return 0;
}

void CFace::boundaryVertices (CVertices &ve) const
{
  HCEdges ed;
  boundaryHEdges(ed);
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    ve.push_back((*e)->tail());
}

bool CFace::boundaryEdge (CEdge *e) const
{
  for (HCEdges::iterator f = e->hedges.begin(); f != e->hedges.end(); ++f)
    if ((*f)->f == this)
      return true;
  return false;
}

void CFace::boundaryHEdges (HCEdges &ed) const
{
  for (HCEdges::const_iterator b = boundary.begin(); b != boundary.end(); ++b)
    (*b)->loop(ed);
}

void CFace::sharedBoundaryVertices (CFace *f, CVertices &vfg) const
{
  CVertexSet vs1, vs2;
  edgeVertices(vs1);
  f->edgeVertices(vs2);
  CVertexSet vs = intersection(vs1, vs2);
  vfg.insert(vfg.end(), vs.begin(), vs.end());
}

void CFace::edgeVertices (CVertexSet &vs) const
{
  HCEdges ed;
  boundaryHEdges(ed);
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e) {
    CVertices ve;
    (*e)->e->edgeVertices(ve);
    vs.insert(ve.begin(), ve.end());
  }
}

bool CFace::sameBoundary (CVertex *a, CVertex *b) const
{
  HCEdges ed;
  boundaryHEdges(ed);
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e) {
    CVertices ve;
    (*e)->e->edgeVertices(ve);
    if (find(ve.begin(), ve.end(), a) != ve.end() &&
	find(ve.begin(), ve.end(), b) != ve.end())
      return true;
  }
  return false;
}

bool CFace::sharedBoundaryVertex (CFace *f, CFace *g) const
{
  CVertexSet vs1, vs2,
  edgeVertices(vs1);
  f->edgeVertices(vs2);
  CVertexSet vs12 = intersection(vs1, vs2);
  if (vs12.empty())
    return false;
  CVertexSet vs3;
  g->edgeVertices(vs3);
  CVertexSet vs123 = intersection(vs3, vs12);
  return !vs123.empty();
}

bool CFace::bboxOverlap (CFace *f, CFace *g) const
{
  double bb[6];
  for (int i = 0; i < 3; ++i) {
    bb[2*i] = max(bbox[2*i], f->bbox[2*i]);
    bb[2*i+1] = min(bbox[2*i+1], f->bbox[2*i+1]);
  }
  return ::bboxOverlap(bb, g->bbox);
}

bool CFace::angleOverlap (CFace *f, CFace *g) const
{
  Angle *s = startAngle(), *e = endAngle(), *fs = f->startAngle(),
    *fe = f->endAngle(), *gs = g->startAngle(), *ge = g->endAngle();
  if (angleOrder(s, fs))
    s = fs;
  if (angleOrder(s, gs))
    s = gs;
  if (angleOrder(fe, e))
    e = fe;
  if (angleOrder(ge, e))
    e = ge;
  return angleOrder(s, e);
}

bool CFace::inInterval (Angle *a) const
{
  return angleOrder(startAngle(), a) && angleOrder(a, endAngle());
}

// patch
bool CFace::contains (CPoint *p) const
{
  double bbp[6];
  p->getBBox(bbp);
  if (!::bboxOverlap(bbox, bbp))
    return false;
  CEdge *l = leftEdge(), *r = rightEdge();
  if (!inInterval(p->a))
    return false;
  SpiralPoint pl(p->a, l->getSpiral()), pr(p->a, r->getSpiral());
  return onEdge(p, &pl, &pr);
}

// general face
bool CFace::contains2 (CPoint *p) const
{
  HCEdges ed;
  boundaryHEdges(ed);
  CPoints pts;
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    if ((*e)->e->contains(p->a))
      pts.push_back(new CPointEdge((*e)->e, p->a));
  sort(pts.begin(), pts.end(), CPointOrderXO());
  bool res = false;
  for (int i = 0; !res && i + 1 < pts.size(); i += 2)
    res = onEdge(p, pts[i], pts[i+1]);
  for (CPoints::iterator i = pts.begin(); i != pts.end(); ++i)
    delete *i;
  return res;
}

PV3 CFace::getN (CPoint *p) const
{
  Angle *th = p->getA();
  PV2 d = p->getP();
  if (type == CFaceVF || type == CFaceFV) {
    Vertex *v = (Vertex *) feature1;
    Face *f = (Face *) feature2;
    PV3 vp = v->getP()->getP(), n = f->getP()->getN();
    PV2 v2(vp.x, vp.y), n2(n.x, n.y);
    if (type == CFaceVF)
      return PV3(n.x, n.y, th->rotateDt(v2).dot(n2));
    PV2 dd = - th->rotate(n2);
    Parameter dt = th->rotateDt(n2).dot(v2 - d);
    return PV3(dd.x, dd.y, dt);
  }
  Edge *e = (Edge *) feature1, *f = (Edge *) feature2;
  PV3 a = e->getT()->getP()->getP(), u = e->getU(),
    b = f->getT()->getP()->getP(), v = f->getU(),
    dd = th->rotateZ(u).cross(v), bd(b.x - d.x, b.y - d.y, b.z),
    bdv = bd.cross(v), ua = u.cross(a);
  PV2  dtu = th->rotateDt(PV2(u.x, u.y)), dtua = th->rotateDt(PV2(ua.x, ua.y));
  Parameter dt = dtu.dot(PV2(bdv.x, bdv.y)) - dtua.dot(PV2(v.x, v.y));
  PV3 df(dd.x, dd.y, dt);
  return type == CFaceEEP ? - df : df;
}

// ((th(a1) + a2).(x, y) + k1*s + k2*c + k3 = 0
void CFace::coeffs (PV2 &a1, PV2 &a2, Parameter &k1, Parameter &k2,
		    Parameter &k3) const
{
  if (type == CFaceVF)
    coeffsVF(a1, a2, k1, k2, k3);
  else if (type == CFaceFV)
    coeffsFV(a1, a2, k1, k2, k3);
  else
    coeffsEE(a1, a2, k1, k2, k3);
}

void CFace::coeffsVF (PV2 &a1, PV2 &a2, Parameter &k1, Parameter &k2,
		      Parameter &k3) const
{
  Vertex *v = (Vertex *) feature1;
  Face *f = (Face *) feature2;
  PV3 p = v->getP()->getP(), a = f->getBoundary(0)->tail()->getP()->getP(),
    n = f->getP()->getN();
  a1.x = a1.y = Parameter::constant(0.0);
  a2.x = n.x;
  a2.y = n.y;
  k1 = p.x*n.y - p.y*n.x;
  k2 = p.x*n.x + p.y*n.y;
  k3 = p.z*n.z - a.dot(n);
}

void CFace::coeffsFV (PV2 &a1, PV2 &a2, Parameter &k1, Parameter &k2,
		      Parameter &k3) const
{
  Vertex *v = (Vertex *) feature1;
  Face *f = (Face *) feature2;
  PV3 p = v->getP()->getP(), a = f->getBoundary(0)->tail()->getP()->getP(),
    n = f->getP()->getN();
  a1.x = n.x;
  a1.y = n.y;
  a2.x = a2.y = Parameter::constant(0.0);
  k1 = p.x*n.y - p.y*n.x;
  k2 = - p.x*n.x - p.y*n.y;
  k3 = a.dot(n) - p.z*n.z;
}

void CFace::coeffsEE (PV2 &a1, PV2 &a2, Parameter &k1, Parameter &k2,
		      Parameter &k3) const
{
  Edge *e = (Edge *) feature1, *f = (Edge *) feature2;
  PV3 a = e->getT()->getP()->getP(), u = e->getH()->getP()->getP() - a,
    b = f->getT()->getP()->getP(), v = f->getH()->getP()->getP() - b, w = u.cross(a);
  a1.x = u.y*v.z;
  a1.y = - u.x*v.z;
  a2.x = - u.z*v.y;
  a2.y = u.z*v.x;
  k1 = v.x*w.y - v.y*w.x - v.z*(b.x*u.x + b.y*u.y) + b.z*(u.x*v.x + u.y*v.y);
  k2 = - v.x*w.x - v.y*w.y + v.z*(b.y*u.x - b.x*u.y) - b.z*(u.x*v.y - u.y*v.x);
  k3 = - v.z*w.z + u.z*(b.x*v.y - b.y*v.x);
}

void CFace::line (Angle *a, Parameter &u, Parameter &v, Parameter &w) const
{
  PV2 a1, a2;
  Parameter k1, k2, k3, s = a->sin(), c = a->cos();
  coeffs(a1, a2, k1, k2, k3);
  PV2 n = a->rotate(a1) + a2;
  u = n.x;
  v = n.y;
  w = k1*s + k2*c + k3;
}

PV2 CFace::intersectionXY (CFace *f, Angle *a) const
{
  Parameter a1, b1, c1, a2, b2, c2;
  line(a, a1, b1, c1);
  f->line(a, a2, b2, c2);
  Parameter den = a1*b2 - a2*b1, x = (b1*c2 - b2*c1)/den,
    y = (a2*c1 - a1*c2)/den;
  return PV2(x, y);
}

CFace * CFace::neighbor (HCEdge *e) const
{
  CEdge *ee = e->e;
  HCEdge *f = e->forward ? e->ccw() : e->cw();
  return f->forward == e->forward ? 0 : f->f;
}

void CFace::formSlabs ()
{
  HHCEdges hhe;
  CPoints hpts;
  initSlabs(hhe, hpts);
  Slabs sl;
  for (HHCEdges::iterator e = hhe.begin(); e != hhe.end(); ++e)
    switch (e->type()) {
    case RightEnd:
      updateRE(sl, *e);
      break;
    case LeftEnd:
      updateLE(sl, *e);
      break;
    case LeftStart:
      updateLS(sl, *e);
      break;
    case RightStart:
      updateRS(sl, *e);
    }
  for (CPoints::iterator p = hpts.begin(); p != hpts.end(); ++p)
    for (Slabs::iterator s = slabs.begin(); s != slabs.end(); ++s)
      (*s)->addPoint(*p);
  for (Slabs::iterator s = slabs.begin(); s != slabs.end(); ++s) {
    CPointEdge bl((*s)->l, (*s)->s), br((*s)->r, (*s)->s),
      tl((*s)->l, (*s)->e), tr((*s)->r, (*s)->e);
    sort((*s)->pb.begin(), (*s)->pb.end(), EdgePointOrderL(&bl, &br));
    sort((*s)->pt.begin(), (*s)->pt.end(), EdgePointOrderL(&tl, &tr));
  } 
  for (Slabs::iterator s = sl.begin(); s != sl.end(); ++s)
    if (find(slabs.begin(), slabs.end(), *s) == slabs.end())
      delete *s;
}

void CFace::initSlabs (HHCEdges &hhe, CPoints &hpts) const
{
  HCEdges ed;
  boundary[0]->loop(ed);
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e) {
    hpts.push_back((*e)->tail()->getP());
    if (!(*e)->e->horizontal()) {
      hhe.push_back(HHCEdge(*e, true));
      hhe.push_back(HHCEdge(*e, false));
    }
  }
  sort(hhe.begin(), hhe.end());
}

void CFace::updateRE (Slabs &sl, const HHCEdge &h)
{
  for (Slabs::iterator s = sl.begin(); s != sl.end(); ++s)
    if ((*s)->r == h.e->e) {
      if ((*s)->l) {
	CPoint *p = h.tail()->p;
	slabs.push_back(new Slab((*s)->l, (*s)->r, (*s)->s, p->a, this));
      }
      (*s)->r = 0;
      return;
    } 
}

void CFace::updateLE (Slabs &sl, const HHCEdge &h)
{
  for (Slabs::iterator s = sl.begin(); s != sl.end(); ++s)
    if ((*s)->l == h.e->e) {
      if ((*s)->r) {
	CPoint *p = h.tail()->p;
	slabs.push_back(new Slab((*s)->l, (*s)->r, (*s)->s, p->a, this));
	for (Slabs::iterator t = sl.begin(); t != sl.end(); ++t)
	  if ((*t)->l && !(*t)->r) {
	    (*t)->r = (*s)->r;
	    (*t)->s = p->a;
	    (*t)->pb = (*s)->pb;
	    (*t)->pt = (*s)->pt;
	    (*s)->r = 0;
	    break;
	  }
      }
      (*s)->l = 0;
      return;
    }
}

void CFace::updateLS (Slabs &sl, const HHCEdge &h)
{
  CPoint *p = h.tail()->p;
  Angle *a = p->a;
  for (Slabs::iterator s = sl.begin(); s != sl.end(); ++s)
    if (!(*s)->l && (*s)->r) {
      (*s)->l = h.e->e;
      (*s)->s = a;
      return;
    }
  for (Slabs::iterator s = sl.begin(); s != sl.end(); ++s)
    if ((*s)->l && (*s)->r && InSlab(p, *s) == 1) {
      slabs.push_back(new Slab((*s)->l, (*s)->r, (*s)->s, a, this));
      CEdge *osl = (*s)->l;
      (*s)->l = h.e->e;
      (*s)->s = a;
      sl.push_back(new Slab(osl, 0, a, 0, this));
      return;
    }
  sl.push_back(new Slab(h.e->e, 0, a, 0, this));
}

void CFace::updateRS (Slabs &sl, const HHCEdge &h)
{
  CPoint *p = h.tail()->p;
  Angle *a = p->a;
  for (Slabs::iterator s = sl.begin(); s != sl.end(); ++s)
    if ((*s)->l && (!(*s)->r || InSlab(p, *s) == 1)) {
      if ((*s)->r)
	slabs.push_back(new Slab((*s)->l, (*s)->r, (*s)->s, a, this));
      (*s)->r = h.e->e;
      (*s)->s = a;
      return;
    }
}

void CFace::splitSlab (Slab *s, const Angles &an, Slabs &sl)
{
  Angle *as = s->s;
  for (int i = 0; i <= an.size(); ++i)
    if (i < an.size()) {
      sl.push_back(new Slab(s->l, s->r, as, an[i], this));
      as = an[i];
    }
    else
      sl.push_back(new Slab(s->l, s->r, as, s->e, this));
  sl[0]->pb = s->pb;
  (*sl.rbegin())->pt = s->pt;
  Slabs::iterator iter = remove(slabs.begin(), slabs.end(), s);
  slabs.erase(iter, slabs.end());
  delete s;
  slabs.insert(slabs.end(), sl.begin(), sl.end());
}

int IndependentFaces::sign ()
{
  Parameter a1, b1, c1, a2, b2, c2;
  f->line(a, a1, b1, c1);
  g->line(a, a2, b2, c2);
  return (a1*b2 - a2*b1).sign();
}

FFFData AnglesFFF::calculate ()
{
  Parameter k1, k2, k3, k4, k5, k6;
  coeffs(k1, k2, k3, k4, k5, k6);
  FFFData fd;
  n = fd.n = sinCosRoots(k1, k2, k3, k4, k5, k6, fd.p);
  return fd;
}

void AnglesFFF::coeffs (Parameter &k1, Parameter &k2, Parameter &k3, Parameter &k4,
			Parameter &k5, Parameter &k6)
{
  PV2 a11, a21, a12, a22, a13, a23;
  Parameter k11, k21, k31, k12, k22, k32, k13, k23, k33;
  f1->coeffs(a11, a21, k11, k21, k31);
  f2->coeffs(a12, a22, k12, k22, k32);
  f3->coeffs(a13, a23, k13, k23, k33);
  k1 = k2 = k3 = k4 = k5 = k6 = Parameter::constant(0.0);
  coeffs(k11, k21, k31, a12, a22, a13, a23, k1, k2, k3, k4, k5, k6);
  coeffs(k12, k22, k32, a13, a23, a11, a21, k1, k2, k3, k4, k5, k6);
  coeffs(k13, k23, k33, a11, a21, a12, a22, k1, k2, k3, k4, k5, k6);
}

void AnglesFFF::coeffs
(Parameter &k1, Parameter &k2, Parameter &k3, PV2 &a1, PV2 &b1, PV2 &a2,
 PV2 &b2, Parameter &x1, Parameter &x2, Parameter &x3, Parameter &x4,
 Parameter &x5, Parameter &x6)
{
  Parameter k4 = a2.dot(b1) - a1.dot(b2), k5 = a1.cross(b2) - a2.cross(b1),
    k6 = a1.cross(a2) + b1.cross(b2);
  x1 = x1 + k1*k4;
  x2 = x2 + k2*k5;
  x3 = x3 + k1*k5 + k2*k4;
  x4 = x4 + k1*k6 + k3*k4;
  x5 = x5 + k2*k6 + k3*k5;
  x6 = x6 + k3*k6;
}

Angle * AnglesFFF::getA (int i)
{
  if (!a[0])
    for (int j = 0; j < n; ++j)
      a[j] = new AngleFFF(this, j);
  return a[i];
}

// k1*s^2 + k2*c^2 + k3*s*c + k4*s + k5*c + k6 = 0
int sinCosRoots (Parameter k1, Parameter k2, Parameter k3, Parameter k4,
		 Parameter k5, Parameter k6, Parameter *x)
{
  Parameter a = k2 + k6 - k5, b = 2.0*(k4 - k3),
    c = 4.0*k1 - 2.0*k2 + 2.0*k6, d = 2.0*(k3 + k4), e = k2 + k5 + k6;
  if (a.sign() == 0) {
    assert(b.sign() == 0);
    return quadraticRoots(c, d, e, x);
  }
  return quarticRoots(a, b, c, d, e, x);
}

int CPointOrderX::sign ()
{
  return (b->getP().x - a->getP().x).sign();
}

void CShell::init ()
{
  CVertexSet vs;
  for (CFaces::iterator f = faces.begin(); f != faces.end(); ++f) {
    (*f)->s = this;
    CVertices vf;
    (*f)->boundaryVertices(vf);
    vs.insert(vf.begin(), vf.end());
  }
  CVertexSet::iterator v = vs.begin();
  vm = *v;
  ++v;
  while (v != vs.end()) {
    if (CPointOrderX((*v)->getP(), vm->getP()) == 1)
      vm = *v;
    ++v;
  }
}

CCell::~CCell ()
{
  for (CShells::iterator s = shells.begin(); s != shells.end(); ++s)
    delete *s;
}

Cspace::~Cspace ()
{
  if (par) {
    for (CVertices::iterator v = vertices.begin(); v != vertices.end(); ++v)
      (*v)->p = 0;
    delete par;
  }
  for (CVertices::iterator v = vertices.begin(); v != vertices.end(); ++v)
    delete *v;
  for (CEdges::iterator e = edges.begin(); e != edges.end(); ++e)
    delete *e;
  for (CFaces::iterator f = faces.begin(); f != faces.end(); ++f)
    delete *f;
  for (CCells::iterator c = cells.begin(); c != cells.end(); ++c)
    delete *c;
  for (vector<AnglesFFF *>::iterator i = fff.begin(); i != fff.end(); ++i)
    delete *i;
  for (EFAMap::iterator i = efamap.begin(); i != efamap.end(); ++i)
    delete i->second;
  for (SpiralSet::iterator s = ss.begin(); s != ss.end(); ++s)
    delete *s;
  for (SFAMap::iterator i = sfamap.begin(); i != sfamap.end(); ++i)
    delete i->second;
}

AnglesCS * Cspace::angles (Edge *e, Face *f, bool aflag)
{
  EFKey ef(e, f);
  EFAMap::iterator i = efamap.find(ef);
  if (i != efamap.end())
    return i->second;
  AnglesCS *a = new AnglesEF(e, f, aflag);
  efamap.insert(EFAPair(ef, a));
  return a;
}

AnglesCS * Cspace::angles (const EFKey &k, bool aflag)
{
  EFAMap::iterator i = efamap.find(k);
  if (i != efamap.end())
    return i->second;
  AnglesCS *a = new AnglesEF(k, aflag);
  efamap.insert(EFAPair(k, a));
  return a;
}

void Cspace::angleIntervals (Edge *e, Face *f, bool aflag, AngleIntervals &pos,
			     AngleIntervals &neg)
{
  AnglesCS *a = angles(e, f, aflag);
  if (a->roots) {
    Angle *a1 = a->getA(true), *a2 = a->getA(false);
    if (a->min) {
      pos.push_back(AngleInterval(Angle::mpi, a1));
      pos.push_back(AngleInterval(a2, Angle::ppi));
      neg.push_back(AngleInterval(a1, a2));
    }
    else {
      pos.push_back(AngleInterval(a1, a2));
      neg.push_back(AngleInterval(Angle::mpi, a1));
      neg.push_back(AngleInterval(a2, Angle::ppi));
    }
  }
  else if (a->min)
    pos.push_back(AngleInterval(Angle::mpi, Angle::ppi));
  else
    neg.push_back(AngleInterval(Angle::mpi, Angle::ppi));
}

void Cspace::angleIntervals (HEdge *e, Face *f, bool aflag, AngleIntervals &pos,
			     AngleIntervals &neg)
{
  if (e->getForward())
    angleIntervals(e->getE(), f, aflag, pos, neg);
  else
    angleIntervals(e->getE(), f, aflag, neg, pos);
}

Spiral * Cspace::getSpiral (Vertex *v, Edge *e, bool aflag)
{
  Spiral *s = new Spiral(v, e, aflag);
  pair<SpiralSet::iterator, bool> x = ss.insert(s);
  if (x.second)
    return s;
  delete s;
  return *x.first;
}

CVertex * Cspace::getVertex (CPoint *p)
{
 CVertex *v = new CVertex(p);
  if (vertices.empty())
    copyBBox(v->bbox, bbox);
  else
    mergeBBox(v->bbox, bbox);
  vertices.push_back(v);
  return v;
}

CVertex * Cspace::getVertex (Angle *a, Spiral *s)
{
  ASPair as(a, s);
  ASVMap::iterator i = asvmap.find(as);
  if (i != asvmap.end())
    return i->second;
  SpiralPoint *p = new SpiralPoint(a, s);
  CVertex *v = getVertex(p);
  asvmap.insert(ASVPair(as, v));
  return v;
}

HCEdge * Cspace::addHEdge (CVertex *a, CVertex *b, bool bflag, HCEdge *hp)
{
  CEdge *e = getEdge(a, b, bflag, hp);
  bool forward = e->t == a;
  HCEdge *he = e->addHEdge(forward);
  he->hp = hp;
  return he;
}

CEdge * Cspace::getEdge (CVertex *a, CVertex *b, bool bflag, HCEdge *hp)
{
  for (CEdges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    if (((*e)->t == a && (*e)->h == b || (*e)->t == b && (*e)->h == a) &&
	(*e)->bflag == bflag && !(*e)->hedges.empty() &&
	(bflag || !hp ||hp->e == (*e)->hedges[0]->hp->e))
      return *e;
  for (CEdges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    if (((*e)->t == a && (*e)->h == b || (*e)->t == b && (*e)->h == a) &&
	(*e)->hedges.empty()) {
      (*e)->bflag = bflag;
      return *e;
    }
  CEdge *e = new CEdge(a, b, bflag);
  edges.push_back(e);
  a->edges.push_back(e);
  b->edges.push_back(e);
  return e;
}

HCEdge * Cspace::addLoop (const CVertices &ve)
{
 int n = ve.size();
 HCEdge *ep = addHEdge(ve[n-1], ve[0], true), *e0 = ep;
 for (int i = 0; i + 1 < n; ++i) {
   HCEdge *e = addHEdge(ve[i], ve[i+1], true);
   ep->next = e;
   ep = e;
  }
 ep->next = e0;
 return e0->next;
}

void Cspace::patches (Polyhedron *a, Polyhedron *b)
{
  patchesVF(a, b, true);
  patchesVF(b, a, false);
  patchesEE(a, b);
}

void Cspace::patchesVF (Polyhedron *a, Polyhedron *b, bool aflag)
{
  for (Vertices::iterator v = a->vertices.begin(); v != a->vertices.end(); ++v) {
    HEdges ed;
    if (convexCone(*v, ed))
      for (Faces::iterator f = b->faces.begin(); f != b->faces.end(); ++f)
	patchVF(*v, *f, aflag);
  }
}

void Cspace::patchVF (Vertex *v, Face *f, bool aflag)
{
  double *vb = v->getBBox(), *fb = f->getBBox();
  if (vb[5] < fb[4] || fb[5] < vb[4])
    return;
  HEdges he;
  v->outgoingHEdges(he);
  AngleIntervals ais, aisd;
  angleIntervals(he[0], f, aflag, ais, aisd);
  for (int i = 1; !ais.empty() && i < he.size(); ++i) {
    AngleIntervals epos, eneg;
    angleIntervals(he[i], f, aflag, epos, eneg);
    ais = intersect(ais, epos);
  }
  for (AngleIntervals::iterator i = ais.begin(); i != ais.end(); ++i)
    patchVF(v, f, aflag, i->first, i->second);
}

void Cspace::patchVF (Vertex *v, Face *f, bool aflag, Angle *s, Angle *e)
{
  Spiral *l = 0, *r = 0;
  HEdges he;
  f->boundaryHEdges(he);
  for (HEdges::iterator h = he.begin(); !r && h != he.end(); ++h) {
    Edge *hh = (*h)->getE();
    if (inIntervalZ(v, hh)) {
      Spiral *sp = getSpiral(v, hh, aflag);
      if (!l)
	l = sp;
      else
	r = sp;
    }
  }
  if (OrientationVF(aflag, f, s, l, r) == -1) {
    Spiral *temp = l;
    l = r;
    r = temp;
  }
  CFace *p = new CFace(v, f, aflag);
  faces.push_back(p);
  patch(s, e, l, r, p);
}

void Cspace::patch (Angle *s, Angle *e, Spiral *l, Spiral *r, CFace *f)
{
  CVertices ve;
  ve.push_back(getVertex(s, l));
  ve.push_back(getVertex(s, r));
  ve.push_back(getVertex(e, r));
  ve.push_back(getVertex(e, l));
  HCEdge *h = addLoop(ve);
  f->addBoundary(h);
  l->edges.push_back(f->leftEdge());
  r->edges.push_back(f->rightEdge());
}

void Cspace::patchesEE (Polyhedron *a, Polyhedron *b)  
{
  Edges aedges, bedges;
  for (Edges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    if (convexEdge(*e) == 1)
      aedges.push_back(*e);
  for (Edges::iterator e = b->edges.begin(); e != b->edges.end(); ++e)
    if (convexEdge(*e) == 1)
      bedges.push_back(*e);
  for (Edges::iterator e = aedges.begin(); e != aedges.end(); ++e)
    for (Edges::iterator f = bedges.begin(); f != bedges.end(); ++f)
      patchEE(*e, *f);
}

void Cspace::patchEE (Edge *e, Edge *f)
{
  double *eb = e->getBBox(), *fb = f->getBBox();
  if (eb[5] < fb[4] || fb[5] < eb[4])
    return;
  AngleIntervals pos1, neg1, pos2, neg2, pos3, neg3, pos4, neg4;
  angleIntervals(f, e->getHEdge(0)->getF(), false, pos1, neg1);
  angleIntervals(f, e->getHEdge(1)->getF(), false, pos2, neg2);
  angleIntervals(e, f->getHEdge(0)->getF(), true, pos3, neg3);
  angleIntervals(e, f->getHEdge(1)->getF(), true, pos4, neg4);
  AngleIntervals pos = intersect(intersect(intersect(pos1, neg2), pos3), neg4),
    neg = intersect(intersect(intersect(neg1, pos2), neg3), pos4);
  for (AngleIntervals::iterator i = pos.begin(); i != pos.end(); ++i)
    patchEE(e, f, i->first, i->second, true);
  for (AngleIntervals::iterator i = neg.begin(); i != neg.end(); ++i)
    patchEE(e, f, i->first, i->second, false);
}

void Cspace::patchEE (Edge *e1, Edge *e2, Angle *s, Angle *e, bool pflag)
{
  Spiral *l = 0, *r = 0;
  Vertex *vs[] = {e1->getT(), e1->getH(), e2->getT(), e2->getH()};
  Edge *es[] = {e2, e2, e1, e1};
  bool fs[] = {true, true, false, false};
  for (int i = 0; i < 4; ++i)
    if (inIntervalZ(vs[i], es[i])) {
      Spiral *sp = getSpiral(vs[i], es[i], fs[i]);
      if (!l)
	l = sp;
      else
	r = sp;
    }
  if (OrientationEE(pflag, e1, e2, s, l, r) == -1) {
    Spiral *temp = l;
    l = r;
    r = temp;
  }
  CFace *p = new CFace(e1, e2, pflag);
  faces.push_back(p);
  patch(s, e, l, r, p);
}

void Cspace::intersectFF ()
{
  Octree<CFace *> *octree = Octree<CFace *>::octree(faces, bbox);
  CFaces fa, fb;
  octree->pairs(fa, fb);
  CEEPairSet eps;
  for (int i = 0; i < fa.size(); ++i) {
    intersectEE(fa[i], fb[i], eps);
    intersectFF(fa[i], fb[i]);
  }
  delete octree;
}

void Cspace::intersectEE (CFace *f, CFace *g, CEEPairSet &eps)
{
  HCEdges ef, eg;
  f->boundaryHEdges(ef);
  g->boundaryHEdges(eg);
  for (HCEdges::iterator i = ef.begin(); i != ef.end(); ++i) {
    CEdge *fi = (*i)->e;
    for (HCEdges::iterator j = eg.begin(); j != eg.end(); ++j) {
      CEdge *gj = (*j)->e;
      if (bboxOverlap(fi->bbox, gj->bbox) &&
	  eps.insert(CEEPair(fi < gj ? fi : gj, fi < gj ? gj : fi)).second)
	intersectEE(fi, gj);
    }
  }
}

void Cspace::intersectEE (CEdge *e, CEdge *f)
{
  Spiral *se = e->getSpiral(), *sf = f->getSpiral();
  if (se) {
    if (se == sf)
      intersectSS(e, f);
    else if (!sf)
      intersectSL(e, f);
  }
  else if (sf)
    intersectSL(f, e);
  else
    intersectLL(e, f);
}

void Cspace::intersectSS (CEdge *e, CEdge *f) const
{
  if (e->contains(f->t->p->a))
    e->addVertex(f->t);
  if (e->contains(f->h->p->a))
    e->addVertex(f->h);
  if (f->contains(e->t->p->a))
    f->addVertex(e->t);
  if (f->contains(e->h->p->a))
    f->addVertex(e->h);
}

void Cspace::intersectSL (CEdge *e, CEdge *f)
{
  if (f->piEdge())
    return;
  Angle *eta = e->t->p->a, *eha = e->h->p->a, *fa = f->t->p->a;
  if (fa == eta && f->contains(e->t->p))
    f->addVertex(e->t);
  else if (fa == eha && f->contains(e->h->p))
    f->addVertex(e->h);
  else if (inInterval(fa, eta, eha) &&
	   e->getSpiral()->sharedVertices(fa) == 3) {
    CVertex *v = getVertex(fa, e->getSpiral());
    if (f->contains(v->p)) {
      e->addVertex(v);
      f->addVertex(v);
    }
  }
}

void Cspace::intersectLL (CEdge *e, CEdge *f)
{
  Angle *a = e->t->p->a;
  if (a != f->t->p->a)
    return;
  if (a->piAngle()) {
    CFace *ef = e->hedges[0]->f, *ff = f->hedges[0]->f;
    if (IndependentFaces(ef, ff, a) == 0)
      return;
    CPoint *p = new FFFPoint(a, ef, ff, 0);
    if (e->contains(p) && f->contains(p)) {
      CVertex *v = getVertex(p);
      e->addVertex(v);
      f->addVertex(v);
    }
    else
      delete p;    
    return;
  }
  if (e->contains(f->t->p))
    e->addVertex(f->t);
  if (e->contains(f->h->p))
    e->addVertex(f->h);
  if (f->contains(e->t->p))
    f->addVertex(e->t);
  if (f->contains(e->h->p))
    f->addVertex(e->h);
}

void Cspace::intersectFF (CFace *f, CFace *g)
{
  CVertices vfg;
  intersectFF(f, g, vfg);
  intersectFF(g, f, vfg);
  formFF(f, g, vfg);
}

void Cspace::intersectFF (CFace *f, CFace *g, CVertices &vfg)
{
  HCEdges eg;
  g->boundaryHEdges(eg);
  for (HCEdges::iterator e = eg.begin(); e != eg.end(); ++e)
    intersectEF(*e, f, vfg);
}

void Cspace::intersectEF (HCEdge *e, CFace *f, CVertices &vfg)
{
  CEdge *ee = e->e;
  if (!bboxOverlap(ee->bbox, f->bbox))
    return;
  AngleInterval eai = ee->angleInterval();
  if (!intervalOverlap(eai.first, eai.second, f->startAngle(), f->endAngle()))
    return;
  CEFPair ef(ee, f);
  CEFVMap::iterator i = efvmap.find(ef);
  if (i != efvmap.end()) {
    vfg.insert(vfg.end(), i->second.begin(), i->second.end());
    return;
  }
  CVertices ve;
  intersectEF(ee, f, ve);
  efvmap.insert(CEFVPair(ef, ve));
  vfg.insert(vfg.end(), ve.begin(), ve.end());
}

void Cspace::intersectEF (CEdge *e, CFace *f, CVertices &ve)
{
  Spiral *s = e->getSpiral();
  if (!s) {
    CVertex *v = intersectEFLine(e, f);
    if (v)
      ve.push_back(v);
    return;
  }
  if (s == f->leftEdge()->s || s == f->rightEdge()->s)
    return;
  AnglesCS *asf = intersectSF(s, f);
  if (asf->roots) {
    AngleInterval eai = e->angleInterval();
    Angle *fs = f->startAngle(), *fe = f->endAngle(),
      *as = angleOrder(eai.first, fs) ? fs : eai.first,
      *ae = angleOrder(eai.second, fe) ? eai.second : fe;
    for (int i = 0; i < 2; ++i) {
      Angle *a = asf->getA(i == 0);
      if (inInterval(a, as, ae)) {
	SpiralPoint p(a, s);
	if (f->contains(&p)) {
	  CVertex *v = getVertex(a, s);
	  e->vertices.push_back(v);
	  ve.push_back(v);
	}
      }
    }
  }
}

CVertex * Cspace::intersectEFLine (CEdge *e, CFace *f)
{
  Angle *a = e->t->p->a;
  AFPair af(a, f);
  CVertex *v = 0;
  AFVMap::iterator i = afvmap.find(af);
  if (i != afvmap.end())
    v = i->second;
  else if (f->sharedVertices(a) < 3) {
    CPoint *p = new LinePatchPoint(e, f);
    if (f->contains(p))
      v = getVertex(p);
    else
      delete p;
    afvmap.insert(AFVPair(af, v));
  }
  if (v && e->contains(v->p)) {
    e->addVertex(v);
    return v;
  }
  return 0;
}

AnglesCS * Cspace::intersectSF (Spiral *s, CFace *f)
{
  AnglesCS *al = intersectSFLine(s, f);
  if (al)
    return al;
  SFPair sf(s, f);
  SFAMap::iterator i = sfamap.find(sf);
  if (i != sfamap.end())
    return i->second;
  AnglesCS *a = new AnglesSF(s, f);
  sfamap.insert(SFAPair(sf, a));
  return a;
}

AnglesCS * Cspace::intersectSFLine (Spiral *s, CFace *f)
{
  VertexSet vsa, vsb;;
  s->vertices(vsa, vsb);
  f->vertices(vsa, vsb);
  if (vsa.size() == 2 && vsb.size() == 3) {
    VertexSet::iterator i = vsb.begin();
    Vertex *v1 = *vsa.begin(), *v2 = *vsa.rbegin(), *v3 = *i, *v5 = *vsb.rbegin();
    ++i;
    EFKey k(v1, v2, v3, *i, v5);
    return angles(k, true);
  }
  if (vsa.size() == 3 && vsb.size() == 2) {
    VertexSet::iterator i = vsa.begin();
    Vertex *v1 = *vsb.begin(), *v2 = *vsb.rbegin(), *v3 = *i, *v5 = *vsa.rbegin();
    ++i;
    EFKey k(v1, v2, v3, *i, v5);
    return angles(k, false);
  }
  return 0;
}

void Cspace::formFF (CFace *f, CFace *g, CVertices &vfg)
{
  if (intersectsFFLine(f, g))
    formFFLine(f, g, vfg);
  else
    formFFCurve(f, g, vfg);
}

bool Cspace::intersectsFFLine (CFace *f, CFace *g) const
{
  VertexSet vsa, vsb;
  f->vertices(vsa, vsb);
  g->vertices(vsa, vsb);
  return vsa.size() == 2 && vsb.size() == 3 || vsa.size() == 3 && vsb.size() == 2;
}

void Cspace::formFFLine (CFace *f, CFace *g, CVertices &vfg)
{
  sort(vfg.begin(), vfg.end(), VertexAngleOrder());
  int i = 0, n = vfg.size();
  while (i < n)
    if (i + 1 < n && vfg[i]->p->a == vfg[i+1]->p->a) {
      formFF(f, g, vfg[i], vfg[i+1]);
      i += 2;
    }
    else {
      Spiral *s = f->sharedSpiral(g);
      formFF(f, g, vfg[i], spiralVertex(vfg[i]->p->a, s));
      ++i;
    }
}

CVertex * Cspace::spiralVertex (Angle *a, Spiral *s)
{
  int nv = vertices.size();
  CVertex *v = getVertex(a, s);
  if (nv < vertices.size())
    for (CEdges::iterator e = s->edges.begin(); e != s->edges.end(); ++e)
      if ((*e)->contains(v->p))
	(*e)->addVertex(v);
  return v;
}

void Cspace::formFFCurve (CFace *f, CFace *g, CVertices &vfg)
{
  int n = vfg.size();
  if (n != 2 && n != 4)
    f->sharedBoundaryVertices(g, vfg);
  sort(vfg.begin(), vfg.end(), VertexAngleOrder());
  for (int i = 0; i + 1 < vfg.size(); i += 2)
    if (!(n == 0 && f->sameBoundary(vfg[i], vfg[i+1]) &&
	  g->sameBoundary(vfg[i], vfg[i+1])))
      formFF(f, g, vfg[i], vfg[i+1]);
}

void Cspace::formFF (CFace *f, CFace *g, CVertex *v, CVertex *w)
{
  CEdge *e = CFFOrder(f, g, v, w) == 1 ? getEdge(v, w, false)
    : getEdge(w, v, false);
  f->edges.push_back(e);
  e->addHEdge(true)->f = f;
  g->edges.push_back(e);
  e->addHEdge(false)->f = g;
} 

void Cspace::intersectFFF ()
{
  CFFEMap ffemap;
  formFFEMap(ffemap);
  for (CFaces::iterator f = faces.begin(); f != faces.end(); ++f)
    intersectFFF(*f, ffemap);
}

void Cspace::formFFEMap (CFFEMap &ffemap) const
{

  CEdges::const_iterator e = edges.begin();
  while (e != edges.end() && (*e)->bflag)
    ++e;
  while (e != edges.end()) {
    CFace *f = (*e)->hedges[0]->f, *g = (*e)->hedges[1]->f;
    CFFPair fg(f < g ? f : g, f < g ? g : f);
    CFFEMap::iterator i = ffemap.find(fg);
    if (i == ffemap.end()) {
      CEdges ed;
      ed.push_back(*e);
      ffemap.insert(CFFEPair(fg, ed));
    }
    else
      i->second.push_back(*e);
    ++e;
  }
}

void Cspace::intersectFFF (CFace *f, const CFFEMap &ffemap)
{
  int n = f->edges.size();
  for (int i = 0; i + 1 < n; ++i) {
    CEdge *ei = f->edges[i];
    CFace *fi = ei->otherFace(f);
    if (f < fi)
      for (int j = i + 1; j < n; ++j) {
	CEdge *ej = f->edges[j];
	CFace *fj = ej->otherFace(f);
	if (f < fj && fi != fj)
	  intersectFFF(f, fi, fj, ei, ej, ffemap);
      }
  }
}

void Cspace::intersectFFF (CFace *f1, CFace *f2, CFace *f3, CEdge *e12, CEdge *e13,
			   const CFFEMap &ffemap)
{
  if (!(f1->bboxOverlap(f2, f3) && f1->angleOverlap(f2, f3)))
    return;
  CFFPair f23(f2 < f3 ? f2 : f3, f2 < f3 ? f3 : f2);
  CFFEMap::const_iterator iter = ffemap.find(f23);
  if (iter == ffemap.end())
    return;
  const CEdges ed = iter->second;
  if (e12->horizontal())
      for (CEdges::const_iterator e23 = ed.begin(); e23 != ed.end(); ++e23)
	intersectFFFH(f1, f2, f3, e12, e13, *e23);
  else if (e13->horizontal())
    for (CEdges::const_iterator e23 = ed.begin(); e23 != ed.end(); ++e23)
      intersectFFFH(f1, f3, f2, e13, e12, *e23);
  else if (ed[0]->horizontal())
    for (CEdges::const_iterator e23 = ed.begin(); e23 != ed.end(); ++e23)
      intersectFFFH(f2, f3, f1, *e23, e12, e13);
  else
    intersectFFFG(f1, f2, f3, e12, e13, ed);
}

void Cspace::intersectFFFH (CFace *f1, CFace *f2, CFace *f3, CEdge *eh,
			    CEdge *e1, CEdge *e2)
{
  Angle *a = eh->t->p->a;
  if (e1->contains(a) && e2->contains(a)) {
    CVertex *v = intersectEFLine(eh, f3);
    if (v && onEdge(v->p, eh->t->p, eh->h->p)) {
      eh->addVertex(v);
      e1->addVertex(v);
      e2->addVertex(v);
    }
  }
}

void Cspace::intersectFFFG (CFace *f1, CFace *f2, CFace *f3, CEdge *e12,
			    CEdge *e13, const CEdges &ed)
{
  AnglesFFF *a = new AnglesFFF(f1, f2, f3);
  if (a->n == 0) {
    delete a;
    return;
  }
  fff.push_back(a);
  for (int i = 0; i < a->n; ++i) {
    Angle *ai = a->getA(i);
    if (e12->contains(ai) && e13->contains(ai))
      for (CEdges::const_iterator e23 = ed.begin(); e23 != ed.end(); ++e23)
	if ((*e23)->contains(ai)) {
	  CPoint *p = new FFFPoint(ai, f1, f2, f3);
	  CVertex *v = getVertex(p);
	  e12->addVertex(v);
	  e13->addVertex(v);
	  (*e23)->addVertex(v);
	  break;
	}
  }
}

void Cspace::sortVertices ()
{
  for (CEdges::const_iterator e = edges.begin(); e != edges.end(); ++e)
    if (!(*e)->vertices.empty())
      if ((*e)->horizontal())
	sort((*e)->vertices.begin(), (*e)->vertices.end(), EdgeVertexOrderL(*e));
      else
	sort((*e)->vertices.begin(), (*e)->vertices.end(), EdgeVertexOrderS(*e));
}

Cspace * Cspace::subfaces ()
{
  Cspace *a = new Cspace(this);
  CVVMap vvmap;
  for (CFaces::const_iterator f = faces.begin(); f != faces.end(); ++f)
    subfaces(*f, a, vvmap);
  return a;
}

void Cspace::subfaces (CFace *f, Cspace *a, CVVMap &vvmap) const
{
  HCEdges he, outer, bad;
  subedges(f, a, vvmap, he);
  for (HCEdges::iterator e = he.begin(); e != he.end(); ++e)
    if (!(*e)->flag)
      if ((*e)->onLoop())
	outer.push_back(*e);
      else
	bad.push_back(*e);
  CFaces nfaces;
  for (HCEdges::iterator e = outer.begin(); e != outer.end(); ++e)
    nfaces.push_back(a->addFace(*e, f));
  a->removeBad(bad);
}

void Cspace::subedges (CFace *f, Cspace *a, CVVMap &vvmap, HCEdges &he) const
{
  HCEdges ed;
  f->boundaryHEdges(ed);
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    subedges(*e, a, vvmap, he);
  for (CEdges::iterator e = f->edges.begin(); e != f->edges.end(); ++e) {
    HCEdge *h = (*e)->hedges[0]->f == f ? (*e)->hedges[0] : (*e)->hedges[1];
    subedges(h, a, vvmap, he);
  }
  setNext(f, he);
}

void Cspace::subedges (HCEdge *e, Cspace *a, CVVMap &vvmap, HCEdges &he) const
{
  CVertices ve;
  e->e->edgeVertices(ve);
  CVertex *t = a->getVertex(ve[0], vvmap);
  for (int i = 1; i < ve.size(); ++i) {
    CVertex *h = a->getVertex(ve[i], vvmap);
    if (e->forward)
      he.push_back(a->addHEdge(t, h, e->e->bflag, e));
    else
      he.push_back(a->addHEdge(h, t, e->e->bflag, e));
    t = h;
  }
}

CVertex * Cspace::getVertex (CVertex *v, CVVMap &vmap)
{
  CVVMap::iterator iter = vmap.find(v);
  if (iter != vmap.end())
    return iter->second;
  CVertex *w = getVertex(v->p);
  vmap.insert(CVVPair(v, w));
  return w;
}

void Cspace::setNext (CFace *f, const HCEdges &he) const
{
  HHCEdges hhe;
  for (HCEdges::const_iterator e = he.begin(); e != he.end(); ++e) {
    hhe.push_back(HHCEdge(*e, true));
    hhe.push_back(HHCEdge(*e, false));
  }
  sort(hhe.begin(), hhe.end(), HHCEdgeOrder(f));
  int i = 0, n = hhe.size();
  while (i < n) {
    int m = 1;
    while (i + m < n && hhe[i].tail() == hhe[i+m].tail())
      ++m;
    int j = 0;
    while (j < m) {
      int k = i + (j + 1)%m;
      if (!hhe[i+j].f && hhe[k].f)
	hhe[i+j].e->next = hhe[k].e;
      ++j;
    }
    i += m;
  }
}

CFace * Cspace::addFace (HCEdge *e, CFace *f)
{
  CFace *g = new CFace(f->type, f->feature1, f->feature2);
  faces.push_back(g);
  g->addBoundary(e);
  return g;
}

void Cspace::removeBad (const HCEdges &ed)
{
  HCEdgeSet es;
  for (HCEdges::const_iterator e = ed.begin(); e != ed.end(); ++e) {
    HCEdge *f = *e;
    while (f) {
      es.insert(f);
      f = f->next;
      if (f == *e)
	break;
    }
  }
  for (HCEdgeSet::iterator e = es.begin(); e != es.end(); ++e)
    (*e)->e->removeHEdge(*e);
}

void Cspace::formCells (Polyhedron *a, Polyhedron *b)
{
  CShells sh;
  formShells(a, b, sh);
  sort(sh.begin(), sh.end(), CShellOrder());
  Octree<CFace *> *octree = Octree<CFace *>::octree(faces, bbox);
  cells.push_back(new CCell);
  for (int i = 0; i < sh.size(); ++i) {
    CShell *s = sh[i];
    CShell *t = i == 0 ? 0 : enclosingShell(s, octree);
    if (!t) {
      s->outer = false;
      cells[0]->shells.push_back(s);
    }
    else if (t->outer) {
      s->outer = false;
      t->c->shells.push_back(s);
    }
    else {
      s->outer = true;
      CCell *c = new CCell;
      c->shells.push_back(s);
      cells.push_back(c);
    }
  }
  delete octree;
}

void Cspace::formShells (Polyhedron *a, Polyhedron *b, CShells &sh)
{
  CShells sh0;
  formShells(sh0);
  for (CShells::iterator s = sh0.begin(); s != sh0.end(); ++s) {
    CEdge *e = maxAngleEdge(*s);
    double eam = 0.5*(e->t->p->a->theta() + e->h->p->a->theta());
    Angle *th = new InputAngle(eam);
    CPointEdge c(e, th);
    Polyhedron *ac = transform(a, &c);
    bool bad = ac->intersects(b);
    delete ac;
    if (!bad) {
      (*s)->init();
      sh.push_back(*s);
    }
    else
      delete *s;
  }
  int i = 0;
  while (i < faces.size())
    if (faces[i]->s)
      ++i;
    else {
      removeFace(faces[i]);
      faces[i] = *faces.rbegin();
      faces.pop_back();
    }
}

void Cspace::formShells (CShells &sh) const
{
  for (CEdges::const_iterator e = edges.begin(); e != edges.end(); ++e)
    (*e)->sortHEdges();
  set<CFace *> done;
  for (CFaces::const_iterator f = faces.begin(); f != faces.end(); ++f)
    if (done.insert(*f).second) {
      CShell *s = new CShell;
      CFaces st;
      st.push_back(*f);
      bool flag = true;
      while (!st.empty()) {
	CFace *g = *st.rbegin();
	st.pop_back();
	s->faces.push_back(g);
	HCEdges ed;
	g->boundaryHEdges(ed);
	for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e)
	  if (!(*e)->e->piEdge()) {
	    CFace *h = g->neighbor(*e);
	    if (!h)
	      flag = false;
	    else if (done.insert(h).second)
	      st.push_back(h);
	  }
      }
      if (flag)
	sh.push_back(s);
      else
	delete s;
    }
}

CEdge * Cspace::maxAngleEdge (CShell *s) const
{
  CEdgeSet es;
  double dtmax = 0.0;
  CEdge *emax = 0;
  for (CFaces::iterator f = s->faces.begin(); f != s->faces.end(); ++f) {
    HCEdges ed;
    (*f)->boundaryHEdges(ed);
    for (HCEdges::iterator h = ed.begin(); h != ed.end(); ++h) {
      CEdge *e = (*h)->e;
      if (!e->horizontal() && es.insert(e).second) {
	double dt = fabs(e->t->p->a->theta() - e->h->p->a->theta());
	if (dt > dtmax) {
	  dtmax = dt;
	  emax = e;
	}
      }
    }
  }
  return emax;
}

CShell * Cspace::enclosingShell (CShell *s, Octree<CFace *> *octree)
{
  CPoint *a = s->vm->p;
  double bb[6];
  a->getBBox(bb);
  bb[0] = bbox[0];
  CFaces fa;
  octree->find(bb, fa);
  CFace *f = 0;
  CPoint *p = 0;
  for (CFaces::iterator g = fa.begin(); g != fa.end(); ++g)
    if ((*g)->s != s) {
      CPoint *q = new RayPatchPointX(a, *g);
      if (CPointOrderX(q, a) == 1 && (*g)->contains2(q) &&
	  (!p || CPointOrderX(q, p) == 1)) {
	if (p)
	  delete p;
	p = q;
	f = *g;
      }
      else
	delete q;
    }
  delete p;
  return f ? f->s : 0;
}

void Cspace::removeFace (CFace *f) const
{
  HCEdges ed;
  f->boundaryHEdges(ed);
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    (*e)->e->removeHEdge(*e);
  delete f;
}

void Cspace::describe () const
{
  for (int i = 0; i < cells.size(); ++i) {
    CCell *c = cells[i];
    cerr << "cell " << i << ": " << c->shells.size() << " shells; ";
    for (CShells::iterator s = c->shells.begin(); s != c->shells.end(); ++s)
      cerr << (*s)->faces.size() << " hfaces ";
    cerr << endl;
  }
}

int PointOrderZ::sign ()
{
  return (b->getP().z - a->getP().z).sign();
}

bool inIntervalZ (Vertex *v, Edge *e)
{
  Point *p = v->getP(), *t = e->getT()->getP(), *h = e->getH()->getP();
  if (p == t || p == h)
    return false;
  int s = PointOrderZ(t, h);
  return PointOrderZ(t, p) == s && PointOrderZ(p, h) == s;
}

int CFFOrder::sign()
{
  CPoint *vp = v->getP(), *wp = w->getP();
  Angle *va = vp->getA(), *wa = wp->getA();
  PV3 u = f->getN(vp).cross(g->getN(vp));
  if (va == wa) {
    PV2 vw = vp->getP() - wp->getP(), u2(u.x, u.y);
    return u2.dot(vw).sign();
  }
   return - u.z.sign();
}

int EdgeOrderL::sign ()
{
  PV2 u = h->getP() - t->getP(), vw = w->getP() - v->getP();
  return u.dot(vw).sign();
}

int CVertexHHEdgeOrder::sign ()
{
  bool efor = v == e->tail(), ffor = v == f->tail();
  CPoint *p = v->getP();
  PV3 eu = efor ? e->getU(p) : - e->getU(p),
    fu = ffor ? f->getU(p) : - f->getU(p),
    r = Rdir.getP();
  int es = eu.dot(r).sign(), fs = fu.dot(r).sign();
  return es == fs ? g->getN(p).tripleProduct(fu, eu).sign() : es;
}

Cspace * cspace (Polyhedron *a, Polyhedron *b)
{
  Cspace *c = new Cspace;
  c->patches(a, b);
  c->intersectFF();
  c->intersectFFF();
  c->sortVertices();
  Cspace *d = c->subfaces();
  d->formCells(a, b);
  return d;
}

Polyhedron * transform (Polyhedron *p, CPoint *c)
{
  Polyhedron *q = new Polyhedron;
  VVMap vvmap;
  for (Vertices::const_iterator v = p->vertices.begin(); v != p->vertices.end(); ++v)
    vvmap.insert(VVPair(*v, q->getVertex(new TransformedPoint(c, (*v)->getP()))));
  for (Faces::const_iterator f = p->faces.begin(); f != p->faces.end(); ++f) {
    HEdge *e = (*f)->getBoundary(0);
    Vertex *u = vvmap.find(e->tail())->second,
      *v = vvmap.find(e->head())->second,
      *w = vvmap.find(e->getNext()->head())->second;
    q->addTriangle(u, v, w);
  }
  return q;
}

Polyhedron * discretize (Cspace *c, double d)
{
  for (CFaces::iterator f = c->faces.begin(); f != c->faces.end(); ++f)
    (*f)->formSlabs();
  CPVMap pvmap;
  Polyhedron *a = discretize(c->faces, d, pvmap);
  facesPI(c, a, pvmap, true);
  facesPI(c, a, pvmap, false);
  return a;
}

Polyhedron * discretize (const CFaces &fa, double d, CPVMap &pvmap)
{
  EPMap epmap = discretizeEdges(fa, d);
  delentil(fa, epmap);
  Polyhedron *a = new Polyhedron;
  FSMap fsmap;
  SFSMap sfmap;
  for (CFaces::const_iterator f = fa.begin(); f != fa.end(); ++f) {
    const Slabs &sl = (*f)->getSlabs();
    for (Slabs::const_iterator s = sl.begin(); s != sl.end(); ++s)
      discretize(*s, pvmap, epmap, a, fsmap, sfmap);
  }
  removeIntersections(a, pvmap, epmap, fsmap, sfmap);
  return a;
}

EPMap discretizeEdges (const CFaces &fa, double d)
{
  EPMap epmap;
  set<CEdge *> es;
  for (CFaces::const_iterator f = fa.begin(); f != fa.end(); ++f) {
    HCEdges ed;
    (*f)->boundaryHEdges(ed);
    for (HCEdges::iterator h = ed.begin(); h != ed.end(); ++h) {
      CEdge *e = (*h)->getE();
      if (!e->piEdge())
	es.insert(e);
    }
  }
  for (set<CEdge *>::iterator e = es.begin(); e != es.end(); ++e)
    epmap.insert(EPPair(*e, discretize(*e, d)));
  return epmap;
}

CPoints discretize (CEdge *e, double d)
{
  CVertices ve;
  for (int i = 0; i < e->HEdgesN(); ++i) {
    CVertices vei;
    e->getHEdge(i)->getF()->boundaryVertices(vei);
    for (CVertices::iterator v = vei.begin(); v != vei.end(); ++v)
      if (e->contains((*v)->getP()->getA()))
	ve.push_back(*v);
  }
  ve.push_back(e->getT());
  ve.push_back(e->getH());
  sort(ve.begin(), ve.end(), VertexAngleOrder());
  int n = ve.size();
  CPoints pts;
  pts.push_back(ve[0]->getP());
  for (int i = 0; i + 1 < n; ++i) {
    Angle *a1 = ve[i]->getP()->getA(), *a2 = ve[i+1]->getP()->getA();
    if (a1 != a2) {
      discretize(e, a1, a2, d, pts);
      if (i + 2 < n)
	pts.push_back(new CPointEdge(e, a2));
    }
  }
  pts.push_back(ve[n-1]->getP());
  return pts;
}

void discretize (CEdge *e, Angle *as, Angle *ae, double d, CPoints &pts)
{
  CPointEdge ps(e, as), pe(e, ae);
  CPoints st, res;
  res.push_back(&ps);
  st.push_back(&pe);
  while (!st.empty()) {
    CPoint *p = *res.rbegin(), *q = *st.rbegin();
    Angle *a = new InputAngle(0.5*(p->getA()->theta() + q->getA()->theta()), false);
    CPoint *m = new CPointEdge(e, a);
    if (close(m, p, q, d)) {
      delete m;
      res.push_back(q);
      st.pop_back();
    }
    else
      st.push_back(m);
  }
  for (int i = 1; i + 1 < res.size(); ++i)
    pts.push_back(res[i]);
}

bool close (CPoint *a, CPoint *t, CPoint *h, double d)
{
  PointCPoint pa(a), pt(t), ph(h);
  return DistancePL(&pa, &pt, &ph, d) == 1;
}

void delentil (const CFaces &fa, EPMap &epmap)
{
  for (CFaces::const_iterator f = fa.begin(); f != fa.end(); ++f) {
    const Slabs &sl = (*f)->getSlabs();
    for (Slabs::const_iterator s = sl.begin(); s != sl.end(); ++s)
      delentil(*s, epmap);
  }
}

void delentil (Slab *s, EPMap &epmap)
{
  CPoints l = getPoints(s->l, s->s, s->e, epmap),
    r = getPoints(s->r, s->s, s->e, epmap);
  if (l.size() == 2 && r.size() == 2 && l[0] == r[0] && l[1] == r[1]) {
    Angle *a = new InputAngle(0.5*(s->s->theta() + s->e->theta()));
    addPoint(s->l, a, epmap);
  }
}

void discretize (Slab *s, CPVMap &pvmap, const EPMap &epmap, Polyhedron *a,
		 FSMap &fsmap, SFSMap &sfmap)
{
  CTriangles tr = discretize(s, epmap);
  Faces fa;
  for (CTriangles::iterator t = tr.begin(); t != tr.end(); ++t) {
    Face *g = a->addTriangle(getVertex(t->a, pvmap, a),
			     getVertex(t->b, pvmap, a),
			     getVertex(t->c, pvmap, a));
    fsmap.insert(FSPair(g, s));
    fa.push_back(g);
  }
  sfmap.insert(SFSPair(s, fa));
}

CTriangles discretize (Slab *s, const EPMap &epmap)
{
  CPoints l = getPoints(s->l, s->s, s->e, epmap),
    r = getPoints(s->r, s->s, s->e, epmap);
  CTriangles trs = discretize(l, r, s), tr;
  int is = s->pb.empty() ? 0 : 1, n = trs.size(),
    ie = s->pt.empty() ? n : n - 1;
  if (is == 1)
    splitB(trs[0], s->pb, trs[0].a == l[1], tr);
  for (int i = is; i < ie; ++i)
    tr.push_back(trs[i]);
  if (ie < n)
    splitT(trs[n-1], s->pt, trs[n-1].a == *l.rbegin(), tr);
  return tr;
}

CPoints getPoints (CEdge *f, Angle *s, Angle *e, const EPMap &epmap)
{
  const CPoints &fpts = epmap.find(f)->second;
  CPoints pts;
  CPoints::const_iterator p = fpts.begin();
  while ((*p)->getA() != s)
    ++p;
  while (true) {
    pts.push_back(*p);
    if ((*p)->getA() == e)
      break;
    ++p;
  }
  return pts;
}

CTriangles discretize (const CPoints &l, const CPoints &r, Slab *s)
{
  CTriangles tr;
  int il = 0, ir = 0, nl = l.size(), nr = r.size();
  while (ir + 1 < nr)
    if (angleOrder(r[ir]->getA(), l[il]->getA())) {
      int i2 = il + 1;
      while (i2 < nl && !angleOrder(r[ir+1]->getA(), l[i2]->getA()))
	++i2;
      double mid = 0.5*(r[ir]->getA()->theta() + r[ir+1]->getA()->theta());
      int m = il;
      for (int j = il + 1; j < i2; ++j)
        if (fabs(l[j]->getA()->theta() - mid) < fabs(l[m]->getA()->theta() - mid))
          m = j;
      for (int j = il + 1; j <= m; ++j)
	ctriangle(l[j], l[j-1], r[ir], tr);
      ctriangle(r[ir], r[ir+1], l[m], tr);
      for (int j = m + 1; j < i2; ++j)
	ctriangle(l[j], l[j-1], r[ir+1], tr);
      il = i2 - 1;
      ++ir;
    }
    else {
      int i2 = ir + 1;
      while (angleOrder(r[i2]->getA(), l[il+1]->getA()))
	++i2;
      double mid = 0.5*(l[il]->getA()->theta() + l[il+1]->getA()->theta());
      int m = ir;
      for (int j = ir + 1; j < i2; ++j)
        if (fabs(r[j]->getA()->theta() - mid) < fabs(r[m]->getA()->theta() - mid))
          m = j;
      for (int j = ir + 1; j <= m; ++j)
	ctriangle(r[j-1], r[j], l[il], tr);
      ctriangle(l[il+1], l[il], r[m], tr);
      for (int j = m + 1; j < i2; ++j)
	ctriangle(r[j-1], r[j], l[il+1], tr);
      ir = i2 - 1;
      ++il;
    }
  return tr;
}

void ctriangle (CPoint *a, CPoint *b, CPoint *c, CTriangles &tr)
{
  if (a == b || a == c || b == c)
    return;
  tr.push_back(CTriangle(a, b, c));
}

void splitB (const CTriangle &t, const CPoints &pb, bool lflag, CTriangles &tr)
{
  int n = pb.size();
  if (lflag) {
    ctriangle(t.a, t.b, pb[0], tr);
    for (int i = 0; i < n; ++i)
      ctriangle(t.a, pb[i], i + 1 < n ? pb[i+1] : t.c, tr);
  }
  else {
    ctriangle(t.a, t.b, pb[n-1], tr);
    for (int i = 0; i < n; ++i)
      ctriangle(t.b, i == 0 ? t.c : pb[i-1], pb[i], tr);
  }
}

void splitT (const CTriangle &t, const CPoints &pt, bool lflag, CTriangles &tr)
{
  int n = pt.size();
  if (lflag) {
    ctriangle(t.a, t.b, pt[0], tr);
    for (int i = 0; i < n; ++i)
      ctriangle(t.b, i + 1 < n ? pt[i+1] : t.c, pt[i], tr);
  }
  else {
    ctriangle(t.a, t.b, pt[n-1], tr);
    for (int i = 0; i < n; ++i)
      ctriangle(t.a, pt[i], i == 0 ? t.c : pt[i-1], tr);
  }
}

Vertex * getVertex (CPoint *p, CPVMap &pvmap, Polyhedron *a)
{
  CPVMap::iterator i = pvmap.find(p);
  if (i != pvmap.end())
    return i->second;
  bool flag = dynamic_cast<CPointEdge *>(p);
  Vertex *v = a->getVertex(new PointCPoint(p, flag));
  pvmap.insert(CPVPair(p, v));
  return v;
}

void removeIntersections (Polyhedron *a, CPVMap &pvmap, EPMap &epmap,
			  FSMap &fsmap, SFSMap &sfmap)
{
  FOctree *octree = a->faceOctree();
  bool flag = false;
  int i = 0;
  while (i < a->faces.size()) {
    Face *f = a->faces[i];
    ++i;
    if (f->getBoundary().empty())
      continue;
    CFace *fc = fsmap.find(f)->second->f;
    Faces fa;
    octree->find(f->getBBox(), fa);
    for (Faces::iterator g = fa.begin(); g != fa.end(); ++g)
      if (!(*g)->getBoundary().empty() && fsmap.find(*g)->second->f != fc &&
	  f->intersects(*g)) {
	flag = true;
	splitSlabs(f, *g, pvmap, epmap, a, fsmap, sfmap, octree);
	break;
      }
  }
  if (flag)
    a->removeNullFaces();
  delete octree;
}

void splitSlabs (Face *f, Face *g, CPVMap &pvmap, EPMap &epmap, Polyhedron *a,
		 FSMap &fsmap, SFSMap &sfmap, FOctree *octree)
{
  Slab *sf = fsmap.find(f)->second, *sg = fsmap.find(g)->second;
  Angles af, ag;
  if (sf->s == sg->s && sf->e == sg->e) {
    Angle *m = new InputAngle(0.5*(sf->s->theta() + sf->e->theta()));
    af.push_back(m);
    ag.push_back(m);
  }
  else {
    af.push_back(sg->s);
    af.push_back(sg->e);
    ag.push_back(sf->s);
    ag.push_back(sf->e);
  }
  splitSlab(sf, af, pvmap, epmap, a, fsmap, sfmap, octree);
  if (sf != sg)
    splitSlab(sg, ag, pvmap, epmap, a, fsmap, sfmap, octree);
}

void splitSlab (Slab *s, const Angles &as, CPVMap &pvmap, EPMap &epmap,
		Polyhedron *a, FSMap &fsmap, SFSMap &sfmap, FOctree *octree)
{
  Slabs sl, ss;
  neighbors(s, s->l, sl);
  neighbors(s, s->r, sl);
  Angles ain;
  for (Angles::const_iterator i = as.begin(); i != as.end(); ++i)
    if (inInterval(*i, s->s, s->e)) {
      addPoint(s->l, *i, epmap);
      addPoint(s->r, *i, epmap);
      ain.push_back(*i);
    }
  if (ain.empty())
    return;
  remove(s, a, fsmap, sfmap);
  for (Slabs::iterator t = sl.begin(); t != sl.end(); ++t)
    remove(*t, a, fsmap, sfmap);
  s->f->splitSlab(s, ain, ss);
  sl.insert(sl.end(), ss.begin(), ss.end());
  for (Slabs::iterator t = sl.begin(); t != sl.end(); ++t)
    rediscretize(*t, pvmap, epmap, a, fsmap, sfmap, octree);
}

void neighbors (Slab *s, CEdge *e, Slabs &sl)
{
  CFace *f = e->getHEdge(0)->getF();
  if (f == s->f)
    f = e->getHEdge(1)->getF();
  for (Slabs::const_iterator t = f->getSlabs().begin(); t != f->getSlabs().end(); ++t)
    if (((*t)->l == e || (*t)->r == e) &&
	intervalOverlap(s->s, s->e, (*t)->s, (*t)->e))
      sl.push_back(*t);
}

void remove (Slab *s, Polyhedron *a, FSMap &fsmap, SFSMap &sfmap)
{
  Faces &fa = sfmap.find(s)->second;
  for (Faces::iterator f = fa.begin(); f != fa.end(); ++f) {
    a->removeLoop((*f)->getBoundary(0));
    fsmap.erase(*f);
  }
  sfmap.erase(s);
}

void addPoint (CEdge *e, Angle *a, EPMap &epmap)
{
  CPoints &pts = epmap.find(e)->second;
  for (CPoints::iterator i = pts.begin(); i + 1 != pts.end(); ++i)
    if (inInterval(a, (*i)->getA(), (*(i+1))->getA())) {
      pts.insert(i+1, new CPointEdge(e, a));
      break;
    }
}

void rediscretize (Slab *s, CPVMap &pvmap, const EPMap &epmap, Polyhedron *a,
		   FSMap &fsmap, SFSMap &sfmap, FOctree *octree)
{
  int nf = a->faces.size();
  discretize(s, pvmap, epmap, a, fsmap, sfmap);
  for (int i = nf; i < a->faces.size(); ++i)
    octree->insert(a->faces[i]);
}

void facesPI (Cspace *c, Polyhedron *a, CPVMap &pvmap, bool flag)
{
  CEdgeSet es;
  VVertices reg;
  for (CEdges::const_iterator e = c->edges.begin(); e != c->edges.end(); ++e)
    if ((*e)->HEdgesN() == 1 && (*e)->piEdge() &&
	(*e)->getT()->getP()->getA() == (flag ? Angle::ppi : Angle::mpi) &&
	es.find(*e) == es.end())
      reg.push_back(loopPI(a, pvmap, es, *e));
  Triangles tr;
  triangulate(reg, flag ? 3 : -3, tr);
  for (Triangles::iterator t = tr.begin(); t != tr.end(); ++t)
    a->addTriangle(t->a, t->b, t->c);
  deleteRegion(reg);
}

Vertices * loopPI (Polyhedron *a, CPVMap &pvmap, CEdgeSet &es, CEdge *e)
{
  Vertices *ve = new Vertices;
  CEdge *e0 = e;
  do {
    es.insert(e);
    ve->push_back(getVertex(e->getH()->getP(), pvmap, a));
    CVertex *t = e->getT();
    for (int i = 0; i < t->EdgesN(); ++i) {
      CEdge *f = t->getEdge(i);
      if (e != f && f->HEdgesN() == 1) {
	e = f;
	break;
      }
    }
  }
  while (e != e0);
  return ve;
}

// debug

PointCPoint * pcpt (Point *p)
{
  return dynamic_cast<PointCPoint *>(p);
}

InputAngle * ai (Angle *a)
{
  return dynamic_cast<InputAngle *>(a);
}

AngleCS * acs (Angle *a)
{
  return dynamic_cast<AngleCS *>(a);
}

AnglesEF * asef (AnglesCS *acs)
{
  return dynamic_cast<AnglesEF *>(acs);
}

AnglesSF * asff (AnglesCS *acs)
{
  return dynamic_cast<AnglesSF *>(acs);
}

SpiralPoint * spt (CPoint *p)
{
  return dynamic_cast<SpiralPoint *>(p);
}

LinePatchPoint * lppt (CPoint *p)
{
  return dynamic_cast<LinePatchPoint *>(p);
}

FFFPoint * fffp (CPoint *p)
{
  return dynamic_cast<FFFPoint *>(p);
}

CPointEdge * cpe (CPoint *p)
{
  return dynamic_cast<CPointEdge *>(p);
}

void pp1 (PV2 p)
{
  cerr << setprecision(16);
  cerr << "(" << p.x.mid() << " " << p.y.mid() << ")";
}

void pp (PV2 p)
{
  pp1(p);
  cerr << endl;
}

void pa1 (Angle *a)
{
  if (a == Angle::mpi)
    cerr << "mpi";
  else if (a == Angle::ppi)
    cerr << "ppi";
  else
    cerr << a->theta();
}

void pa (Angle *a)
{
  pa1(a);
  cerr << endl;
}

void pai1 (AngleInterval &ai)
{
  cerr << "(";
  pa1(ai.first);
  cerr << " ";
  pa1(ai.second);
  cerr << ")";
}

void pai (AngleInterval &ai)
{
  pai1(ai);
  cerr << endl;
}

void pais (AngleIntervals &ais)
{
  cerr << "(";
  for (AngleIntervals::iterator i = ais.begin(); i != ais.end(); ++i)
    pai1(*i);
  cerr << ")" << endl;
}

void pp1 (CPoint *p)
{
  PV2 xy = p->get();
  cerr << "(" << xy.x.mid() << " " << xy.y.mid() << " ";
  pa1(p->getA());
  cerr << ")";
}

void pp (CPoint *p)
{
  cerr << setprecision(16);
  pp1(p);
  cerr << endl;
}

void pv (CVertex *v)
{
  pp(v->getP());
}

void pvs (const CVertices &ve)
{
  cerr << "(" << endl;
  for (CVertices::const_iterator v = ve.begin(); v != ve.end(); ++v)
    pv(*v);
  cerr << ")" << endl;
}

void pe (CEdge *e)
{
  cerr << "(";
  pp1(e->getT()->getP());
  cerr << " ";
  pp1(e->getH()->getP());
  cerr << ")" << endl;
}

void pes (const CEdges &ed)
{
  cerr << "(";
  for (CEdges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e);
  cerr << ")" << endl;
}

void pe (HCEdge *e)
{
  cerr << "(";
  pp1(e->tail()->getP());
  cerr << " ";
  pp1(e->head()->getP());
  cerr << ")" << endl;
}

void pes (const HCEdges &ed)
{
  cerr << "(";
  for (HCEdges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e);
  cerr << ")" << endl;
}

void pedge (CEdge *e, double d, PV3s &pts)
{
  CPoints cpts;
  cpts.push_back(e->getT()->getP());
  if (!e->horizontal())
    discretize(e, e->getT()->getP()->getA(), e->getH()->getP()->getA(), d, cpts);
  cpts.push_back(e->getH()->getP());
  for (CPoints::iterator i = cpts.begin(); i != cpts.end(); ++i) {
    PointCPoint p(*i);
    pts.push_back(p.getP());
  }
}

void pp1 (PV3);

void ppts (const PV3s &pts)
{
  cerr << "(";
  for (PV3s::const_iterator p = pts.begin();  p != pts.end(); ++p)
    pp1(*p);
  cerr << ")" << endl;
}

void pe (CEdge *e, double d)
{
  PV3s pts;
  pedge(e, d, pts);
  ppts(pts);
}

void pes (const CEdges &ed, double d)
{
  cerr << "(";
  for (CEdges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e, d);
  cerr << ")" << endl;
}

void pe (CEdge *e, Angle *as, Angle *ae, double d)
{
  CPoints cpts;
  discretize(e, as, ae, d, cpts);
  PV3s pts;
  for (CPoints::iterator i = cpts.begin(); i != cpts.end(); ++i) {
    PointCPoint p(*i);
    pts.push_back(p.getP());
  }
  ppts(pts);
}

void pe (HCEdge *e, double d)
{
  PV3s pts;
  pedge(e->getE(), d, pts);
  if (!e->getForward())
    reverse(pts.begin(), pts.end());
  ppts(pts);
}

void pes (const HCEdges &ed, double d)
{
  cerr << "(";
  for (HCEdges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e, d);
  cerr << ")" << endl;
}

void pl (HCEdge *e)
{
  HCEdges ed;
  e->loop(ed);
  CVertices ve;
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    ve.push_back((*e)->head());
  pvs(ve);
}

void pl (HCEdge *e, double d)
{
  HCEdges ed;
  e->loop(ed);
  cerr << "(";
  for (HCEdges::iterator f = ed.begin(); f != ed.end(); ++f)
    pe(*f, d);
  cerr << ")" << endl;
}

void pf (CFace *f, Angle *a)
{
  PV2 p = f->leftEdge()->xy(a), q = f->rightEdge()->xy(a);
  cerr << "(";
  pp1(p);
  cerr << " ";
  pp1(q);
  cerr << ")" << endl;
}

HCEdge * nth (HCEdge *e, int n)
{
  for (int i = 0; i < n; ++i)
    e = e->getNext();
  return e;
}

bool find (Spiral *s, const SpiralSet &ss)
{
  return ss.find(s) != ss.end();
}

int find (CPoint *p, const CPoints &pts)
{
  for (int i = 0; i < pts.size(); ++i)
    if (p == pts[i])
      return i;
  return -1;
}

int find (CVertex *v, const CVertices &vertices)
{
  for (int i = 0; i < vertices.size(); ++i)
    if (v == vertices[i])
      return i;
  return -1;
}

int find (CPoint *p, const CVertices &vertices)
{
  for (int i = 0; i < vertices.size(); ++i)
    if (p == vertices[i]->getP())
      return i;
  return -1;
}

int find (CEdge *e, const CEdges &edges)
{
  for (int i = 0; i < edges.size(); ++i)
    if (e == edges[i])
      return i;
  return -1;
}

int find (CFace *f, const CFaces &faces)
{
  for (int i = 0; i < faces.size(); ++i)
    if (f == faces[i])
      return i;
  return -1;
}

int find (Slab *s, const Slabs &sl)
{
  for (int i = 0; i < sl.size(); ++i)
    if (s == sl[i])
      return i;
  return -1;
}

void plines (const vector<PV3s> &lines, int i);;

void pedge (CEdge *e, double d, int i)
{
  PV3s pts;
  pedge(e, d, pts);
  vector<PV3s> lines;
  lines.push_back(pts);
  plines(lines, i);
}

void pedges (const CEdges &edges, double d, int i)
{
  vector<PV3s> lines;
  for (CEdges::const_iterator e = edges.begin(); e != edges.end(); ++e) {
    PV3s pts;
    pedge(*e, d, pts);
    lines.push_back(pts);
  }
  plines(lines, i);
}

void pedges (const HCEdges &edges, double d, int i)
{
  vector<PV3s> lines;
  for (HCEdges::const_iterator e = edges.begin(); e != edges.end(); ++e) {
    PV3s pts;
    pedge((*e)->getE(), d, pts);
    lines.push_back(pts);
  }
  plines(lines, i);
}

void pboundary (CFace *f, double d, int i)
{
  HCEdges ed;
  f->boundaryHEdges(ed);
  pedges(ed, d, i);
}

void pfaces (const Faces &fa, int i);

void pfaces (const CFaces &fa, double d, int i)
{
  CPVMap pvmap;
  Polyhedron *a = discretize(fa, d, pvmap);
  pfaces(a->faces, i);
  delete a;
}

void pfaces (Cspace *c, double d, int i)
{
  pfaces(c->faces, d, i);
}

void pface (CFace *f, double d, int i)
{
  CFaces fa;
  fa.push_back(f);
  pfaces(fa, d, i);
}

void pfaces (const CFaces &fa, double d, int i, int is, int ie)
{
  CFaces fad;
  for (int j = is; j < ie; ++j)
    fad.push_back(fa[j]);
  pfaces(fad, d, i);
}

void pslab (Slab *s, const SFSMap &sfmap, int i)
{
  const Faces &fa = sfmap.find(s)->second;
  pfaces(fa, i);
}

void pfs (const Faces &fa);

void ps (Slab *s, const SFSMap &sfmap)
{
  const Faces &fa = sfmap.find(s)->second;
  pfs(fa);
}

double distance (CPoint *a, CPoint *b)
{
  PV2 u = a->getP() - b->getP();
  Parameter uu = u.dot(u);
  double da = fabs(a->getA()->theta() - b->getA()->theta());
  return uu.ub() + da;
}

double distance (CPoint *p, CFace *f)
{
  Angle *pa = p->getA();
  Spiral *l = f->leftEdge()->getSpiral(), *r = f->rightEdge()->getSpiral();
  PV2 a = p->getP(), t = l->xy(pa), h = r->xy(pa), u = (h - t).unit(), v = a - t;
  return u.cross(v).mid();
}

void findSpiral (Cspace *a, Spiral *s)
{
  for (int i = 0; i < a->faces.size(); ++i)
    if (a->faces[i]->leftEdge()->getSpiral() == s ||
	a->faces[i]->rightEdge()->getSpiral() == s)
      cerr << i << " ";
  cerr << endl;
}

void findAngle (Cspace *a, Angle *x)
{
  for (int i = 0; i < a->faces.size(); ++i)
    if (a->faces[i]->startAngle() == x || a->faces[i]->endAngle() == x)
      cerr << i << " ";
  cerr << endl;
}

class InputCPoint : public CPoint {
 public:
  InputCPoint (double x, double y, double th) : CPoint(new InputAngle(th)) {
    set(PV2::input(x, y));
  }
  ~InputCPoint () { delete a; }
};

Polyhedron * transform (Polyhedron *a, double x, double y, double th)
{
  InputCPoint p(x, y, th);
  return transform(a, &p);
}

void describe (Slab *s, EPMap &epmap)
{
  CPoints l = getPoints(s->l, s->s, s->e, epmap),
    r = getPoints(s->r, s->s, s->e, epmap);
  CTriangles tr = discretize(s, epmap);
  bool dummy = true;
}

Slab * getSlab (Face *f, const FSMap &fsmap)
{
  return fsmap.find(f)->second;
}

Primitive5(Intersects, CEdge *, el, CEdge *, er, CEdge *, fl, CEdge *, fr, Angle *, s);

int Intersects::sign ()
{
  PV2 a = el->xy(s), b = er->xy(s), c = fl->xy(s), d = fr->xy(s), u = b - a, v = d - c;
  int s1 = u.cross(c - a).sign(), s2 = u.cross(d - a).sign(),
    s3 = v.cross(a - c).sign(), s4 = v.cross(b - c).sign();
  return s1*s2 == -1 && s3*s4 == -1 ? 1 : -1;
}

bool intersects (CEdge *el, CEdge *er, CEdge *fl, CEdge *fr, double a)
{
  PTR<Angle> th = new InputAngle(a, false);
  return Intersects(el, er, fl, fr, th) == 1;
}
  
