#include "polyhedron.h"

double getTime ()
{
  timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1e-6*tv.tv_usec;
}

Parameter cross (const PV3 &a, const PV3 &b, int coord)
{
  int s = coord > 0 ? 1 : -1;
  switch (s*coord) {
  case 1:
    return s*(a.getY()*b.getZ() - a.getZ()*b.getY());
  case 2:
    return s*(a.getZ()*b.getX() - a.getX()*b.getZ());
  case 3:
    return s*(a.getX()*b.getY() - a.getY()*b.getX());
  default:
    assert(0);
  }
}

set<ID> intersection (const set<ID> &a, const set<ID> &b)
{
  set<ID> ab;
  set_intersection(a.begin(), a.end(), b.begin(), b.end(),
		   inserter(ab, ab.begin()));
  return ab;
}

void Point::getBBox (double *bbox)
{
  PV3 p = getApprox(1.0);
  bbox[0] = p.x.lb();
  bbox[1] = p.x.ub();
  bbox[2] = p.y.lb();
  bbox[3] = p.y.ub();
  bbox[4] = p.z.lb();
  bbox[5] = p.z.ub();
}

bool Point::identical (Point *a)
{
  return order(a) == 0;
}

int Point::order (Point *a)
{
  if (this == a)
    return 0;
  SumPoint *s1 = dynamic_cast<SumPoint *>(this);
  if (s1) {
    SumPoint *s2 = dynamic_cast<SumPoint *>(a);
    if (s2 &&
	(s1->a->identicalI(s2->a) && s1->b->identicalI(s2->b) ||
	 s1->a->identicalI(s2->b) && s1->b->identicalI(s2->a)))
      return 0;
  }
  return PointOrderR(this, a);
}

bool Point::identicalI (Point *a)
{
  if (!(input() && a->input()))
    return false;
  PV3 p = getApprox(1.0), q = a->getApprox(1.0);
  return p.x.lb() == q.x.lb() && p.y.lb() == q.y.lb() &&
    p.z.lb() == q.z.lb();
}

bool Point::onLine (Point *a, Point *b)
{
  return this == a || this == b || OnLine(this, a, b) == 0;
}

int Point::side (Plane *a)
{
  if (ps.find(a->id) != ps.end())
    return 0;
  int s = Side(a, this);
  if (s == 0) 
    ps.insert(a->id);
  return s;
}

PTR<Point> Rdir = new Point(0.8401877171547095, 0.394382926819093, 0.7830992237586059);

bool onEdge (Point *a, Point *t, Point *h, bool strict)
{
  if (a->identical(t) || a->identical(h))
    return !strict;
  return Order(a, t, h) == 1 && Order(a, h, t) == 1;
}

bool closerPair (Point *a, Point *b, Point *c, Point *d)
{
  if (a == c && b == d || a == d && b == c)
    return false;
  return CloserPair(a, b, c, d) == 1;
}

ID Plane::planeid = 1u;

double bboxSize (double *bb)
{
  return max(bb[1] - bb[0], max(bb[3] - bb[2], bb[5] - bb[4]));
}

void copyBBox (const double *bbf, double *bbt)
{
  for (int i = 0; i < 6; ++i)
    bbt[i] = bbf[i];
}

void mergeBBox (const double *bbf, double *bbt)
{
  for (int i = 0; i < 3; ++i) {
    bbt[2*i] = min(bbt[2*i], bbf[2*i]);
    bbt[2*i+1] = max(bbt[2*i+1], bbf[2*i+1]);
  }
}

bool bboxOverlap (const double *a, const double *b, double s)
{
  for (int i = 0; i < 3; ++i)
    if (a[2*i+1] + s < b[2*i] || b[2*i+1] + s < a[2*i])
      return false;
  return true;
}

bool bboxOverlap (Point *a, const double *bbox)
{
  double abox[6];
  a->getBBox(abox);
  return bboxOverlap(abox, bbox);
}

Vertex::Vertex (Point *p, bool perturbed) : p(p), node(0)
{
  p->getBBox(bbox);
  if (!perturbed) {
    Parameter k = Rdir->getApprox(1.0).dot(p->getApprox(1.0));
    rint[0] = k.lb();
    rint[1] = k.ub();
  }
}

HEdges Vertex::outgoingHEdges () const
{
  HEdges ed;
  for (Edges::const_iterator e = edges.begin(); e != edges.end(); ++e)
    for (HEdges::iterator f = (*e)->hedges.begin(); f != (*e)->hedges.end(); ++f)
      if ((*f)->tail() == this)
	ed.push_back(*f);
  return ed;
}

Faces Vertex::incidentFaces () const
{
  Faces fa;
  for (Edges::const_iterator e = edges.begin(); e != edges.end(); ++e)
    for (HEdges::iterator f = (*e)->hedges.begin(); f != (*e)->hedges.end(); ++f)
      if ((*f)->tail() == this)
	fa.push_back((*f)->getF());
  return fa;
}

HEdge * Vertex::connected (Vertex *a) const
{
  HEdges ed = outgoingHEdges();
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    if ((*e)->head() == a)
      return *e;
  return 0;
}

int Vertex::order (Vertex *v) const
{
  if (rint[1] < v->rint[0]) return 1;
  if (v->rint[1] < rint[0]) return -1;
  return p->order(v->p);
}

Vertex * HEdge::tail () const
{
  return forward ? e->t : e->h;
}

Vertex * HEdge::head () const
{
  return forward ? e->h: e->t;
}

PV3 HEdge::getU ()
{
  return forward ? e->getU() : - e->getU();
}

PV3 HEdge::getN ()
{
  return forward ? f->p->get().n : - f->p->get().n; 
}

HEdge * HEdge::cw () const
{
  int n = e->hedges.size();
  for (int i = 0; i < n; ++i)
    if (e->hedges[i] == this)
      return e->hedges[(i+1)%n];
  return 0;
}

HEdge * HEdge::ccw () const
{
  int n = e->hedges.size();
  for (int i = 0; i < n; ++i)
    if (e->hedges[i] == this)
      return e->hedges[i == 0 ? n - 1 : i - 1];
  return 0;
}

Vertices HEdge::loop () const
{
  Vertices ve;
  const HEdge *e = this;
  do {
    ve.push_back(e->tail());
    e = e->next;
  }
  while (e != this);
  return ve;
}

Points HEdge::pointLoop () const
{
  Vertices ve = loop();
  Points pts;
  for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
    pts.push_back((*v)->p);
  return pts;
}

HEdges HEdge::edgeLoop ()
{
  HEdges ed;
  HEdge *e = this;
  do {
    ed.push_back(e);
    e = e->next;
  }
  while (e != this);
  return ed;
}

Edge::Edge (Vertex *t, Vertex *h) : t(t), h(h)
{ 
  hedges.reserve(2);
  setBBox();
}

void Edge::setBBox ()
{
  copyBBox(t->bbox, bbox);
  mergeBBox(h->bbox, bbox);
}

Edge::~Edge ()
{
  for (HEdges::iterator h = hedges.begin(); h != hedges.end(); ++h)
    delete *h;
}

HEdge * Edge::addHEdge (bool forward)
{
  HEdge *e = new HEdge(this, forward);
  hedges.push_back(e);
  return e;
}

void Edge::removeHEdge (HEdge *e)
{
  HEdges::iterator j = remove(hedges.begin(), hedges.end(), e);
  hedges.erase(j, hedges.end());
  delete e;
}

void Edge::sortHEdges ()
{
  if (hedges.size() > 2)
    sort(hedges.begin(), hedges.end(), EdgeOrder(this));
}

bool HFace::pos () const 
{ 
  return this == f->hfaces; 
}
  
HFace * HFace::twin () const 
{
  return this == f->hfaces ? f->hfaces + 1 : f->hfaces;
}

PV3 HFace::getN ()
{
  return pos() ? f->p->get().n : - f->p->get().n;
}

HFaces HFace::neighbors () const
{
  HEdges ed = f->boundaryHEdges();
  HFaces hf;
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    hf.push_back(neighbor(*e));
  return hf;
}

HFace * HFace::neighbor (HEdge *e) const
{
  bool ccw = pos() == e->forward;
  HEdge *f = ccw ? e->ccw() : e->cw();
  bool flag = e->forward == f->forward ? pos() : !pos();
  Face *g = f->f;
  return flag ? g->hfaces + 1 : g->hfaces;
}

Face::Face (HEdge *h, int pc) : h(h), pc(pc)
{
  p = new TrianglePlane(h->tail()->p, h->head()->p, h->next->head()->p);
  hfaces[0].f = hfaces[1].f = this;
  addLoop(h, true);
}

Face::Face (Point *a, Point *b, Point *c) : h(0), pc(0)
{
  p = new TrianglePlane(a, b, c);
  hfaces[0].f = hfaces[1].f = this;
  a->getBBox(bbox);
  double bb[6];
  b->getBBox(bb);
  mergeBBox(bb, bbox);
  c->getBBox(bb);
  mergeBBox(bb, bbox);
}

Face::Face (HEdge *h, TrianglePlane *p, int pc) : h(h), p(p), pc(pc)
{
  hfaces[0].f = hfaces[1].f = this;
}

void Face::addLoop (HEdge *h, bool flag)
{
  HEdges ed = h->edgeLoop();
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e) {
    (*e)->f = this;
    (*e)->tail()->p->ps.insert(p->id);
  }
  if (flag) {
    copyBBox(ed[0]->tail()->bbox, bbox);
    for (int i = 1; i < ed.size(); ++i)
      mergeBBox(ed[i]->tail()->bbox, bbox);
  }
}

void Face::update ()
{
  Vertices ve = h->loop();
  p = new TrianglePlane(ve[0]->p, ve[1]->p, ve[2]->p);
  copyBBox(ve[0]->bbox, bbox);
  for (int i = 1; i < ve.size(); ++i)
    mergeBBox(ve[i]->bbox, bbox);
  pc = 0;
}

int Face::getPC () {
  if (pc == 0)
    pc = ProjectionCoordinate(p);
  return pc;
}

bool Face::boundaryVertex (Point *a) const
{
  Vertices ve = h->loop();
  for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
    if ((*v)->p->identical(a))
      return true;
  return false;
}

bool Face::boundaryVertex (Vertex *v) const
{
  Vertices ve = h->loop();
  return find(ve.begin(), ve.end(), v) != ve.end();
}

bool Face::boundaryEdge (Edge *e) const
{
  for (HEdges::iterator f = e->hedges.begin(); f != e->hedges.end(); ++f)
    if ((*f)->f == this)
      return true;
  return false;
}

bool Face::sharedEdge (Face *f) const
{
  if (h && f->h) {
    HEdges ed = h->edgeLoop();
    for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
      for (HEdges::iterator i = (*e)->e->hedges.begin(); 
	   i != (*e)->e->hedges.end(); ++i)
	if ((*i)->f == f)
	  return true;
    return false;
  }
  Points p1 = boundaryPoints(), p2 = f->boundaryPoints();
  int n = 0;
  for (Points::iterator p = p1.begin(); p != p1.end(); ++p)
    if (find(p2.begin(), p2.end(), *p) != p2.end())
      ++n;
  return n > 1;
}

PTR<Point> Face::sharedVertex (Face *f) const
{
  Vertices ve = h->loop(), vef = f->h->loop();
  for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
    if (find(vef.begin(), vef.end(), *v) != vef.end())
      return (*v)->p;
  return 0;
}

Points Face::boundaryPoints () const
{
  Points pts;
  if (h) {
    Vertices ve = h->loop();
    for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
      pts.push_back((*v)->p);
  }
  else {
    pts.push_back(p->getA());
    pts.push_back(p->getB());
    pts.push_back(p->getC());
  }
  return pts;
}

bool Face::coplanar (Face *f)
{
  if (p->id == f->p->id)
    return true;
  return p->getA()->side(f->p) == 0 && p->getB()->side(f->p) == 0 &&
    p->getC()->side(f->p) == 0;
};

bool Face::intersects (Face *g, bool strict)
{
  if (coplanar(g))
    return !strict && intersectsFP(g);
  if (sharedEdge(g))
    return false;
  int sf[4], sg[4];
  if (intersectsFP(g, sg) || g->intersectsFP(this, sf) ||
      verifyFP(g, sg) || g->verifyFP(this, sf))
    return false;
  return intersectsFE(g, sg, strict) || g->intersectsFE(this, sf, strict);  
}

bool Face::intersectsFP (Face *g, int *sg)
{
  PV3 nf = getP()->getApprox(1.0).n;
  Parameter kf = getP()->getApprox(1.0).k;
  Points pf = boundaryPoints(), pg = g->boundaryPoints();
  int n = pg.size();
  for (int i = 0; i < n; ++i)
    if (find(pf.begin(), pf.end(), pg[i]) != pf.end())
      sg[i] = 2;
    else {
      int s = (nf.dot(pg[i]->getApprox(1.0)) + kf).sign(false);
      sg[i] = s ? s : 3;
    }
  return checkFP(sg, n);
}

bool Face::checkFP (int *s, int n) const
{
  bool pflag = false, nflag = false;
  int nv = 0;
  for (int i = 0; i < n; ++i)
    if (s[i] == 0 || s[i] == 3)
      return false;
    else if (s[i] == 1)
      pflag = true;
    else if (s[i] == -1)
      nflag = true;
    else
      ++nv;
  return nv < 2 && !(pflag && nflag);
}

bool Face::verifyFP (Face *g, int *sg)
{
  Points pg = g->boundaryPoints();
  int n = pg.size();
  for (int i = 0; i < n; ++i)
    if (sg[i] == 3)
      sg[i] = pg[i]->side(getP());
  return checkFP(sg, n);
}

bool Face::intersectsFE (Face *g, int *sg, bool strict)
{
  Points pg = g->boundaryPoints();
  int n = pg.size();
  if (!strict) {
    for (int i = 0; i < n; ++i)
      if ((sg[i] == 0 || sg[i] == 2) && (sg[(i+1)%n] == 0 || sg[(i+1)%n] == 2))
	return intersectsFEP(pg[i], pg[(i+1)%n], strict);
    for (int i = 0; i < n; ++i)
      if (sg[i] == 0 && contains(pg[i], strict))
	return true;
  }
  for (int i = 0; i < n; ++i)
    if (sg[i] && sg[i] == - sg[(i+1)%n]) {
      PTR<Point> p = new EPPoint(pg[i], pg[(i+1)%n], getP());
      if (bboxOverlap(p, getBBox()) && contains(p, strict))
	return true;
    }
  return false;
}

bool Face::intersectsFEP (Point *et, Point *eh, bool strict)
{
  if (contains(et, strict) || contains(eh, strict))
    return true;
  Points pf = boundaryPoints();
  int n = pf.size();
  for (int i = 0; i < n; ++i)
    if (intersectsEE(pf[i], pf[(i+1)%n], et, eh, strict))
      return true;
  return false;
}

bool Face::intersectsEE (Point *et, Point *eh, Point *ft, Point *fh, bool strict)
{
  if (et == ft || et == fh || eh == ft || eh == fh)
    return false;
  int c = getPC(), tp1 = LeftTurn(et, ft, fh, c);
  if (!strict && tp1 == 0 && onEdge(et, ft, fh, true))
    return true;
  int tp2 = LeftTurn(eh, ft, fh, c);
  if (!strict && tp2 == 0 && onEdge(eh, ft, fh, true))
    return true;
  if (tp1*tp2 > -1)
    return false;
  int tp3 = LeftTurn(ft, et, eh, c);
  if (!strict && tp3 == 0 && onEdge(ft, et, eh, true))
    return true;
  int tp4 = LeftTurn(fh, et, eh, c);
  return !strict && tp4 == 0 && onEdge(fh, et, eh, true) ||
    tp3*tp4 == -1;
}

bool Face::intersectsFP (Face *f)
{
  Points pf = f->boundaryPoints();
  int n = pf.size();
  for (int i = 0; i < n; ++i)
    if (intersectsFEP(pf[i], pf[(i+1)%n], false))
      return true;
  Points pts = boundaryPoints();
  int m = pts.size();
  for (int i = 0; i < m; ++i)
    if (f->intersectsFEP(pts[i], pts[(i+1)%m], false))
      return true;
  return false;
}

bool Face::contains (Point *a, bool strict, int *ie)
{
  return bboxOverlap(a, bbox) &&
    ::contains(boundaryPoints(), getPC(), a, strict, ie);
}

PTR<Point> Face::centroid () const
{
  PTR<Point> a = h->tail()->p, b = h->head()->p, c = h->next->head()->p;
  return new CentroidPoint(a, b, c);
}

PTR<Point> Face::rayIntersection (Point *a, Point *r)
{
  if (a->side(p) == 0)
    return contains(a, false) ? a : 0;
  if (PlaneRayAlignment(p, r) == 0)
    return 0;
  PTR<Point> q = new RayPlanePoint(a, r, p);
  if (contains(q, false) && PointOrder(a, q, r) == 1)
    return q;
  return 0;
}

void Face::triangulate (Triangles &tr)
{
  Vertices ve = h->loop();
  int n = ve.size();
  for (int i = 1; i + 1 < n; ++i)
    tr.push_back(Triangle(ve[0]->getP(), ve[i]->getP(), ve[i+1]->getP()));
}

bool contains (const Points &pts, int c, Point *a, bool strict, int *ie)
{
  int n = pts.size();
  Point *et = pts.back();
  for (int i = 0; i < n; ++i) {
    Point *eh = pts[i];
    if (a == et || a == eh)
      return false;
    int s = LeftTurn(a, et, eh, c);
    if (s == -1)
      return false;
    if (s == 0)
      if (strict)
	return false;
      else if (onEdge(a, et, eh, true)) {
	if (ie)
	  *ie = i == 0 ? n - 1 : i - 1;
      	return true;
      }
      else
	return false;
    et = eh;
  }
  return true;
}

Shell::~Shell ()
{
  if (octreef)
    delete octreef;
}

Shell::Shell (const HFaces &hf) : hfaces(hf), c(0)
{
  for (HFaces::iterator h = hfaces.begin(); h != hfaces.end(); ++h)
    (*h)->s = this;
  setBBox();
  setOctree();
}

void Shell::setBBox ()
{
  bbox[0] = bbox[2] = bbox[4] = 1e20;
  bbox[1] = bbox[3] = bbox[5] = -1e20;
  for (HFaces::iterator f = hfaces.begin(); f != hfaces.end(); ++f)
    mergeBBox((*f)->f->bbox, bbox);
}

void Shell::setOctree ()
{
  Faces fa;
  for (HFaces::iterator f = hfaces.begin(); f != hfaces.end(); ++f)
    fa.push_back((*f)->f);
  octreef = Octree<Face *>::octree(fa, bbox);
}

bool Shell::outer () const
{
  PTR<Point> r = new Point(0.0, 0.0, 1.0);
  Vertex *vm = vmax(r);
  HEdge *em = 0;
  HFace *fm = 0;
  for (Edges::iterator e = vm->edges.begin(); e != vm->edges.end(); ++e)
    for (HEdges::iterator h = (*e)->hedges.begin(); h != (*e)->hedges.end(); ++h)
      for (int i = 0; i < 2; ++i)
	if ((*h)->f->hfaces[i].s == this) {
	  if (!em || em->e != *e && SlopeOrder(*e, em->e, r) == 1) {
	    HEdge *nem = *h;
	    HFace *nfm = nem->f->hfaces + i;
	    if (!nfm->getF()->coplanar(nfm->neighbor(nem)->getF())) {
	      em = nem;
	      fm = nfm;
	      }
	  }
	  break;
	}
  return Convex(em, fm) == 1;
}

Vertex * Shell::vmax (Point *r) const
{
  Vertex *v = hfaces[0]->f->h->tail();
  while (true) {
    Vertex *w = v;
    for (Edges::iterator e = v->edges.begin(); e != v->edges.end(); ++e)
      for (HEdges::iterator h = (*e)->hedges.begin(); h != (*e)->hedges.end(); ++h)
	if ((*h)->f->hfaces[0].s == this || (*h)->f->hfaces[1].s == this) {
	  Vertex *nw = vmax((*h)->f, r);
	  if (w != nw && PointOrder(w->p, nw->p, r) == 1)
	    w = nw;
	}
    if (v == w) {
      double rb[6];
      rayBBox(w->p, r, rb);
      Faces fa;
      octreef->find(rb, fa);
      for (Faces::iterator f = fa.begin(); v == w && f != fa.end(); ++f)
	if ((*f)->rayIntersection(w->p, r) != 0)
	  w = vmax(*f, r);
    }
    if (v == w)
      break;
    v = w;
  }
  return v;
}

Vertex * Shell::vmax (Face *f, Point *r) const
{
  Vertices ve = f->h->loop();
  Vertex *v = ve[0];
  for (int i = 1; i < ve.size(); ++i)
    if (PointOrder(v->p, ve[i]->p, r) == 1)
      v = ve[i];
  return v;
}

bool Shell::contains (Shell *s) const
{
  if (!bboxOverlap(bbox, s->bbox))
    return false;
  for (HFaces::const_iterator h = s->hfaces.begin(); h != s->hfaces.end(); ++h) {
    Vertices ve = (*h)->f->h->loop();
    for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v) {
      int i = contains((*v)->p);
      if (i == 1)
	return true;
      if (i == -1)
	return false;
    }
  }
  return false;
}

int Shell::contains (Point *a) const
{
  if (!bboxOverlap(a, bbox))
    return -1;
  PTR<Point> r = new Point(0.0, 0.0, 1.0);
  double rb[6];
  rayBBox(a, r, rb);
  Faces fa;
  octreef->find(rb, fa);
  bool res = false;
  for (Faces::iterator f = fa.begin(); f != fa.end(); ++f)
    if ((*f)->boundaryVertex(a) ||
	a->side((*f)->getP()) == 0 && (*f)->contains(a, false))
      return 0;
    else if ((*f)->rayIntersection(a, r) != 0)
      res = !res;
  return res ? 1 : -1;
}

void Shell::rayBBox (Point *a, Point *r, double *rb) const
{
  a->getBBox(rb);
  RayZPlanePoint q(a, r, bbox[5]);
  double qb[6];
  q.getBBox(qb);
  mergeBBox(qb, rb);
}

int Shell::euler () const
{
  set<Vertex *> vs;
  set<Edge *> es;
  set<Face *> fs;
  for (HFaces::const_iterator f = hfaces.begin(); f != hfaces.end(); ++f)
    fs.insert((*f)->f);
  for (set<Face *>::iterator f = fs.begin(); f != fs.end(); ++f) {
    HEdge *e = (*f)->h;
    do {
      vs.insert(e->tail());
      es.insert(e->e);
      e = e->next;
    }
    while (e != (*f)->h);
  }
  int nv = vs.size(), ne = es.size(), nf = fs.size();
  return nv - ne + nf;
}

bool Shell::pos () const
{
  for (HFaces::const_iterator h = hfaces.begin(); h != hfaces.end(); ++h)
    if (!(*h)->pos())
      return false;
  return true;
}

void deleteShells (const Shells &sh)
{
  for (Shells::const_iterator s = sh.begin(); s != sh.end(); ++s)
    delete *s;
}

bool Cell::contains (Point *p) const
{
  if (outer && outer->contains(p) < 1)
    return false;
  for (Shells::const_iterator s = inner.begin(); s != inner.end(); ++s)
    if ((*s)->contains(p) > -1)
      return false;
  return true;
}

PTR<Point> Cell::interiorPoint () const
{
  Shell *s = getShell(0);
  for (HFaces::iterator hf = s->hfaces.begin(); hf != s->hfaces.end(); ++hf) {
    Face *f = (*hf)->f;
    PTR<Point> p = f->centroid(), n = new HFaceNormal(*hf), qmin = 0;
    for (int i = 0; i < nShells(); ++i) {
      Shell *s = getShell(i);
      for (HFaces::iterator h = s->hfaces.begin(); h != s->hfaces.end(); ++h) {
	Face *g = (*h)->f;
	if (g != f) {
	  PTR<Point> q = g->rayIntersection(p, n);
	  if (q && p != q && (!qmin || CloserPair(p, q, p, qmin) == 1))
	    qmin = q;
	}
      }
    }
    if (qmin)
      return new CentroidPoint(p, qmin);
  }
  return 0;
}

Polyhedron::~Polyhedron ()
{
  for (Vertices::iterator v = vertices.begin(); v != vertices.end(); ++v)
    delete *v;
  for (Edges::iterator e = edges.begin(); e != edges.end(); ++e)
    delete *e;
  for (Faces::iterator f = faces.begin(); f != faces.end(); ++f)
    delete *f;
  for (Cells::iterator c = cells.begin(); c != cells.end(); ++c)
    delete *c;
}

bool Polyhedron::findPoint (Point *p) const
{
  if (perturbed)
    return false;
  Vertex v(p, false);
  return vtree.find(&v);
}

Vertex * Polyhedron::getVertex (Point *p)
{
  Vertex *v = new Vertex(p, perturbed);
  if (!perturbed) {
    RBTree<Vertex *>::Node *n = vtree.insert(v);
    if (n->v != v) {
      delete v;
      return n->v;
    }
    v->node = n;
  }
  if (vertices.empty())
    copyBBox(v->bbox, bbox);
  else
    mergeBBox(v->bbox, bbox);
  vertices.push_back(v);
  return v;
}

Vertex * Polyhedron::getVertex (Point *p, PVMap &pvmap)
{
  PVMap::iterator iter = pvmap.find(p);
  if (iter != pvmap.end())
    return iter->second;
  Vertex *w = getVertex(p);
  pvmap.insert(PVPair(p, w));
  return w;
}

Edge * Polyhedron::getEdge (Vertex *a, Vertex *b)
{
  for (Edges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    if ((*e)->t == a && (*e)->h == b || (*e)->t == b && (*e)->h == a)
      return *e;
  Edge *e = new Edge(a, b);
  edges.push_back(e);
  a->edges.push_back(e);
  b->edges.push_back(e);
  return e;
}

HEdge * Polyhedron::addHEdge (Vertex *a, Vertex *b)
{
  Edge *e = getEdge(a, b);
  bool forward = e->t == a;
  return e->addHEdge(forward);
}

HEdge * Polyhedron::getHEdge (Vertex *a, Vertex *b)
{
  Edge *e = getEdge(a, b);
  bool forward = e->t == a;
  for (HEdges::iterator f = e->hedges.begin(); f != e->hedges.end(); ++f)
    if ((*f)->forward == forward && !(*f)->f)
      return *f;
  return e->addHEdge(forward);
}

Face * Polyhedron::addTriangle (Vertex *a, Vertex *b, Vertex *c, int pc)
{
  HEdge *u = addHEdge(a, b), *v = addHEdge(b, c), *w = addHEdge(c, a);
  u->setNext(v);
  v->setNext(w);
  w->setNext(u);
  Face *f = new Face(u, pc);
  faces.push_back(f);
  return f;
}

Face * Polyhedron::addRectangle (Vertex *a, Vertex *b, Vertex *c, Vertex *d)
{
  HEdge *u = addHEdge(a, b), *v = addHEdge(b, c), *w = addHEdge(c, d), 
    *x = addHEdge(d, a);
  u->setNext(v);
  v->setNext(w);
  w->setNext(x);
  x->setNext(u);
  Face *f = new Face(u, 0);
  faces.push_back(f);
  return f;
}

Face * Polyhedron::addFace (HEdge *h, int pc)
{
  Face *f = new Face(h, pc);
  faces.push_back(f);
  return f;
}

void Polyhedron::formCells () {
  Shells shells;
  formShells(shells);
  formCellsAux(shells);
}

void Polyhedron::formCellsAux (const Shells &shells)
{
  Shells inner;
  cells.push_back(new Cell(0));
  for (Shells::const_iterator s = shells.begin(); s != shells.end(); ++s) {
    (*s)->setBBox();
    (*s)->setOctree();
    if ((*s)->outer())
      cells.push_back(new Cell(*s));
    else
      inner.push_back(*s);
  }
  Octree<Cell *> *octreec = cellOctree();
  for (Shells::iterator s = inner.begin(); s != inner.end(); ++s)
    enclosingCell(*s, octreec)->addInner(*s);
  delete octreec;
}

void Polyhedron::formShells (Shells &shells)
{
  for (Edges::iterator e = edges.begin(); e != edges.end(); ++e)
    (*e)->sortHEdges();
  for (Faces::iterator f = faces.begin(); f != faces.end(); ++f)
    for (int i = 0; i < 2; ++i)
      if (!(*f)->hfaces[i].s)
	shells.push_back(formShell((*f)->hfaces + i));
}

Shell * Polyhedron::formShell (HFace *f) const
{
  Shell *s = new Shell;
  HFaces st;
  st.push_back(f);
  while (!st.empty()) {
    HFace *g = *(st.end()-1);
    st.pop_back();
    if (!g->s) {
      g->s = s;
      s->hfaces.push_back(g);
      HFaces gn = g->neighbors();
      st.insert(st.end(), gn.begin(), gn.end());
    }
  }
  return s;
}

Cell * Polyhedron::enclosingCell (Shell *s, Octree<Cell *> *octreec) const
{
  Cells ce;
  octreec->find(s->hfaces[0]->f->h->tail()->getBBox(), ce);
  Cell *c = 0;
  for (Cells::iterator d = ce.begin(); d != ce.end(); ++d)
    if ((*d)->outer->contains(s) &&
	(!c || c->outer->contains((*d)->outer)))
      c = *d;
  return c ? c : cells[0];
}

void Polyhedron::clearCells ()
{
  for (Faces::iterator f = faces.begin(); f != faces.end(); ++f)
    for (int i = 0; i < 2; ++i)
      (*f)->hfaces[i].s = 0;
  for (Cells::iterator c = cells.begin(); c != cells.end(); ++c)
    delete *c;
  cells.clear();
}

Face * Polyhedron::addTriangle (PTR<Point> a, PTR<Point> b, PTR<Point> c,
				PVMap &pvmap, int pc)
{
  Vertex *u = getVertex(a, pvmap), *v = getVertex(b, pvmap),
    *w = getVertex(c, pvmap);
  return addTriangle(u, v, w, pc);
}

Face * Polyhedron::addTriangle (Face *f, PVMap &pvmap)
{
  Points pf = f->h->pointLoop();
  int pc = f->getPC();
  return addTriangle(pf[0], pf[1], pf[2], pvmap, pc);
}

Polyhedron * Polyhedron::copy () const
{
  Polyhedron *a = new Polyhedron;
  PVMap pvmap;
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f)
    a->addTriangle(*f, pvmap);
  return a;
}

Polyhedron * Polyhedron::scale (double unit) const
{
  Polyhedron *a = new Polyhedron;
  PVMap pvmap;
  for (Vertices::const_iterator v = vertices.begin(); v != vertices.end(); ++v)
    pvmap.insert(PVPair((*v)->p, a->getVertex(new ScalePoint((*v)->p, unit))));
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f)
    a->addTriangle(*f, pvmap);
  return a;
}

Polyhedron * Polyhedron::negative () const
{
  Polyhedron *a = new Polyhedron(perturbed);
  PVMap pvmap;
  for (Vertices::const_iterator v = vertices.begin(); v != vertices.end(); ++v)
    pvmap.insert(PVPair((*v)->p, a->getVertex(new NegPoint((*v)->p))));
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f) {
    Points pf = (*f)->h->pointLoop();
    a->addTriangle(pf[2], pf[1], pf[0], pvmap);
  }
  return a;
}

Polyhedron * Polyhedron::translate (Point *t) const
{
  Polyhedron *a = new Polyhedron(perturbed);
  PVMap pvmap;
  for (Vertices::const_iterator v = vertices.begin(); v != vertices.end(); ++v)
    pvmap.insert(PVPair((*v)->p, a->getVertex(new SumPoint(t, (*v)->p))));
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f)
    a->addTriangle(*f, pvmap);
  return a;
}

Polyhedron * Polyhedron::negativeTranslate (Point *t) const
{
  Polyhedron *a = new Polyhedron(perturbed);
  PVMap pvmap;
  for (Vertices::const_iterator v = vertices.begin(); v != vertices.end(); ++v)
    pvmap.insert(PVPair((*v)->p, a->getVertex(new DiffPoint(t, (*v)->p))));
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f) {
    Points pf = (*f)->h->pointLoop();
    a->addTriangle(pf[2], pf[1], pf[0], pvmap);
  }
  return a;
}

bool Polyhedron::intersects (Polyhedron *a, bool strict) const
{
  return contains(a->vertices[0]->p) || a->contains(vertices[0]->p) ||
    intersectsEdges(a, strict) || a->intersectsEdges(this, strict);
}

bool Polyhedron::contains (Point *p) const
{
  if (!bboxOverlap(p, bbox))
    return false;
  for (int i = 1; i < cells.size(); ++i)
    if (cells[i]->wn == 1 && cells[i]->contains(p))
      return true;
  return false;
}

int Polyhedron::containingCell (Point *p) const
{
  for (int i = 0; i< cells.size(); i++)
    if (cells[i]->contains(p))
      return i;
  return -1;
}

bool Polyhedron::intersectsEdges (const Polyhedron *a, bool strict) const
{
  Octree<Face *> *octree = a->faceOctree();
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f) {
    Faces fa;
    octree->find((*f)->bbox, fa);
    for (Faces::iterator g = fa.begin(); g != fa.end(); ++g)
      if ((*f)->intersects(*g, strict)) {
	delete octree;
	return true;
      }
  }
  delete octree;
  return false;
}

Polyhedron * Polyhedron::boolean (Polyhedron *a, SetOp op)
{
  Polyhedron *poly[] = {this, a},
    *c = overlay(poly, 2), *d = new Polyhedron(false);
  PVMap pvmap;
  computeWindingNumbers();
  a->computeWindingNumbers();
  c->formCells();
  set<Cell *> cin;
  for (int i = 1; i < c->cells.size(); ++i) {
    Cell *ci = c->cells[i];
    PTR<Point> p = ci->interiorPoint();
    bool ina = contains(p), inb = a->contains(p);
    if (inSet(ina, inb, op))
      cin.insert(ci);
  }
  for (Faces::iterator f = c->faces.begin(); f != c->faces.end(); ++f) {
    bool in1 = cin.find((*f)->hfaces[0].s->c) != cin.end(),
      in2 = cin.find((*f)->hfaces[1].s->c) != cin.end();
    if (in1 != in2) {
      Points pf = (*f)->h->pointLoop();
      int pc = (*f)->getPC();
      if (in2)
	d->addTriangle(pf[0], pf[1], pf[2], pvmap, pc);
      else
	d->addTriangle(pf[2], pf[1], pf[0], pvmap, - pc);
    }
  }
  delete c;
  return d;
}

// used in simplify from here on

Polyhedron * Polyhedron::cellPolyhedron (int i) const
{
  Cell *c = cells[i];
  Polyhedron *a = new Polyhedron(perturbed);
  PVMap pvmap;
  if (c->outer)
    a->addHFaces(c->outer->hfaces, pvmap);
  for (Shells::const_iterator s = c->inner.begin(); s != c->inner.end(); ++s)
    a->addHFaces((*s)->hfaces, pvmap);
  a->formCells();
  return a;  
}

void Polyhedron::addHFaces (const HFaces &hf, PVMap &pvmap)
{
  for (HFaces::const_iterator h = hf.begin(); h != hf.end(); ++h)
    addTriangle((*h)->f, pvmap);
}

void Polyhedron::replaceVertex (Face *f, Vertex *v, Vertex *w)
{
  HEdge *e = f->h;
  while (e->next->head() != v)
    e = e->next;
  removeHEdge(e->next->next);
  removeHEdge(e->next);
  Vertex *t = e->tail(), *h = e->head();
  HEdge *en = getHEdge(h, w), *enn = getHEdge(w, t);
  en->f = enn->f = f;
  e->next = en;
  en->next = enn;
  enn->next = e;
  f->h = e;
  f->p = new TrianglePlane(t->p, h->p, w->p);
  f->pc = 0;
  copyBBox(e->e->bbox, f->bbox);
  mergeBBox(en->e->bbox, f->bbox);
  mergeBBox(enn->e->bbox, f->bbox);
}

void Polyhedron::removeLoop (HEdge *e)
{
  if (e->f)
    e->f->h = 0;
  HEdges ed = e->edgeLoop();
  for (HEdges::iterator h = ed.begin(); h != ed.end(); ++h)
    removeHEdge(*h);
}

void Polyhedron::removeHEdge (HEdge *h)
{
  Edge *e = h->e;
  e->removeHEdge(h);
  if (e->hedges.empty()) {
    Edges::iterator k = remove(e->t->edges.begin(), e->t->edges.end(), e);
    e->t->edges.erase(k, e->t->edges.end());
    k = remove(e->h->edges.begin(), e->h->edges.end(), e);
    e->h->edges.erase(k, e->h->edges.end());
  }
}

void Polyhedron::moveVertex (Vertex *v, Point *p)
{
  v->p = p;
  if (!perturbed) {
    vtree.remove(v->node);
    v->node = vtree.insert(v);
  }
  p->getBBox(v->bbox);
  mergeBBox(v->bbox, bbox);
  HEdges ed = v->outgoingHEdges();
  for (HEdges::iterator h = ed.begin(); h != ed.end(); ++h) {
    (*h)->e->setBBox();
    if ((*h)->f)
      (*h)->f->update();
    else
      addFace(*h);
  }
}

void Polyhedron::removeNullFaces ()
{
  int i = 0;
  while (i < vertices.size())
    if (vertices[i]->edges.empty()) {
      if (vertices[i]->node)
	vtree.remove(vertices[i]->node);
      delete vertices[i];
      vertices[i] = vertices.back();
      vertices.pop_back();
    }
    else
      ++i;
  i = 0;
  while (i < edges.size())
    if (edges[i]->hedges.empty()) {
      delete edges[i];
      edges[i] = edges.back();
      edges.pop_back();
    }
    else
      ++i;
  i = 0;
  while (i < faces.size())
    if (!faces[i]->h) {
      delete faces[i];
      faces[i] = faces.back();
      faces.pop_back();
    }
    else
      ++i;
}

Octree<Face *> * Polyhedron::faceOctree (double s) const
{
  return Octree<Face *>::octree(faces, bbox, s);
}

Octree<Cell *> * Polyhedron::cellOctree () const
{
  Cells ce;
  for (int i = 1; i < cells.size(); ++i)
    ce.push_back(cells[i]);
  return Octree<Cell *>::octree(ce, bbox);
}

void Polyhedron::computeWindingNumbers ()
{
  if (cells.empty())
    formCells();
  cells[0]->wn = 0;
  HFaces st;
  updateWN(cells[0], st);
  set<Cell *> done;
  while (!st.empty()) {
    HFace *f = st.back();
    st.pop_back();
    if (done.insert(f->s->c).second) {
      f->s->c->wn = f->twin()->s->c->wn + (f->pos() ? -1 : 1);
      updateWN(f->s->c, st);
    }
  }
}

void Polyhedron::updateWN (Cell *c, HFaces &st) const
{
  for (int i = 0; i < c->nShells(); ++i) {
    Shell *s = c->getShell(i);
    for (HFaces::iterator h = s->hfaces.begin(); h != s->hfaces.end(); ++h)
      st.push_back((*h)->twin());
  }
}

void Polyhedron::describe () const
{
  if (cells.empty())
    return;
  bool wflag = false;
  for (int i = 1; !wflag && i < cells.size(); ++i)
    wflag = cells[i]->wn;
  Cell *c = cells[0];
  cerr << "unbounded cell: " << c->inner.size() << " shells: ";
  for (Shells::iterator s = c->inner.begin(); s != c->inner.end(); ++s)
    cerr << (*s)->hfaces.size() << " hfaces, euler = " << (*s)->euler() << "; ";
  cerr << endl;
  for (int i = 1; i < cells.size(); ++i) {
    Cell *c = cells[i];
    cerr << "bounded cell " << i << ": ";
    if (wflag)
      cerr << "wn = " << c->wn << "; ";
    cerr << c->nShells() << " shells; ";
    for (int  j = 0; j < c->nShells(); ++j) {
      Shell *s = c->getShell(j);
      cerr << s->hfaces.size() << " hfaces, euler = " << s->euler() << "; ";
    }
    cerr << endl;
  }
}

Face * faceVertices (Vertex *a, Vertex *b, Vertex *c)
{
  HEdges ea = a->outgoingHEdges();
  for (HEdges::iterator u = ea.begin(); u != ea.end(); ++u)
    if ((*u)->head() == b) {
      Face *f = (*u)->getF();
      if (f) {
	HEdge *v = (*u)->getNext();
	if (v->head() == c && v->getF() == f) {
	  HEdge *w = v->getNext();
	  if (w->head() == a && w->getF() == f)
	    return f;
	}
      }
    }
  return 0;
}

Face * faceVertices (Vertex *a, Vertex *b, Vertex *c, Vertex *d)
{
  HEdges ea = a->outgoingHEdges();
  for (HEdges::iterator u = ea.begin(); u != ea.end(); ++u)
    if ((*u)->head() == b) {
      Face *f = (*u)->getF();
      if (f) {
	HEdge *v = (*u)->getNext();
	if (v->head() == c && v->getF() == f) {
	  HEdge *w = v->getNext();
	  if (w->head() == d && w->getF() == f) {
	    HEdge *x = w->getNext();
	    if (x->head() == a && x->getF() == f)
	      return f;
	  }
	}
      }
    }
  return 0;
}

Polyhedron * subdivide (Polyhedron *a, bool oneway)
{
  map<Edge *, Points *> epsmap;
  map<pair<Face *, Face *>, FFE> ffemap;
  vector<pair<Face *, Edge *>> fe;
  intersectFF(a, epsmap, ffemap, fe);
  map<Face *, vector<FFE>> femap = feMap(epsmap, ffemap, fe);
  intersectFFF(ffemap, femap);
  set<Point *> bpts;
  subedgesE(epsmap, oneway, bpts);
  subedgesFF(ffemap, oneway, bpts);
  Polyhedron *b = subfaces(a, oneway, epsmap, femap);
  deleteMaps(epsmap, ffemap);
  return b;
}

void intersectFF (Polyhedron *a, map<Edge *, Points *> &epsmap,
		  map<pair<Face *, Face *>, FFE> &ffemap,
		  vector<pair<Face *, Edge *>> &fe)
{
  Octree<Face *> *octree = a->faceOctree();
  vector<pair<Face *, Face *>> ff1, ff;
  octree->pairs(ff1);
  delete octree;
  map<pair<Face *, Edge *>, PTR<Point>> fepmap;
  map<Face *, Points *> fpsmap;
  intersectFE(ff1, fepmap, epsmap, fpsmap, ff, fe);
  intersectFF(ff, a->perturbed, fepmap, epsmap, fpsmap, ffemap);
  for (map<Face *, Points *>::iterator i = fpsmap.begin(); i != fpsmap.end(); ++i)
    delete i->second;
}

void intersectFE (const vector<pair<Face *, Face *>> &ff1,
		  map<pair<Face *, Edge *>, PTR<Point>> &fepmap,
		  map<Edge *, Points *> &epsmap,
		  map<Face *, Points *> &fpsmap,
		  vector<pair<Face *, Face *>> &ff,
		  vector<pair<Face *, Edge *>> &fe)
{
  const unsigned int n = 8;
  unsigned int k = ff1.size(), m = k/n, is = 0;
  FEData fed[n];
  for (int i = 0; i < n; ++i) {
    fed[i].i = i;
    fed[i].is = is;
    is = i + 1 == n ? k : is + m;
    fed[i].ie = is;
    fed[i].ff1 = &ff1;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) intersectFET, (void *) (fed + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  for (int i = 0; i < n; ++i) {
    for (vector<pair<Edge *, PTR<Point>>>::iterator x = fed[i].ep.begin();
	 x != fed[i].ep.end(); ++x)
      update(x->first, x->second, epsmap);
    for (vector<pair<pair<Face *, Edge *>, PTR<Point>>>::iterator
	   x = fed[i].fep.begin(); x != fed[i].fep.end(); ++x)
      fepmap.insert(*x);
    for (vector<pair<Face *, PTR<Point>>>::iterator x = fed[i].fp.begin();
	 x != fed[i].fp.end(); ++x)
      update(x->first, x->second, fpsmap);
    ff.insert(ff.end(), fed[i].ff.begin(), fed[i].ff.end());
    fe.insert(fe.end(), fed[i].fe.begin(), fed[i].fe.end());
  }
}

void intersectFET (void *ptr)
{
  FEData *fed = (FEData *) ptr;
  BaseObject::addThread(fed->i);
  for (int i = fed->is; i < fed->ie; ++i) {
    const pair<Face *, Face *> &ff = fed->ff1->at(i);
    intersectFE(ff.first, ff.second, fed->ep, fed->fep, fed->fp, fed->ff, fed->fe);
  }
}

void intersectFE (Face *f, Face *g,
		  vector<pair<Edge *, PTR<Point>>> &ep,
		  vector<pair<pair<Face *, Edge *>, PTR<Point>>> &fep,
		  vector<pair<Face *, PTR<Point>>> &fp,
		  vector<pair<Face *, Face *>> &ff,
		  vector<pair<Face *, Edge *>> &fe)
{
  int sf[] = {0, 0, 0, 0}, sg[] = {0, 0, 0, 0};
  if (f->intersectsFP(g, sg) || g->intersectsFP(f, sf) ||
      f->verifyFP(g, sg) || g->verifyFP(f, sf))
    return;
  intersectFE(f, g, sg, ep, fep, fe, fp);
  intersectFE(g, f, sf, ep, fep, fe,  fp);
  if (signChange(sf) && signChange(sg))
    ff.push_back(pair<Face *, Face *>(f, g));
}

void intersectFEP (Face *f, Face *g,
		   vector<pair<Edge *, PTR<Point>>> &ep,
		   vector<pair<Face *, Edge *>> &fe)
{
  HEdges eg = g->getBoundary()->edgeLoop();
  for (HEdges::iterator e = eg.begin(); e != eg.end(); ++e)
    intersectFEP(f, *e, ep, fe);
}

void intersectFEP (Face *f, HEdge *h,
		   vector<pair<Edge *, PTR<Point>>> &ep,
		   vector<pair<Face *, Edge *>> &fe)
{
  if (!first(h))
    return;
  Edge *e = h->getE();
  if (f->boundaryEdge(e) || !bboxOverlap(f->getBBox(), e->getBBox()))
    return;
  HEdges ef = f->getBoundary()->edgeLoop();
  for (HEdges::iterator i = ef.begin(); i != ef.end(); ++i)
    if (first(*i))
      intersectEE(e, (*i)->getE(), f->getPC(), ep);
  fe.push_back(pair<Face *, Edge *>(f, e));
}

void intersectEE (Edge *e, Edge *f, int c,
		  vector<pair<Edge *, PTR<Point>>> &ep)
{
  Points pe, pf;
  intersectEE(e->getT()->getP(), e->getH()->getP(), f->getT()->getP(),
	      f->getH()->getP(), c, pe, pf);
  for (Points::iterator p = pe.begin(); p != pe.end(); ++p)
    ep.push_back(pair<Edge *, PTR<Point>>(e, *p));
  for (Points::iterator p = pf.begin(); p != pf.end(); ++p)
    ep.push_back(pair<Edge *, PTR<Point>>(f, *p));
}

void intersectEE (Point *et, Point *eh, Point *ft, Point *fh,
		  int c, Points &pe, Points &pf)
{
  int tp1 = et == ft || et == fh ? 0 : LeftTurn(et, ft, fh, c),
     tp2 = eh == ft || eh == fh ? 0 :LeftTurn(eh, ft, fh, c);
  if (tp1*tp2 == 1)
    return;
  if (tp1 == 0 && onEdge(et, ft, fh, true))
    pf.push_back(et);
  if (tp2 == 0 && onEdge(eh, ft, fh, true))
    pf.push_back(eh);
  int tp3 = ft == et || ft == eh ? 0 : LeftTurn(ft, et, eh, c),
    tp4 = fh == et || fh == eh ? 0 : LeftTurn(fh, et, eh, c);
  if (tp3 == 0 && onEdge(ft, et, eh, true))
    pe.push_back(ft);
  if (tp4 == 0 && onEdge(fh, et, eh, true))
    pe.push_back(fh);
  if (tp1*tp2 == -1 && tp3*tp4 == -1) {
    PTR<Point> p = new EEPoint(et, eh, ft, fh, c);
    pe.push_back(p);
    pf.push_back(p);
    }
}

void intersectFE (Face *f, Face *g, int *sg,
		  vector<pair<Edge *, PTR<Point>>> &ep,
		  vector<pair<pair<Face *, Edge *>, PTR<Point>>> &fep,
		  vector<pair<Face *, Edge *>> &fe,
		  vector<pair<Face *, PTR<Point>>> &fp)
{
  HEdges eg = g->getBoundary()->edgeLoop();
  int n = eg.size();
  bool flag = false;
  for (int i = 0; i < n; ++i)
    if ((sg[i] == 0 || sg[i] == 2) && (sg[(i+1)%n] == 0 || sg[(i+1)%n] == 2)) {
      intersectFEP(f, eg[i], ep, fe);
      flag = true;
    }
  if (flag || f->sharedEdge(g))
    return;
  for (int i = 0; i < n; ++i)
    if (sg[i] == 0)
      intersectFV(f, eg[i], ep, fp);
    else if (sg[i] == - sg[(i+1)%n])
      intersectFEG(f, eg[i], ep, fep);
}

void intersectFV (Face *f, HEdge *h,
		  vector<pair<Edge *, PTR<Point>>> &ep,
		  vector<pair<Face *, PTR<Point>>> &fp)
{
  int ie = -1;
  PTR<Point> p = h->tail()->getP();
  if (f->contains(p, false, &ie)) {
    if (ie == -1)
      fp.push_back(pair<Face *, PTR<Point>>(f, p));
    else {
      HEdges ef = f->getBoundary()->edgeLoop();
      ep.push_back(pair<Edge *, PTR<Point>>(ef[ie]->getE(), p));
    }
  }
}

void intersectFEG (Face *f, HEdge *h,
		   vector<pair<Edge *, PTR<Point>>> &ep,
		   vector<pair<pair<Face *, Edge *>, PTR<Point>>> &fep)
{
  if (!first(h))
    return;
  Edge *e = h->getE();
  if (!bboxOverlap(f->getBBox(), e->getBBox()))
    return;
  PTR<Point> p = new EPPoint(e->getT()->getP(), e->getH()->getP(), f->getP());
  int ie = -1;
  if (f->contains(p, false, &ie)) {
    ep.push_back(pair<Edge *, PTR<Point>>(e, p));
    if (ie == -1) {
      pair<Face *, Edge *> fe(f, e);
      fep.push_back(pair<pair<Face *, Edge *>, PTR<Point>>(fe, p));
    }
    else {
      HEdges ef = f->getBoundary()->edgeLoop();
      ep.push_back(pair<Edge *, PTR<Point>>(ef[ie]->getE(), p));
    }
  }
}

bool first (HEdge *h)
{
  return h == h->getE()->getHEdge(0);
}

bool signChange (int *s)
{
  bool pflag = false, nflag = false;
  for (int i = 0; i < 4; ++i)
    if (s[i] == 1)
      pflag = true;
    else if (s[i] == -1)
      nflag = true;
  return pflag && nflag;
}

void update (Edge *e, PTR<Point> p, map<Edge *, Points *> &epsmap)
{
  map<Edge *, Points *>::iterator i = epsmap.find(e);
  if (i == epsmap.end()) {
    Points *pts = new Points;
    pts->push_back(p);
    epsmap.insert(pair<Edge *, Points *>(e, pts));
  }
  else
    i->second->push_back(p);
}

void update (Face *f, PTR<Point> p, map<Face *, Points *> &fpsmap)
{
  map<Face *, Points *>::iterator i = fpsmap.find(f);
  if (i == fpsmap.end()) {
    Points *pts = new Points;
    pts->push_back(p);
    fpsmap.insert(pair<Face *, Points *>(f, pts));
  }
  else
    i->second->push_back(p);
}

void intersectFF (const vector<pair<Face *, Face *>> &ff, bool perturbed,
		  const map<pair<Face *, Edge *>, PTR<Point>> &fepmap,
		  const map<Edge *, Points *> &epsmap,
		  const map<Face *, Points *> &fpsmap,
		  map<pair<Face *, Face *>, FFE> &ffemap)
{
  const unsigned int n = 8;
  unsigned int k = ff.size(), m = k/n, is = 0;
  FFData ffd[n];
  for (int i = 0; i < n; ++i) {
    ffd[i].i = i;
    ffd[i].is = is;
    is = i + 1 == n ? k : is + m;
    ffd[i].ie = is;
    ffd[i].ff = &ff;
    ffd[i].perturbed = perturbed;
    ffd[i].fepmap = &fepmap;
    ffd[i].epsmap = &epsmap;
    ffd[i].fpsmap = &fpsmap;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) intersectFFT, (void *) (ffd + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  vector<pair< pair<Face *, Face *>, pair<PTR<Point>, PTR<Point>>>> ffpp;
  for (int i = 0; i < n; ++i)
    for (vector<pair< pair<Face *, Face *>, pair<PTR<Point>, PTR<Point>>>>::iterator
	 j = ffd[i].ffpp.begin(); j != ffd[i].ffpp.end(); ++j)
    update(j->first.first, j->first.second, j->second.first, j->second.second,
	   ffemap);
}

void intersectFFT (void *ptr)
{
  FFData *ffd = (FFData *) ptr;
  BaseObject::addThread(ffd->i);
  for (int i = ffd->is; i < ffd->ie; ++i) {
    const pair<Face *, Face *> &ff = ffd->ff->at(i);
    intersectFF(ff.first, ff.second, ffd->perturbed, *ffd->fepmap, *ffd->epsmap,
		*ffd->fpsmap, ffd->ffpp);
  }
}

void intersectFF (Face *f, Face *g, bool perturbed,
		  const map<pair<Face *, Edge *>, PTR<Point>> &fepmap,
		  const map<Edge *, Points *> &epsmap,
		  const map<Face *, Points *> &fpsmap,
		  vector<pair< pair<Face *, Face *>, pair<PTR<Point>, PTR<Point>>>> &ffpp)
{
  Points pts;
  findFE(f, g, fepmap, pts);
  findFE(g, f, fepmap, pts);
  if (pts.size() == 1) {
    Point *p = f->sharedVertex(g);
    if (p)
      pts.push_back(p);
  }
  if (!perturbed && pts.size() < 2)
    sharedVertices(f, g, epsmap, fpsmap, pts);
  if (pts.size() != 2)
    return;
  Face *ff = f < g ? f : g, *gg = f < g ? g : f;
  bool flag = FFOrder(ff, gg, pts[0], pts[1]) == 1;
  pair<Face *, Face *> fg(ff, gg);
  pair<PTR<Point>, PTR<Point>> pq(pts[flag ? 0 : 1], pts[flag ? 1 : 0]);
  ffpp.push_back(pair<pair<Face *, Face *>, pair<PTR<Point>, PTR<Point>>>(fg, pq));
}

void findFE (Face *f, Face *g,
	     const map<pair<Face *, Edge *>, PTR<Point>> &fepmap,
	     Points &pts)
{
  HEdges eg = g->getBoundary()->edgeLoop();
  for (HEdges::iterator h = eg.begin(); h != eg.end(); ++h) {
    pair<Face *, Edge *> fe(f, (*h)->getE());
    map<pair<Face *, Edge *>, PTR<Point>>::const_iterator i = fepmap.find(fe);
    if (i != fepmap.end())
      pts.push_back(i->second);
  }
}

void sharedVertices (Face *f, Face *g,
		     const map<Edge *, Points *> &epsmap,
		     const map<Face *, Points *> &fpsmap,
		     Points &pts)
{
  set<Point *> pf = vertices(f, epsmap, fpsmap), pg = vertices(g, epsmap, fpsmap),
    pfg;
  set_intersection(pf.begin(), pf.end(), pg.begin(), pg.end(),
		   inserter(pfg, pfg.begin()));
  pts.insert(pts.end(), pfg.begin(), pfg.end());
  removeDuplicates(pts);
}

set<Point *> vertices (Face *f, const map<Edge *, Points *> &epsmap,
		       const map<Face *, Points *> &fpsmap)
{
  set<Point *> res = boundaryVertices(f, epsmap);
  map<Face *, Points *>::const_iterator i = fpsmap.find(f);
  if (i != fpsmap.end())
    res.insert(i->second->begin(), i->second->end());
  return res;
}

set<Point *> boundaryVertices (Face *f, const map<Edge *, Points *> &epsmap)
{
  set<Point *> res;
  HEdges ed = f->getBoundary()->edgeLoop();
  for (HEdges::iterator h = ed.begin(); h != ed.end(); ++h) {
    res.insert((*h)->tail()->getP());
    map<Edge *, Points *>::const_iterator i = epsmap.find((*h)->getE());
    if (i != epsmap.end())
      res.insert(i->second->begin(), i->second->end());
  }
  return res;
}

void removeDuplicates (Points &pts)
{
  for (int i = 0; i + 1 < pts.size(); ++i) {
    int j = i + 1;
    while (j < pts.size()) {
      if (pts[i]->identical(pts[j])) {
	pts[j] = pts.back();
	pts.pop_back();
      }
      else
	++j;
    }
  }
}

void update (Face *f, Face *g, PTR<Point> a, PTR<Point> b,
	     map<pair<Face *, Face *>, FFE> &ffemap)
{
  pair<Face *, Face *> fg(f, g);
  FFE ffe(f, g, a, b); 
  ffemap.insert(pair<pair<Face *, Face *>, FFE>(fg, ffe));
}

map<Face *, vector<FFE>> feMap (map<Edge *, Points *> &epsmap,
				map<pair<Face *, Face *>, FFE> &ffemap,
				const vector<pair<Face *, Edge *>> &fe)
{
  map<Face *, vector<FFE>> femap;
  for (map<pair<Face *, Face *>, FFE>::iterator i = ffemap.begin();
       i != ffemap.end(); ++i) {
    update(i->first.first, i->second, femap);
    update(i->first.second, i->second, femap);
  }
  for (vector<pair<Face *, Edge *>>::const_iterator i = fe.begin(); i != fe.end(); ++i)
    feMapFE(epsmap, ffemap, i->first, i->second, femap);
  return femap;
}

void feMapFE (map<Edge *, Points *> &epsmap,
	      map<pair<Face *, Face *>, FFE> &ffemap,
	      Face *f, Edge *e, map<Face *, vector<FFE>> &femap)
{
  HEdges ef = f->getBoundary()->edgeLoop();
  for (HEdges::iterator h = ef.begin(); h != ef.end(); ++h)
    if ((*h)->tail()->getP()->onLine(e->getT()->getP(), e->getH()->getP()) &&
	(*h)->head()->getP()->onLine(e->getT()->getP(), e->getH()->getP()))
      return;
  Points *pts;
  map<Edge *, Points *>::const_iterator i = epsmap.find(e);
  if (i == epsmap.end()) {
    pts = new Points;
    epsmap.insert(pair<Edge *, Points *>(e, pts));
  }
  else
    pts = i->second;
  pts->insert(pts->begin(), (e->getT()->getP()));
  pts->push_back(e->getH()->getP());
  FFE ffe(f, 0, pts);
  update(f, ffe, femap);
  for (int i = 0; i < e->HEdgesN(); ++i) {
    Face *g = e->getHEdge(i)->getF();
    pair<Face *, Face *> fg(f < g ? f : g, f < g ? g : f);
    ffemap.insert(pair<pair<Face *, Face *>, FFE>(fg, FFE()));
  }
}

void update (Face *f, const FFE &ffe, map<Face *, vector<FFE>> &femap)
{
  map<Face *, vector<FFE>>::iterator i = femap.find(f);
  if (i == femap.end()) {
    vector<FFE> ffes;
    ffes.push_back(ffe);
    femap.insert(pair<Face *, vector<FFE>>(f, ffes));
  }
  else
    i->second.push_back(ffe);
}

void intersectFFF (const map<pair<Face *, Face *>, FFE> &ffemap,
		   const map<Face *, vector<FFE>> &femap)
{
  vector<pair<Points *, PTR<Point>>> psps;
  intersectFFFAux(ffemap, femap, psps);
  for (vector<pair<Points *, PTR<Point>>>::iterator i = psps.begin();
       i != psps.end(); ++i)
    i->first->push_back(i->second);
}

void intersectFFFAux (const map<pair<Face *, Face *>, FFE> &ffemap,
		      const map<Face *, vector<FFE>> &femap,
		      vector<pair<Points *, PTR<Point>>> &psps)
{
  const unsigned int n = 8;
  unsigned int k = femap.size(), m = k/n, is = 0;
  FFFData fd[n];
  for (int i = 0; i < n; ++i) {
    fd[i].i = i;
    fd[i].is = is;
    is = i + 1 == n ? k : is + m;
    fd[i].ie = is;
    fd[i].ffemap = &ffemap;
    fd[i].femap = &femap;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) intersectFFFT, (void *) (fd + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  for (int i = 0; i < n; ++i)
    psps.insert(psps.end(), fd[i].psps.begin(), fd[i].psps.end());
}

void intersectFFFT (void *ptr)
{
  FFFData *fd = (FFFData *) ptr;
  BaseObject::addThread(fd->i);
  map<Face *, vector<FFE>>::const_iterator x = fd->femap->begin();
  for (int i = 0; i < fd->is; ++i, ++x)
    ;
  for (int i = fd->is; i < fd->ie; ++i, ++x)
    intersectFFF(x->first, x->second, *fd->ffemap, fd->psps);
}

void intersectFFF (Face *f, const vector<FFE> &ed,
		   const map<pair<Face *, Face *>, FFE> &ffemap,
		   vector<pair<Points *, PTR<Point>>> &psps)
{
  int n = ed.size(), pc = f->getPC();
  for (int i = 0; i + 1 < n; ++i) {
    const FFE &ei = ed[i];
    Face *fi = ei.f == f ? ei.g : ei.f;
    for (int j = i + 1; j < n; ++j) {
      const FFE &ej = ed[j];
      Face *fj = ej.f == f ? ej.g : ej.f;
      if ((fi || fj) && bboxOverlap(ei.bbox, ej.bbox)) {
	bool flag = !fi || !fj;
	map<pair<Face *, Face *>, FFE>::const_iterator y = ffemap.end();
	if (!flag) {
	  pair<Face *, Face *> fifj(fi < fj ? fi : fj, fi < fj ? fj : fi);
	  y = ffemap.find(fifj);
	  if (y != ffemap.end())
	    flag = !y->second.pts || f < fi && f < fj;
	}
	if (flag) {
	  Points pi, pj;
	  intersectEE(ei.pts->at(0), ei.pts->at(1), ej.pts->at(0), ej.pts->at(1),
		      pc, pi, pj);
	  if (fi)
	    for (Points::iterator p = pi.begin(); p != pi.end(); ++p)
	      psps.push_back(pair<Points *, PTR<Point>>(ei.pts, *p));
	  if (fj)
	    for (Points::iterator p = pj.begin(); p != pj.end(); ++p)
	      psps.push_back(pair<Points *, PTR<Point>>(ej.pts, *p));
	  if (pi.size() == 1 && pj.size() == 1 && y != ffemap.end() && y->second.pts) {
	    psps.push_back(pair<Points *, PTR<Point>>(y->second.pts, pi[0]));
	    EEPoint *pee = dynamic_cast<EEPoint *>((Point *) pi[0]);
	    if (pee) {
	      pee->addFace(f);
	      pee->addFace(fi);
	      pee->addFace(fj);
	    }
	  }
	}
      }
    }
  }
}

void subedgesE (map<Edge *, Points *> &epsmap, bool oneway, set<Point *> &bpts)
{
  const unsigned int n = 8;
  unsigned int k = epsmap.size(), m = k/n, is = 0;
  SEEData see[n];
  for (int i = 0; i < n; ++i) {
    see[i].i = i;
    see[i].is = is;
    is = i + 1 == n ? k : is + m;
    see[i].ie = is;
    see[i].epsmap = &epsmap;
    see[i].oneway = oneway;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) subedgesET, (void *) (see + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  for (int i = 0; i < n; ++i)
    bpts.insert(see[i].bpts.begin(), see[i].bpts.end());
}

void subedgesET (void *ptr)
{
  SEEData *see = (SEEData *) ptr;
  BaseObject::addThread(see->i);
  map<Edge *, Points *>::const_iterator x = see->epsmap->begin();
  for (int i = 0; i < see->is; ++i, ++x)
    ;
  for (int i = see->is; i < see->ie; ++i, ++x)
    subedgesE(x->first, *x->second, see->oneway, see->bpts);
}

void subedgesE (Edge *e, Points &pts, bool oneway, set<Point *> &bpts)
{
  PTR<Point> t = e->getT()->getP(), h = e->getH()->getP();
  if (pts.size() > 1)
    sort(pts.begin(), pts.end(), PointOrderPPP(t, h));
  pts.insert(pts.begin(), t);
  pts.push_back(h);
  subedgesE(pts, oneway, bpts);
}

void subedgesE (Points &pts, bool oneway, set<Point *> &bpts)
{
  if (oneway)
    subedgesE1(pts, bpts);
  else
    subedgesP(pts);
}

void subedgesE1 (Points &pts, set<Point *> &bpts)
{
  Points pp;
  PTR<Point> u = new DiffPoint(pts.back(), pts[0]);
  int i = 0, n = pts.size();
  bool lfree = true, pblocked = false;
  while (i + 1 < n) {
    int j = i + 1;
    while (j + 1 < n && pts[j]->identical(pts[j+1]))
      ++j;
    int s = 0;
    if (i + 1 == j) {
      EPPoint *epp = dynamic_cast<EPPoint *>((Point *) pts[j]);
      if (epp)
	s = PlaneRayAlignment(epp->getP(), u);
    }
    if (lfree && s <= 0) {
      pp.push_back(pts[i]);
      pp.push_back(pts[j]);
      pblocked = false;
    }
    else {
      if (pblocked)
	bpts.insert(pts[i]);
      pblocked = true;
    }
    i = j;
    lfree = s >= 0;
  }
  pts = pp;
}

void subedgesP (Points &pts)
{
  Points pp;
  for (int i = 0; i + 1 < pts.size(); ++i) {
    pp.push_back(pts[i]);
    pp.push_back(pts[i+1]);
  }
  pts = pp;
}

void subedgesFF (map<pair<Face *, Face *>, FFE> &ffemap, bool oneway,
		 const set<Point *> &bpts)
{
  const unsigned int n = 8;
  unsigned int k = ffemap.size(), m = k/n, is = 0;
  SEFData sef[n];
  for (int i = 0; i < n; ++i) {
    sef[i].i = i;
    sef[i].is = is;
    is = i + 1 == n ? k : is + m;
    sef[i].ie = is;
    sef[i].ffemap = &ffemap;
    sef[i].oneway = oneway;
    sef[i].bpts = &bpts;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) subedgesFFT, (void *) (sef + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
}

void subedgesFFT (void *ptr)
{
  SEFData *sef = (SEFData *) ptr;
  BaseObject::addThread(sef->i);
  map<pair<Face *, Face *>, FFE>::iterator x = sef->ffemap->begin();
  for (int i = 0; i < sef->is; ++i, ++x)
    ;
  for (int i = sef->is; i < sef->ie; ++i, ++x)
    subedgesFF(x->second, sef->oneway, *sef->bpts);
}

void subedgesFF (FFE &ffe, bool oneway, const set<Point *> &bpts)
{
  if (ffe.pts) {
    if (ffe.pts->size() > 1)
      sort(ffe.pts->begin(), ffe.pts->end(), 
	   PointOrderPPP(ffe.pts->at(0), ffe.pts->at(1)));
    if (oneway)
      subedgesFF(ffe, bpts);
    else
      subedgesP(*ffe.pts);
  }
}

void subedgesFF (FFE &ffe, const set<Point *> &bpts)
{
  Points &pts = *ffe.pts, pp;
  PTR<Point> u = new DiffPoint(pts[1], pts[0]);
  int i = 0, n = pts.size();
  bool lfree = true;
  while (i + 1 < n) {
    int j = i + 1;
    while (j + 1 < n && pts[j]->identical(pts[j+1]))
      ++j;
    int s = 0;
    if (i + 1 == j && i + 2 < n) {
      EEPoint *eh = dynamic_cast<EEPoint *>((Point *) pts[i+1]);
      if (eh) {
	Face *g = otherFace(ffe.f, ffe.g, eh->getFaces());
	if (g)
	  s = PlaneRayAlignment(g->getP(), u);
      }
    }
    if ((i > 0 || bpts.find(pts[i]) == bpts.end()) &&
	(i + 2 < n || bpts.find(pts[i+1]) == bpts.end()) &&
	lfree && s <= 0) {
      pp.push_back(pts[i]);
      pp.push_back(pts[j]);
    }
    lfree = s >= 0;
    i = j;
  }
  pts = pp;
}

Face * otherFace (Face *f, Face *g, const Faces &fa)
{
  for (Faces::const_iterator h = fa.begin(); h != fa.end(); ++h)
    if (*h != f && *h != g)
      return *h;
  return 0;
}

Points SEdge::loop () const
{
  Points pts;
  const SEdge *e = this;
  do {
    pts.push_back(e->tail);
    e = e->next();
  }
  while (e != this);
  return pts;
}

Polyhedron * subfaces (Polyhedron *a, bool oneway,
		       const map<Edge *, Points *> &epsmap,
		       const map<Face *, vector<FFE>> &femap)
{
  const unsigned int n = 16;
  unsigned int k = a->faces.size(), m = k/n, is = 0;
  SUData sud[n];
  for (int i = 0; i < n; ++i) {
    sud[i].i = i;
    sud[i].is = is;
    is = i + 1 == n ? k : is + m;
    sud[i].ie = is;
    sud[i].a = a;
    sud[i].oneway = oneway;
    sud[i].epsmap = &epsmap;
    sud[i].femap = &femap;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) subfacesT, (void *) (sud + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  Polyhedron *b = new Polyhedron(a->perturbed);
  PVMap pvmap;
  for (int i = 0; i < n; ++i)
    for (SFaces::iterator s = sud[i].sf.begin(); s != sud[i].sf.end(); ++s)
      addFace(*s, b, pvmap);
  return b;
}

void subfacesT (void *ptr)
{
  SUData *sud = (SUData *) ptr;
  BaseObject::addThread(sud->i);
  for (int i = sud->is; i < sud->ie; ++i)
    subfaces(sud->a->faces[i], sud->oneway, *sud->epsmap,
	     *sud->femap, sud->a->perturbed, sud->sf);
}

void subfaces (Face *f, bool oneway, const map<Edge *, Points *> &epsmap,
	       const map<Face *, vector<FFE>> &femap, bool perturbed,
	       SFaces &sf)
{
  int c = f->getPC();
  SEdges ed = subedges(f, oneway, epsmap, femap, perturbed), outer, inner;
  for (SEdges::iterator e = ed.begin(); e != ed.end(); ++e) {
    SEdge *l = findLoop(*e);
    if (l)
      if (LeftTurn(l->tail, l->head(), l->next()->head(), c) == 1)
	outer.push_back(l);
      else
	inner.push_back(l);
  }
  sort(outer.begin(), outer.end(), SEdgeHeadOrder());
  sort(inner.begin(), inner.end(), SEdgeHeadOrder());
  SFaces nf;
  for (SEdges::iterator e = outer.begin(); e != outer.end(); ++e)
    nf.push_back(SFace(f, (*e)->loop()));
  for (SEdges::iterator e = inner.begin(); e != inner.end(); ++e)
    addInner(*e, oneway, c, nf);
  sf.insert(sf.end(), nf.begin(), nf.end());
  for (SEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    delete *e;
}

SEdges subedges (Face *f, bool oneway,
		 const map<Edge *, Points *> &epsmap,
		 const map<Face *, vector<FFE>> &femap, bool perturbed)
{
  vector<pair<PTR<Point>, PTR<Point>>> pp;
  HEdges ef = f->getBoundary()->edgeLoop();
  for (HEdges::iterator h = ef.begin(); h != ef.end(); ++h)
    subedgesH(*h, epsmap, pp);
  map<Face *, vector<FFE>>::const_iterator i = femap.find(f);
  if (i != femap.end())
    for (vector<FFE>::const_iterator j = i->second.begin();
	 j != i->second.end(); ++j)
      subedgesFFE(f, *j, oneway, pp);
  return subedgesPP(pp, f->getPC(), perturbed);
}

void subedgesH (HEdge *h, const map<Edge *, Points *> &epsmap,
		vector<pair<PTR<Point>, PTR<Point>>> &pp)
{
  map<Edge *, Points *>::const_iterator i = epsmap.find(h->getE());
  if (i == epsmap.end()) {
    Points pts;
    pts.push_back(h->tail()->getP());
    pts.push_back(h->head()->getP());
    subedges(pts, true, false, pp);
  }
  else
    subedges(*i->second, h->getForward(), !h->getForward(), pp);
}

void subedges (const Points &pts, bool fflag, bool bflag,  
	       vector<pair<PTR<Point>, PTR<Point>>> &pp)
{
  for (int i = 0; i + 1 < pts.size(); i += 2) {
    if (fflag)
      pp.push_back(pair<PTR<Point>, PTR<Point>>(pts[i], pts[i+1]));
    if (bflag)
      pp.push_back(pair<PTR<Point>, PTR<Point>>(pts[i+1], pts[i]));
  }
}

void subedgesFFE (Face *f, const FFE &ffe, bool oneway,
		  vector<pair<PTR<Point>, PTR<Point>>> &pp)
{
  Points &pts = *ffe.pts;
  if (ffe.g) {
    bool fflag = !oneway || ffe.f == f, gflag = !oneway || ffe.f != f;
    subedges(pts, fflag, gflag, pp);
  }
  else
    for (int i = 0; i + 1 < pts.size(); i += 2) {
      CentroidPoint c(pts[i], pts[i+1]);
      if (f->contains(&c, true)) {
	pp.push_back(pair<PTR<Point>, PTR<Point>>(pts[i], pts[i+1]));
	pp.push_back(pair<PTR<Point>, PTR<Point>>(pts[i+1], pts[i]));
      }
    }
}

SEdges subedgesPP (const vector<pair<PTR<Point>, PTR<Point>>> &pp, int c,
		   bool perturbed)
{
  set<PTR<Point>, PointOrderRP> ps;
  map<PTR<Point>, PTR<Point>> ppmap;
  set<pair<PTR<Point>, PTR<Point>>> pps;
  for (vector<pair<PTR<Point>, PTR<Point>>>::const_iterator i = pp.begin();
       i != pp.end(); ++i) {
    PTR<Point> t = perturbed ? i->first : getPoint(i->first, ps, ppmap),
      h = perturbed ? i->second : getPoint(i->second, ps, ppmap);
    if (t != h)
      pps.insert(pair<PTR<Point>, PTR<Point>>(t, h));
  }
  SEdges ed;
  for (set<pair<PTR<Point>, PTR<Point>>>::const_iterator i = pps.begin();
       i != pps.end(); ++i)
    sedges(i->first, i->second, ed);
  sort(ed.begin(), ed.end(), SEdgeCWOrder(c));
  int i = 0, n = ed.size();
  while (i < n) {
    int j = i + 1;
    while (j < n && ed[i]->tail == ed[j]->tail) {
      ed[j-1]->cw = ed[j];
      ++j;
    }
    ed[j-1]->cw = ed[i];
    i = j;
  }
  return ed;
}

PTR<Point> getPoint (PTR<Point> p, set<PTR<Point>, PointOrderRP> &ps,
		     map<PTR<Point>, PTR<Point>> &ppmap)
{
  map<PTR<Point>, PTR<Point>>::iterator i = ppmap.find(p);
  if (i != ppmap.end())
    return i->second;
  PTR<Point> q = *ps.insert(p).first;
  ppmap.insert(pair<PTR<Point>, PTR<Point>>(p, q));
  return q;
}

void sedges (PTR<Point> t, PTR<Point> h, SEdges &ed)
{
  SEdge *f = new SEdge(t, true), *b = new SEdge(h, false, f);
  f->twin = b;
  ed.push_back(f);
  ed.push_back(b);
}

SEdge * findLoop (SEdge *e)
{
  if (e->flag || !e->forward)
    return 0;
  SEdge *l = e, *f = e;
  do {
    f->flag = true;
    if (SEdgeHeadOrder()(f, l))
      l = f;
    f = f->next();
    if (!(f && f->forward))
      return 0;
  }
  while (f != e);
  return l;
}

void addInner (SEdge *e, bool oneway, int c, SFaces &sf)
{
  Point *p = e->head();
  for (int i = sf.size() - 1; i >= 0; --i)
    if (contains(sf[i].b[0], c, p, true, 0)) {
      if (oneway)
	for (int j = 1; j < sf[i].b.size(); ++j)
	  if (contains(sf[i].b[j], c, p, true, 0))
	    return;
      sf[i].b.push_back(e->loop());
      return;
    }
}

HEdges MFace::boundaryHEdges () const
{
  HEdges ed = h->edgeLoop();
  for (HEdges::const_iterator e = inner.begin(); e != inner.end(); ++e) {
    HEdges eed = (*e)->edgeLoop();
    ed.insert(ed.end(), eed.begin(), eed.end());
  }
  return ed;
}

void MFace::triangulate (Triangles &tr)
{
  vector<Points> pp;
  pp.push_back(h->pointLoop());
  for (HEdges::iterator e = inner.begin(); e != inner.end(); ++e)
    pp.push_back((*e)->pointLoop());
  if (pp.size() == 1 && pp[0].size() == 3)
    tr.push_back(Triangle(pp[0][0], pp[0][1], pp[0][2]));
  else {
    vector<Points *> reg;
    for (vector<Points>::iterator i = pp.begin(); i != pp.end(); ++i)
      reg.push_back(&*i);
    ::triangulate(reg, getPC(), tr);
  }
}

void addFace (const SFace &sf, Polyhedron *a, PVMap &pvmap)
{
  Vertices ve = loop(sf.b[0], a, pvmap);
  if (!a->perturbed && currentFace(ve, a))
    return;
  MFace *f = new MFace(addLoop(ve, a), sf.f->getP(), sf.f->getPC());
  for (int i = 1; i < sf.b.size(); ++i)
    f->addInner(addLoop(loop(sf.b[i], a, pvmap), a));
  a->faces.push_back(f);
}

Vertices loop (const Points &pts, Polyhedron *a, PVMap &pvmap)
{
  Vertices ve;
  for (Points::const_iterator p = pts.begin(); p != pts.end(); ++p)
    ve.push_back(a->getVertex(*p, pvmap));
  return ve;
}

bool currentFace (const Vertices &ve, Polyhedron *a)
{
  int n = ve.size();
  Edge *e = a->getEdge(ve[0], ve[1]);
  for (int i = 0; i < e->HEdgesN(); ++i) {
    HEdge *h = e->getHEdge(i);
    Vertices vh = h->loop();
    if (vh.size() == n) {
      bool flag = true;
      if (h->tail() == ve[0])
      	for (int j = 0; flag && j < n; ++j)
	  flag = ve[j] == vh[j];
      else
	for (int j = 0; flag && j < n; ++j)
	  flag = ve[j] == vh[(n+1-j)%n];
      if (flag)
	return true;
    }
  }
  return false;
}

HEdge * addLoop (const Vertices &ve, Polyhedron *a)
{
  Vertex *t = ve.back();
  HEdges ed;
  for (Vertices::const_iterator h = ve.begin(); h != ve.end(); ++h) {
    HEdge *f = a->getHEdge(t, *h);
    ed.push_back(f);
    t = *h;
  }
  for (int i = 0; i + 1 < ed.size(); ++i)
    ed[i]->setNext(ed[i+1]);
  ed.back()->setNext(ed[0]);
  return ed[0];
}

void deleteMaps (map<Edge *, Points *> &epsmap,
		 map<pair<Face *, Face *>, FFE> &ffemap)
{
  for (map<Edge *, Points *>::iterator i = epsmap.begin(); i != epsmap.end(); ++i)
    delete i->second;
  for (map<pair<Face *, Face *>, FFE>::iterator i = ffemap.begin();
       i != ffemap.end(); ++i)
    delete i->second.pts;
}

Polyhedron * triangulate (Polyhedron *a)
{
  Triangles tr = triangulate(a->faces);
  Polyhedron *b = new Polyhedron(a->perturbed);
  PVMap pvmap;
  for (Triangles::iterator t = tr.begin(); t != tr.end(); ++t) 
    b->addTriangle(t->a, t->b, t->c, pvmap);
  return b;
}

Triangles triangulate (const Faces &fa)
{
  const unsigned int n = 1; // debug 8
  unsigned int k = fa.size(), m = k/n, is = 0;
  TRData trd[n];
  for (int i = 0; i < n; ++i) {
    trd[i].i = i;
    trd[i].is = is;
    is = i + 1 == n ? k : is + m;
    trd[i].ie = is;
    trd[i].fa = &fa;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) triangulateT, (void *) (trd + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  Triangles tr;
  for (int i = 0; i < n; ++i)
    tr.insert(tr.end(), trd[i].tr.begin(), trd[i].tr.end());
  return tr;
}

void triangulateT (void *ptr)
{
  TRData *trd = (TRData *) ptr;
  BaseObject::addThread(trd->i);
  for (int i = trd->is; i < trd->ie; ++i)
    trd->fa->at(i)->triangulate(trd->tr);
}

bool inSet (bool ina, bool inb, SetOp op)
{
  switch (op) {
  case Union:
    return ina || inb;
  case Intersection:
    return ina && inb;
  case Complement:
    return ina && !inb;
  }
}

Polyhedron * overlay (Polyhedron **poly, int n)
{
  Polyhedron *c = new Polyhedron(false);
  PVMap pvmap;
  for (int i = 0; i < n; ++i) {
    const Faces &fa = poly[i]->faces;
    for (Faces::const_iterator f = fa.begin(); f != fa.end(); ++f) {
      Points pf = (*f)->getBoundary()->pointLoop();
      c->addTriangle(pf[0], pf[1], pf[2], pvmap);
    }
  }
  Polyhedron *d = subdivide(c, false), *e = triangulate(d);
  delete c;
  delete d;
  return e;
}

Polyhedron * multiUnion (Polyhedron **poly, int n)
{
  Polyhedron *c = overlay(poly, n), *d = new Polyhedron(false);
  for (int i = 0; i < n; ++i)
    poly[i]->computeWindingNumbers();
  c->formCells();
  set<Cell *> cin;
  for (int i = 1; i < c->cells.size(); ++i) {
    Cell *ci = c->cells[i];
    PTR<Point> p = ci->interiorPoint();
    bool flag = false;
    for (int j = 0; !flag && j < n; ++j)
      flag = poly[j]->contains(p);
    if (flag)
      cin.insert(ci);
  }
  PVMap pvmap;
  for (Faces::iterator f = c->faces.begin(); f != c->faces.end(); ++f) {
    bool in1 = cin.find((*f)->getHFace(0)->getS()->getC()) != cin.end(),
      in2 = cin.find((*f)->getHFace(1)->getS()->getC()) != cin.end();
    if (in1 != in2) {
      Points pf = (*f)->getBoundary()->pointLoop();
      int pc = (*f)->getPC();
      if (in2)
	d->addTriangle(pf[0], pf[1], pf[2], pvmap, pc);
      else
	d->addTriangle(pf[2], pf[1], pf[0], pvmap, - pc);
    }
  }
  delete c;
  return d;
}

Polyhedron * coalesce (Polyhedron *a)
{
  Polyhedron *b = coalesceFaces(a),
    *c = coalesceEdges(b);
  delete b;
  return c;
}

Polyhedron * coalesceFaces (Polyhedron *a)
{
  vector<set<Face *>> fg = groupFaces(a);
  Polyhedron *b = new Polyhedron(a->perturbed);
  PVMap pvmap;
  for (vector<set<Face *>>::iterator f = fg.begin(); f != fg.end(); ++f)
    coalesceFace(b, pvmap, *f);
  return b;
}

vector<set<Face *>> groupFaces (Polyhedron *a)
{
  vector<set<Face *>> res;
  set<Face *> done;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f)
    if (done.insert(*f).second) {
      Faces st;
      set<Face *> fg;
      st.push_back(*f);
      while (!st.empty()) {
	Face *g = st.back();
	st.pop_back();
	if (g->coplanar(*f) && fg.insert(g).second) {
	  HEdges ed = g->getBoundary()->edgeLoop();
	  for (HEdges::iterator h = ed.begin(); h != ed.end(); ++h)
	    st.push_back((*h)->getF());
	}
      }
      res.push_back(fg);
    }
  return res;
}

void coalesceFace (Polyhedron *a, PVMap &pvmap, const set<Face *> &fs)
{
  /* to do
  HEdges ed, inner;
  for (set<Face *>::const_iterator f = fs.begin(); f != fs.end(); ++f) {
    HEdges ef = (*f)->getBoundary()->edgeLoop();
    for (HEdges::iterator h = ef.begin(); h != ef.end(); ++h) {
      bool flag = true;
      for (HEdge *i = (*h)->ccw(); flag && i != *h; i = i->ccw())
	flag = fs.find(i->getF()) == fs.end();
      if (flag)
	ed.push_back(a->addHEdge(a->getVertex((*h)->tail(), pvmap),
				 a->getVertex((*h)->head(), pvmap)));
    }
  }
  int pc = (*fs.begin())->getPC();
  a->setNext(ed, pc);
  HEdge *outer = 0;
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    if (!(*e)->getFlag()) {
      Vertices ve;
      HEdge *l = (*e)->findLoop(ve);
      if (outerLoop(ve, pc))
	outer = l;
      else 
	inner.push_back(l);
    }
  Face *f = a->addFace(outer, pc);
  for (HEdges::iterator e = inner.begin(); e != inner.end(); ++e)
    f->addBoundary(*e);
  */
}

Polyhedron * coalesceEdges (Polyhedron *a)
{
  Polyhedron *b = new Polyhedron(a->perturbed);
  PVMap pvmap;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f)
    coalesceEdges(b, pvmap, *f);
  return b;
}

void coalesceEdges (Polyhedron *a, PVMap &pvmap, Face *f)
{
  /* to do
  VVertices reg;
  const HEdges &fe = f->getBoundary();
  for (int i = 0; i < fe.size(); ++i) {
    HEdges ed;
    fe[i]->loop(ed);
    int j0 = 0, n = ed.size();
    while (ed[n-1]->tail()->getP()->onLine(ed[n-1]->head()->getP(),
					   ed[j0]->head()->getP()))
      ++j0;
    Vertices *ve = new Vertices;
    ve->push_back(a->getVertex(ed[j0]->tail(), pvmap));
    int j = 0;
    while (j + 1 < n) {
      int k = j + 1;
      while (ed[j0]->ccw()->getF() == ed[j0+k]->ccw()->getF() &&
	     ve->back()->getP()->onLine(ed[(j0+j)%n]->head()->getP(),
					ed[(j0+k)%n]->head()->getP()))
	++k;
      ve->push_back(a->getVertex(ed[(j0+k)%n]->tail(), pvmap));
      j = k;
    }
    reg.push_back(ve);
  }
  a->addFace(reg, f->getPC());
  deleteRegion(reg);
  */
}  

// shapes

Polyhedron * box (double *b, bool perturb)
{
  Polyhedron *p = new Polyhedron(perturb);
  p->getVertex(b[0], b[2], b[4]);
  p->getVertex(b[1], b[2], b[4]);
  p->getVertex(b[1], b[3], b[4]);
  p->getVertex(b[0], b[3], b[4]);
  p->getVertex(b[0], b[2], b[5]);
  p->getVertex(b[1], b[2], b[5]);
  p->getVertex(b[1], b[3], b[5]);
  p->getVertex(b[0], b[3], b[5]);
  p->addTriangle(p->vertices[0], p->vertices[2], p->vertices[1]);
  p->addTriangle(p->vertices[0], p->vertices[3], p->vertices[2]);
  p->addTriangle(p->vertices[0], p->vertices[1], p->vertices[5]);
  p->addTriangle(p->vertices[0], p->vertices[5], p->vertices[4]);
  p->addTriangle(p->vertices[1], p->vertices[2], p->vertices[6]);
  p->addTriangle(p->vertices[1], p->vertices[6], p->vertices[5]);
  p->addTriangle(p->vertices[2], p->vertices[3], p->vertices[7]);
  p->addTriangle(p->vertices[2], p->vertices[7], p->vertices[6]);
  p->addTriangle(p->vertices[3], p->vertices[0], p->vertices[4]);
  p->addTriangle(p->vertices[3], p->vertices[4], p->vertices[7]);
  p->addTriangle(p->vertices[4], p->vertices[5], p->vertices[6]);
  p->addTriangle(p->vertices[4], p->vertices[6], p->vertices[7]);
  return p;
}

Polyhedron * sphere (double ox, double oy, double oz, double r, double err)
{
  Polyhedron *a = sphere(err/r);
  Polyhedron *b = new Polyhedron;
  PVMap pvmap;
  PV3 o = PV3::constant(ox, oy, oz);
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f) {
    Vertices va = (*f)->getBoundary()->loop(), vb;
    for (Vertices::iterator v = va.begin(); v != va.end(); ++v) {
      PVMap::iterator i = pvmap.find((*v)->getP());
      if (i == pvmap.end()) {
	PV3 p = o + r*(*v)->getP()->getApprox(1.0);
	Vertex *w = b->getVertex(p.x.mid(), p.y.mid(), p.z.mid());
	pvmap.insert(PVPair((*v)->getP(), w));
	vb.push_back(w);
      }
      else
	vb.push_back(i->second);
    }
    b->addTriangle(vb[0], vb[1], vb[2]);
  }
  delete a;
  return b;
}

Polyhedron * sphere (double err)
{
  Polyhedron *a = octohedron();
  while (sphereError(a) > err) {
    Polyhedron *b = sphereRefine(a);
    delete a;
    a = b;
  }
  return a;
}

Polyhedron * octohedron ()
{
  Polyhedron *a = new Polyhedron;
  a->getVertex(1.0, 0.0, 0.0);
  a->getVertex(-1.0, 0.0, 0.0);
  a->getVertex(0.0, 1.0, 0.0);
  a->getVertex(0.0, -1.0, 0.0);
  a->getVertex(0.0, 0.0, 1.0);
  a->getVertex(0.0, 0.0, -1.0);
  a->addTriangle(a->vertices[2], a->vertices[4], a->vertices[0]);
  a->addTriangle(a->vertices[1], a->vertices[4], a->vertices[2]);
  a->addTriangle(a->vertices[3], a->vertices[4], a->vertices[1]);
  a->addTriangle(a->vertices[0], a->vertices[4], a->vertices[3]);
  a->addTriangle(a->vertices[5], a->vertices[2], a->vertices[0]);
  a->addTriangle(a->vertices[5], a->vertices[1], a->vertices[2]);
  a->addTriangle(a->vertices[5], a->vertices[3], a->vertices[1]);
  a->addTriangle(a->vertices[5], a->vertices[0], a->vertices[3]);
  return a;
}

double sphereError (Polyhedron *a)
{
  double d = 1.0;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f) {
    Vertices ve = (*f)->getBoundary()->loop();
    PV3 p = (ve[0]->getP()->getApprox(1.0) + ve[1]->getP()->getApprox(1.0) 
			+ ve[2]->getP()->getApprox(1.0))/3.0;
    d = min(d, p.dot(p).mid());
  }
  return 1.0 - sqrt(d);
}

Polyhedron * sphereRefine (Polyhedron *a)
{
  Polyhedron *b = new Polyhedron;
  map<Edge *, Vertex *> evmap;
  for (Edges::iterator e = a->edges.begin(); e != a->edges.end(); ++e) {
    PV3 pm = (0.5*((*e)->getT()->getP()->getApprox(1.0) +
		   (*e)->getH()->getP()->getApprox(1.0))).unit();
    Vertex *vm = b->getVertex(pm.x.mid(), pm.y.mid(), pm.z.mid());
    evmap.insert(pair<Edge *, Vertex *>(*e, vm));
  }
  PVMap pvmap;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f) {
    HEdge *e = (*f)->getBoundary();
    Vertices ve = e->loop();
    Vertex *u = b->getVertex(ve[0], pvmap), *v = b->getVertex(ve[1], pvmap),
      *w = b->getVertex(ve[2], pvmap), *uv = evmap.find(e->getE())->second,
      *vw = evmap.find(e->getNext()->getE())->second,
      *wu = evmap.find(e->getNext()->getNext()->getE())->second;
    b->addTriangle(u, uv, wu);
    b->addTriangle(uv, v, vw);
    b->addTriangle(vw, w, wu);
    b->addTriangle(uv, vw, wu);
  }
  return b;
}

// n tetrahedra (t, t + a, t + b, t + c) with the coordinates of t in [0, u]
// and the coordinates of a, b, and c in [0, v]
Polyhedron * randomTets (int n, double u, double v)
{
  Polyhedron *a = new Polyhedron;
  for (int i = 0; i < n; ++i)
    randomTet(a, u, v);
  Polyhedron *b = subdivide(a, false), *c = triangulate(b);
  delete b;
  delete a;
  return c;
}

void randomTet (Polyhedron *a, double u, double v)
{
  PTR<Point> p0 = randomPoint(u),
    p[] = {p0, new SumPoint(p0, randomPoint(v)),
	   new SumPoint(p0, randomPoint(v)), new SumPoint(p0, randomPoint(v))};
  if (Orientation(p[0], p[1], p[2], p[3]) == -1)
    reverse(p, p + 4);
  Vertex *v1 = a->getVertex(p[0]), *v2 = a->getVertex(p[1]),
    *v3 = a->getVertex(p[2]), *v4 = a->getVertex(p[3]);
  a->addTriangle(v1, v2, v4);
  a->addTriangle(v2, v3, v4);
  a->addTriangle(v3, v1, v4);
  a->addTriangle(v3, v2, v1);
}

Point * randomPoint (double d)
{
  return new Point(randomNumber(0.0, d), randomNumber(0.0, d),
		   randomNumber(0.0, d), false);
}

// debug

void find (Point *p, const Vertices &vertices)
{
  for (int i = 0; i < vertices.size(); ++i)
    if (p == vertices[i]->getP())
      cerr << i << " ";
  cerr << endl;
}

int find (Vertex *v, const Vertices &vertices)
{
  for (int i = 0; i < vertices.size(); ++i)
    if (v == vertices[i])
      return i;
  return -1;
}

int find (Edge *e, const Edges &edges)
{
  for (int i = 0; i < edges.size(); ++i)
    if (e == edges[i])
      return i;
  return -1;
}

int find (Face *f, const Faces &faces)
{
  for (int i = 0; i < faces.size(); ++i)
    if (f == faces[i])
      return i;
  return -1;
}

int find (HFace *f, const HFaces &faces)
{
  for (int i = 0; i < faces.size(); ++i)
    if (f == faces[i])
      return i;
  return -1;
}

void find (ID id, const Faces &faces)
{
  for (int i = 0; i < faces.size(); ++i)
    if (faces[i]->getP()->getid() == id)
      cerr << i << " ";
  cerr << endl;
}

void find (ID id, const Faces &faces, Faces &res)
{
 for (int i = 0; i < faces.size(); ++i)
   if (faces[i]->getP()->getid() == id)
      res.push_back(faces[i]);
}

int find (Shell *s, const Shells &shells)
{
  for (int i = 0; i < shells.size(); ++i)
    if (s == shells[i])
      return i;
  return -1;
}

int find (Cell *c, const Cells &cells)
{
  for (int i = 0; i < cells.size(); ++i)
    if (c == cells[i])
      return i;
  return -1;
}

void pp1 (PV3 p)
{
  cerr << setprecision(16);
  cerr << "(" << p.x.mid() << " " << p.y.mid() << " " << p.z.mid() << ")";
}

void pp (PV3 p)
{
  pp1(p);
  cerr << endl;
}

void ppla (Plane *p)
{
  PV3 n = p->getApprox().n;
  Parameter k = p->getApprox().k;
  cerr << "(" << n.x.mid() << " " << n.y.mid() << " " << n.z.mid() << " "
       << k.mid() << ")" << endl;
}

void pp (Point *p)
{
  pp(p->getApprox());
}

void pps (const Points &pts)
{
  cerr << "(" << endl;
  for (Points::const_iterator p = pts.begin(); p != pts.end(); ++p)
    pp(*p);
  cerr << ")" << endl;
}

void pv (Vertex *v)
{
  pp(v->getP());
}

void pvs (const Vertices &ve)
{
  cerr << "(" << endl;
  for (Vertices::const_iterator v = ve.begin(); v != ve.end(); ++v)
    pv(*v);
  cerr << ")" << endl;
}

void ptr (const Triangle &t)
{
  cerr << "(";
  pp(t.a);
  pp(t.b);
  pp(t.c);
  cerr << ")";
}

void ptrs (const Triangles &tr)
{
  cerr << "(";
  for (Triangles::const_iterator t = tr.begin(); t != tr.end(); ++t)
    ptr(*t);
  cerr << ")" << endl;
}

void pids (const set<ID> &ids)
{
  for (set<ID>::const_iterator i = ids.begin(); i != ids.end(); ++i)
    cerr << *i << " ";
  cerr << endl;
}

void pe (Edge *e)
{
  cerr << "(";
  pp1(e->getT()->getP()->getApprox());
  pp1(e->getH()->getP()->getApprox());
  cerr << ")" << endl;
}

void pes (const Edges &ed)
{
  cerr << "(";
  for (Edges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e);
  cerr << ")" << endl;
}

void pe (HEdge *e)
{
  cerr << "(";
  pp1(e->tail()->getP()->getApprox());
  pp1(e->head()->getP()->getApprox());
  cerr << ")" << endl;
}

void pes (const HEdges &ed)
{
  cerr << "(";
  for (HEdges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e);
  cerr << ")" << endl;
}

void pl (HEdge *e)
{
  Vertices ve = e->loop();
  pvs(ve);
}

void pe (SEdge *e)
{
  cerr << "(";
  pp1(e->tail->getApprox());
  pp1(e->head()->getApprox());
  cerr << ")" << endl;
}

void pes (const SEdges &ed)
{
  cerr << "(";
  for (SEdges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e);
  cerr << ")" << endl;
}

void pl (SEdge *e)
{
  Points pts = e->loop();
  pps(pts);
}

void pf (Face *f)
{
  pvs(f->getBoundary()->loop());
}

void pfs (const Faces &fa)
{
  cerr << "(";
  for (Faces::const_iterator f = fa.begin(); f != fa.end(); ++f)
    pf(*f);
  cerr << ")" << endl;
}

void pfs (const HFaces &fa)
{
  cerr << "(";
  for (HFaces::const_iterator f = fa.begin(); f != fa.end(); ++f)
    pf((*f)->getF());
  cerr << ")" << endl;
}

TrianglePlane * tpl (Plane *p)
{
  return dynamic_cast<TrianglePlane *>(p);
}

SumPoint * spt (Point *v)
{
  return dynamic_cast<SumPoint *>(v);
}

EPPoint * eppt (Point *v)
{
  return dynamic_cast<EPPoint *>(v);
}		      

EEPoint * eept (Point *v)
{
  return dynamic_cast<EEPoint *>(v);
}

double distance (Point *v, Point *w)
{
  PV3 u = v->getApprox() - w->getApprox();
  return sqrt(fabs(u.dot(u).mid()));
}

double distance (Point *a, Point *t, Point *h)
{
  PV3 tp = t->getApprox(), u = h->getApprox() - tp,
    p = tp + ((a->getApprox() - tp).dot(u)/u.dot(u))*u, w = a->getApprox() - p;
  Parameter ww = w.dot(w);
  return sqrt(fabs(ww.mid()));
}

double distance (Point *v, Plane *p)
{
  PV3 n = p->getApprox().n;
  double k = sqrt(fabs(n.dot(n).mid()));
  Parameter d = n.dot(v->getApprox()) + p->getApprox().k;
  return d.mid()/k;
}

double distance (Point *v, Point *a, Point *b, Point *c)
{
  PV3 n = (c->getApprox() - b->getApprox()).cross(a->getApprox() - b->getApprox());
  double k = sqrt(fabs(n.dot(n).mid()));
  Parameter d = n.dot(v->getApprox() - b->getApprox());
  return d.mid()/k;
}

double distanceEE (Point *a, Point *b, Point *c, Point *d)
{
  PV3 u = b->getApprox() - a->getApprox(), v = d->getApprox() - c->getApprox(),
    w = u.cross(v).unit();
  return (c->getApprox() - a->getApprox()).dot(w).mid();
}

double distance (Edge *e, Edge *f)
{
  return distanceEE(e->getT()->getP(), e->getH()->getP(),
		    f->getT()->getP(), f->getH()->getP());
}

void edgesNM (Shell *s)
{
  set<Edge *> es;
  const HFaces &hf = s->getHFaces();
  for (HFaces::const_iterator f = hf.begin(); f != hf.end(); ++f) {
    HEdges he = (*f)->getF()->getBoundary()->edgeLoop();
    for (HEdges::iterator e = he.begin(); e != he.end(); ++e)
      es.insert((*e)->getE());
  }
  for (set<Edge *>::iterator e = es.begin(); e != es.end(); ++e) {
    int n = 0;
    for (int i = 0; i < (*e)->HEdgesN(); ++i) {
      Face *f = (*e)->getHEdge(i)->getF();
      for (int j = 0; j < 2; ++j)
	if (f->getHFace(j)->getS() == s)
	  ++n;
    }
    if (n != 2)
      cerr << "non-manifold edge " << n << endl;
  }
}

void edgesNM (Polyhedron *a, Edges &ed)
{
  for (int i = 0; i < a->edges.size(); ++i) {
    Edge *e = a->edges[i];
    int n = e->HEdgesN();
    if (n == 0) {
      cerr << "no hedges " << i << endl;
      ed.push_back(e);
    }
    else if (n == 1) {
      cerr << "dangling edge " << i << endl;
      ed.push_back(e);
    }
    else if (n > 2) {
      cerr << "non-manifold edge " << i << " hedges " << n << endl;
      ed.push_back(e);
    }
    else if (e->getHEdge(0)->getForward() == e->getHEdge(1)->getForward()) {
      cerr << "same sign hedges" << endl;
      ed.push_back(e);
    }
  }
}

void edgesNM (Polyhedron *a)
{
  Edges ed;
  edgesNM(a, ed);
}

double * rotationMatrix (double t, double *u)
{
  double c = cos(t), d = 1.0 - c, s = sin(t), ux = u[0], uy = u[1], uz = u[2],
    *m = new double[9];
  m[0] = ux*ux*d + c;
  m[1] = ux*uy*d + uz*s;
  m[2] = ux*uz*d - uy*s;
  m[3] = ux*uy*d - uz*s;
  m[4] = uy*uy*d + c;
  m[5] = uy*uz*d + ux*s;
  m[6] = ux*uz*d + uy*s;
  m[7] = uy*uz*d - ux*s;
  m[8] = uz*uz*d + c;
  return m;
}

Point * rotate (Point *p, double *m)
{
  PV3 q = p->getApprox();
  double x = m[0]*q.x.mid() + m[3]*q.y.mid() + m[6]*q.z.mid(),
    y = m[1]*q.x.mid() + m[4]*q.y.mid() + m[7]*q.z.mid(),
    z = m[2]*q.z.mid() + m[5]*q.y.mid() + m[8]*q.z.mid();
  return new Point(x, y, z, false);
}

Polyhedron * rotate (Polyhedron *a, double t,  double *u)
{
  Polyhedron *b = new Polyhedron(false);
  PVMap pvmap;
  double *m = rotationMatrix(t, u);
  for (Vertices::const_iterator v = a->vertices.begin(); v != a->vertices.end(); ++v)
    pvmap.insert(PVPair((*v)->getP(), b->getVertex(rotate((*v)->getP(), m))));
  for (Faces::const_iterator f = a->faces.begin(); f != a->faces.end(); ++f) {
    Points pf = (*f)->getBoundary()->pointLoop();
    b->addTriangle(pf[0], pf[1], pf[2], pvmap);
  }
  delete [] m;
  return b;
}

int side (Plane *p, Point *a)
{
  return Side(p, a);
}
