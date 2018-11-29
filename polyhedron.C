#include "polyhedron.h"

double getTime ()
{
  timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1e-6*tv.tv_usec;
}

Parameter cross (const PV3 &a, const PV3 &b, const int coord)
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

ID Plane::planeid = 1u;

int Plane::projectionCoordinate ()
{
  PV3 n = getApprox(1.0).n;
  double x = n.x.mid(), y = n.y.mid(), z = n.z.mid(),
    ax = fabs(x), ay = fabs(y), az = fabs(z);
  if (ax >= ay && ax >= az)
    return x > 0.0 ? 1 : -1;
  if (ay >= ax && ay >= az)
    return y > 0.0 ? 2 : -2;
  return z > 0.0 ? 3 : -3;
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
  HFaces hf;
  HEdge *e = f->h;
  do {
    hf.push_back(neighbor(e));
    e = e->next;
  }
  while (e != f->h);
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

Face::Face (HEdge *h, int pc)
  : h(h), pc(pc)
{
  p = new TrianglePlane(h->tail()->p, h->head()->p, h->next->head()->p);
  hfaces[0].f = hfaces[1].f = this;
  HEdges ed = h->edgeLoop();
  copyBBox(ed[0]->tail()->bbox, bbox);
  for (int i = 1; i < ed.size(); ++i)
    mergeBBox(ed[i]->tail()->bbox, bbox);
  for (int i = 0; i < ed.size(); ++i) {
    ed[i]->f = this;
    ed[i]->tail()->p->ps.insert(p->id);
  }
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
    pc = p->projectionCoordinate();
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
  return SamePlane(p, f->p) == 0;
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

/*
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
	    em = *h;
	    fm = (*h)->f->hfaces + i;
	  }
	  break;
	}
  return !fm->getF()->coplanar(fm->neighbor(em)->getF())
    && Convex(em, fm) == 1;
}
*/

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
  HFace *hf = getShell(0)->hfaces[0];
  Face *f = hf->f;
  PTR<Point> p = f->centroid(), n = new HFaceNormal(hf), qmin = 0;
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
  return new CentroidPoint(p, qmin);
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

Face * Polyhedron::addTriangle (Vertex *ta, Vertex *tb, Vertex *tc,
				PVMap &pvmap, int pc)
{
  Vertex *a = getVertex(ta, pvmap), *b = getVertex(tb, pvmap),
    *c = getVertex(tc, pvmap);
  return addTriangle(a, b, c, pc);
}

Face * Polyhedron::addTriangle (Face *f, PVMap &pvmap)
{
  Vertices vf = f->h->loop();
  int pc = f->getPC();
  return addTriangle(vf[0], vf[1], vf[2], pvmap, pc);
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
    Vertices ve = (*f)->h->loop();
    a->addTriangle(ve[2], ve[1], ve[0], pvmap);
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
    Vertices ve = (*f)->h->loop();
    a->addTriangle(ve[2], ve[1], ve[0], pvmap);
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
      Vertices vf = (*f)->h->loop();
      int pc = (*f)->getPC();
      if (in2)
	d->addTriangle(vf[0], vf[1], vf[2], pvmap, pc);
      else
	d->addTriangle(vf[2], vf[1], vf[0], pvmap, - pc);
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
  map<Face *, Edges> fesmap;
  intersectFF(a, epsmap, ffemap, fesmap);
  double t0 = getTime();
  for (map<Edge *, Points *>::iterator i = epsmap.begin(); i != epsmap.end(); ++i)
    sortPoints(i->first, i->second);
  map<Face *, vector<FFE> > femap = feMap(epsmap, ffemap, fesmap);
  intersectFFF(ffemap, femap);
  double t1 = getTime() - t0;
  t0 = getTime();
  for (map<pair<Face *, Face *>, FFE>::iterator i = ffemap.begin();
       i != ffemap.end(); ++i)
    sortPoints(i->second.pts);
  Polyhedron *b = subfaces(a, oneway, epsmap, femap);
  for (map<Edge *, Points *>::iterator i = epsmap.begin(); i != epsmap.end(); ++i)
    delete i->second;
  deleteFEmap(femap);
  double t2 = getTime() - t0;
  cerr << "; fff = " << t1 << "; sub " << t2 << endl;
  return b;
}

void intersectFF (Polyhedron *a, map<Edge *, Points *> &epsmap,
		  map<pair<Face *, Face *>, FFE> &ffemap,
		  map<Face *, Edges> &fesmap)
{
  double t0 = getTime();
  Octree<Face *> *octree = a->faceOctree();
  vector<pair<Face *, Face *> > ff1, ff;
  octree->pairs(ff1);
  delete octree;
  map<pair<Face *, Edge *>, PTR<Point> > fepmap;
  map<Face *, Points *> fpsmap;
  intersectFE(ff1, fepmap, epsmap, fesmap, fpsmap, ff);
  double t1 = getTime() - t0;
  t0 = getTime();
  intersectFF(ff, a->perturbed, fepmap, epsmap, fpsmap, ffemap);
  for (map<Face *, Points *>::iterator i = fpsmap.begin(); i != fpsmap.end(); ++i)
    delete i->second;
  double t2 = getTime() - t0;
  cerr << "fe = " << t1 << "; ff = " << t2;
}

void intersectFE (const vector<pair<Face *, Face *> > &ff1,
		  map<pair<Face *, Edge *>, PTR<Point> > &fepmap,
		  map<Edge *, Points *> &epsmap,
		  map<Face *, Edges> &fesmap,
		  map<Face *, Points *> &fpsmap,
		  vector<pair<Face *, Face *> > &ff)
{
  static unsigned int n = 8;
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
    for (vector<pair<pair<Face *, Edge *>, PTR<Point> > >::iterator
	   x = fed[i].fep.begin(); x != fed[i].fep.end(); ++x)
      fepmap.insert(*x);
    for (vector<pair<Edge *, PTR<Point> > >::iterator x = fed[i].ep.begin();
	 x != fed[i].ep.end(); ++x)
      update(x->first, x->second, epsmap);
    for (vector<pair<Face *, Edge *> >::iterator x = fed[i].fe.begin();
       x != fed[i].fe.end(); ++x)
      update(x->first, x->second, fesmap);
    for (vector<pair<Face *, PTR<Point> > >::iterator x = fed[i].fp.begin();
	 x != fed[i].fp.end(); ++x)
      update(x->first, x->second, fpsmap);
    ff.insert(ff.end(), fed[i].ff.begin(), fed[i].ff.end());
  }
}

void intersectFET (void *ptr)
{
  FEData *fed = (FEData *) ptr;
  BaseObject::addThread(fed->i);
  for (int i = fed->is; i < fed->ie; ++i) {
    const pair<Face *, Face *> &ff = fed->ff1->at(i);
    intersectFE(ff.first, ff.second, fed->ep, fed->fep, fed->fe, fed->fp, fed->ff);
  }
}

void intersectFE (Face *f, Face *g,
		  vector<pair<Edge *, PTR<Point> > > &ep,
		  vector<pair<pair<Face *, Edge *>, PTR<Point> > > &fep,
		  vector<pair<Face *, Edge *> > &fe,
		  vector<pair<Face *, PTR<Point> > > &fp,
		  vector<pair<Face *, Face *> > &ff)
{
  if (f->coplanar(g)) {
    intersectFEP(f, g, ep, fe);
    intersectFEP(g, f, ep, fe);
    return;
  }
  if (f->sharedEdge(g))
    return;
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
		   vector<pair<Edge *, PTR<Point> > > &ep,
		   vector<pair<Face *, Edge *> > &fe)
{
  HEdges eg = g->getBoundary()->edgeLoop();
  for (HEdges::iterator e = eg.begin(); e != eg.end(); ++e)
    intersectFEP(f, *e, ep, fe);
}

void intersectFEP (Face *f, HEdge *h,
		   vector<pair<Edge *, PTR<Point> > > &ep,
		   vector<pair<Face *, Edge *> > &fe)
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
		  vector<pair<Edge *, PTR<Point> > > &ep)
{
  Points pe, pf;
  intersectEE(e->getT()->getP(), e->getH()->getP(), f->getT()->getP(),
	      f->getH()->getP(), c, pe, pf);
  for (Points::iterator p = pe.begin(); p != pe.end(); ++p)
    ep.push_back(pair<Edge *, PTR<Point> >(e, *p));
  for (Points::iterator p = pf.begin(); p != pf.end(); ++p)
    ep.push_back(pair<Edge *, PTR<Point> >(f, *p));
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
  int tp3 = LeftTurn(ft, et, eh, c), tp4 = LeftTurn(fh, et, eh, c);
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
		  vector<pair<Edge *, PTR<Point> > > &ep,
		  vector<pair<pair<Face *, Edge *>, PTR<Point> > > &fep,
		  vector<pair<Face *, Edge *> > &fe,
		  vector<pair<Face *, PTR<Point> > > &fp)
{
  HEdges eg = g->getBoundary()->edgeLoop();
  int n = eg.size();
  for (int i = 0; i < n; ++i)
    if ((sg[i] == 0 || sg[i] == 2) && (sg[(i+1)%n] == 0 || sg[(i+1)%n] == 2)) {
      intersectFEP(f, eg[i], ep, fe);
      return;
    }
  for (int i = 0; i < n; ++i)
    if (sg[i] == 0)
      intersectFV(f, eg[i], ep, fp);
    else if (sg[i] == - sg[(i+1)%n])
      intersectFEG(f, eg[i], ep, fep);
}

void intersectFV (Face *f, HEdge *h,
		  vector<pair<Edge *, PTR<Point> > > &ep,
		  vector<pair<Face *, PTR<Point> > > &fp)
{
  int ie = -1;
  PTR<Point> p = h->tail()->getP();
  if (f->contains(p, false, &ie)) {
    if (ie == -1)
      fp.push_back(pair<Face *, PTR<Point> >(f, p));
    else {
      HEdges ef = f->getBoundary()->edgeLoop();
      ep.push_back(pair<Edge *, PTR<Point> >(ef[ie]->getE(), p));
    }
  }
}

void intersectFEG (Face *f, HEdge *h,
		   vector<pair<Edge *, PTR<Point> > > &ep,
		   vector<pair<pair<Face *, Edge *>, PTR<Point> > > &fep)
{
  if (!first(h))
    return;
  Edge *e = h->getE();
  if (!bboxOverlap(f->getBBox(), e->getBBox()))
    return;
  PTR<Point> p = new EPPoint(e->getT()->getP(), e->getH()->getP(), f->getP());
  int ie = -1;
  if (f->contains(p, false, &ie)) {
    ep.push_back(pair<Edge *, PTR<Point> >(e, p));
    if (ie == -1) {
      pair<Face *, Edge *> fe(f, e);
      fep.push_back(pair<pair<Face *, Edge *>, PTR<Point> >(fe, p));
    }
    else {
      HEdges ef = f->getBoundary()->edgeLoop();
      ep.push_back(pair<Edge *, PTR<Point> >(ef[ie]->getE(), p));
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

void update (Face *f, Edge *e, map<Face *, Edges> &fesmap)
{
  map<Face *, Edges>::iterator i = fesmap.find(f);
  if (i == fesmap.end()) {
    Edges ed;
    ed.push_back(e);
    fesmap.insert(pair<Face *, Edges>(f, ed));
  }
  else
    i->second.push_back(e);
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

void intersectFF (const vector<pair<Face *, Face *> > &ff, bool perturbed,
		  const map<pair<Face *, Edge *>, PTR<Point> > &fepmap,
		  const map<Edge *, Points *> &epsmap,
		  const map<Face *, Points *> &fpsmap,
		  map<pair<Face *, Face *>, FFE> &ffemap)
{
  vector<pair< pair<Face *, Face *>, pair<PTR<Point>, PTR<Point> > > > ffpp;
  for (vector<pair<Face *, Face *> >::const_iterator i = ff.begin();
       i != ff.end(); ++i)
    intersectFF(i->first, i->second, perturbed, fepmap, epsmap, fpsmap, ffpp);
  for (vector<pair< pair<Face *, Face *>, pair<PTR<Point>, PTR<Point> > > >::iterator
	 i = ffpp.begin(); i != ffpp.end(); ++i)
    update(i->first.first, i->first.second, i->second.first, i->second.second,
	   ffemap);
}

void intersectFF (Face *f, Face *g, bool perturbed,
		  const map<pair<Face *, Edge *>, PTR<Point> > &fepmap,
		  const map<Edge *, Points *> &epsmap,
		  const map<Face *, Points *> &fpsmap,
		  vector<pair< pair<Face *, Face *>, pair<PTR<Point>, PTR<Point> > > > &ffpp)
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
  pair<PTR<Point>, PTR<Point> > pq(pts[flag ? 0 : 1], pts[flag ? 1 : 0]);
  ffpp.push_back(pair<pair<Face *, Face *>, pair<PTR<Point>, PTR<Point> > >(fg, pq));
}

void findFE (Face *f, Face *g,
	     const map<pair<Face *, Edge *>, PTR<Point> > &fepmap,
	     Points &pts)
{
  HEdges eg = g->getBoundary()->edgeLoop();
  for (HEdges::iterator h = eg.begin(); h != eg.end(); ++h) {
    pair<Face *, Edge *> fe(f, (*h)->getE());
    map<pair<Face *, Edge *>, PTR<Point> >::const_iterator i = fepmap.find(fe);
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

map<Face *, vector<FFE> > feMap (map<Edge *, Points *> &epsmap,
				const map<pair<Face *, Face *>, FFE> &ffemap,
				map<Face *, Edges> &fesmap)
{
  map<Face *, vector<FFE> > femap;
  for (map<pair<Face *, Face *>, FFE>::const_iterator i = ffemap.begin();
       i != ffemap.end(); ++i) {
    update(i->first.first, i->second, femap);
    update(i->first.second, i->second, femap);
  }
  for (map<Face *, Edges>::iterator i = fesmap.begin(); i != fesmap.end(); ++i)
    for (Edges::iterator e = i->second.begin(); e != i->second.end(); ++e)
      feMapFE(epsmap, i->first, *e, femap);
  return femap;
}

void feMapFE (const map<Edge *, Points *> &epsmap,
	      Face *f, Edge *e, map<Face *, vector<FFE> > &femap)
{
  HEdges ef = f->getBoundary()->edgeLoop();
  for (HEdges::iterator h = ef.begin(); h != ef.end(); ++h)
    if ((*h)->tail()->getP()->onLine(e->getT()->getP(), e->getH()->getP()) &&
	(*h)->head()->getP()->onLine(e->getT()->getP(), e->getH()->getP()))
      return;
  set<Point *> vpts1 = boundaryVertices(f, epsmap);
  set<Point *, PointOrderRP> vpts(vpts1.begin(), vpts1.end());
  Points epts;
  epts.push_back(e->getT()->getP());
  map<Edge *, Points *>::const_iterator i = epsmap.find(e);
  if (i != epsmap.end())
    epts.insert(epts.end(), i->second->begin(), i->second->end());
  epts.push_back(e->getH()->getP());
  int is = 0, n = epts.size();
  while (is < n && vpts.find(epts[is]) == vpts.end() &&
	 !f->contains(epts[is], true))
    ++is;
  if (is == n)
    return;
  while (is + 1 < n && epts[is]->identical(epts[is+1]))
    ++is;
  int ie = is + 1;
  while (ie < n && f->contains(epts[ie], true))
    ++ie;
  if (ie < n && vpts.find(epts[ie]) != vpts.end())
    ++ie;
  if (ie == is + 1)
    return;
  FFE ffe(f, 0, epts[is], epts[is+1]);
  for (int i = is + 2; i < ie; ++i)
    ffe.pts->push_back(epts[i]);
  update(f, ffe, femap);
}

void update (Face *f, const FFE &ffe, map<Face *, vector<FFE> > &femap)
{
  map<Face *, vector<FFE> >::iterator i = femap.find(f);
  if (i == femap.end()) {
    vector<FFE> ffes;
    ffes.push_back(ffe);
    femap.insert(pair<Face *, vector<FFE> >(f, ffes));
  }
  else
    i->second.push_back(ffe);
}

void intersectFFF (const map<pair<Face *, Face *>, FFE> &ffemap,
		   const map<Face *, vector<FFE> > &femap)
{
  vector<pair<Points *, PTR<Point> > > psps;
  intersectFFFAux(ffemap, femap, psps);
  for (vector<pair<Points *, PTR<Point> > >::iterator i = psps.begin();
       i != psps.end(); ++i)
    i->first->push_back(i->second);
}

void intersectFFFAux (const map<pair<Face *, Face *>, FFE> &ffemap,
		      const map<Face *, vector<FFE> > &femap,
		      vector<pair<Points *, PTR<Point> > > &psps)
{
  static unsigned int n = 8;
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
  double t0 = getTime();
  FFFData *fd = (FFFData *) ptr;
  BaseObject::addThread(fd->i);
  map<Face *, vector<FFE> >::const_iterator x = fd->femap->begin();
  for (int i = 0; i < fd->is; ++i, ++x)
    ;
  for (int i = fd->is; i < fd->ie; ++i, ++x)
    intersectFFF(x->first, x->second, *fd->ffemap, fd->psps);
}

void intersectFFF (Face *f, const vector<FFE> &ed,
		   const map<pair<Face *, Face *>, FFE> &ffemap,
		   vector<pair<Points *, PTR<Point> > > &psps)
{
  int n = ed.size(), pc = f->getPC();
  for (int i = 0; i + 1 < n; ++i) {
    const FFE &ei = ed[i];
    Face *fi = ei.f == f ? ei.g : ei.f;
    for (int j = i + 1; j < n; ++j) {
      const FFE &ej = ed[j];
      Face *fj = ej.f == f ? ej.g : ej.f;
      if (bboxOverlap(ei.bbox, ej.bbox)) {
	pair<Face *, Face *> fifj(fi < fj ? fi : fj, fi < fj ? fj : fi);
	map<pair<Face *, Face *>, FFE>::const_iterator y = ffemap.find(fifj);
	if ((fi == 0 || f < fi) && (fj == 0 || f < fj) || y == ffemap.end()) {
	  Points pi, pj;
	  intersectEE(ei.pts->at(0), ei.pts->at(1), ej.pts->at(0), ej.pts->at(1),
		      pc, pi, pj);
	  for (Points::iterator p = pi.begin(); p != pi.end(); ++p)
	    psps.push_back(pair<Points *, PTR<Point> >(ei.pts, *p));
	  for (Points::iterator p = pj.begin(); p != pj.end(); ++p)
	    psps.push_back(pair<Points *, PTR<Point> >(ej.pts, *p));
	  if (!pi.empty() && y != ffemap.end())
	    psps.push_back(pair<Points *, PTR<Point> >(y->second.pts, pi[0]));
	}
      }
    }
  }
}

void sortPoints (Edge *e, Points *pts)
{
  if (pts->size() > 1)
    sort(pts->begin(), pts->end(),
	 PointOrderPPP(e->getT()->getP(), e->getH()->getP()));
}

void sortPoints (Points *pts)
{
  if (pts->size() > 2) {
    Point *t = pts->at(0), *h = pts->at(1);
    sort(pts->begin(), pts->end(), PointOrderPPP(t, h));
  }
}

Polyhedron * subfaces (Polyhedron *a, bool oneway,
		       const map<Edge *, Points *> &epsmap,
		       const map<Face *, vector<FFE> > &femap)
{
  PTriangles tr;
  subfacesAux(a, oneway, epsmap, femap, tr);
  Polyhedron *b = new Polyhedron(a->perturbed);
  PVMap pvmap;
  for (PTriangles::iterator t = tr.begin(); t != tr.end(); ++t) {
    Vertex *u = b->getVertex(t->a, pvmap), *v = b->getVertex(t->b, pvmap),
      *w = b->getVertex(t->c, pvmap);
    if (!faceVertices(u, v, w) && !faceVertices(u, w, v))
      b->addTriangle(u, v, w);
  }
  return b;
}

void subfacesAux (Polyhedron *a, bool oneway,
		  const map<Edge *, Points *> &epsmap,
		  const map<Face *, vector<FFE> > &femap,
		  PTriangles &tr)
{
  static unsigned int n = 16;
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
  for (int i = 0; i < n; ++i)
    tr.insert(tr.end(), sud[i].tr.begin(), sud[i].tr.end());
}

void subfacesT (void *ptr)
{
  SUData *sud = (SUData *) ptr;
  BaseObject::addThread(sud->i);
  for (int i = sud->is; i < sud->ie; ++i)
    subfaces(sud->a->faces[i], sud->oneway, *sud->epsmap, *sud->femap,
	     sud->a->perturbed, sud->tr);
}

void subfaces (Face *f, bool oneway, const map<Edge *, Points *> &epsmap,
	       const map<Face *, vector<FFE> > &femap,
	       bool perturbed, PTriangles &tr)
{
  Polyhedron *a = new Polyhedron(perturbed);
  PVMap pvmap;
  int c = f->getPC();
  HEdges he = subedges(f, oneway, epsmap, femap, a, pvmap), outer, inner;
  for (HEdges::iterator e = he.begin(); e != he.end(); ++e)
    if (!(*e)->getFlag()) {
      HEdge *l = findLoop(*e);
      if (l)
	if (LeftTurn(l->tail()->getP(), l->head()->getP(),
		     l->getNext()->head()->getP(), c) == 1)
	  outer.push_back(l);
	else 
	  inner.push_back(l);
    }
  sort(outer.begin(), outer.end(), HEdgeOrder());
  sort(inner.begin(), inner.end(), HEdgeOrder());
  vector<HEdges> nfaces;
  for (HEdges::iterator e = outer.begin(); e != outer.end(); ++e) {
    HEdges he;
    he.push_back(*e);
    nfaces.push_back(he);
  }
  for (HEdges::iterator e = inner.begin(); e != inner.end(); ++e)
    addInner(*e, oneway, c, nfaces);
  for (vector<HEdges>::iterator i = nfaces.begin(); i != nfaces.end(); ++i)
    triangulate(*i, c, tr);
  delete a;
}

HEdges subedges (Face *f, bool oneway,
		 const map<Edge *, Points *> &epsmap,
		 const map<Face *, vector<FFE> > &femap,
		 Polyhedron *a, PVMap &pvmap)
{
  set<pair<Vertex *, Vertex *> > vv;
  HEdges ef = f->getBoundary()->edgeLoop(), ed;
  for (HEdges::iterator h = ef.begin(); h != ef.end(); ++h)
    subedges(*h, epsmap, a, pvmap, vv);
  map<Face *, vector<FFE> >::const_iterator i = femap.find(f);
  if (i != femap.end())
    for (vector<FFE>::const_iterator j = i->second.begin();
	 j != i->second.end(); ++j) {
      bool oflag = !oneway || !j->g,
	fflag = oflag || j->f == f, gflag = oflag || j->f != f;
      subedges(*j->pts, fflag, gflag, a, pvmap, vv);
    }
  for (set<pair<Vertex *, Vertex *> >::iterator v = vv.begin();
       v != vv.end(); ++v)
    ed.push_back(a->addHEdge(v->first, v->second));
  setNext(ed, f->getPC());
  return ed;
}

void subedges (HEdge *h, const map<Edge *, Points *> &epsmap,
	       Polyhedron *a, PVMap &pvmap, set<pair<Vertex *, Vertex *> > &vv)
{
  Edge *e = h->getE();
  Points pts;
  pts.push_back(e->getT()->getP());
  map<Edge *, Points *>::const_iterator i = epsmap.find(e);
  if (i != epsmap.end())
    pts.insert(pts.end(), i->second->begin(), i->second->end());
  pts.push_back(e->getH()->getP());
  subedges(pts, h->getForward(), !h->getForward(), a, pvmap, vv);
}

void subedges (const Points &pts, bool fflag, bool bflag, Polyhedron *a,
	       PVMap &pvmap, set<pair<Vertex *, Vertex *> > &vv)
{
  Vertex *v = a->getVertex(pts[0], pvmap);
  for (int i = 1; i < pts.size(); ++i) {
    Vertex *w = a->getVertex(pts[i], pvmap);
    if (v != w) {
      if (fflag)
	vv.insert(pair<Vertex *, Vertex *>(v, w));
      if (bflag)
	vv.insert(pair<Vertex *, Vertex *>(w, v));
    }
    v = w;
  }
}

void setNext (const HEdges &he, int c)
{
  HHEdges hhe;
  for (HEdges::const_iterator e = he.begin(); e != he.end(); ++e) {
    hhe.push_back(HHEdge(*e, true));
    hhe.push_back(HHEdge(*e, false));
  }
  sort(hhe.begin(), hhe.end(), HHEdgeOrder(c));
  int i = 0, n = hhe.size();
  while (i < n) {
    int m = 1;
    while (i + m < n && hhe[i].tail() == hhe[i+m].tail())
      ++m;
    int j = 0;
    while (j < m) {
      int k = i + (j + 1)%m;
      if (!hhe[i+j].f && hhe[k].f)
	hhe[i+j].e->setNext(hhe[k].e);
      ++j;
    }
    i += m;
  }
}

HEdge * findLoop (HEdge *h)
{
  HEdge *l = h, *h0 = h;
  do {
    h->setFlag(true);
    if (HEdgeOrder()(h, l))
      l = h;
    h = h->getNext();
    if (!h)
      return 0;
  }
  while (h != h0);
  return l;
}

void addInner (HEdge *e, bool oneway, int c, vector<HEdges> &fa)
{
  Point *p = e->head()->getP();
  for (int i = fa.size() - 1; i >= 0; --i) {
    HEdges &f = fa[i];
    if (contains(f[0], c, p)) {
      if (oneway)
	for (int j = 1; j < f.size(); ++j)
	  if (contains(f[j], c, p))
	    return;
      f.push_back(e);
      return;
    }
  }
}

bool contains (HEdge *e, int c, Point *p)
{
  Vertices ve = e->loop();
  Points pts;
  for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
    pts.push_back((*v)->getP());
  return contains(pts, c, p, true, 0);
}
    
void triangulate (const HEdges &fa, int c, PTriangles &ptr)
{
  Triangles tr;
  vector<Vertices> vv;
  for (HEdges::const_iterator h = fa.begin(); h != fa.end(); ++h)
    vv.push_back((*h)->loop());
  if (vv.size() == 1 && vv[0].size() == 3)
    tr.push_back(Triangle(vv[0][0], vv[0][1], vv[0][2]));
  else {
    VVertices reg;
    for (vector<Vertices>::iterator v = vv.begin(); v != vv.end(); ++v)
      reg.push_back(&*v);
    triangulate(reg, c, tr);
  }
  for (Triangles::iterator t = tr.begin(); t != tr.end(); ++t)
    ptr.push_back(PTriangle(t->a->getP(), t->b->getP(), t->c->getP()));
}

void deleteFEmap (map<Face *, vector<FFE> > &femap)
{
  for (map<Face *, vector<FFE> >::iterator i = femap.begin(); i != femap.end(); ++i)
    for (vector<FFE>::iterator j = i->second.begin(); j != i->second.end(); ++j)
      if (j->f == i->first)
	delete j->pts;
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
  int m = 0;
  for (int i = 0; i < n; ++i) {
    m += poly[i]->faces.size();
    for (int j = 0; j < poly[i]->faces.size(); ++j)
      if (poly[i]->faces[j]->getBoundary()->loop().size() > 3)
	cerr << "bad" << endl;
  }      

  Polyhedron *c = new Polyhedron(false);
  PVMap pvmap;
  for (int i = 0; i < n; ++i) {
    const Faces &fa = poly[i]->faces;
    for (Faces::const_iterator f = fa.begin(); f != fa.end(); ++f) {
      Vertices ve = (*f)->getBoundary()->loop();
      c->addTriangle(ve[0], ve[1], ve[2], pvmap);
    }
  }
  Polyhedron *d = subdivide(c, false);
  delete c;
  return d;
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
    for (int i = 0; !flag && i < n; ++i)
      flag = poly[i]->contains(p);
    if (flag)
      cin.insert(ci);
  }
  PVMap pvmap;
  for (Faces::iterator f = c->faces.begin(); f != c->faces.end(); ++f) {
    bool in1 = cin.find((*f)->getHFace(0)->getS()->getC()) != cin.end(),
      in2 = cin.find((*f)->getHFace(1)->getS()->getC()) != cin.end();
    if (in1 != in2) {
      Vertices vf = (*f)->getBoundary()->loop();
      int pc = (*f)->getPC();
      if (in2)
	d->addTriangle(vf[0], vf[1], vf[2], pvmap, pc);
      else
	d->addTriangle(vf[2], vf[1], vf[0], pvmap, - pc);
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
  vector<set<Face *> > fg = groupFaces(a);
  Polyhedron *b = new Polyhedron(a->perturbed);
  PVMap pvmap;
  for (vector<set<Face *> >::iterator f = fg.begin(); f != fg.end(); ++f)
    coalesceFace(b, pvmap, *f);
  return b;
}

vector<set<Face *> > groupFaces (Polyhedron *a)
{
  vector<set<Face *> > res;
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
  Polyhedron *b = subdivide(a, false);
  delete a;
  return b;
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
  pv(t.a);
  pv(t.b);
  pv(t.c);
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

void pes (const HHEdges &ed)
{
  cerr << "(";
  for (HHEdges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(e->e);
  cerr << ")" << endl;
}

void pl (HEdge *e)
{
  Vertices ve = e->loop();
  pvs(ve);
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
    Vertices ve = (*f)->getBoundary()->loop();
    b->addTriangle(ve[0], ve[1], ve[2], pvmap);
  }
  delete [] m;
  return b;
}

int po (Point *a, Point *b, Point *r)
{
  return PointOrder(a, b, r);
}
