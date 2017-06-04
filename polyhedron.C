#include "polyhedron.h"

bool inputPerturbed = true;

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
  }
}

IDSet intersection (const IDSet &a, const IDSet &b)
{
  IDSet ab;
  set_intersection(a.begin(), a.end(), b.begin(), b.end(),
		   inserter(ab, ab.begin()));
  return ab;
}

ID Plane::planeid = 1u;

int SamePlane::sign ()
{
  PV3 np = p->getN(), nq = q->getN();
  Parameter kp = p->getK(), kq = q->getK();
  int kps = kp.sign(), kqs = kq.sign();
  if (kps == 0 || kqs == 0)
    return kps == 0 && kqs == 0 && Rdir.getP().tripleProduct(np, nq).sign() ? 1 : -1;
  return Rdir.getP().dot(kp*nq - kq*np).sign() == 0 ? 1 : -1;
}

int projectionCoordinate (Plane *p)
{
  return ProjectionCoordinate(p).getPC();
}

void Point::decref ()
{
  if (ref > 0u)
    --ref; 
  if (ref > 0u)
    return;
  delete this;
}

void Point::getBBox (double *bbox)
{
  PV3 p = getP();
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
  if (inputPerturbed) 
    return PointOrderR(this, a);
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
  PV3 p = getP(), q = a->getP();
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

int TripleProduct::sign ()
{
  return a->getP().tripleProduct(b->getP(), c->getP()).sign();
}

InputPoint Rdir(0.8401877171547095, 0.394382926819093, 0.7830992237586059);

int OnLine::sign ()
{
  PV3 u = a->getP() - t->getP(), v = h->getP() - t->getP();
  return Rdir.getP().tripleProduct(u, v).sign();
}

int Order::sign ()
{
  return (a->getP() - t->getP()).dot(h->getP() - t->getP()).sign();
}

bool onEdge (Point *a, Point *t, Point *h, bool strict)
{
  if (a->identical(t) || a->identical(h))
    return !strict;
  return Order(a, t, h) == 1 && Order(a, h, t) == 1;
}

int Side::sign ()
{
  return (a->getP().dot(p->getN()) + p->getK()).sign();
}

int PointOrder::sign () 
{
  return r->getP().dot(b->getP() - a->getP()).sign();
}

int PointOrderR::sign () 
{
  return Rdir.getP().dot(b->getP() - a->getP()).sign();
}

int TripleProductR::sign ()
{
  return Rdir.getP().tripleProduct(a->getP() - t->getP(),
				   h->getP() - t->getP()).sign();
}

int Orientation::sign ()
{
  
  PV3 pb = b->getP();
  PV3 pc = c->getP();
  PV3 pd = d->getP();
  PV3 u = d->getP() - a->getP(), v = b->getP() - a->getP(), 
    w = c->getP() - a->getP();

  return u.tripleProduct(v, w).sign();
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

Vertex::Vertex (Point *p) : p(p), node(0)
{
  p->incref();
  p->getBBox(bbox);
  if (!inputPerturbed) {
    Parameter k = Rdir.getP().dot(p->getP());
    rint[0] = k.lb();
    rint[1] = k.ub();
  }
}

void Vertex::outgoingHEdges (HEdges &ed) const
{
  for (Edges::const_iterator e = edges.begin(); e != edges.end(); ++e)
    for (HEdges::iterator f = (*e)->hedges.begin(); f != (*e)->hedges.end(); ++f)
      if ((*f)->tail() == this)
	ed.push_back(*f);
}

void Vertex::incidentFaces (FaceSet &fs) const
{
  for (Edges::const_iterator e = edges.begin(); e != edges.end(); ++e)
    for (HEdges::iterator f = (*e)->hedges.begin(); f != (*e)->hedges.end(); ++f)
      if ((*f)->tail() == this)
	fs.insert((*f)->getF());
}

HEdge * Vertex::connected (Vertex *a) const
{
  HEdges ed;
  outgoingHEdges(ed);
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

Edge::Edge (Vertex *t, Vertex *h) : t(t), h(h)
{ 
  setBBox();
}

Edge::~Edge ()
{
  for (HEdges::iterator e = hedges.begin(); e != hedges.end(); ++e)
    if (*e != hei && *e != hei + 1)
      delete *e;
}

void Edge::setBBox ()
{
  copyBBox(t->bbox, bbox);
  mergeBBox(h->bbox, bbox);
}

HEdge * Edge::addHEdge (bool forward)
{
  for (int i = 0; i < 2; ++i)
    if (!hei[i].iflag) {
      hei[i].e = this;
      hei[i].forward = forward;
      hei[i].iflag = true;
      hei[i].f = 0;
      hedges.push_back(hei + i);
      return hei + i;
    }
  HEdge *e = new HEdge(this, forward);
  hedges.push_back(e);
  return e;
}

void Edge::removeHEdge (HEdge *e)
{
  HEdges::iterator j = remove(hedges.begin(), hedges.end(), e);
  hedges.erase(j, hedges.end());
  if (e != hei && e != hei + 1)
    delete e;
}

void Edge::edgeVertices (Vertices &ve) const
{
  ve.push_back(t);
  ve.insert(ve.end(), vertices.begin(), vertices.end());
  ve.push_back(h);
}

Face * Edge::otherFace (Face *f) const
{
  if (hedges.size() != 2 || hedges[0]->f == hedges[1]->f)      
    return 0;
  if (hedges[0]->f == f)
    return hedges[1]->f;
  return hedges[0]->f;
}

void Edge::sortVertices ()
{
  sort(vertices.begin(), vertices.end(), EdgeVertexOrder(this));
}

void Edge::sortHEdges ()
{
  if (hedges.size() > 2)
    sort(hedges.begin(), hedges.end(), EdgeOrder(this));
}

int EdgeVertexOrderP::sign ()
{
  return e->getU().dot(w->getP()->getP() - v->getP()->getP()).sign();
}

int EdgeOrderP::sign ()
{
  PV3 nf = f->getN(), ng = g->getN(), r = Rdir.getP();
  int yf = nf.dot(r).sign(), yg = ng.dot(r).sign();
  if (yf != yg)
    return yf;
  // for degenerate mink
  if (f->getF()->coplanar(g->getF()))
    return f < g ? -1 : 1;
  PV3 u = e->getU();
  return u.tripleProduct(ng, nf).sign();
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
  return forward ? f->p.getN() : - f->p.getN(); 
}

HEdge * HEdge::cw () const
{
  int n = e->hedges.size();
  for (int i = 0; i < n; ++i)
    if (e->hedges[i] == this)
      return e->hedges[(i+1)%n];
}

HEdge * HEdge::ccw () const
{
  int n = e->hedges.size();
  for (int i = 0; i < n; ++i)
    if (e->hedges[i] == this)
      return e->hedges[i == 0 ? n - 1 : i - 1];
}

void HEdge::loop (Vertices &ve)
{
  HEdge *e = this;
  do {
    ve.push_back(e->tail());
    e = e->next;
  }
  while (e != this);
}

HEdge * HEdge::findLoop (Vertices &ve)
{
  HEdge *h = this, *l = h, *n = h->next;
  Face *f = 0;
  if (n)
    for (HEdges::iterator i = e->hedges.begin(); !f && i != e->hedges.end(); ++i)
      if ((*i)->forward == forward && (*i)->f)
	for (HEdges::iterator j = n->e->hedges.begin();
	     !f && j != n->e->hedges.end(); ++j)
	  if ((*i)->f == (*j)->f)
	    f = (*i)->f;
  do {
    ve.push_back(h->tail());
    h->flag = true;
    if (HEdgeOrder()(h, l))
      l = h;
    h = h->next;
    if (!h)
      return 0;
    if (f) {
      bool flag = false;
      for (HEdges::iterator i = h->e->hedges.begin();
	   !flag && i != h->e->hedges.end(); ++i)
	flag = (*i)->f == f;
      if (!flag)
	f = 0;
    }
  }
  while (h != this);
  return f ? 0 : l;
}

void HEdge::loop (HEdges &ed)
{
  HEdge *e = this;
  do {
    ed.push_back(e);
    e = e->next;
  }
  while (e != this);
}

int VertexHHEdgeOrder::sign ()
{
  PV3 eu = v == e->tail() ? e->getU() : - e->getU(),
    fu = v == f->tail() ? f->getU() : - f->getU(),
    r = Rdir.getP();
  int es = eu.dot(r).sign(), fs = fu.dot(r).sign();
  return es == fs ? cross(fu, eu, c).sign() : es;
}

bool HFace::pos () const 
{ 
  return this == f->hfaces; 
}
  
HFace * HFace::twin () const 
{
  return this == f->hfaces ? f->hfaces + 1 : f->hfaces;
}

PV3 HFace::getN () const
{
  return pos() ? f->p.getN() : - f->p.getN();
}

void HFace::neighbors (HFaces &hfaces) const
{
  for (HEdges::iterator b = f->boundary.begin(); b != f->boundary.end(); ++b) {
    HEdge *e = *b;
    do {
      hfaces.push_back(neighbor(e));
      e = e->next;
    }
    while (e != *b);
  }
}

HFace * HFace::neighbor (HEdge *e) const
{
  bool ccw = pos() == e->forward;
  HEdge *f = ccw ? e->ccw() : e->cw();
  bool flag = e->forward == f->forward ? pos() : !pos();
  Face *g = f->f;
  return flag ? g->hfaces + 1 : g->hfaces;
}

void HFace::triangulate (VVMap &vvmap, Polyhedron *a) const
{
  Triangles tr;
  f->triangulate(tr);
  int pc = f->getPC();
  if (pos())
    for (Triangles::iterator t = tr.begin(); t != tr.end(); ++t)
      a->addTriangle(a->getVertex(t->a, vvmap), a->getVertex(t->b, vvmap),
		     a->getVertex(t->c, vvmap), pc);
  else
    for (Triangles::iterator t = tr.begin(); t != tr.end(); ++t)
      a->addTriangle(a->getVertex(t->c, vvmap), a->getVertex(t->b, vvmap),
		     a->getVertex(t->a, vvmap), pc);
}

Face::Face (const Vertices &b, Vertex *u, Vertex *v, Vertex *w, int pc, bool flag)
  : p(u->p, v->p, w->p), pc(pc), flag(flag)
{
  hfaces[0].f = hfaces[1].f = this;
  copyBBox(b[0]->bbox, bbox);
  for (int i = 1; i < b.size(); ++i)
    mergeBBox(b[i]->bbox, bbox);
}

void Face::update ()
{
  Vertices ve;
  boundaryVertices(ve);
  Vertex *u, *v, *w;
  extremalR(ve, u, v, w);
  p = TrianglePlane(u->p, v->p, w->p);
  copyBBox(ve[0]->bbox, bbox);
  for (int i = 1; i < ve.size(); ++i)
    mergeBBox(ve[i]->bbox, bbox);
  pc = 0;
}

void Face::addBoundary (HEdge *e)
{
  boundary.push_back(e);
  while (!e->f) {
    e->f = this;
    e->tail()->p->ps.insert(p.id);
    e = e->next;
  }
}

int Face::getPC () {
  if (pc == 0)
    pc = projectionCoordinate(getP());
  return pc;
}

bool Face::boundaryVertex (Point *a) const
{
  Vertices ve;
  boundaryVertices(ve);
  for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
    if ((*v)->p->identical(a))
      return true;
  return false;
}

bool Face::boundaryVertex (Vertex *v) const
{
  Vertices ve;
  boundaryVertices(ve);
  return find(ve.begin(), ve.end(), v) != ve.end();
}

void Face::boundaryVertices (Vertices &ve) const
{
 for (HEdges::const_iterator b = boundary.begin(); b != boundary.end(); ++b)
    (*b)->loop(ve);
}

bool Face::boundaryEdge (Edge *e) const
{
  for (HEdges::iterator f = e->hedges.begin(); f != e->hedges.end(); ++f)
    if ((*f)->f == this)
      return true;
  return false;
}

void Face::boundaryHEdges (HEdges &ed) const
{
  for (HEdges::const_iterator b = boundary.begin(); b != boundary.end(); ++b)
    (*b)->loop(ed);
}

void Face::sharedBoundaryVertices (Face *f, Vertices &vfg) const
{
  sharedBoundaryVertices1(f, vfg);
  f->sharedBoundaryVertices1(this, vfg);
}

void Face::sharedBoundaryVertices1 (const Face *f, Vertices &vfg) const
{
  VertexSet vs;
  edgeVertices(vs);
  Vertices vf;
  f->boundaryVertices(vf);
  for (Vertices::iterator v = vf.begin(); v != vf.end(); ++v)
    if (vs.find(*v) != vs.end() &&
	find(vfg.begin(), vfg.end(), *v) == vfg.end())
      vfg.push_back(*v);
}

void Face::containedBoundaryVertices (Face *f, Vertices &vfg)
{
  containedBoundaryVertices1(f, vfg);
  f->containedBoundaryVertices1(this, vfg);
}

void Face::containedBoundaryVertices1 (const Face *f, Vertices &vfg)
{
  Vertices vf;
  f->boundaryVertices(vf);
  for (Vertices::iterator v = vf.begin(); v != vf.end(); ++v)
    if ((*v)->p->side(&p) == 0 && contains((*v)->p) &&
	find(vfg.begin(), vfg.end(), *v) == vfg.end())
      vfg.push_back(*v);
}

bool Face::sharedBoundaryEdge (Face *f) const
{
  HEdges ed;
  boundaryHEdges(ed);
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    for (HEdges::iterator h = (*e)->e->hedges.begin(); 
	 h != (*e)->e->hedges.end(); ++h)
      if ((*h)->f == f)
	return true;
  return false;
}

void Face::edgeVertices (VertexSet &vs) const
{
  HEdges ed;
  boundaryHEdges(ed);
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e) {
    Vertices ve;
    (*e)->e->edgeVertices(ve);
    vs.insert(ve.begin(), ve.end());
  }
}

bool Face::coplanar (Face *f)
{
  if (p.id == f->p.id)
    return true;
  return !inputPerturbed && SamePlane(&p, &f->p) == 1;
};

bool Face::intersectRay (Point *a, Point *r)
{
  if (a->side(&p) == 0)
    return false;
  RayPlanePoint q(a, r, &p);
  return contains(&q) && PointOrderR(a, &q) == 1;
}

bool Face::contains (Point *a)
{
  if (!bboxOverlap(a, bbox))
    return false;
  if (triangle())
    return containsConvex(a);
  if (!boundaryContains(a, 0))
    return false;
  for (int i = 1; i < boundary.size(); ++i)
    if (boundaryContains(a, i))
      return false;
  return true;
}

bool Face::triangle () const
{
  return boundary.size() == 1 &&
    boundary[0] == boundary[0]->next->next->next;
}

bool Face::containsConvex (Point *a, bool strict)
{
  Vertices ve;
  boundary[0]->loop(ve);
  Point *t = (*ve.rbegin())->getP();
  int s = strict ? 1 : 0;
  for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v) {
    Point *h = (*v)->getP();
    if (a == t || a == h || LeftTurn(a, t, h, getPC()) < s)
      return false;
    t = h;
  }
  return true;
}

bool Face::boundaryContains (Point *a, int i)
{
  Vertices ve;
  boundary[i]->loop(ve);
  bool res = false;
  Point *vt = (*ve.rbegin())->p;
  for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v) {
    Point *vh = (*v)->p;
    if (a->onLine(vt, vh)) {
      if (onEdge(a, vt, vh, false))
	return i > 0;
    }
    else if (RayEdgeIntersection(a, &Rdir, vt, vh, getPC()) == 1)
      res = !res;
    vt = vh;
  }
  return res;
}

int RayEdgeIntersection::sign ()
{
  PV3 ap = a->getP(), u = r->getP(), tp = t->getP(), hp = h->getP(), v = hp - tp;
  if (cross(tp - ap, u, c).sign() == cross(hp - ap, u, c).sign())
    return -1;
  return cross(ap - tp, v, c).sign() == cross(u, v, c).sign() ? -1 : 1;
}

void Face::triangulate (Triangles &tr)
{
  if (boundary.size() == 1) {
    Vertices ve;
    boundary[0]->loop(ve);
    if (ve.size() == 3) {
      tr.push_back(Triangle(ve[0], ve[1], ve[2]));
      return;
    }
  }
  VVertices reg;
  for (HEdges::const_iterator b = boundary.begin(); b != boundary.end(); ++b) {
    HEdge *e = *b;
    Vertices *rb = new Vertices;
    reg.push_back(rb);
    e->loop(*rb);
  }
  ::triangulate(reg, getPC(), tr);
  deleteRegion(reg);
}

Face * newFace (const Vertices &ve, int pc, bool flag)
{
  Vertex *u, *v, *w;
  extremalR(ve, u, v, w);
  return new Face(ve, u, v, w, pc, flag);
}

void extremalR (const Vertices &ve, Vertex *&u, Vertex *&v, Vertex *&w)
{
  int im = 0, n = ve.size();
  if (n == 3) {
    u = ve[0];
    v = ve[1];
    w = ve[2];
    return;
  }
  for (int i = 1; i < n; ++i)
    if (ve[im] != ve[i] && PointOrderR(ve[im]->getP(), ve[i]->getP()) == 1)
      im = i;
  u = ve[im == 0 ? n - 1 : im - 1];
  v = ve[im];
  w = ve[(im+1)%n];
}

void deleteRegion (const VVertices &reg)
{
  for (VVertices::const_iterator r = reg.begin(); r != reg.end(); ++r)
    delete(*r);
}

Shell::~Shell ()
{
  if (octreef)
    delete octreef;
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
  InputPoint r(0.0, 0.0, 1.0);
  Vertex *vm = vmax(&r);
  HEdge *em = 0;
  HFace *fm = 0;
  for (Edges::iterator e = vm->edges.begin(); e != vm->edges.end(); ++e)
    for (HEdges::iterator h = (*e)->hedges.begin(); h != (*e)->hedges.end(); ++h)
      for (int i = 0; i < 2; ++i)
	if ((*h)->f->hfaces[i].s == this) {
	  if (!em || em->e != *e && SlopeOrder(*e, em->e, &r) == 1) {
	    em = *h;
	    fm = (*h)->f->hfaces + i;
	  }
	  break;
	}
  return Convex(em, fm) == 1;
}

Vertex * Shell::vmax (Point *r) const
{
  Vertex *v = hfaces[0]->f->boundary[0]->tail();
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
	if ((*f)->intersectRay(w->p, r))
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
  Vertices ve;
  f->boundaryVertices(ve);
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
  HEdge *e = s->hfaces[0]->f->boundary[0];
  int i = contains(e->tail()->p);
  if (i == 0)
    i = contains(e->head()->p);
  return i == 1;
}

int Shell::contains (Point *a) const
{
  InputPoint r(0.0, 0.0, 1.0);
  double rb[6];
  rayBBox(a, &r, rb);
  Faces fa;
  octreef->find(rb, fa);
  bool res = false;
  for (Faces::iterator f = fa.begin(); f != fa.end(); ++f)
    if ((*f)->boundaryVertex(a))
      return 0;
  else if ((*f)->intersectRay(a, &r))
    res = !res;
  return res ? 1 : -1;
}

void Shell::rayBBox (Point *a, Point *r, double *rb) const
{
  a->getBBox(rb);
  PV3 ap = a->getP(), rp = r->getP();
  Parameter k = (bbox[5] - ap.z)/rp.z;
  InputPoint q(ap + k*rp);
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
  for (set<Face *>::iterator f = fs.begin(); f != fs.end(); ++f)
    for (HEdges::iterator b = (*f)->boundary.begin(); 
	 b != (*f)->boundary.end(); ++b) {
      HEdge *e = *b, *e0 = e;
      do {
	vs.insert(e->tail());
	es.insert(e->e);
	e = e->next;
      }
      while (e != e0);
    }
  int nv = vs.size(), ne = es.size(), nf = fs.size();
  return nv - ne + nf;
}

void deleteShells (const Shells &sh)
{
  for (Shells::const_iterator s = sh.begin(); s != sh.end(); ++s)
    delete *s;
}

int SlopeOrder::sign ()
{
  PV3 u = e->getU(), v = f->getU(), r = x->getP();
  Parameter ur = u.dot(r), vr = v.dot(r);
  return (vr*vr*u.dot(u) - ur*ur*v.dot(v)).sign();
}

int Convex::sign ()
{
  HFace *g = f->neighbor(e);
  Face *ff = f->getF(), *gf = g->getF();
  if (ff->getP() == gf->getP())
    return -1;
  PV3 nf = ff->getP()->getN(), ng = g->getN(), u = e->getU();
  return u.tripleProduct(ng, nf).sign();
}

void Cell::addBoundary (Shell *s)
{
  boundary.push_back(s);
  s->c = this;
}

bool Cell::contains (Point *p) const
{
  if (boundary[0]->contains(p) < 1)
    return false;
  for (int i = 1; i < boundary.size(); ++i)
    if (boundary[i]->contains(p) > -1)
      return false;
  return true;
}

FFPair ffpair (Face *f, Face *g)
{
  return FFPair(f < g ? f : g, f < g ? g : f);
}

int FFOrderP::sign ()
{
  PV3 u = f->getP()->getN().cross(g->getP()->getN()), 
    vw = v->getP()->getP() - w->getP()->getP();
  return u.dot(vw).sign();
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

Vertex * Polyhedron::getVertex (Point *p)
{
  Vertex *v = new Vertex(p);
  if (!inputPerturbed) {
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

Vertex * Polyhedron::getVertex (Vertex *v, VVMap &vmap)
{
  VVMap::iterator iter = vmap.find(v);
  if (iter != vmap.end())
    return iter->second;
  Vertex *w = getVertex(v->p);
  vmap.insert(VVPair(v, w));
  return w;
}

HEdge * Polyhedron::addHEdge (Vertex *a, Vertex *b)
{
  Edge *e = getEdge(a, b);
  bool forward = e->t == a;
  return e->addHEdge(forward);
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

HEdge * Polyhedron::getHEdge (Vertex *a, Vertex *b)
{
  Edge *e = getEdge(a, b);
  bool forward = e->t == a;
  for (HEdges::iterator f = e->hedges.begin(); f != e->hedges.end(); ++f)
    if ((*f)->forward == forward && !(*f)->f)
      return *f;
  return e->addHEdge(forward);
}

HEdge * Polyhedron::findHEdge (Vertex *a, Vertex *b)
{
  for (Edges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    if ((*e)->t == a && (*e)->h == b || (*e)->t == b && (*e)->h == a) {
      for (HEdges::iterator h = (*e)->hedges.begin(); h != (*e)->hedges.end(); ++h)
	if ((*h)->tail() == a && (*h)->head() == b)
	  return *h;
      return 0;
    }
  return 0;
}

void Polyhedron::addVertex (Edge *e, Vertex *v)
{
  if (find(e->vertices.begin(), e->vertices.end(), v) == e->vertices.end()) {
    e->vertices.push_back(v);
    iedges.insert(e);
  }
}

Face * Polyhedron::addFace (const VVertices &reg, int pc, bool flag)
{
  Face *f = newFace(*reg[0], pc, flag);
  faces.push_back(f);
  for (VVertices::const_iterator r = reg.begin(); r != reg.end(); ++r)
    f->addBoundary(addLoop(**r));
  return f;
}

HEdge * Polyhedron::addLoop (const Vertices &ve)
{
 int n = ve.size();
 HEdge *ep = addHEdge(ve[n-1], ve[0]), *e0 = ep;
 for (int i = 0; i + 1 < n; ++i) {
   HEdge *e = addHEdge(ve[i], ve[i+1]);
   ep->next = e;
   ep = e;
  }
 ep->next = e0;
 return e0->next;
}

HEdge * Polyhedron::addLoop (Vertex **v, int n)
{
  Vertices ve(v, v + n);
  return addLoop(ve);
}

Face * Polyhedron::addFace (const Vertices &ve, int pc, bool flag)
{
  Face *f = newFace(ve, pc, flag);
  faces.push_back(f);
  f->addBoundary(addLoop(ve));
  return f;
}

Face * Polyhedron::addFace (HEdge *e, int pc, bool flag)
{
  Vertices ve;
  e->loop(ve);
  Face *f = newFace(ve, pc, flag);
  faces.push_back(f);
  f->addBoundary(e);
  return f;
}

Face * Polyhedron::addTriangle (Vertex *a, Vertex *b, Vertex *c, int pc, bool flag)
{
  Vertices ve;
  ve.push_back(a);
  ve.push_back(b);
  ve.push_back(c);
  return addFace(ve, pc, flag);
}

Face * Polyhedron::addRectangle (Vertex *a, Vertex *b, Vertex *c, Vertex *d)
{
  Vertices ve;
  ve.push_back(a);
  ve.push_back(b);
  ve.push_back(c);
  ve.push_back(d);
  return addFace(ve);
}

void Polyhedron::formCells () {
  if (!cells.empty()) {
    for (Faces::iterator f = faces.begin(); f != faces.end(); ++f)
      for (int i = 0; i < 2; ++i)
	(*f)->hfaces[i].s = 0;
    for (Cells::iterator c = cells.begin(); c != cells.end(); ++c)
      delete *c;
    cells.clear();
  }
  Shells shells;
  formShells(shells);
  formCellsAux(shells);
}

void Polyhedron::formCellsAux (const Shells &shells)
{
  Shells inner;
  cells.push_back(new Cell);
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
    enclosingCell(*s, octreec)->addBoundary(*s);
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
      g->neighbors(st);
    }
  }
  return s;
}

Cell * Polyhedron::enclosingCell (Shell *s, Octree<Cell *> *octreec) const
{
  Cells ce;
  octreec->find(s->hfaces[0]->f->boundary[0]->tail()->getBBox(), ce);
  Cell *c = 0;
  for (Cells::iterator d = ce.begin(); d != ce.end(); ++d)
    if ((*d)->boundary[0]->contains(s) &&
	(!c || c->boundary[0]->contains((*d)->boundary[0])))
      c = *d;
  return c ? c : cells[0];
}

Polyhedron * Polyhedron::triangulate () const
{
  Polyhedron *a = new Polyhedron;
  VVMap vvmap;
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f) {
    Triangles tr;
    (*f)->triangulate(tr);
    for (Triangles::iterator t = tr.begin(); t != tr.end(); ++t)
      a->addTriangleUnique(a->getVertex(t->a, vvmap), a->getVertex(t->b, vvmap),
			   a->getVertex(t->c, vvmap), (*f)->getPC(), (*f)->flag);
  }
  return a;
}

Polyhedron * Polyhedron::subdivide ()
{
  intersectFF();
  return subdivideAux();
}

void Polyhedron::intersectFF ()
{
  Octree<Face *> *octree = faceOctree();
  Faces fa, fb;
  octree->pairs(fa, fb);
  EFVMap efvmap;
  for (int i = 0; i < fa.size(); ++i)
    if (!fa[i]->boundary.empty() && !fb[i]->boundary.empty())
      intersectFF(fa[i], fb[i], efvmap);
  delete octree;
}

void Polyhedron::intersectFF (Face *f, Face *g, EFVMap &efvmap)
{
  if (f->coplanar(g)) {
    intersectFFP(f, g, efvmap);
    intersectFFP(g, f, efvmap);
  }
  else if (!f->sharedBoundaryEdge(g)) {
    Edges iedges;
    Faces ifaces;
    if (intersectPE(f, g, efvmap, iedges, ifaces) &&
	intersectPE(g, f, efvmap, iedges, ifaces)) {
      Vertices vfg;
      for (int i = 0; i < iedges.size(); ++i)
	intersectEF(iedges[i], ifaces[i], efvmap, vfg);
      formFF(f, g, vfg);
    }
  }
}

void Polyhedron::intersectFFP (Face *f, Face *g, EFVMap &efvmap)
{
  HEdges ed;
  g->boundaryHEdges(ed);
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e) {
    intersectFEP(f, (*e)->e, efvmap);
  }
}

void Polyhedron::intersectFEP (Face *f, Edge *e, EFVMap &efvmap)
{
  if (f->boundaryEdge(e) || !bboxOverlap(f->bbox, e->bbox) ||
      !efvmap.insert(EFVPair(EFPair(e, f), 0)).second)
    return;
  HEdges ed;
  f->boundaryHEdges(ed);
  for (HEdges::iterator h = ed.begin(); h != ed.end(); ++h)
    if (e->t->p->onLine((*h)->tail()->p, (*h)->head()->p) &&
	e->h->p->onLine((*h)->tail()->p, (*h)->head()->p)) {
      e->dedges.push_back((*h)->e);
      (*h)->e->dedges.push_back(e);
      return;
    }
  VertexSet vs;
  for (HEdges::iterator h = ed.begin(); h != ed.end(); ++h) {
    Vertex *v = intersectEE(e, (*h)->e, f->getPC());
    if (v)
      vs.insert(v);
  }
  if (f->contains(e->t->p))
    vs.insert(e->t);
  if (f->contains(e->h->p))
    vs.insert(e->h);
  Vertices ve;
  e->edgeVertices(ve);
  for (HEdges::iterator h = ed.begin(); vs.size() < 2 && h != ed.end(); ++h)
    if (find(ve.begin(), ve.end(), (*h)->tail()) != ve.end())
      vs.insert((*h)->tail());
    else {
      Vertices vh;
      (*h)->e->edgeVertices(vh);
      if (find(vh.begin(), vh.end(), e->t) != vh.end())
	vs.insert(e->t);
      else if (find(vh.begin(), vh.end(), e->h) != vh.end())
	vs.insert(e->h);
    }
  if (vs.size() < 2)
    return;
  Edge *fg = getEdge(*vs.begin(), *vs.rbegin());
  fg->addHEdge(true)->f = f;
  fg->addHEdge(false)->f = f;
  f->edges.push_back(fg);
  if (fg != e) {
    fg->dedges.push_back(e);
    e->dedges.push_back(fg);
  }  
}

Vertex * Polyhedron::intersectEE (Edge *e, Edge *f, int pc)
{
  if (!bboxOverlap(e->bbox, f->bbox))
    return 0;
  EEPair ef(e < f ? e : f, e < f ? f : e);
  EEVMap::iterator i = eevmap.find(ef);
  if (i != eevmap.end())
    return i->second;
  bool f1 = intersectEV(e, f->t), f2 = intersectEV(e, f->h),
    f3 = intersectEV(f, e->t), f4 = intersectEV(f, e->h);
  if (f1 + f2 + f3 + f4 > 1) {
    intersectEEL(e, f);
    eevmap.insert(EEVPair(ef, 0));
    return 0;
  }
  if (f1 || f2 || f3 || f4) {
    eevmap.insert(EEVPair(ef, 0));
    return 0;
  }
  int tp1 = TripleProductR(e->t->p, f->t->p, f->h->p),
    tp2 = TripleProductR(e->h->p, f->t->p, f->h->p);
  if (tp1 == 0 || tp2 == 0 || tp1 == tp2) {
    eevmap.insert(EEVPair(ef, 0));
    return 0;
  }
  int tp3 = TripleProductR(f->t->p, e->t->p, e->h->p),
    tp4 = TripleProductR(f->h->p, e->t->p, e->h->p);
  if (tp3 == 0 || tp4 == 0 || tp3 == tp4) {
    eevmap.insert(EEVPair(ef, 0));
    return 0;
  }
  Vertex *v = getVertex(new EEPoint(e, f, pc));
  eevmap.insert(EEVPair(ef, v));
  addVertex(e, v);
  addVertex(f, v);
  return v;
}

bool Polyhedron::intersectEV (Edge *e, Vertex *v)
{
  if (v == e->t || v == e->h)
    return true;
  if (v->p->onLine(e->t->p, e->h->p)) {
    if (bboxOverlap(v->getBBox(), e->getBBox()) &&
	onEdge(v->p, e->t->p, e->h->p, true))
      addVertex(e, v);
    return true;
  }
  return false;
}

void Polyhedron::intersectEEL (Edge *e, Edge *f)
{
  if (find(e->vertices.begin(), e->vertices.end(), f->t) == e->vertices.end() &&
      find(e->vertices.begin(), e->vertices.end(), f->h) == e->vertices.end() &&
      find(f->vertices.begin(), f->vertices.end(), e->t) == f->vertices.end() &&
      find(f->vertices.begin(), f->vertices.end(), e->h) == f->vertices.end())
    return;
  e->dedges.push_back(f);
  f->dedges.push_back(e);
}

bool Polyhedron::intersectPE (Face *f, Face *g, EFVMap &efvmap, Edges &iedges,
			      Faces &ifaces)
{
  vector<int> s;
  HEdges ed;
  g->boundaryHEdges(ed);
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e) {
    Vertex *t = (*e)->tail();
    s.push_back(t->p->side(&f->p));
  }
  int n = ed.size();
  bool sp = false, sm = false;
  for (int i = 0; i < n; ++i) {
    int si = s[i], sj = s[(i+1)%n];
    if (si == 0 && sj == 0)
      intersectFEP(f, ed[i]->e, efvmap);
    else if (si == 0)
      intersectFV(f, ed[i]->tail());
    else {
      if (si == 1)
	sp = true;
      else
	sm = true;
      if (si*sj == -1) {
	ifaces.push_back(f);
	iedges.push_back(ed[i]->e);
      }
    }
  }
  return sp && sm;
}

void Polyhedron::intersectFV (Face *f, Vertex *v)
{
  if (f->boundaryVertex(v))
    return;
  HEdges ed;
  f->boundaryHEdges(ed);
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    if (bboxOverlap(v->p, (*e)->e->getBBox()) && 
	intersectEV((*e)->e, v))
      return;
}

void Polyhedron::intersectEF (Edge *e, Face *f, EFVMap &efvmap, Vertices &vfg)
{
  if (!bboxOverlap(e->bbox, f->bbox))
    return;
  EFPair ef(e, f);
  EFVMap::iterator i = efvmap.find(ef);
  if (i != efvmap.end()) {
    if (i->second && find(vfg.begin(), vfg.end(), i->second) == vfg.end())
      vfg.push_back(i->second);
    return;
  }
  Point *p = new EPPoint(e->t->p, e->h->p, &f->p);
  if (bboxOverlap(p, f->bbox)) {
    HEdges ed;
    f->boundaryHEdges(ed);
    int ie = -1, n = ed.size();
    bool flag = true;
    for (int i = 0; ie == -1 && flag && i < n; ++i) {
      Point *t = ed[i]->tail()->p, *h = ed[i]->head()->p;
      int s = LeftTurn(p, t, h, f->getPC());
      if (s == -1)
	flag = false;
      else if (s == 0)
	if (onEdge(p, t, h, true))
	  ie = i;
	else
	  flag = false;
    }
    if (flag) {
      Vertex *v = getVertex(p);
      addVertex(e, v);
      if (ie != -1)
	addVertex(ed[ie]->e, v);
      if (find(vfg.begin(), vfg.end(), v) == vfg.end())
	vfg.push_back(v);
      efvmap.insert(pair<EFPair, Vertex *>(ef, v));
      return;
    }
  }
  p->decref();
  efvmap.insert(pair<EFPair, Vertex *>(ef, 0));
}

void Polyhedron::formFF (Face *f, Face *g, Vertices &vfg)
{
  if (vfg.size() < 2)
    f->sharedBoundaryVertices(g, vfg);
  if (vfg.size() < 2)
    f->containedBoundaryVertices(g, vfg);
  if (vfg.size() < 2)
    return;
  sort(vfg.begin(), vfg.end(), FFOrder(f, g));
  Edge *e = getEdge(vfg[0], vfg[1]);
  f->edges.push_back(e);
  g->edges.push_back(e);
  e->t->p->ps.insert(f->p.id);
  e->t->p->ps.insert(g->p.id);
  e->h->p->ps.insert(f->p.id);
  e->h->p->ps.insert(g->p.id);
  e->addHEdge(true)->f = f;
  e->addHEdge(false)->f = g;
}

Polyhedron * Polyhedron::subdivideAux ()
{
  FFFVMap fffvmap;
  for (Faces::iterator f = faces.begin(); f != faces.end(); ++f)
    intersectFFF(*f, fffvmap);
  for (Edges::const_iterator e = edges.begin(); e != edges.end(); ++e)
    (*e)->sortVertices();
  Polyhedron *a = new Polyhedron;
  VVMap vvmap;
  for (Faces::iterator f = faces.begin(); f != faces.end(); ++f)
    subfaces(*f, a, vvmap);
  return a;
}

void Polyhedron::intersectFFF (Face *f, FFFVMap &fffvmap)
{
  int n = f->edges.size(), pc = f->getPC();
  for (int i = 0; i + 1 < n; ++i) {
    Edge *ei = f->edges[i];
    Face *fi = ei->otherFace(f);
    for (int j = i + 1; j < n; ++j) {
      Edge *ej = f->edges[j];
      Face *fj = ej->otherFace(f);
      if (fi && fj && fi != fj) {
	FFF fff(f, fi, fj);
	FFFVMap::iterator iter = fffvmap.find(fff);
	if (iter == fffvmap.end()) {
	  Vertex *v = intersectEE(ei, ej, pc);
	  if (v)
	    fffvmap.insert(FFFVPair(fff, v));
	}
	else {
	  Vertex *v = iter->second;
	  if (find(ei->vertices.begin(), ei->vertices.end(), v) == ei->vertices.end())
	    addVertex(ei, v);
	  else
	    addVertex(ej, v);
	}
      }
      else
	intersectEE(ei, ej, pc);
    }
  }
}

void Polyhedron::subfaces (Face *f, Polyhedron *a, VVMap &vvmap)
{
  HEdges he, outer, inner, bad;
  subedges(f, a, vvmap, he);
  for (HEdges::iterator e = he.begin(); e != he.end(); ++e)
    if (!(*e)->flag) {
      Vertices ve;
      HEdge *l = (*e)->findLoop(ve);
      if (!l)
	bad.push_back(*e);
      else if (outerLoop(ve, f->getPC()))
	outer.push_back(l);
      else 
	inner.push_back(l);
    }
  sort(outer.begin(), outer.end(), HEdgeOrder());
  if (oneWayInt)
    sort(inner.begin(), inner.end(), HEdgeOrder());
  Faces nfaces;
  for (HEdges::iterator e = outer.begin(); e != outer.end(); ++e)
    nfaces.push_back(a->addFace(*e, f->getPC(), f->flag));
  for (HEdges::iterator e = inner.begin(); e != inner.end(); ++e) {
    Face *g = enclosingFace(nfaces, (*e)->tail()->p);
    if (g)
      g->addBoundary(*e);
    else
      bad.push_back(*e);
  }
  a->removeLoops(bad);
}

void Polyhedron::subedges (Face *f, Polyhedron *a, VVMap &vvmap, HEdges &he)
{
  HEdges ed;
  f->boundaryHEdges(ed);
  VVPairSet vvps;
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    subedges(*e, vvps);
  for (Edges::iterator e = f->edges.begin(); e != f->edges.end(); ++e)
    for (HEdges::iterator h = (*e)->hedges.begin(); h != (*e)->hedges.end(); ++h)
      if (!oneWayInt || (*h)->f == f)
	subedges(*h, vvps);
  for (VVPairSet::iterator i = vvps.begin(); i != vvps.end(); ++i)
    he.push_back(a->addHEdge(a->getVertex(i->first, vvmap), 
			     a->getVertex(i->second, vvmap)));
  setNext(he, f->getPC());
}

void Polyhedron::subedges (HEdge *e, VVPairSet &vvps)
{
  Vertices ve;
  e->e->edgeVertices(ve);
  Vertex *t = ve[0];
  for (int i = 1; i < ve.size(); ++i) {
    Vertex *h = ve[i];
    if (e->forward)
      subedge(e->e, t, h, vvps);
    else
      subedge(e->e, h, t, vvps);
    t = h;
  }
}

void Polyhedron::subedge (Edge *e, Vertex *t, Vertex *h, VVPairSet &vvps)
{
  vvps.insert(VVPair(t, h));
  PPPair pp(t->p < h->p ? t->p : h->p, t->p < h->p ? h->p : t->p);
  PPEMap::iterator i = ppemap.find(pp);
  if (i == ppemap.end()) {
    EdgeSet es;
    es.insert(e);
    es.insert(e->dedges.begin(), e->dedges.end());
    ppemap.insert(PPEPair(pp, es));
  }
  else {
    i->second.insert(e);
    i->second.insert(e->dedges.begin(), e->dedges.end());
  }
}

void Polyhedron::setNext (const HEdges &he, int c) const
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
	hhe[i+j].e->next = hhe[k].e;
      ++j;
    }
    i += m;
  }
}

Face * Polyhedron::enclosingFace (const Faces &fa, Point *p) const
{
  for (Faces::const_iterator f = fa.begin(); f != fa.end(); ++f)
    if ((*f)->boundaryContains(p, 0)) {
      if (oneWayInt)
	for (int i = 1; i < (*f)->boundary.size(); ++i)
	  if ((*f)->boundaryContains(p, i))
	    return 0;
      return *f;
    }
  return 0;
}

void Polyhedron::removeLoops (const HEdges &ed)
{
  HEdgeSet es;
  for (HEdges::const_iterator e = ed.begin(); e != ed.end(); ++e) {
    HEdge *f = *e;
    while (f) {
      es.insert(f);
      f = f->next;
      if (f == *e)
	break;
    }
  }
  for (HEdgeSet::iterator e = es.begin(); e != es.end(); ++e)
    (*e)->e->removeHEdge(*e);
}

// used in pack and MinkowskiSumFull from here on

Polyhedron * Polyhedron::negative () const
{
  Polyhedron *a = new Polyhedron;
  VVMap vvmap;
  for (Vertices::const_iterator v = vertices.begin(); v != vertices.end(); ++v)
    vvmap.insert(VVPair(*v, a->getVertex(new NegPoint((*v)->p))));
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f) {
    HEdge *e = (*f)->boundary[0];
    Vertex *u = vvmap.find(e->tail())->second,
      *v = vvmap.find(e->head())->second,
      *w = vvmap.find(e->next->head())->second;
    a->addTriangle(w, v, u);
  }
  return a;
}

Polyhedron * Polyhedron::translate (Point *t) const
{
  Polyhedron *a = new Polyhedron;
  VVMap vvmap;
  for (Vertices::const_iterator v = vertices.begin(); v != vertices.end(); ++v)
    vvmap.insert(VVPair(*v, a->getVertex(new SumPoint(t, (*v)->p))));
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f) {
    HEdge *e = (*f)->boundary[0];
    Vertex *u = vvmap.find(e->tail())->second,
      *v = vvmap.find(e->head())->second,
      *w = vvmap.find(e->next->head())->second;
    a->addTriangle(u, v, w);
  }
  return a;
}

Polyhedron * Polyhedron::negativeTranslate (Point *t) const
{
  Polyhedron *a = new Polyhedron;
  VVMap vvmap;
  for (Vertices::const_iterator v = vertices.begin(); v != vertices.end(); ++v)
    vvmap.insert(VVPair(*v, a->getVertex(new DiffPoint(t, (*v)->p))));
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f) {
    HEdge *e = (*f)->boundary[0];
    Vertex *u = vvmap.find(e->tail())->second,
      *v = vvmap.find(e->head())->second,
      *w = vvmap.find(e->next->head())->second;
    a->addTriangle(w, v, u);
  }
  return a;
}

bool Polyhedron::intersects (Polyhedron *a) const
{
  return contains(a->vertices[0]->p) || a->contains(vertices[0]->p) ||
    intersectsEdges(a) || a->intersectsEdges(this);
}

bool Polyhedron::contains (Point *p) const
{
  for (int i = 1; i < cells.size(); ++i)
    if (cells[i]->contains(p))
      return true;
  return false;
}

bool Polyhedron::intersectsEdges (const Polyhedron *a) const
{
  Octree<Face *> *octree = a->faceOctree();
  for (Edges::const_iterator e = edges.begin(); e != edges.end(); ++e) {
    Faces fa;
    octree->find((*e)->bbox, fa);
    for (Faces::iterator f = fa.begin(); f != fa.end(); ++f)
      if (intersectsEF(*e, *f)) {
	delete octree;
	return true;
      }
  }
  delete octree;
  return false;
}

bool Polyhedron::intersectsEF (Edge *e, Face *f) const
{
  if (e->t->p->side(&f->p)*e->h->p->side(&f->p) != -1)
    return false;
  EPPoint p(e->t->p, e->h->p, &f->p);
  if (!bboxOverlap(&p, f->bbox))
    return false;
  HEdges ed;
  f->boundaryHEdges(ed);
  for (HEdges::iterator g = ed.begin(); g != ed.end(); ++g) {
    Point *t = (*g)->tail()->p, *h = (*g)->head()->p;
    if (LeftTurn(&p, t, h, f->getPC()) < 1)
      return false;
  }
  return true;
}

Polyhedron * Polyhedron::overlay (Polyhedron *a) const
{
  Polyhedron *c = new Polyhedron;
  VVMap vvmap;
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f)
    c->copyFace(*f, vvmap);
  for (Faces::iterator f = c->faces.begin(); f != c->faces.end(); ++f)
    (*f)->flag = true;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f)
    c->copyFace(*f, vvmap);
  Polyhedron *d = c->subdivide(), *e = d->triangulate();
  delete c;
  delete d;
  return e;
}

Face * Polyhedron::copyFace (Face *f, VVMap &vvmap)
{
  VVertices reg;
  for (HEdges::const_iterator e = f->boundary.begin(); e != f->boundary.end(); ++e) {
    Vertices ve, *nve = new Vertices;
    (*e)->loop(ve);
    for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
      nve->push_back(getVertex(*v, vvmap));
    reg.push_back(nve);
  }
  Face *g = addFace(reg);
  deleteRegion(reg);
  return g;
}

Polyhedron * Polyhedron::boolean (Polyhedron *a, SetOp op) const
{
  Polyhedron *c = overlay(a), *d = new Polyhedron;
  VVMap vvmap;
  c->formCells();
  CellSet cin = c->boolean(op);
  for (Faces::iterator f = c->faces.begin(); f != c->faces.end(); ++f) {
    bool in1 = cin.find((*f)->hfaces[0].s->c) != cin.end(),
      in2 = cin.find((*f)->hfaces[1].s->c) != cin.end();
    if (!in1 && in2)
      (*f)->hfaces[0].triangulate(vvmap, d);
    else if (in1 && !in2)
      (*f)->hfaces[1].triangulate(vvmap, d);
  }
  delete c;
  return d;
}

CellSet Polyhedron::boolean (SetOp op) const
{
  CellSet done, res;
  boolean(op, cells[0], false, false, done, res);
  return res;
}

void Polyhedron::boolean (SetOp op, Cell *c, bool ina, bool inb, CellSet &done,
			  CellSet &res) const
{
  if (!done.insert(c).second)
    return;
  if (inSet(ina, inb, op))
    res.insert(c);
  for (Shells::iterator s = c->boundary.begin(); s != c->boundary.end(); ++s)
    for (HFaces::iterator f = (*s)->hfaces.begin(); f != (*s)->hfaces.end(); ++f) {
      Cell *d = (*f)->twin()->s->c;
      bool inad = (*f)->f->flag ? !ina : ina, inbd = (*f)->f->flag ? inb : !inb;
      boolean(op, d, inad, inbd, done, res);
    }
}

// used in simplify from here on

void Polyhedron::replaceVertex (Face *f, Vertex *v, Vertex *w)
{
  HEdge *e = f->boundary[0];
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
  f->boundary[0] = e;
  f->p = TrianglePlane(t->p, h->p, w->p);
  f->pc = 0;
  copyBBox(e->e->bbox, f->bbox);
  mergeBBox(en->e->bbox, f->bbox);
  mergeBBox(enn->e->bbox, f->bbox);
}

void Polyhedron::removeLoop (HEdge *e)
{
  if (e->f)
    e->f->boundary.clear();
  HEdges ed;
  e->loop(ed);
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
  v->p->decref();
  p->incref();
  v->p = p;
  if (!inputPerturbed) {
    vtree.remove(v->node);
    v->node = vtree.insert(v);
  }
  p->getBBox(v->bbox);
  mergeBBox(v->bbox, bbox);
  HEdges ed;
  v->outgoingHEdges(ed);
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
      vertices[i] = *vertices.rbegin();
      vertices.pop_back();
    }
    else
      ++i;
  i = 0;
  while (i < edges.size())
    if (edges[i]->hedges.empty()) {
      delete edges[i];
      edges[i] = *edges.rbegin();
      edges.pop_back();
    }
    else
      ++i;
  i = 0;
  while (i < faces.size())
    if (faces[i]->boundary.empty()) {
      delete faces[i];
      faces[i] = *faces.rbegin();
      faces.pop_back();
    }
    else
      ++i;
}

void Polyhedron::addTriangleUnique (Vertex *a, Vertex *b, Vertex *c, int pc, bool flag)
{
  for (Edges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    for (HEdges::iterator h = (*e)->hedges.begin(); h != (*e)->hedges.end(); ++h)
      if ((*h)->tail() == a &&
	  ((*h)->head() == b && (*h)->next->head() == c ||
	   (*h)->head() == c && (*h)->next->head() == b))
	return;
  addTriangle(a, b, c, pc, flag);
}

Octree<Face *> * Polyhedron::faceOctree () const
{
  return Octree<Face *>::octree(faces, bbox);
}

Octree<Cell *> * Polyhedron::cellOctree () const
{
  Cells ce;
  for (int i = 1; i < cells.size(); ++i)
    ce.push_back(cells[i]);
  return Octree<Cell *>::octree(ce, bbox);
}

void Polyhedron::describe (int i0)
{
  if (cells.empty())
    formCells();
  for (int i = i0; i < cells.size(); ++i) {
    Cell *c = cells[i];
    cerr << "cell " << i << ": " << c->boundary.size() << " shells; ";
    for (Shells::iterator s = c->boundary.begin(); s != c->boundary.end(); ++s)
      cerr << (*s)->hfaces.size() << " hfaces; euler = " << (*s)->euler() << " ";
    cerr << endl;
  }
}

bool outerLoop (const Vertices &ve, int pc)
{
  Vertex *u, *v, *w;
  extremalR(ve, u, v, w);
  return u != v && u != w && v != w &&
    LeftTurn(u->getP(), v->getP(), w->getP(), pc) == 1;
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

Polyhedron * complement (Polyhedron *a, Polyhedron *b)
{
  return a->boolean(b, Complement);
}

Polyhedron * intersection (Polyhedron *a, Polyhedron *b)
{
  return a->boolean(b, Intersection);
}

Polyhedron * box (double *b)
{
  Polyhedron *p = new Polyhedron;
  p->getVertex(new InputPoint(b[0], b[2], b[4]));
  p->getVertex(new InputPoint(b[1], b[2], b[4]));
  p->getVertex(new InputPoint(b[1], b[3], b[4]));
  p->getVertex(new InputPoint(b[0], b[3], b[4]));
  p->getVertex(new InputPoint(b[0], b[2], b[5]));
  p->getVertex(new InputPoint(b[1], b[2], b[5]));
  p->getVertex(new InputPoint(b[1], b[3], b[5]));
  p->getVertex(new InputPoint(b[0], b[3], b[5]));
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
  VVMap vvmap;
  PV3 o = PV3::constant(ox, oy, oz);
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f) {
    Vertices va, vb;
    (*f)->getBoundary(0)->loop(va);
    for (Vertices::iterator v = va.begin(); v != va.end(); ++v) {
      VVMap::iterator i = vvmap.find(*v);
      if (i == vvmap.end()) {
	PV3 p = o + r*(*v)->getP()->getP();
	Point *q = new InputPoint(p.x.mid(), p.y.mid(), p.z.mid());
	Vertex *w = b->getVertex(q);
	vvmap.insert(VVPair(*v, w));
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
    Vertices ve;
    (*f)->getBoundary(0)->loop(ve);
    PV3 p = (ve[0]->getP()->getP() + ve[1]->getP()->getP() 
	     + ve[2]->getP()->getP())/3.0;
    d = min(d, p.dot(p).mid());
  }
  return 1.0 - sqrt(d);
}

Polyhedron * sphereRefine (Polyhedron *a)
{
  Polyhedron *b = new Polyhedron;
  map<Edge *, Vertex *> evmap;
  for (Edges::iterator e = a->edges.begin(); e != a->edges.end(); ++e) {
    PV3 pm = (0.5*((*e)->getT()->getP()->getP() +
		   (*e)->getH()->getP()->getP())).unit();
    Point *ppm = new InputPoint(pm.x.mid(), pm.y.mid(), pm.z.mid());
    Vertex *vm = b->getVertex(ppm);
    evmap.insert(pair<Edge *, Vertex *>(*e, vm));
  }
  VVMap vvmap;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f) {
    Vertices ve;
    HEdge *e = (*f)->getBoundary(0);
    e->loop(ve);
    Vertex *u = b->getVertex(ve[0], vvmap), *v = b->getVertex(ve[1], vvmap),
      *w = b->getVertex(ve[2], vvmap), *uv = evmap.find(e->getE())->second,
      *vw = evmap.find(e->getNext()->getE())->second,
      *wu = evmap.find(e->getNext()->getNext()->getE())->second;
    b->addTriangle(u, uv, wu);
    b->addTriangle(uv, v, vw);
    b->addTriangle(vw, w, wu);
    b->addTriangle(uv, vw, wu);
  }
  return b;
}

Polyhedron * lbox (double x1, double x2, double y1, double y2, double z1, double z2)
{
  Polyhedron *a = new Polyhedron;
  Vertex *p[] = {a->getVertex(0, 0, z1), a->getVertex(x2, 0, z1),
		 a->getVertex(x2, y1, z1), a->getVertex(x1, y1, z1),
		 a->getVertex(x1, y2, z1), a->getVertex(0, y2, z1)};
  Vertex *q[] = {a->getVertex(0, 0, z2), a->getVertex(x2, 0, z2),
		 a->getVertex(x2, y1, z2), a->getVertex(x1, y1, z2),
		 a->getVertex(x1, y2, z2), a->getVertex(0, y2, z2)};
  for (int i = 1; i < 5; ++i) {
    a->addTriangle(p[0], p[i+1], p[i]);
    a->addTriangle(q[0], q[i], q[i+1]);
  }
  for (int i = 0; i < 6; ++i) {
    a->addTriangle(p[i], p[(i+1)%6], q[i]);
    a->addTriangle(q[i], p[(i+1)%6], q[(i+1)%6]);
  }
  return a;
}

Polyhedron * room (double x1, double x2, double x3, double y1, double y2,
		   double y3, double y4, double y5, double z1, double z2)
{
  Polyhedron *a = new Polyhedron;
  Vertex *p[] = {a->getVertex(0, 0, z1), a->getVertex(x3, 0, z1),
		 a->getVertex(x3, y2, z1), a->getVertex(x2, y2, z1),
		 a->getVertex(x2, y1, z1), a->getVertex(x1, y1, z1),
		 a->getVertex(x1, y4, z1), a->getVertex(x2, y4, z1),
		 a->getVertex(x2, y3, z1), a->getVertex(x3, y3, z1),
		 a->getVertex(x3, y5, z1), a->getVertex(0, y5, z1)};
  Vertex *q[] = {a->getVertex(0, 0, z2), a->getVertex(x3, 0, z2),
		 a->getVertex(x3, y2, z2), a->getVertex(x2, y2, z2),
		 a->getVertex(x2, y1, z2), a->getVertex(x1, y1, z2),
		 a->getVertex(x1, y4, z2), a->getVertex(x2, y4, z2),
		 a->getVertex(x2, y3, z2), a->getVertex(x3, y3, z2),
		 a->getVertex(x3, y5, z2), a->getVertex(0, y5, z2)};
  a->addTriangle(p[5], p[1], p[0]);
  a->addTriangle(p[5], p[4], p[1]);
  a->addTriangle(p[4], p[2], p[1]);
  a->addTriangle(p[4], p[3], p[2]);
  a->addTriangle(p[6], p[5], p[0]);
  a->addTriangle(p[11], p[6], p[0]);
  a->addTriangle(p[10], p[7], p[6]);
  a->addTriangle(p[11], p[10], p[6]);
  a->addTriangle(p[9], p[8], p[7]);
  a->addTriangle(p[10], p[9], p[7]);
  
  a->addTriangle(q[0], q[1], q[5]);
  a->addTriangle(q[1], q[4], q[5]);
  a->addTriangle(q[1], q[2], q[4]);
  a->addTriangle(q[2], q[3], q[4]);
  a->addTriangle(q[0], q[5], q[6]);
  a->addTriangle(q[0], q[6], q[11]);
  a->addTriangle(q[6], q[7], q[10]);
  a->addTriangle(q[6], q[10], q[11]);
  a->addTriangle(q[7], q[8], q[9]);
  a->addTriangle(q[7], q[9], q[10]);
  for (int i = 0; i < 12; ++i) {
    a->addTriangle(p[i], p[(i+1)%12], q[i]);
    a->addTriangle(q[i], p[(i+1)%12], q[(i+1)%12]);
  }
  return a;
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

int find (ID id, const Faces &faces)
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
  PV3 n = p->getN();
  Parameter k = p->getK();
  cerr << "(" << n.x.mid() << " " << n.y.mid() << " " << n.z.mid() << " "
       << k.mid() << ")" << endl;
}

void pp (Point *p)
{
  pp(p->getP());
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

void pids (const IDSet &ids)
{
  for (IDSet::const_iterator i = ids.begin(); i != ids.end(); ++i)
    cerr << *i << " ";
  cerr << endl;
}

void pe (Edge *e)
{
  cerr << "(";
  pp1(e->getT()->getP()->getP());
  pp1(e->getH()->getP()->getP());
  cerr << ")" << endl;
}

void pes (const Edges &ed)
{
  cerr << "(";
  for (Edges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e);
  cerr << ")" << endl;
}

void pes (const EdgeSet &ed)
{
  cerr << "(";
  for (EdgeSet::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e);
  cerr << ")" << endl;
}

void pe (HEdge *e)
{
  cerr << "(";
  pp1(e->tail()->getP()->getP());
  pp1(e->head()->getP()->getP());
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
  Vertices ve;
  e->loop(ve);
  pvs(ve);
}

void pf (Face *f)
{
  cerr << "(";
  const HEdges &fb = f->getBoundary();
  for (HEdges::const_iterator b = fb.begin(); b != fb.end(); ++b) {
    Vertices ve;
    (*b)->loop(ve);
    pvs(ve);
  }
  cerr << ")" << endl;
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

InputPoint * ipt (Point *v)
{
  return dynamic_cast<InputPoint *>(v);
}

SumPoint * spt (Point *v)
{
  return dynamic_cast<SumPoint *>(v);
}

CentroidPoint * cpt (Point *p)
{
  return dynamic_cast<CentroidPoint *>(p);
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
  PV3 u = v->getP() - w->getP();
  return sqrt(fabs(u.dot(u).mid()));
}

double distance (Point *a, Point *t, Point *h)
{
  PV3 tp = t->getP(), u = h->getP() - tp,
    p = tp + ((a->getP() - tp).dot(u)/u.dot(u))*u, w = a->getP() - p;
  Parameter ww = w.dot(w);
  return sqrt(fabs(ww.mid()));
}

double distance (Point *v, Plane *p)
{
  PV3 n = p->getN();
  double k = sqrt(fabs(n.dot(n).mid()));
  Parameter d = n.dot(v->getP()) + p->getK();
  return d.mid()/k;
}

double distance (Point *v, Point *a, Point *b, Point *c)
{
  PV3 n = (c->getP() - b->getP()).cross(a->getP() - b->getP());
  double k = sqrt(fabs(n.dot(n).mid()));
  Parameter d = n.dot(v->getP() - b->getP());
  return d.mid()/k;
}

double distanceEE (Point *a, Point *b, Point *c, Point *d)
{
  PV3 u = b->getP() - a->getP(), v = d->getP() - c->getP(),
    w = u.cross(v).unit();
  return (c->getP() - a->getP()).dot(w).mid();
}

double distance (Edge *e, Edge *f)
{
  return distanceEE(e->getT()->getP(), e->getH()->getP(),
		    f->getT()->getP(), f->getH()->getP());
}

void facesBB (const Faces &fa, double *bb, Faces &res)
{
  for (Faces::const_iterator f = fa.begin(); f != fa.end(); ++f)
    if (bboxOverlap((*f)->getBBox(), bb))
      res.push_back(*f);
}

void pfacesBB (const Faces &fa, double *bb, int i)
{
  Faces ff;
  facesBB(fa, bb, ff);
  void pfaces (const Faces &, int i);
  pfaces(ff, i);
}

void edgesNM (Shell *s)
{
  EdgeSet es;
  const HFaces &hf = s->getHFaces();
  for (HFaces::const_iterator f = hf.begin(); f != hf.end(); ++f) {
    HEdges he;
    (*f)->getF()->boundaryHEdges(he);
    for (HEdges::iterator e = he.begin(); e != he.end(); ++e)
      es.insert((*e)->getE());
  }
  for (EdgeSet::iterator e = es.begin(); e != es.end(); ++e) {
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
  int n1 = 0, n2 = 0;
  for (int i = 0; i < a->edges.size(); ++i) {
    Edge *e = a->edges[i];
    int n = e->HEdgesN();
    if (n != 0)
      if (n == 1) {
	ed.push_back(e);
	++n1;
      }
      else if (n != 2) {
	cerr << "non-manifold edge " << i << " hedges " << n << endl;
	++n2;
      }
  }
  if (n1 || n2)
    cerr << "dangling " << n1 << " non-manifold " << n2 << endl;
}

void edgesNM (Polyhedron *a)
{
  Edges ed;
  edgesNM(a, ed);
}
