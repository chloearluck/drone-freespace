#include "mink.h"

Polyhedron * minkowskiSum (Polyhedron *a, Polyhedron *b)
{
  Convolution *con = convolution(a, b);
  Octree<Face *> *octree = con->faceOctree();
  Polyhedron *c = new Polyhedron;
  VVMap vvmap;
  FaceSet fdone, fs;
  Face *f0 = con->minkowskiInit(fdone, octree, c, vvmap);
  Faces st;
  st.push_back(f0);
  while (!st.empty()) {
    Face *f = *st.rbegin();
    st.pop_back();
    if (fs.insert(f).second)
      con->expand(fdone, octree, c, vvmap, f, st);
  }
  Polyhedron *d = triangulate(fs);
  delete octree;
  delete c;
  delete con;
  return d;
}

class HullFace {
 public:
  HullFace (HEdge *e, HullFace *prev, HullFace *next) : e(e), prev(prev), next(next) {}
  bool inCset (HEdge *f) const
  { return find(cset.begin(), cset.end(), f) != cset.end(); }
  void updateCset (HullFace *h, HEdge *f);
  bool conflict (HEdge *f) const;
  void cone (HEdges &hedges) const;

  HEdges cset;
  HEdge *e;
  HullFace *prev, *next;
};
HullFace * initHull (HEdges &hedges);
HullFace * updateHull (HullFace *hull, HEdge *e, bool &flag);
HullFace * updateHullAux (HullFace *fs, HullFace *fe, HEdge *e);
void deleteHull (HullFace *hull);

Face * Convolution::minkowskiInit (FaceSet &fdone, Octree<Face *> *octree,
				   Polyhedron *a, VVMap &vvmap)
{
  Vertex *v;
  Face *f = rmaxFace(v);
  subdivide(f, fdone, octree, a, vvmap);
  Vertex *w = a->getVertex(v, vvmap);
  for (Faces::iterator g = a->faces.begin(); g != a->faces.end(); ++g)
    if ((*g)->boundaryVertex(w))
      return *g;
  return 0;
}

Face * Convolution::rmaxFace (Vertex *&v)
{
  v = vertices[0];
  for (int i = 1; i < vertices.size(); ++i)
    if (PointOrderR(v->p, vertices[i]->p) == 1)
      v = vertices[i];
  Face *f = 0;
  for (Edges::const_iterator e = v->edges.begin(); e != v->edges.end(); ++e)
    for (HEdges::iterator h = (*e)->hedges.begin(); h != (*e)->hedges.end(); ++h)
      if ((*h)->tail() == v) {
	Face *g = (*h)->f;
	if (!f || f->p.id != g->p.id && NormalOrderR(f->getP(), g->getP()) == 1)
	  f = g;
      }
  return f;
}

void Convolution::subdivide (Face *f, FaceSet &fdone, Octree<Face *> *octree,
			     Polyhedron *a, VVMap &vvmap)
{
  if (!fdone.insert(f).second)
    return;
  iedges.clear();
  intersectFF(f, fdone, octree);
  intersectFFF(f, fffvmap);
  for (EdgeSet::iterator e = iedges.begin(); e != iedges.end(); ++e)
    (*e)->sortVertices();
  int nf = a->faces.size();
  subfaces(f, a, vvmap);
  sortHEdges(a, nf);
}

void Convolution::intersectFF (Face *f, const FaceSet &fdone, 
			       Octree<Face *> *octree)
{
  Faces ff;
  octree->find(f->bbox, ff);
  for (Faces::iterator g = ff.begin(); g != ff.end(); ++g)
    if (fdone.find(*g) == fdone.end())
      Polyhedron::intersectFF(f, *g, efvmap);
}

void Convolution::expand (FaceSet &fdone, Octree<Face *> *octree,
			  Polyhedron *a, VVMap &vvmap, Face *f, Faces &st)
{
  HEdges he;
  f->boundaryHEdges(he);
  for (HEdges::iterator e = he.begin(); e != he.end(); ++e) {
    Faces fa;
    neighborFaces(*e, fa);
    for (Faces::iterator g = fa.begin(); g != fa.end(); ++g)
      subdivide(*g, fdone, octree, a, vvmap);
    (*e)->e->sortHEdges();
    Face *h = f->hfaces[0].neighbor(*e)->f;
    st.push_back(h);
  }
}

void Convolution::neighborFaces (HEdge *e, Faces &fa) const
{
  Vertex *t = e->tail(), *h = e->head();
  PPPair pp(t->p < h->p ? t->p : h->p, t->p < h->p ? h->p : t->p);
  const EdgeSet &es = ppemap.find(pp)->second;
  for (EdgeSet::const_iterator f = es.begin(); f != es.end(); ++f) {
    const HEdges &he = (*f)->hedges;
    for (HEdges::const_iterator g = he.begin(); g != he.end(); ++g)
      fa.push_back((*g)->f);
  }
}

int NormalOrderR::sign ()
{
  PV3 m = p->getN(), n = q->getN(), r = Rdir.getP();
  Parameter mm = m.dot(m), mr = m.dot(r), nn = n.dot(n), nr = n.dot(r),
    k = nr*nr*mm - mr*mr*nn;
  return k.sign();
}

void sortHEdges (Polyhedron *a, int nf)
{
  EdgeSet es;
  for (int i = nf; i < a->faces.size(); ++i) {
    HEdges ed;
    a->faces[i]->boundaryHEdges(ed);
    for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
      es.insert((*e)->getE());
  }
  for (EdgeSet::iterator e = es.begin(); e != es.end(); ++e)
    (*e)->sortHEdges();
}

void sphereBBox (Edge *e, double *bbox)
{
  FNormal en1(e->getHEdge(0)->getF()), en2(e->getHEdge(1)->getF());
  sphereBBox(&en1, &en2, bbox);
}

void sphereBBox (Point *t, Point *h, double *bbox)
{
  GreatCircleN n(t, h);
  GreatCircleMinX px(&n);
  GreatCircleMinY py(&n);
  GreatCircleMinZ pz(&n);
  int tx = TripleProduct(&n, t, &px), hx = TripleProduct(&n, &px, h),
    ty = TripleProduct(&n, t, &py), hy = TripleProduct(&n, &py, h),
    tz = TripleProduct(&n, t, &pz), hz = TripleProduct(&n, &pz, h);
  UnitVector tu(t), hu(h);
  PV3 tp = tu.getP(), hp = hu.getP();
  bbox[0] = min(tp.x.lb(), hp.x.lb());
  bbox[1] = max(tp.x.ub(), hp.x.ub());
  bbox[2] = min(tp.y.lb(), hp.y.lb());
  bbox[3] = max(tp.y.ub(), hp.y.ub());
  bbox[4] = min(tp.z.lb(), hp.z.lb());
  bbox[5] = max(tp.z.ub(), hp.z.ub());
  if (tx == 1 && hx == 1) {
    UnitVector u(&px);
    bbox[0] = min(bbox[0], u.getP().x.lb());
  }
  else if (tx == -1 && hx == -1) {
    UnitVector u(&px);
    bbox[1] = max(bbox[1], - u.getP().x.lb());
  }
  if (ty == 1 && hy == 1) {
    UnitVector u(&py);
    bbox[2] = min(bbox[2], u.getP().y.lb());
  }
  else if (ty == -1 && hy == -1) {
    UnitVector u(&py);
    bbox[3] = max(bbox[3], - u.getP().y.lb());
  }
  if (tz == 1 && hz == 1) {
    UnitVector u(&pz);
    bbox[4] = min(bbox[4], u.getP().z.lb());
  }
  else if (tz == -1 && hz == -1) {
    UnitVector u(&pz);
    bbox[5] = max(bbox[5], - u.getP().z.lb());
  }
}

void sphereBBox (const HEdges &ed, double *bbox)
{
  bbox[0] = bbox[2] = bbox[4] = 1.0;
  bbox[1] = bbox[3] = bbox[5] = -1.0;
  bool flags[6] = {false, false, false, false, false, false};
  int n = ed.size();
  Points pts;
  for (int i = 0; i < n; ++i)
    pts.push_back(new EENormal(ed[i], ed[(i+1)%n]));
  for (int i = 0; i < n; ++i) {
    Point *t = pts[i], *h = pts[(i+1)%n];
    double eb[6];
    sphereBBox(t, h, eb);
    merge(eb, bbox); 
    GreatCircleN n(t, h);
    for (int j = 0; j < 3; ++j) {
      int s = coordinate(&n, j);
      if (s == 0) {
	bbox[2*j] = -1.0;
	bbox[2*j+1] = 1.0;
      }
      else if (s == -1)
	flags[2*j] = true;
      else
	flags[2*j+1] = true;
    }
  }
  for (int i = 0; i < 3; ++i) {
    if (flags[2*i] && !flags[2*i+1])
      bbox[2*i] = -1.0;
    if (!flags[2*i] && flags[2*i+1])
      bbox[2*i+1] = 1.0;
  }
}

int coordinate (Point *a, int c)
{
  Parameter p = a->getP()[c];
  if (p.lb() > 0.0)
    return 1;
  if (p.ub() < 0.0)
    return -1;
  return 0;
}

void merge (double *bbox1, double *bbox2)
{
  for (int i = 0; i < 3; ++i) {
    bbox2[2*i] = min(bbox1[2*i], bbox2[2*i]);
    bbox2[2*i+1] = max(bbox1[2*i+1], bbox2[2*i+1]);
  }
}

BSPElt::BSPElt (Vertex *v, const HEdges &ed) : l(0u), ed(ed)
{
  d.v = v;
  int m = ed.size();
  sphereBBox(ed, bbox);
}

BSPElt::BSPElt (Edge *e) : l(0)
{
  d.e = e;
  sphereBBox(e, bbox);
}

BSPElt::BSPElt (Face *f) : l(0)
{
  d.f = f;
  FNormal n(f);
  UnitVector u(&n);
  u.getBBox(bbox);
}

BSPElt::BSPElt (const BSPElt &x)
{
  d.v = x.d.v;
  ed = x.ed;
  l = x.l;
  for (int i = 0; i < 6; ++i)
    bbox[i] = x.bbox[i];
}

int BSPElt::side (Point *r) const
{
  PV3 p(Parameter::interval(bbox[0], bbox[1]),
	Parameter::interval(bbox[2], bbox[3]),
	Parameter::interval(bbox[4], bbox[5]));
  Parameter k = p.dot(r->getP());
  return k.sign(false);
}

/* acp0 version
int BSPElt::side (Point *r) const
{
  PV3 p = r->getP();
  IParameter x(bbox[0], bbox[1]), y(bbox[2], bbox[3]), z(bbox[4], bbox[5]),
    k = x*p.x + y*p.y + z*p.z;
  return k.sign();
}
*/

void BSPTree (BSPElts &aelts, BSPElts &belts, BSPElts &ea, BSPElts &eb, 
	      int nmax, int dmax, ID c)
{
  if (dmax == 0 || aelts.size() < nmax && belts.size() < nmax)
    BSPLeaf(aelts, belts, ea, eb);
  else {
    InputPoint r(randomNumber(-1.0, 1.0), randomNumber(-1.0, 1.0),
		 randomNumber(-1.0, 1.0));
    BSPElts aelts1, aelts2, belts1, belts2;
    BSPPartition(aelts, &r, c, aelts1, aelts2);
    BSPPartition(belts, &r, c, belts1, belts2);
    if (aelts1.size() == aelts.size() || aelts2.size() == aelts.size() ||
	belts1.size() == belts.size() || belts2.size() == belts.size())
      BSPLeaf(aelts, belts, ea, eb);
    else {
      BSPTree(aelts1, belts1, ea, eb, nmax, dmax - 1, c + 1u);
      BSPTree(aelts2, belts2, ea, eb, nmax, dmax - 1, c + 1u);
    }
  }
}

void BSPPartition (BSPElts &elts, Point *r, ID c, BSPElts &elts1, BSPElts &elts2)
{
  for (BSPElts::iterator e = elts.begin(); e != elts.end(); ++e) {
    int s = e->side(r);
    if (s == -1)
      elts1.push_back(*e);
    else if (s == 1)
      elts2.push_back(*e);
    else {
      elts1.push_back(*e);
      BSPElt f(*e);
      f.l += 1u << c;
      elts2.push_back(f);
    }
  }
}

void BSPLeaf (const BSPElts &aelts, const BSPElts &belts, BSPElts &ea, 
	      BSPElts &eb)
{
  for (BSPElts::const_iterator e = aelts.begin(); e != aelts.end(); ++e)    
    for (BSPElts::const_iterator f = belts.begin(); f != belts.end(); ++f)
      if (e->compatible(*f) && bboxOverlap(e->bbox, f->bbox)) {
	ea.push_back(*e);
	eb.push_back(*f);
      }
}

bool HullFace::conflict (HEdge *f) const
{
  int s = circulationEEE(f, e, next->e);
  if (s != 0) return s == 1;
  Point *a = e->tail()->getP(), *b = e->head()->getP(),
    *c = next->e->head()->getP(), *d = f->head()->getP();
  return d->onLine(a, b) || d->onLine(a, c) ||
    DegenerateConflict(a, b, c, d) == 1;
}

void HullFace::updateCset (HullFace *h, HEdge *f)
{
  for (HEdges::iterator g = h->cset.begin(); g != h->cset.end(); ++g)
    if (*g != f && conflict(*g))
      cset.push_back(*g);
}

void HullFace::cone (HEdges &hedges) const
{
  const HullFace *h = this;
  do {
    hedges.push_back(h->e);
    h = h->next;
  }
  while (h != this);
}

int DegenerateConflict::sign ()
{
  PV3 u = b->getP() - a->getP(), v = c->getP() - a->getP(), w = d->getP() - a->getP(),
    x = u.cross(v);
  if (x.tripleProduct(v, w).sign() == 1)
    return 1;
  return x.tripleProduct(w, u).sign();
}

bool convexCone (Vertex *v, HEdges &hedges)
{
  HEdges he;
  for (int i = 0; i < v->EdgesN(); ++i) {
    HEdge *e = v->getEdge(i)->getHEdge(0);
    he.push_back(e->tail() == v ? e : e->ccw());
  }
  HullFace *hull = initHull(he);
  if (!hull)
    return false;
  for (int i = 3; i < he.size(); ++i) {
    bool flag;
    hull = updateHull(hull, he[i], flag);
    if (!flag) {
      deleteHull(hull);
      return false;
    }
  }
  hull->cone(hedges);
  deleteHull(hull);
  return convexOrder(hedges);
}

HullFace * initHull (HEdges &hedges)
{
  HEdge *e = hedges[0], *f = 0, *g = 0;
  for (int i = 2; i < hedges.size(); ++i) {
    int s = circulationEEE(e, hedges[1], hedges[i]);
    if (s) {
      if (s == -1) {
	f = hedges[1];
	g = hedges[i];
	hedges[i] = hedges[2];
	hedges[2] = g;
      }
      else {
	f = hedges[i];
	g = hedges[1];
	hedges[i] = hedges[2];
	hedges[2] = f;
      }
      break;
    }
  }
  if (!f)
    return 0;
  HullFace *he = new HullFace(e, 0, 0), *hf = new HullFace(f, he, 0), 
    *hg = new HullFace(g, hf, he);
  he->next = hf;
  he->prev = hg;
  hf->next = hg;
  HullFace *faces[3] = {he, hf, hg};
  for (int i = 0; i < 3; ++i)
    for (int j = 3; j < hedges.size(); ++j)
      if (faces[i]->conflict(hedges[j]))
	faces[i]->cset.push_back(hedges[j]);
  return he;
}

int circulationEEE (HEdge *e, HEdge *f, HEdge *g)
{
  return Orientation(e->tail()->getP(), e->head()->getP(), 
		     f->head()->getP(), g->head()->getP());
}

HullFace * updateHull (HullFace *hull, HEdge *e, bool &flag)
{
  flag = true;
  HullFace *fs = hull;
  if (fs->inCset(e)) {
    while (fs->prev != hull && fs->prev->inCset(e))
      fs = fs->prev;
    if (fs->prev == hull) {
      flag = false;
      return hull;
    }
  }
  else {
    do
      fs = fs->next;
    while (fs != hull && !fs->inCset(e));
    if (fs == hull)
      return hull;
  }
  HullFace *fe = fs;
  while (fe->next->inCset(e))
    fe = fe->next;
  return updateHullAux(fs, fe, e);
}

HullFace * updateHullAux (HullFace *fs, HullFace *fe, HEdge *e)
{
  HullFace *f1 = new HullFace(fs->e, fs->prev, 0), *f2 = new HullFace(e, f1, fe->next);
  f1->next = f2;
  f1->updateCset(fs->prev, e);
  f1->updateCset(fs, e);
  f2->updateCset(fe, e);
  f2->updateCset(fe->next, e);
  f1->prev->next = f1;
  f2->next->prev = f2;
  fe->next = 0;
  while (fs) {
    HullFace *ptr = fs->next;
    delete fs;
    fs = ptr;
  }
  return f1;
}

void deleteHull (HullFace *hull)
{
  hull->prev->next = 0;
  while (hull) {
    HullFace *nhull = hull->next;
    delete hull;
    hull = nhull;
  }
}

bool convexOrder (const HEdges &hedges)
{
  HEdge *e1 = hedges[0], *e2 = hedges[1], *e3 = hedges[2];
  do {
    e1 = e1->getNext()->getNext()->ccw();
    if (e1 == e2) 
      return true;
    if (e1 == e3)
      return false;
  }
  while (e1 != hedges[0]);
  return false;
}

Polyhedron * triangulate (const FaceSet &fs)
{
  Polyhedron *a = new Polyhedron;
  VVMap vvmap;
  for (FaceSet::const_iterator f = fs.begin(); f != fs.end(); ++f) {
    Triangles tr;
    (*f)->triangulate(tr);
    int pc = (*f)->getPC();
    for (Triangles::iterator t = tr.begin(); t != tr.end(); ++t)
      a->getTriangle(t->a, t->b, t->c, vvmap, pc);
  }
  return a;
}

Polyhedron * minkowskiSumFull (Polyhedron *a, Polyhedron *b)
{
  Polyhedron *con = convolution(a, b);
  Polyhedron *c = con->subdivide(), *d = new Polyhedron;
  Shells sh;
  c->formShells(sh);
  VVMap vvmap;
  for (Shells::iterator s = sh.begin(); s != sh.end(); ++s)
    if (minkowskiShell(a, b, *s)) {
      const HFaces &hf = (*s)->getHFaces();
      for (HFaces::const_iterator f = hf.begin(); f != hf.end(); ++f)
	d->getTriangle((*f)->getF(), vvmap);
    }
  delete c;
  delete con;
  deleteShells(sh);
  return d;
}

bool minkowskiShell (Polyhedron *a, Polyhedron *b, Shell *s)
{
  const HFaces &hf = s->getHFaces();
  for (HFaces::const_iterator f = hf.begin(); f != hf.end(); ++f)
    if (!(*f)->pos())
      return false;
  Point *p = hf[0]->getF()->getBoundary(0)->tail()->getP();
  Polyhedron *c = a->negativeTranslate(p);
  bool res = !c->intersects(b);
  delete c;
  return res;
}

Convolution * convolution (Polyhedron *a, Polyhedron *b)
{
  Convolution *c = new Convolution;
  VVPairVMap vmap;
  FaceDescSet fds;
  sumVF(a, b, true, fds, vmap, c);
  sumVF(b, a, false, fds, vmap, c);
  sumEE(a, b, fds, vmap, c);
  for (Edges::iterator e = c->edges.begin(); e != c->edges.end(); ++e)
    (*e)->sortHEdges();
  return c;
}

void sumVF (Polyhedron *a, Polyhedron *b, bool avflag, FaceDescSet &fds,
	    VVPairVMap &vmap, Polyhedron *con)
{
  BSPElts aelts, belts, ea, eb;
  convexVertices(a, aelts);
  for (Faces::iterator f = b->faces.begin(); f != b->faces.end(); ++f)
    belts.push_back(BSPElt(*f));
  BSPTree(aelts, belts, ea, eb);
  for (int i = 0; i < ea.size(); ++i)
    if (compatibleVF(ea[i].ed, eb[i].d.f))
      sumVF(ea[i].d.v, eb[i].d.f, avflag, fds, vmap, con);
}

void convexVertices (Polyhedron *a, BSPElts &elts)
{
  for (Vertices::iterator v = a->vertices.begin(); v != a->vertices.end(); ++v) {
    HEdges ed;
    if (convexCone(*v, ed))
      elts.push_back(BSPElt(*v, ed));
  }
}

bool compatibleVF (HEdges &ed, Face *f)
{
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    if (InnerProductEF(*e, f) == 1)
      return false;
  return true;
}

int InnerProductEF::sign ()
{
  return e->getU().dot(f->getP()->getN()).sign();
}

void sumVF (Vertex *v, Face *f, bool avflag, FaceDescSet &fds, 
	    VVPairVMap &vmap, Polyhedron *con)
{
  Vertices ve, vo;
  f->getBoundary(0)->loop(ve);
  for (Vertices::iterator w = ve.begin(); w != ve.end(); ++w)
    vo.push_back(sumVV(v, *w, avflag, vmap, con));
  if (!inputPerturbed && !fds.insert(vo).second)
    return;
  con->addTriangle(vo[0], vo[1], vo[2]);
}

Vertex * sumVV (Vertex *a, Vertex *b, bool aflag, VVPairVMap &vmap, 
		Polyhedron *con)
{
  Vertex *aa = aflag ? a : b, *bb = aflag ? b : a;
  VVPair vp(aa, bb);
  VVPairVMap::iterator iter = vmap.find(vp);
  if (iter != vmap.end())
    return iter->second;
  Vertex *v = con->getVertex(new SumPoint(aa->getP(), bb->getP()));
  vmap.insert(pair<VVPair, Vertex *>(vp, v));
  return v;
}

void sumEE (Polyhedron *a, Polyhedron *b, FaceDescSet &fds,
	    VVPairVMap &vmap, Polyhedron *con)
{
  BSPElts aelts, belts, ea, eb;
  convexEdges(a, aelts);
  convexEdges(b, belts);
  BSPTree(aelts, belts, ea, eb);
  for (int i = 0; i < ea.size(); ++i) {
    Edge *eai = ea[i].d.e, *ebi = eb[i].d.e;
    bool aflag;
    if (compatibleEE(eai, ebi, aflag))
      sumEE(eai, ebi, aflag, fds, vmap, con);
  }
}

void convexEdges (Polyhedron *a, BSPElts &elts)
{
  for (Edges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    if (convexEdge(*e) == 1)
      elts.push_back(BSPElt(*e));
}

int convexEdge (Edge *e)
{
  HEdge *e1 = e->getHEdge(0), *e2 = e->getHEdge(1);
  if (e1->getF()->getP()->getid() == e2->getF()->getP()->getid())
    return 0;
  return ConvexEdge(e1, e2);
}

int ConvexEdge::sign ()
{
  return e1->getU().tripleProduct(e1->getF()->getP()->getN(),
				  e2->getF()->getP()->getN()).sign();
}

bool compatibleEE (Edge *e, Edge *f, bool &aflag)
{
  int s1 = CompatibleEdge(f, e->getHEdge(0)->getF()), 
    s2 = CompatibleEdge(f, e->getHEdge(1)->getF());
  aflag = s1 == 1 || s1 == 0 && s2 == -1;
  if (s1 == s2)
    return false;
  int s3 = CompatibleEdge(e, f->getHEdge(0)->getF());
  if (s1 != 0 && s1 == s3 || s2 != 0 && s2 == - s3)
    return false;
  int s4 = CompatibleEdge(e, f->getHEdge(1)->getF());
  if (s3 == s4)
    return false;
  return s4 == 0 || (s1 == 0 ? s2 == - s4 : s1 == s4);
}

int CompatibleEdge::sign ()
{
  return e->getU().dot(f->getP()->getN()).sign();
}

void sumEE (Edge *e, Edge *f, bool aflag, FaceDescSet &fds,
	    VVPairVMap &vmap, Polyhedron *con)
{
  Edge *ee = aflag ? e : f, *ff = aflag ? f : e;
  Vertex *a = sumVV(ee->getT(), ff->getT(), aflag, vmap, con),
    *b = sumVV(ee->getH(), ff->getT(), aflag, vmap, con),
    *c = sumVV(ee->getH(), ff->getH(), aflag, vmap, con),
    *d = sumVV( ee->getT(), ff->getH(), aflag, vmap, con);
  if (inputPerturbed)
    con->addRectangle(a, b, c, d);
  else if (!(a->getP()->onLine(b->getP(), c->getP()))) {
    Vertex *vo[] = {a, b, c, d};
    if (fds.insert(FaceDesc(vo, 4)).second)
      con->addRectangle(a, b, c, d);
  }
}

// debug


void neighborFaces (HEdge *e, Octree<Face *> *octree, Faces &fa)
{
  Vertex *t = e->tail(), *h = e->head();
  Faces fcand;
  octree->find(t->getBBox(), fcand);
    for (Faces::iterator f = fcand.begin(); f != fcand.end(); ++f) {
      Plane *p = (*f)->getP();
      if (bboxOverlap((*f)->getBBox(), h->getBBox()) &&
	  t->getP()->side(p) == 0 && h->getP()->side(p) == 0)
	fa.push_back(*f);
    }
}
