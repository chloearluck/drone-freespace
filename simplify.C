#include "simplify.h"

void simplify (Polyhedron *a, double d, bool opt2)
{
  simplify1(a, d);
  simplify2(a, d, opt2);
  a->removeNullFaces();
}

void simplify1 (Polyhedron *a, double d)
{
  double t = getTime();
  IFeatureQ fq, nfq;
  for (Edges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    addCollapseFlips(*e, d, fq);
  Octree<Face *> *octree = faceOctree(a);
  VPMap vpmap;
  int n1 = 0, n2 = 0;
  while (true) {
    int n = n1 + n2;
    simplify1(a, d, fq, nfq, octree, vpmap, n1, n2);
    if (n == n1 + n2 || nfq.empty())
      break;
    fq = nfq;    
  }
  delete octree;
  t = getTime() - t;
  describe1(t, n1, n2, vpmap);
}  

bool closeVE (Vertex *v, Vertex *w, Vertex *x, double d)
{
  Point *a = v->getP(), *t = w->getP(), *h = x->getP();
  return DistancePL(a, t, h, d) == 1 && onEdge(a, t, h, true);
}

void addCollapseFlips (Edge *e, double d, IFeatureQ &fq)
{
  if (e->HEdgesN() != 2)
    return;
  if (DistancePP(e->getT()->getP(), e->getH()->getP(), d) == 1)
    fq.push(IFeature(e));
  else
    addFlips(e, d, fq);
}

void addFlips (Edge *e, double d, IFeatureQ &fq)
{
  if (e->HEdgesN() == 2)
    for (int i = 0; i < 2; ++i)
      if (closeVE(e->getHEdge(i)->getNext()->head(),
		  e->getT(), e->getH(), d))
	fq.push(IFeature(e->getHEdge(i)));
}

Octree<Face *> * faceOctree (Polyhedron *a, double s)
{
  Faces fa;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f)
    if ((*f)->getBoundary())
      fa.push_back(*f);
  return Octree<Face *>::octree(fa, a->bbox, s);
}

void simplify1 (Polyhedron *a, double d, IFeatureQ &fq, IFeatureQ &nfq,
		Octree<Face *> *octree, VPMap &vpmap, int &n1, int &n2)
{
  while (!fq.empty()) {
    IFeature f = fq.top();
    fq.pop();
    if (f.valid())
      if (!f.v) {
	if (collapse(a, d, f.e, fq, nfq, octree, vpmap))
	  ++n1;
      }
      else if (flip(a, d, f.getH(), fq, nfq, octree, vpmap))
	++n2;
  }
}

bool collapse (Polyhedron *a, double d, Edge *e, IFeatureQ &fq, IFeatureQ &nfq,
	       Octree<Face *> *octree, VPMap &vpmap)
{
  if (!collapsible(e))
    return false;
  HEdge *e1 = e->getHEdge(0), *e2 = e->getHEdge(1);
  Vertex *t = e1->tail(), *h = e1->head(), *v = e1->getNext()->head(),
    *w = e2->getNext()->head();
  PTR<Point> tp = t->getP(), hp = h->getP(), m = new CentroidPoint(tp, hp);
  PV3 p = m->getApprox();
  PTR<Point> q = new Point(p.x.mid(), p.y.mid(), p.z.mid(), false),
    cp = a->findPoint(q) ? m : q;
  a->removeLoop(e1);
  a->removeLoop(e2);
  vpmap.insert(VPPair(h, h->getP()));
  a->moveVertex(h, cp);
  Faces ft = t->incidentFaces();
  for (Faces::iterator f = ft.begin(); f != ft.end(); ++f)
    a->replaceVertex(*f, t, h);
  if (badCollapse(h, octree)) {
    for (Faces::iterator f = ft.begin(); f != ft.end(); ++f)
      a->replaceVertex(*f, h, t);
    a->moveVertex(h, hp);
    octree->insert(a->addTriangle(t, h, v));
    octree->insert(a->addTriangle(h, t, w));
    addStar(t, d, nfq);
    addStar(h, d, nfq);
    return false;
  }
  addStar(h, d, fq);
  Faces fh = h->incidentFaces();
  for (Faces::iterator f = fh.begin(); f != fh.end(); ++f)
    octree->insert(*f);
  return true;
}

bool collapsible (Edge *e)
{
  int n = e->HEdgesN();
  if (n == 1)
    return true;
  if (n > 2)
    return false;
  Vertex *t = e->getT(), *h = e->getH(), *v = e->getHEdge(0)->getNext()->head(),
    *w = e->getHEdge(1)->getNext()->head();
  HEdges et = t->outgoingHEdges(), eh = h->outgoingHEdges();
  set<Vertex *> vs;
  for (HEdges::iterator f = et.begin(); f != et.end(); ++f) {
    Vertex *x = (*f)->head();
    if (x != h && x != v && x != w)
      vs.insert(x);
  }
  for (HEdges::iterator f = eh.begin(); f != eh.end(); ++f)
    if (vs.find((*f)->head()) != vs.end())
      return false;
  return true;
}

bool badCollapse (Vertex *v, Octree<Face *> *octree)
{
  Faces fv = v->incidentFaces();
  for (Faces::iterator f = fv.begin(); f != fv.end(); ++f) {
    Vertices ve = (*f)->getBoundary()->loop();
    if (ve[0]->getP()->onLine(ve[1]->getP(), ve[2]->getP()))
      return true;
  }
  for (Faces::iterator f = fv.begin(); f != fv.end(); ++f) {
    Faces fa;
    octree->find((*f)->getBBox(), fa);
    for (Faces::iterator g = fa.begin(); g != fa.end(); ++g)
      if ((*g)->getBoundary() &&
	  find(fv.begin(), fv.end(), *g) == fv.end() &&
	  (*f)->intersects(*g, false))
	return true;
  }
  return false;
}

void addStar (Vertex *v, double d, IFeatureQ &fq)
{
  HEdges vout = v->outgoingHEdges();
  for (HEdges::iterator e = vout.begin(); e != vout.end(); ++e) {
    if (DistancePP((*e)->tail()->getP(), (*e)->head()->getP(), d) == 1)
      fq.push(IFeature((*e)->getE()));
    HEdges ev = (*e)->edgeLoop();
    for (HEdges::iterator h = ev.begin(); h != ev.end(); ++h)
      if (closeVE((*h)->getNext()->head(), (*h)->tail(), (*h)->head(), d))
	fq.push(IFeature(*h));
  }
}

bool flip (Polyhedron *a, double d, HEdge *e, IFeatureQ &fq, IFeatureQ &nfq,
	   Octree<Face *> *octree, VPMap &vpmap)
{
  Vertex *t = e->tail(), *h = e->head(), *v = e->getNext()->head();
  if (closerPair(t->getP(), h->getP(), v->getP(), t->getP()))
    return collapse(a, d, e->getNext()->getE(), fq, nfq, octree, vpmap);
  if (closerPair(t->getP(), h->getP(), h->getP(), v->getP()))
    return collapse(a, d, e->getNext()->getNext()->getE(), fq, nfq, octree, vpmap);
  if (!flippable(e, d))
    return false;
  HEdge *vw = flip(a, e);
  octree->insert(vw->getF());
  octree->insert(vw->ccw()->getF());
  if (badFlip(vw, octree)) {
    vw = flip(a, vw);
    octree->insert(vw->getF());
    octree->insert(vw->ccw()->getF());
    addFlips(vw->getE(), d, nfq);
    return false;
  }
  HEdge *wh = vw->getNext(), *hv = wh->getNext(), *vt = vw->ccw()->getNext(),
    *th = vt->getNext();
  addCollapseFlips(vw->getE(), d, fq);
  addFlips(wh->getE(), d, fq);
  addFlips(hv->getE(), d, fq);
  addFlips(vt->getE(), d, fq);
  addFlips(th->getE(), d, fq);
  return true;
}

bool flippable (HEdge *e, double d)
{
  if (e->getE()->HEdgesN() != 2)
    return false;
  HEdge *en = e->getNext(), *fn = e->ccw()->getNext();
  Vertex *t = e->tail(), *h = e->head(), *v = en->head(), *w = fn->head();
  if (v->connected(w) || t->getP()->onLine(v->getP(), w->getP()) ||
      h->getP()->onLine(v->getP(), w->getP()))
    return false;
  if (!(closeVE(t, v, w, d) || closeVE(h, v, w, d)))
    return true;
  return closerPair(v->getP(), w->getP(), t->getP(), h->getP());
}

HEdge * flip (Polyhedron *a, HEdge *e)
{
  Vertex *t = e->tail(), *h = e->head(), *v = e->getNext()->head(),
    *w = e->ccw()->getNext()->head();
  a->removeLoop(e->ccw());
  a->removeLoop(e);
  Face *f = a->addTriangle(v, w, h), *g = a->addTriangle(w, v, t);
  return f->getBoundary();
}

bool badFlip (HEdge *e, Octree<Face *> *octree)
{
  for (int i = 0; i < 2; ++i, e = e->ccw()) {
    Face *f = e->getF();
    Faces fa;
    octree->find(f->getBBox(), fa);
    for (Faces::iterator g = fa.begin(); g != fa.end(); ++g)
      if ((*g)->getBoundary() && f->intersects(*g, false))
	return true;
  }
  return false;
}

void describe1 (double t, int n1, int n2, const VPMap &vpmap)
{
  cerr << "simplify 1: time = " << t;
  if (n1 > 0 || n2 > 0) {
    int nv = vpmap.size();
    double dmax = 0.0, dr = nv ? displacement(vpmap, dmax)/nv : 0.0;
    if (n1)
      cerr << " collapse = " << n1;
    if (n2)
      cerr << " flip = " << n2;
    if (nv)
      cerr << endl << "            moved = " << nv << " dav = " << dr
	   << " dmax = " << dmax;
  }
  cerr << endl;
}

double displacement (const VPMap &vpmap, double &dmax)
{
  dmax = 0.0;
  if (vpmap.empty())
    return 0.0;
  double dr = 0.0;
  for (VPMap::const_iterator i = vpmap.begin(); i != vpmap.end(); ++i) {
    PV3 u = i->first->getP()->getApprox(1e-6) - i->second->getApprox(1e-6);
    double di = 0.0;
    for (int j = 0; j < 3; ++j)
      di += fabs(u[j].mid());
    dr += di;
    dmax = max(dmax, di);
  }
  return dr;
}

void simplify2 (Polyhedron *a, double d, bool opt2)
{
  double t = getTime();
  int n = 0;
  VPMap vpmap;
  FeatureSet fslp = smallFeatures(a, 3.465*d, 0), fs = smallFeatures(fslp, d);
  while (separation(fslp) < 0.999*d) {
    ++n;
    simplify2s(a, d, fslp, vpmap);
  }
  if (opt2 && n > 0) {
    describe2(n, fs, getTime() - t, vpmap);
    double dmax, dprev = displacement(vpmap, dmax), bound = 1.0;
    VPMap opos = savePos(fslp), ovpmap = vpmap;
    int rcount = 0;
    while (bound > 1e-5) {
      ++n;
      simplify2s(a, d, fslp, vpmap, false, bound);
      double dbest = displacement(vpmap, dmax);
      if (dbest + 1e-8 > dprev)
        break;
      while (separation(fslp) < 0.999*d) {
	++n;
	simplify2s(a, d, fslp, vpmap);
      }
      double dcurr = displacement(vpmap, dmax);
      if (dcurr > dprev) {
        restorePos(a, opos);
	vpmap = ovpmap;
        dcurr = dprev;
      }
      else {
	opos = savePos(fslp);
	ovpmap = vpmap;
      }
      double ratio = (dprev - dcurr)/(dprev - dbest);
      if (ratio < 0.25) {
	bound *= 0.5;
        rcount = 0;
      }
      else if (ratio > 0.5) {
	++rcount;
	if (rcount == 3) {
	  bound = min(2.0*bound, 1.0);
	  rcount = 0;
	}
      }
      dprev = dcurr;
    }
  }
  describe2(n, fs, getTime() - t, vpmap);
}

FeatureSet smallFeatures (Polyhedron *a, double d, VPMap *vpmap)
{
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f)
    if ((*f)->getBoundary())
      (*f)->getPC();
  Octree<Face *> *octree = faceOctree(a, d);
  vector<pair<Face *, Face *>> ff;
  octree->pairs(ff, d);  
  delete octree;
  static unsigned int n = 8;
  unsigned int k = ff.size(), m = k/n, is = 0;
  SFData sfd[n];
  for (int i = 0; i < n; ++i) {
    sfd[i].i = i;
    sfd[i].is = is;
    is = i + 1 == n ? k : is + m;
    sfd[i].ie = is;
    sfd[i].ff = &ff;
    sfd[i].d = d;
    sfd[i].vpmap = vpmap;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) smallFeaturesT, (void *) (sfd + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  FeatureSet fs;
  for (int i = 0; i < n; ++i)
    fs.insert(sfd[i].fs.begin(), sfd[i].fs.end());
  return fs;
}

void smallFeaturesT (void *ptr)
{
  SFData *sfd = (SFData *) ptr;
  BaseObject::addThread(sfd->i);
  VPMap *vpmap = sfd->vpmap;
  FeatureSet fs;
  vector<pair<Face *, Face *>>::iterator x = sfd->ff->begin() + sfd->is;
  for (int i = sfd->is; i < sfd->ie; ++i, ++x)
    if (x->first->getBoundary() && x->second->getBoundary() &&
	(!vpmap || moved(x->first, *vpmap) || moved(x->second, *vpmap)))
      smallCandidates(x->first, x->second, sfd->d, fs);
  for (FeatureSet::iterator f = fs.begin(); f != fs.end(); ++f)
    if (small(*f, sfd->d))
      sfd->fs.insert(*f);
}

bool moved (Face *f, const VPMap &vpmap)
{
  Vertices ve = f->getBoundary()->loop();
  for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
    if (vpmap.find(*v) != vpmap.end())
      return true;
  return false;
}

void smallCandidates (Face *f, Face *g, double d, FeatureSet &fs)
{
  HEdge *fb = f->getBoundary(), *gb = g->getBoundary();
  Vertex *vf[] = {fb->tail(), fb->head(), fb->getNext()->head()};
  Vertex *vg[] = {gb->tail(), gb->head(), gb->getNext()->head()};
  HEdge *ef[] = {fb, fb->getNext(), fb->getNext()->getNext()};
  HEdge *eg[] = {gb, gb->getNext(), gb->getNext()->getNext()};
  for (int i = 0; i < 3; ++i)
    if (firstFace(vf[i], f) &&
	bboxOverlap(vf[i]->getBBox(), g->getBBox(), d) &&
	closeVFT(vf[i], g, d))
      fs.insert(Feature(vf[i], g));
  for (int i = 0; i < 3; ++i)
    if (firstFace(vg[i], g) &&
	bboxOverlap(vg[i]->getBBox(), f->getBBox(), d) &&
	closeVFT(vg[i], f, d))
      fs.insert(Feature(vg[i], f));
  for (int i = 0; i < 3; ++i) {
    Edge *e1 = ef[i]->getE();
    if (e1->getHEdge(0) == ef[i])
      for (int j = 0; j < 3; ++j) {
	Edge *e2 = eg[j]->getE();
	if (e1 != e2 && e2->getHEdge(0) == eg[j] &&
	    bboxOverlap(e1->getBBox(), e2->getBBox(), d) &&
	    closeEET(e1, e2, d))
	  fs.insert(Feature(e1, e2));
      }
  }
}

bool firstFace (Vertex *v, Face *f)
{
  Edge *e = v->getEdge(0);
  for (int i = 0; i < e->HEdgesN(); ++i) {
    HEdge *h = e->getHEdge(i);
    if (h->tail() == v && h->getF() == f)
      return true;
  }
  return false;
}

bool closeVVT (Vertex *v, Vertex *w, double d)
{
  PV3 u = v->getP()->getApprox(1.0) - w->getP()->getApprox(1.0);
  Parameter k = d*d - u.dot(u);
  return k.ub() > 0.0;
}

bool closeVET (Vertex *v, Edge *e, double d)
{
  Vertex *et = e->getT(), *eh = e->getH();
  if (v == et || v == eh)
    return false;
  PV3 a = v->getP()->getApprox(1.0), t = et->getP()->getApprox(1.0),
    h = eh->getP()->getApprox(1.0), u = h - t;
  Parameter k = u.dot(a - t), uu = u.dot(u);
  if (k.ub() <= 0.0 || k.lb() >= uu.ub())
    return false;
  PV3 w = uu*(t - a) + k*u;
  Parameter s = uu*uu*d*d - w.dot(w); 
  return s.sign(false) > -1;
}

bool closeVFT (Vertex *v, Face *f, double d)
{
  HEdge *fb = f->getBoundary(), *fbn = fb->getNext(), *fbnn = fbn->getNext();
  Vertex *ve[] = {fb->tail(), fb->head(), fbn->head()};
  if (v == ve[0] || v == ve[1] || v == ve[2])
    return false;
  PV3 p = v->getP()->getApprox(1.0), n = f->getP()->getApprox(1.0).n;
  Parameter k1 = n.dot(p) + f->getP()->getApprox(1.0).k, k2 = d*d*n.dot(n) - k1*k1;
  if (k2.sign(false) == -1)
    return false;
  if (closeVET(v, fb->getE(), d) || closeVET(v, fbn->getE(), d) ||
      closeVET(v, fbnn->getE(), d) || closeVVT(v, ve[0], d) ||
      closeVVT(v, ve[1], d) || closeVVT(v, ve[2], d))
    return true;
  PV3 pts[3];
  for (int i = 0; i < 3; ++i)
    pts[i] = ve[i]->getP()->getApprox(1.0);
  int c = f->getPC();
  for (int i = 0; i < 3; ++i)
    if (cross(pts[(i+1)%3] - pts[i], p - pts[i], c).sign(false) == -1)
      return false;
  return true;
}

bool closeEET (Edge *e, Edge *f, double d)
{
  Vertex *et = e->getT(), *eh = e->getH(), *ft = f->getT(), *fh = f->getH();
  Point *a = et->getP(), *b = ft->getP(), *p = eh->getP(), *q = fh->getP();
  if (a == b || a == q || p == b || p == q)
    return false;
  if (closeVET(et, f, d) || closeVET(eh, f, d) || closeVET(ft, e, d) ||
      closeVET(fh, e, d) || closeVVT(et, ft, d) || closeVVT(et, fh, d) ||
      closeVVT(eh, ft, d) || closeVVT(eh, fh, d))
    return true;
  PV3 aa = a->getApprox(1.0), u = p->getApprox(1.0) - aa,
    bb = b->getApprox(1.0), v = q->getApprox(1.0) - bb, w = bb - aa;
  Parameter uu = u.dot(u), uv = u.dot(v), vv = v.dot(v), uw = u.dot(w),
    vw = v.dot(w), den = uu*vv - uv*uv, s = uw*vv - uv*vw, t = uv*uw - uu*vw;
  double dub = den.ub();
  if (s.ub() < 0.0 || s.lb() > dub || t.ub() < 0.0 || t.lb() > dub)
    return false;
  PV3 n = u.cross(v);
  Parameter wn = w.dot(n), nn = n.dot(n), k = d*d*nn - wn*wn;
  return k.sign(false) > -1;
}

FeatureSet smallFeatures (const FeatureSet &fs, double d)
{
  FeatureSet res;
  for (FeatureSet::const_iterator f = fs.begin(); f != fs.end(); ++f) {
    PTR<Point> p, q;
    Feature g = minimalFeature(*f, p, q);
    if (DistancePP(p, q, d) == 1)
      res.insert(g);
  }
  return res;
}

double separation (const FeatureSet &fs)
{
  double d = 1e10;
  for (FeatureSet::const_iterator f = fs.begin(); f != fs.end(); ++f) {
    PTR<Point> p, q;
    minimalFeature(*f, p, q);
    PV3 u = DiffPoint(p, q).getApprox(1e-8);
    Parameter k = u.dot(u);
    d = min(d, k.ub());
  }
  return sqrt(d);
}

void simplify2s (Polyhedron *a, double d, FeatureSet &fs, VPMap &vpmap,
		 bool vobj, double vbound)
{
  set<Vertex *> vs;
  Expander2 *exp = solveLP(fs, d, vpmap, vobj, vbound, vs);
  moveVertices(a, d, vpmap, vs, exp);
  delete exp;
  fs = smallFeatures(a, 3.465*d, fs, &vpmap);
}

Expander2 * solveLP (const FeatureSet &fs, double d, const VPMap &vpmap,
		     bool vobj, double vbound, set<Vertex *> &vs)
{
  Expander2 *exp = new Expander2(d);
  setupLP(fs, d, vs, exp);
  for (set<Vertex *>::iterator v = vs.begin(); v != vs.end(); ++v) {
    double dv[3];
    getDV(*v, d, vpmap, dv);
    exp->addDisplacement(*v, Expander2::Point(dv));
  }
  Parameter::disable();
  exp->expandV(1.0, vobj, vbound);
  Parameter::enable();
  return exp;
}

void setupLP (const FeatureSet &fs, double d, set<Vertex *> &vs, Expander2 *exp)
{
  for (FeatureSet::const_iterator f = fs.begin(); f != fs.end(); ++f) {
    Vertices ve[2];
    if (f->v) {
      ve[0].push_back(f->v);
      ve[1] = f->fa->getBoundary()->loop();
    }
    else {
      ve[0].push_back(f->e->getT());
      ve[0].push_back(f->e->getH());
      ve[1].push_back(f->f->getT());
      ve[1].push_back(f->f->getH());
    }
    Separator s(*f, d);
    SeparatorData data = s.getApprox(1e-7);
    PV3 uu = data.u.unit(), vu = data.v.unit(), wu = data.w.unit();
    Expander2::Point u(uu.x.mid(), uu.y.mid(), uu.z.mid());
    Expander2::Point v(vu.x.mid(), vu.y.mid(), vu.z.mid());
    Expander2::Point w(wu.x.mid(), wu.y.mid(), wu.z.mid());
    Expander2::Pair pair(u, v, w, data.pq.sqrt().mid()/d);
    for (int ab = 0; ab < 2; ++ab)
      for (Vertices::iterator x = ve[ab].begin(); x != ve[ab].end(); ++x) {
	double vwr[3];
	s.getVWR((*x)->getP(), ab == 0, vwr);
	pair.addConstraint(ab, Expander2::Constraint(*x, vwr[0], vwr[1], vwr[2]));
	vs.insert(*x);
      }
    exp->addPair(pair);
  }
}

Feature minimalFeature (const Feature &f, PTR<Point> &p, PTR<Point> &q)
{
  return f.type == VF ? minimalFeature(f.v, f.fa, p, q)
    : minimalFeature(f.e, f.f, p, q);
}

Feature minimalFeature (Vertex *v, Face *f, PTR<Point> &p, PTR<Point> &q)
{
  p = v->getP();
  if (p->side(f->getP()) != 0) {
    q = new PlanePoint(f->getP(), p);
    if (f->contains(q, true))
      return Feature(v, f);
  }
  HEdges ed = f->getBoundary()->edgeLoop();
  q = 0;
  Edge *e = 0;
  for (int i = 0; i < 3; ++i) {
    Edge *ei = ed[i]->getE();
    Point *t = ei->getT()->getP(), *h = ei->getH()->getP();
    if (onEdge(p, t, h, true)) {
      PTR<Point> qi = new LinePoint(p, t, h);
      if (!q || closerPair(p, qi, p, q)) {
	q = qi;
	e = ei;
      }
    }
  }
  Vertex *w = 0;
  for (int i = 0; i < 3; ++i)
    if (!q || closerPair(p, ed[i]->tail()->getP(), p, q)) {
      w = ed[i]->tail();
      e = 0;
    }
  if (e)
    return Feature(v, e);
  q = w->getP();
  return Feature(v, w);
}

Feature minimalFeature (Edge *e, Edge *f, PTR<Point> &p, PTR<Point> &q)
{
  Vertex *vs[] = {e->getT(), e->getH(), f->getT(), f->getH()};
  Point *et = vs[0]->getP(), *eh = vs[1]->getP(), *ft = vs[2]->getP(),
    *fh = vs[3]->getP();
  if (Orientation(et, eh, ft, fh) != 0) {
    p = new LineLinePoint(et, eh, ft, fh);
    q = new LineLinePoint(ft, fh, et, eh);
    if (onEdge(p, et, eh, true) && onEdge(q, ft, fh, true))
      return Feature(e, f);
  }
  Edge *es[] = {f, f, e, e};
  PTR<Point> r = 0;
  int fi = -1;
  for (int i = 0; i < 4; ++i) {
    Point *a = vs[i]->getP(), *t = es[i]->getT()->getP(),
      *h = es[i]->getH()->getP();
    if (onEdge(a, t, h, true)) {
      PTR<Point> ri = new LinePoint(a, t, h);
      if (!r || closerPair(a, ri, vs[fi]->getP(), r)) {
	r = ri;
	fi = i;
      }
    }
  }
  Vertex *v = 0, *w = 0;
  for (int i = 0; i < 2; ++i)
    for (int j = 2; j < 4; ++j)
      if (!v || closerPair(vs[i]->getP(), vs[j]->getP(), v->getP(), w->getP())) {
	v = vs[i];
	w = vs[j];
      }
  if (r && closerPair(vs[fi]->getP(), r, v->getP(), w->getP())) {
    p = vs[fi]->getP();
    q = r;
    return Feature(vs[fi], es[fi], fi < 2);
  }
  p = v->getP();
  q = w->getP();
  return Feature(v, w, v == vs[0] || v == vs[1]);
}

bool small (const Feature &f, double d)
{
  PTR<Point> p, q;
  minimalFeature(f, p, q);
  return DistancePP(p, q, d) == 1;
}

PV3 orthogonal (const PV3 &u)
{
  double ux = u.x.sign(false) ? u.x.abs().mid() : min(- u.x.lb(), u.x.ub()),
    uy = u.y.sign(false) ? u.y.abs().mid() : min(- u.y.lb(), u.y.ub()),
    uz = u.z.sign(false) ? u.z.abs().mid() : min(- u.z.lb(), u.z.ub());
  if (ux <= uy && ux <= uz)
    return PV3(Parameter::constant(0.0), - u.z, u.y);
  if (uy <= ux && uy <= uz)
    return PV3(- u.z, Parameter::constant(0.0), u.x);
  return PV3(- u.y, u.x, Parameter::constant(0.0));
}

void getDV (Vertex *v, double d, const VPMap &vpmap, double *dv)
{
  VPMap::const_iterator i = vpmap.find(v);
  if (i == vpmap.end()) {
    dv[0] = dv[1] = dv[2] = 0.0;
    return;
  }
  PV3 u = (v->getP()->getApprox(1e-7) - i->second->getApprox(1e-7))/d;
  for (int i = 0; i < 3; ++i)
    dv[i] = u[i].mid();
}

void moveVertices (Polyhedron *a, double d, VPMap &vpmap, const set<Vertex *> &vs,
		   Expander2 *exp)
{
  for (set<Vertex *>::const_iterator v = vs.begin(); v != vs.end(); ++v) {
    Expander2::Point dv = exp->motion(*v);
    if (dv.x[0] != 0.0 || dv.x[1] != 0.0 || dv.x[2] != 0.0) {
      PTR<Point> vp = (*v)->getP(),
	dp = new Point(d*dv.x[0], d*dv.x[1], d*dv.x[2], false),
	q = new SumPoint(vp, dp);
      vpmap.insert(VPPair(*v, vp));
      a->moveVertex(*v, q);
    }
  }
}

FeatureSet smallFeatures (Polyhedron *a, double d, const FeatureSet &fs,
			  VPMap *vpmap)
{
  FeatureSet nfs = smallFeatures(a, d, vpmap);
  for (FeatureSet::iterator f = fs.begin(); f != fs.end(); ++f)
    if (small(*f, d))
      nfs.insert(*f);
  return nfs;  
}

VPMap savePos (const FeatureSet &fs)
{
  set<Vertex *> vs;
  for (FeatureSet::const_iterator f = fs.begin(); f != fs.end(); ++f)
    f->vertices(vs);
  VPMap vpmap;
  for (set<Vertex *>::iterator v = vs.begin(); v != vs.end(); ++v)
    vpmap.insert(VPPair(*v, (*v)->getP()));
  return vpmap;
}

void restorePos (Polyhedron *a, const VPMap &vpmap)
{
  for (VPMap::const_iterator i = vpmap.begin(); i != vpmap.end(); ++i)
    a->moveVertex(i->first, i->second);
}

void describe2 (int n, const FeatureSet &fs, double t, const VPMap &vpmap)
{
  int nv = vpmap.size(), vv = 0, ve = 0, vf = 0, ee = 0;
  for (FeatureSet::const_iterator f = fs.begin(); f != fs.end(); ++f)
    switch (f->type) {
    case VV: ++vv; break;
    case VE: ++ve; break;
    case VF: ++vf; break;
    case EE: ++ee; break;
    }
  double dmax, dr = displacement(vpmap, dmax)/vpmap.size();
  cerr << "simplify 2: time = " << t;
  if (n) {
    cerr << " n = " << n;
    if (vv)
      cerr << " vv = " << vv;
    if (ve)
      cerr << " ve = " << ve;
    if (vf)
      cerr << " vf = " << vf;
    if (ee)
      cerr << " ee = " << ee;
    if (nv)
    cerr << endl << "            moved = " << nv << " dav = " << dr
	 << " dmax = " << dmax;
  }
  cerr << endl;
}

bool intersects (Expander2 *e, Expander2::Pair *p)
{
  Points pts, dsp;
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < p->constraints[i].size(); ++j) {
      Expander2::Constraint c = p->constraints[i][j];
      pts.push_back(c.i->getP());
      Expander2::Point m = e->motion(c.i);
      PTR<Point> q = new Point(e->d*m.x[0], e->d*m.x[1], e->d*m.x[2], false);
      dsp.push_back(q);
    }
  bool enabled = Parameter::isEnabled();
  Parameter::enable();
  bool res = intersects(pts, dsp, p->constraints[0].size() == 1);
  if (!enabled)
    Parameter::disable();
  return res;
}

bool intersects (const Points &pts, const Points &dsp, bool vf)
{
  PTR<Object<PPoly>> p = new FIPoly(pts, dsp);
  Roots rts = PolySolver(p).getRoots(0, 1);
  for (Roots::iterator r = rts.begin(); r != rts.end(); ++r) {
    Points rpts;
    for (int i = 0; i < 4; ++i)
      rpts.push_back(new DisplacedPoint(pts[i], dsp[i], *r));
    if (vf ? intersectsVF(rpts) : intersectsEE(rpts))
      return true;
  }
  return false;
}

bool intersectsVF (const Points &pts)
{
  Face f(pts[1], pts[2], pts[3]);
  return f.contains(pts[0], false);
}

bool intersectsEE (const Points &pts)
{
  EEPoint p(pts[0], pts[1], pts[2], pts[3], 3);
  return onEdge(&p, pts[0], pts[1], false) &&
    onEdge(&p, pts[2], pts[3], false);
}

Polyhedron * round (Polyhedron *a)
{
  Polyhedron *b = new Polyhedron;
  PVMap pvmap;
  for (Vertices::const_iterator v = a->vertices.begin(); v != a->vertices.end(); ++v) {
    Point *p = (*v)->getP();
    PV3 q = p->getApprox();
    pvmap.insert(PVPair(p, b->getVertex(new Point(q.x.lb(), q.y.lb(), q.z.lb()))));
  }
  for (Faces::const_iterator f = a->faces.begin(); f != a->faces.end(); ++f)
    b->addTriangle(*f, pvmap);
  return b;
}
