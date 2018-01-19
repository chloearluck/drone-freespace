#include "simplify.h"

void simplify (Polyhedron *a, double d, bool opt2)
{
  simplify1(a, d);
  simplify2(a, d, opt2);
  a->updateCells();
  a->removeNullFaces();
}

void simplify1 (Polyhedron *a, double d)
{
  double t = getTime();
  FOctree *octree = faceOctree(a);
  VPMap vpmap;
  int n1 = 0, n2 = 0;
  while (!simplify1(a, d, octree, vpmap, n1, n2))
    ;
  delete octree;
  t = getTime() - t;
  describe1(t, n1, n2, vpmap);
}  

FOctree * faceOctree (Polyhedron *a, double s)
{
  Faces fa;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f)
    if (!(*f)->getBoundary().empty())
      fa.push_back(*f);
  return Octree<Face *>::octree(fa, a->bbox, s);
}

bool simplify1 (Polyhedron *a, double d, FOctree *octree, VPMap &vpmap, int &n1, int &n2)
{
  bool bad = false;
  int n = n1 + n2;
  IFeatureQ fq;
  for (Edges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    addCollapseFlips(*e, d, fq);
  while (!fq.empty()) {
    IFeature f = fq.top();
    fq.pop();
    if (f.valid())
      if (!f.v) {
	if (collapse(a, d, f.e, fq, octree, vpmap, bad))
	  ++n1;
      }
      else if (flip(a, d, f.getH(), fq, octree, vpmap, bad))
	++n2;
  }
  return n1 + n2 == n || !bad;
}

bool closeVE (Vertex *v, Vertex *w, Vertex *x, double d)
{
  Point *a = v->getP(), *t = w->getP(), *h = x->getP();
  return DistancePL1(a, t, h, d) == 1 && onEdge(a, t, h, true);
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

bool collapse (Polyhedron *a, double d, Edge *e, IFeatureQ &fq, FOctree *octree,
	       VPMap &vpmap, bool &bad)
{
  if (!collapsible(e))
    return false;
  HEdge *e1 = e->getHEdge(0), *e2 = e->getHEdge(1);
  Vertex *t = e1->tail(), *h = e1->head(), *v = e1->getNext()->head(),
    *w = e2->getNext()->head();
  PTR<Point> tp = t->getP(), hp = h->getP();
  PTR<Point> cp = new MidPoint(tp, hp);
  a->removeLoop(e1);
  a->removeLoop(e2);
  vpmap.insert(VPPair(h, h->getP()));
  a->moveVertex(h, cp);
  Faces ft = t->incidentFaces();
  for (Faces::iterator f = ft.begin(); f != ft.end(); ++f)
    a->replaceVertex(*f, t, h);
  double dc = distanceUB(tp, hp) + 1.733*Parameter::delta;
  if (badCollapse(h, octree, dc)) {
    bad = true;
    for (Faces::iterator f = ft.begin(); f != ft.end(); ++f)
      a->replaceVertex(*f, h, t);
    a->moveVertex(h, hp);
    octree->insert(a->addTriangle(t, h, v));
    octree->insert(a->addTriangle(h, t, w));
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
  HEdges et, eh;
  t->outgoingHEdges(et);
  h->outgoingHEdges(eh);
  VertexSet vs;
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

double distanceUB (Point *t, Point *h)
{
  PV3 u = h->getApprox(1.0) - t->getApprox(1.0);
  Parameter k = u.dot(u);
  return 0.5*sqrt(k.ub());
}

bool badCollapse (Vertex *v, FOctree *octree, double dc)
{
  Faces fv = v->incidentFaces();
  for (Faces::iterator f = fv.begin(); f != fv.end(); ++f) {
    bool ef[3];
    newEdges(*f, v, ef);
    Faces fa;
    octree->find((*f)->getBBox(), fa, dc);
    for (Faces::iterator g = fa.begin(); g != fa.end(); ++g)
      if (!(*g)->getBoundary().empty() &&
	  find(fv.begin(), fv.end(), *g) == fv.end() &&
	  intersects(*f, ef, *g))
	return true;
  }
  return false;
}

void newEdges (Face *f, Vertex *v, bool *ef)
{
  HEdges ed;
  f->boundaryHEdges(ed);
  for (int i = 0; i < 3; ++i)
    ef[i] = ed[i]->tail() == v || ed[i]->head() == v;
}

bool intersects (Face *f, bool *ef, Face *g)
{
  if (f->sharedBoundaryEdge(g))
    return false;
  bool eg[3] = {false, false, false};
  if (f->coplanar(g))
    return intersectsFFP(f, ef, g, eg) || intersectsFFP(g, eg, f, ef);
  Edges iedges;
  Faces ifaces;
  bool pflag = false, flag = intersectsPE(f, ef, g, eg, iedges, ifaces, pflag);
  if (pflag)
    return true;
  if (!flag)
    return false;
  flag = intersectsPE(g, eg, f, ef, iedges, ifaces, pflag);
  if (pflag)
    return true;
  if (!flag)
    return false;
  for (int i = 0; i < iedges.size(); ++i)
    if (ifaces[i]->intersects(iedges[i]))
      return true;
  return false;
}

bool intersectsFFP (Face *f, bool *ef, Face *g, bool *eg)
{
  HEdges ed;
  g->boundaryHEdges(ed);
  for (int i = 0; i < 3; ++i)
    if (intersectsFEP(f, ef, ed[i]->getE(), eg[i]))
      return true;
  return false;
}

bool intersectsPE (Face *f, bool *ef, Face *g, bool *eg, Edges &iedges,
		   Faces &ifaces, bool &pflag)
{
  int s[3];
  HEdges ed;
  g->boundaryHEdges(ed);
  for (int i = 0; i < 3; ++i)
    s[i] = ed[i]->tail()->getP()->side(f->getP());
  bool sp = false, sm = false, ff = ef[0] || ef[1] || ef[2];
  for (int i = 0; i < 3; ++i) {
    int si = s[i], sj = s[(i+1)%3];
    if (si == 0 && sj == 0) {
      if (intersectsFEP(f, ef, ed[i]->getE(), eg[i])) {
	pflag = true;
	return true;
      }
    }
    if (si == 0) {
      Vertex *vi = ed[i]->tail();
      if (!f->boundaryVertex(vi)) {
	pflag = f->containsConvex(vi->getP(), false);
	if (pflag)
	  return true;
      }
    }
    else {
      if (si == 1)
	sp = true;
      else
	sm = true;
      if (si*sj == -1 && (ff || eg[i])) {
	ifaces.push_back(f);
	iedges.push_back(ed[i]->getE());
      }
    }
  }
  return sp && sm;
}

bool intersectsFEP (Face *f, bool *ef, Edge *e, bool ee)
{
  bool ff = ef[0] || ef[1] || ef[2];
  if (!((ff || ee) && bboxOverlap(f->getBBox(), e->getBBox())))
    return false;
  HEdges ed;
  f->boundaryHEdges(ed);
  for (int i = 0; i < 3; ++i)
    if ((ef[i] || ee) && f->intersectsEE(ed[i]->getE(), e))
      return true;
  return false;
}

void addStar (Vertex *v, double d, IFeatureQ &fq)
{
  HEdges vout;
  v->outgoingHEdges(vout);
  for (HEdges::iterator e = vout.begin(); e != vout.end(); ++e) {
    if (DistancePP((*e)->tail()->getP(), (*e)->head()->getP(), d) == 1)
      fq.push(IFeature((*e)->getE()));
    HEdges ev;
    (*e)->loop(ev);
    for (HEdges::iterator h = ev.begin(); h != ev.end(); ++h)
      if (closeVE((*h)->getNext()->head(), (*h)->tail(), (*h)->head(), d))
	fq.push(IFeature(*h));
  }
}

bool flip (Polyhedron *a, double d, HEdge *e, IFeatureQ &fq, FOctree *octree,
	   VPMap &vpmap, bool &bad)
{
  Vertex *t = e->tail(), *h = e->head(), *v = e->getNext()->head();
  if (CloserPair(t->getP(), h->getP(), v->getP(), t->getP()) == 1)
    return collapse(a, d, e->getNext()->getE(), fq, octree, vpmap, bad);
  if (CloserPair(t->getP(), h->getP(), h->getP(), v->getP()) == 1)
    return collapse(a, d, e->getNext()->getNext()->getE(), fq, octree,
		    vpmap, bad);
  if (!flippable(e, d))
    return false;
  HEdge *vw = flip(a, e);
  octree->insert(vw->getF());
  octree->insert(vw->ccw()->getF());
  if (badFlip(vw, octree)) {
    bad = true;
    vw = flip(a, vw);
    octree->insert(vw->getF());
    octree->insert(vw->ccw()->getF());
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
  return CloserPair(v->getP(), w->getP(), t->getP(), h->getP()) == 1;
}

HEdge * flip (Polyhedron *a, HEdge *e)
{
  Vertex *t = e->tail(), *h = e->head(), *v = e->getNext()->head(),
    *w = e->ccw()->getNext()->head();
  a->removeLoop(e->ccw());
  a->removeLoop(e);
  Vertex *vwh[] = {v, w, h}, *wvt[] = {w, v, t};
  HEdge *vw = a->addLoop(vwh, 3), *wv = a->addLoop(wvt, 3);
  a->addFace(vw);
  a->addFace(wv);
  return vw;
}

bool badFlip (HEdge *e, FOctree *octree)
{
  for (int i = 0; i < 2; ++i, e = e->ccw()) {
    Face *f = e->getF();
    bool ef[3];
    newEdges(f, e->getE(), ef);
    Faces fa;
    octree->find(f->getBBox(), fa);
    for (Faces::iterator g = fa.begin(); g != fa.end(); ++g)
      if (!(*g)->getBoundary().empty() &&
	  intersects(f, ef, *g))
	return true;
  }
  return false;
}

void newEdges (Face *f, Edge *e, bool *ef)
{
  HEdges ed;
  f->boundaryHEdges(ed);
  for (int i = 0; i < 3; ++i)
    ef[i] = ed[i]->getE() == e;
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
  while (separation(fslp) < 0.99999*d) {
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
      if (dbest >= dprev)
        break;
      while (separation(fslp) < 0.99999*d) {
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

bool closeVF (Vertex *v, Face *f, double d)
{
  Point *p = v->getP();
  Plane *pl = f->getP();
  PTR<Point> q = new PlanePoint(pl, p);
  return DistancePP(p, q, d) == 1 && f->contains(q, false);
}

bool closeEE (Edge *e, Edge *f, double d)
{
  Point *et = e->getT()->getP(), *eh = e->getH()->getP(),
    *ft = f->getT()->getP(), *fh = f->getH()->getP();
  if (Orientation(et, eh, ft, fh) == 0)
    return false;
  PTR<Point> p = new LineLinePoint(et, eh, ft, fh),
    q = new LineLinePoint(ft, fh, et, eh);
  return DistancePP(p, q, d) == 1 && onEdge(p, et, eh, true) &&
    onEdge(q, ft, fh, true);
}

FeatureSet smallFeatures (Polyhedron *a, double d, VPMap *vpmap)
{
  FOctree *octree = faceOctree(a, d);
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f)
    if (!(*f)->getBoundary().empty()) {
      (*f)->getPC();
      (*f)->getP()->getApprox(1e-6);
    }
  int m = 16;
  vector<FOctree *> nodes;
  octree->expand(m, nodes);
  int n = nodes.size();
  TData *td = new TData [n];
  for (int i = 0; i < n; ++i) {
    td[i].octree = nodes[i];
    td[i].d = d;
    td[i].vpmap = vpmap;
  }
  pthread_t *threads = new pthread_t [n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) smallCandidatesT, (void *) (td + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  delete [] threads;
  FeatureSet fs;
  for (int i = 0; i < n; ++i)
    for (FeatureSet::iterator f = td[i].fs.begin(); f != td[i].fs.end(); ++f)
      if (f->small(d))
	fs.insert(*f);
  delete [] td;
  delete octree;
  return fs;
}

void smallCandidatesT (void *ptr)
{
  TData *td = (TData *) ptr;
  Faces fa, fb;
  td->octree->pairs(fa, fb, td->d);
  VPMap *vpmap = td->vpmap;
  for (int i = 0; i < fa.size(); ++i)
    if (!fa[i]->getBoundary().empty() && !fb[i]->getBoundary().empty() &&
	(!vpmap || moved(fa[i], *vpmap) || moved(fb[i], *vpmap)))
      smallCandidates(fa[i], fb[i], td->d, td->fs);
}

bool moved (Face *f, const VPMap &vpmap)
{
  Vertices ve;
  f->boundaryVertices(ve);
  for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
    if (vpmap.find(*v) != vpmap.end())
      return true;
  return false;
}

void smallCandidates (Face *f, Face *g, double d, FeatureSet &fs)
{
  HEdge *fb = f->getBoundary(0), *gb = g->getBoundary(0);
  Vertex *vf[] = {fb->tail(), fb->head(), fb->getNext()->head()};
  Vertex *vg[] = {gb->tail(), gb->head(), gb->getNext()->head()};
  HEdge *ef[] = {fb, fb->getNext(), fb->getNext()->getNext()};
  HEdge *eg[] = {gb, gb->getNext(), gb->getNext()->getNext()};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j)
      if (vf[i] != vg[j] &&
      	  bboxOverlap(vf[i]->getBBox(), vg[j]->getBBox(), d) &&
	  closeVVT(vf[i], vg[j], d))
	fs.insert(Feature(vf[i], vg[j]));
    bool flag = firstFace(vf[i], f);
    for (int k = 0; k < 3; ++k) {
      Edge *ee = eg[k]->getE();
      if ((flag || eg[k]->getNext()->head() == vf[i]) &&
	  ee->getHEdge(0) == eg[k] &&
	  bboxOverlap(vf[i]->getBBox(), ee->getBBox(), d) &&
	  closeVET(vf[i], ee, d))
	fs.insert(Feature(vf[i], ee));
    }
    if (flag && bboxOverlap(vf[i]->getBBox(), g->getBBox(), d) &&
	closeVFT(vf[i], g, d))
      fs.insert(Feature(vf[i], g));
  }
  for (int i = 0; i < 3; ++i) {
    bool flag = firstFace(vg[i], g);
    for (int j = 0; j < 3; ++j) {
      Edge *ee = ef[j]->getE();
      if ((flag || ef[j]->getNext()->head() == vg[i]) &&
	  ee->getHEdge(0) == ef[j] &&
	  bboxOverlap(vg[i]->getBBox(), ee->getBBox(), d) &&
	  closeVET(vg[i], ee, d))
	fs.insert(Feature(vg[i], ee));
    }
    if (flag && bboxOverlap(vg[i]->getBBox(), f->getBBox(), d) &&
	closeVFT(vg[i], f, d))
      fs.insert(Feature(vg[i], f));
  }
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
  if (v == e->getT() || v == e->getH())
    return false;
  PV3 a = v->getP()->getApprox(1.0), t = e->getT()->getP()->getApprox(1.0),
    h = e->getH()->getP()->getApprox(1.0), u = h - t;
  Parameter k = u.dot(a - t), uu = u.dot(u);
  if (k.ub() <= 0.0 || k.lb() >= uu.ub())
    return false;
  PV3 w = uu*(t - a) + k*u;
  Parameter s = uu*uu*d*d - w.dot(w); 
  return s.sign(false) > -1;
}

bool closeVFT (Vertex *v, Face *f, double d)
{
  HEdge *fb = f->getBoundary(0);
  Vertex *ve[] = {fb->tail(), fb->head(), fb->getNext()->head()};
  if (v == ve[0] || v == ve[1] || v == ve[2])
    return false;
  PV3 p = v->getP()->getApprox(1.0), n = f->getP()->getApprox(1.0).n;
  Parameter k1 = n.dot(p) + f->getP()->getApprox(1.0).k, k2 = d*d*n.dot(n) - k1*k1;
  if (k2.sign(false) == -1)
    return false;
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
  Point *a = e->getT()->getP(), *b = f->getT()->getP(), 
    *p = e->getH()->getP(), *q = f->getH()->getP();
  if (a == b || a == q || p == b || p == q)
    return false;
  PV3 aa = a->getApprox(1.0), u = p->getApprox(1.0) - aa,
    bb = b->getApprox(1.0), v = q->getApprox(1.0) - bb, w = bb - aa;
  Parameter uu = u.dot(u), uv = u.dot(v), vv = v.dot(v), uw = u.dot(w),
    vw = v.dot(w), den = uu*vv - uv*uv, s = uw*vv - uv*vw, t = uv*uw - uu*vw;
  double dub = den.ub();
  if (s.ub() < 0.0 || s.lb() > dub || t.ub() < 0.0 || t.lb() > dub)
    return false;
  PV3 x = den*w + t*v - s*u;
  Parameter k = d*d*den*den - x.dot(x);
  return k.sign(false) > -1;
}

FeatureSet smallFeatures (const FeatureSet &fs, double d)
{
  FeatureSet res;
  for (FeatureSet::const_iterator f = fs.begin(); f != fs.end(); ++f)
    if (f->small(d))
      res.insert(*f);
  return res;
}

PV3 orthogonal (const PV3 &u)
{
  Parameter ux = u.x.abs(), uy = u.y.abs(), uz = u.z.abs();
  if ((uy - ux).sign() > -1 && (uz - ux).sign() > -1)
    return PV3(Parameter(0.0), - u.z, u.y);
  if ((ux - uy).sign() > -1 && (uz - uy).sign() > -1)
    return PV3(- u.z, Parameter(0.0), u.x);
  return PV3(- u.y, u.x, Parameter(0.0));
}

double separation (const FeatureSet &fs)
{
  double d = 1e10;
  for (FeatureSet::const_iterator f = fs.begin(); f != fs.end(); ++f) {
    PTR<Point> p, q;
    f->closestPoints(p, q);
    PV3 u = p->getApprox(1e-8) - q->getApprox(1e-8);
    Parameter k = u.dot(u);
    d = min(d, k.ub());
  }
  return sqrt(d);
}

void simplify2s (Polyhedron *a, double d, FeatureSet &fs, VPMap &vpmap,
		 bool vobj, double vbound)
{
  VIMap vimap;
  double *v = solveLP(fs, d, vpmap, vobj, vbound, vimap);
  moveVertices(a, vpmap, vimap, v);
  delete [] v;
  fs = smallFeatures(a, 3.465*d, fs, &vpmap);
}

double * solveLP (const FeatureSet &fs, double d, const VPMap &vpmap,
		  bool vobj, double vbound, VIMap &vimap)
{
  Expander2 exp;
  setupLP(fs, d, vimap, exp);
  for (VIMap::iterator i = vimap.begin(); i != vimap.end(); ++i) {
    double dv[3];
    getDV(i->first, d, vpmap, dv);
    exp.addDisplacement(i->second, Expander2::Point(dv));
  }
  exp.expandV(1.0, vobj, vbound);
  double *x = new double [3*vimap.size()];
  for (VIMap::iterator i = vimap.begin(); i != vimap.end(); ++i) {
    int k = i->second;
    Expander2::Point p = exp.motion(k);
    for (int j = 0; j < 3; ++j)
      x[3*k+j] = d*p.x[j];
  }
  return x;
}

void setupLP (const FeatureSet &fs, double d, VIMap &vimap, Expander2 &exp)
{
  FeatureSet fsnm;
  for (FeatureSet::const_iterator f = fs.begin(); f != fs.end(); ++f)
    if (f->type == VV)
      featuresLPVV(f->v, f->w, fsnm);
    else if (f->type == VE)
      featuresLPVE(f->v, f->e, fsnm);
    else
      setupLP(*f, d, vimap, exp, true);
  for (FeatureSet::const_iterator f = fsnm.begin(); f != fsnm.end(); ++f)
    setupLP(*f, d, vimap, exp, false);
}

int InnerProduct4::sign ()
{
  return (a->getP() - b->getP()).dot(c->getP() - d->getP()).sign();
}

void featuresLPVV (Vertex *v, Vertex *w, FeatureSet &fs)
{
  Point *p = v->getP(), *q = w->getP();
  HEdges hv, hw;
  Edges edv, edw;
  v->outgoingHEdges(hv);
  w->outgoingHEdges(hw);
  for (HEdges::iterator h = hv.begin(); h != hv.end(); ++h)
    if (InnerProduct4((*h)->head()->getP(), p, q, p) == -1)
      edv.push_back((*h)->getE());
  for (HEdges::iterator h = hw.begin(); h != hw.end(); ++h)
    if (InnerProduct4((*h)->head()->getP(), q, p, q) == -1)
      edw.push_back((*h)->getE());
  for (Edges::iterator e = edv.begin(); e != edv.end(); ++e)
    for (Edges::iterator f = edw.begin(); f != edw.end(); ++f)
      fs.insert(Feature(*e, *f));
  for (HEdges::iterator h = hv.begin(); h != hv.end(); ++h)
    if (InnerProduct4((*h)->head()->getP(), p, q, p) == -1 &&
	InnerProduct4((*h)->getNext()->head()->getP(), p, q, p) == -1)
      fs.insert(Feature(w, (*h)->getF()));
  for (HEdges::iterator h = hw.begin(); h != hw.end(); ++h)
    if (InnerProduct4((*h)->head()->getP(), q, p, q) == -1 &&
	InnerProduct4((*h)->getNext()->head()->getP(), q, p, q) == -1)
      fs.insert(Feature(v, (*h)->getF()));
}

void featuresLPVE (Vertex *v, Edge *e, FeatureSet &fs)
{
  PTR<Point> p = v->getP(),
    q = new LinePoint(p, e->getT()->getP(), e->getH()->getP());
  HEdges hv;
  v->outgoingHEdges(hv);
  for (HEdges::iterator h = hv.begin(); h != hv.end(); ++h)
    if (InnerProduct4((*h)->head()->getP(), p, q, p) == -1)
      fs.insert(Feature((*h)->getE(), e));
  for (int i = 0; i < 2; ++i) {
    HEdge *h = e->getHEdge(i);
    if (InnerProduct4(h->getNext()->head()->getP(), q, p, q) == -1)
      fs.insert(Feature(v, h->getF()));
  }
}

void setupLP (const Feature &f, double d, VIMap &vimap, Expander2 &exp,
	      bool minimal)
{
  Vertices ve[2];
  if (f.type == VF) {
    ve[0].push_back(f.v);
    f.fa->boundaryVertices(ve[1]);
  }
  else {
    ve[0].push_back(f.e->getT());
    ve[0].push_back(f.e->getH());
    ve[1].push_back(f.f->getT());
    ve[1].push_back(f.f->getH());
  }
  Separator s(minimal ? f : minimalFeature(f), d);
  if (!check(&s, ve[0], ve[1]))
    cerr << "bad feature" << endl;
  SeparatorData data = s.getApprox(1e-7);
#ifdef UNIT_U
  Expander2::Point u(data.u.x.mid(), data.u.y.mid(), data.u.z.mid());
  Expander2::Point v(data.v.x.mid(), data.v.y.mid(), data.v.z.mid());
  Expander2::Point w(data.w.x.mid(), data.w.y.mid(), data.w.z.mid());
  Expander2::Pair pair(u, v, w, data.pq.mid()/d);
#else
  PV3 uu = data.u.unit();
  PV3 vu = data.v.unit();
  PV3 wu = data.w.unit();
  Expander2::Point u(uu.x.mid(), uu.y.mid(), uu.z.mid());
  Expander2::Point v(vu.x.mid(), vu.y.mid(), vu.z.mid());
  Expander2::Point w(wu.x.mid(), wu.y.mid(), wu.z.mid());
  Expander2::Pair pair(u, v, w, data.pq.sqrt().mid()/d);
#endif
  for (int ab = 0; ab < 2; ++ab)
    for (Vertices::iterator x = ve[ab].begin(); x != ve[ab].end(); ++x) {
      int xi = getVertex(*x, vimap);
      double vwr[3];
      s.getVWR((*x)->getP(), ab == 0, vwr);
      pair.addConstraint(ab, Expander2::Constraint(xi, vwr[0], vwr[1], vwr[2]));
    }
  exp.addPair(pair);
}

Feature minimalFeature (const Feature &f)
{
  return f.type == VF ? minimalFeature(f.v, f.fa)
    : minimalFeature(f.e, f.f);
}

Feature minimalFeature (Vertex *v, Face *f)
{
  Point *a = v->getP();
  if (a->side(f->getP()) != 0) {
    PlanePoint pp(f->getP(), a);
    if (f->containsConvex(&pp, false))
      return Feature(v, f);
  }
  HEdges ed;
  f->boundaryHEdges(ed);
  PTR<Point> p = 0;
  Edge *e = 0;
  for (int i = 0; i < 3; ++i) {
    Edge *ei = ed[i]->getE();
    Point *t = ei->getT()->getP(), *h = ei->getH()->getP();
    if (onEdge(a, t, h, true)) {
      PTR<Point> pi = new LinePoint(a, t, h);
      if (!e || CloserPair(a, pi, a, p) == 1) {
	e = ei;
	p = pi;
      }
    }
  }
  Vertex *w = ed[0]->tail();
  for (int i = 1; i < 3; ++i)
    if (CloserPair(a, ed[i]->tail()->getP(), a, w->getP()) == 1)
      w = ed[i]->tail();
  if (p && CloserPair(a, p, a, w->getP()) == 1)
    return Feature(v, e);
  return Feature(v, w, v < w);
}

Feature minimalFeature (Edge *e, Edge *f)
{
  Vertex *vs[] = {e->getT(), e->getH(), f->getT(), f->getH()};
  Point *et = vs[0]->getP(), *eh = vs[1]->getP(), *ft = vs[2]->getP(),
    *fh = vs[3]->getP();
  if (Orientation(et, eh, ft, fh) != 0) {
    LineLinePoint pe(et, eh, ft, fh), pf(ft, fh, et, eh);
    if (onEdge(&pe, et, eh, true) && onEdge(&pf, ft, fh, true))
      return Feature(e, f);
  }
  Edge *es[] = {f, f, e, e};
  PTR<Point> p = 0;
  int fi = -1;
  for (int i = 0; i < 4; ++i) {
    Point *a = vs[i]->getP(), *t = es[i]->getT()->getP(),
      *h = es[i]->getH()->getP();
    if (onEdge(a, t, h, true)) {
      PTR<Point> pi = new LinePoint(a, t, h);
      if (!p || CloserPair(a, pi, vs[fi]->getP(), p) == 1) {
	p = pi;
	fi = i;
      }
    }
  }
  Vertex *v = vs[0], *w = vs[2];
  if (CloserPair(vs[0]->getP(), vs[3]->getP(), v->getP(), w->getP()) == 1)
    w = vs[3];
  if (CloserPair(vs[1]->getP(), vs[2]->getP(), v->getP(), w->getP()) == 1) {
    v = vs[1];
    w = vs[2];
  }
  if (CloserPair(vs[1]->getP(), vs[3]->getP(), v->getP(), w->getP()) == 1) {
    v = vs[1];
    w = vs[3];
  }
  if (p && CloserPair(vs[fi]->getP(), p, v->getP(), w->getP()) == 1)
    return Feature(vs[fi], es[fi], fi < 2);
  bool flag = (v == vs[0] || v == vs[1]) == (v < w);
  return Feature(v, w, flag);
}

int getVertex (Vertex *v, VIMap &vimap)
{
  VIMap::iterator iter = vimap.find(v);
  if (iter != vimap.end())
    return iter->second;
  int k = vimap.size();
  vimap.insert(VIPair(v, k));
  return k;
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

bool check (Separator *s, const Vertices &v, const Vertices &w)
{
  bool res = true;
  for (Vertices::const_iterator vi = v.begin(); vi != v.end(); ++vi)
    for (Vertices::const_iterator wi = w.begin(); wi != w.end(); ++wi)
      if (SeparatorOrder(s, (*vi)->getP(), (*wi)->getP()) < 1)
	res = false;
  return res;
}

int SeparatorOrder::sign ()
{
  return s->getU().dot(b->getP() - a->getP()).sign();
}

void moveVertices (Polyhedron *a, VPMap &vpmap, const VIMap &vimap, double *dvw)
{
  for (VIMap::const_iterator i = vimap.begin(); i != vimap.end(); ++i) {
    double *di = dvw + 3*i->second;
    if (di[0] != 0.0 || di[1] != 0.0 || di[2] != 0.0) {
      Vertex *v = i->first;
      vpmap.insert(VPPair(v, v->getP()));
      PTR<Point> vp = v->getP(),
	dp = new InputPoint(di[0], di[1], di[2], false),
	q = new SumPoint(vp, dp);
      a->moveVertex(v, q);
    }
  }
}

FeatureSet smallFeatures (Polyhedron *a, double d, const FeatureSet &fs,
			  VPMap *vpmap)
{
  FeatureSet nfs = smallFeatures(a, d, vpmap);
  for (FeatureSet::iterator f = fs.begin(); f != fs.end(); ++f)
    if (f->small(d))
      nfs.insert(*f);
  return nfs;  
}

VPMap savePos (const FeatureSet &fs)
{
  VertexSet vs;
  for (FeatureSet::const_iterator f = fs.begin(); f != fs.end(); ++f)
    f->vertices(vs);
  VPMap vpmap;
  for (VertexSet::iterator v = vs.begin(); v != vs.end(); ++v)
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

// debug

int fsn (const FeatureSet &fs)
{
  return fs.size();
}

double amax (int n, double *a)
{
  double x = 0.0;
  for (int i = 0; i < n; ++i)
    x = max(x, fabs(a[i]));
  return x;
}

void pedges (const Edges &, int);

void pedges (const FeatureSet &fs, int i)
{
  EdgeSet es;
  for (FeatureSet::const_iterator f = fs.begin(); f != fs.end(); ++f)
    switch (f->type) {
    case VV: {
      HEdge *e = f->v->connected(f->w);
      if (e)
	es.insert(e->getE());
      break;
    }
    case VE: {
      es.insert(f->e);
      for (int j = 0; j < 2; ++j) {
	HEdge *e = f->e->getHEdge(j);
	if (e->getNext()->head() == f->v) {
	  es.insert(e->getE());
	  es.insert(e->getNext()->getE());
	  break;
	}
      }
      break;
    }
    case VF: {
      HEdges ed;
      f->fa->boundaryHEdges(ed);
      for (int i = 0; i < 3; ++i)
	es.insert(ed[i]->getE());
      break;
    }
    case EE:
      es.insert(f->e);
      es.insert(f->f);
    }
  Edges ed(es.begin(), es.end());
  pedges(ed, i);
}
  
      
