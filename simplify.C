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

void insert (Vertex *v, VPMap &vpmap)
{
  Point *p = v->getP();
  vpmap.insert(VPPair(v, p)).second;
}

bool simplify1 (Polyhedron *a, double d, FOctree *octree, VPMap &vpmap,
		int &n1, int &n2)
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
  if (e->HEdgesN() > 0) {
    if (DistancePP(e->getT()->getP(), e->getH()->getP(), d) == 1)
      fq.push(IFeature(e));
    else
      addFlips(e, d, fq);
  }
}

void addFlips (Edge *e, double d, IFeatureQ &fq)
{
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
  double dc = collapseDistance(e);
  HEdge *e1 = e->getHEdge(0), *e2 = e->getHEdge(1);
  Vertex *t = e1->tail(), *h = e1->head(), *v = e1->getNext()->head(),
    *w = e2->getNext()->head();
  PTR<Point> tp = t->getP(), hp = h->getP();
  PV3 q = 0.5*(tp->getApprox() + hp->getApprox());
  PTR<Point> cp = new InputPoint(PV3::input(q.x.mid(), q.y.mid(), q.z.mid()));
  a->removeLoop(e1);
  a->removeLoop(e2);
  insert(h, vpmap);
  a->moveVertex(h, cp);
  Faces ft = t->incidentFaces();
  for (Faces::iterator f = ft.begin(); f != ft.end(); ++f)
    a->replaceVertex(*f, t, h);
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

double collapseDistance (Edge *e)
{
  PV3 u = e->getU();
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
  Vertices ve;
  f->boundaryVertices(ve);
  for (int i = 0; i < 3; ++i)
    ef[i] = ve[i] == v;
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
    else if (si == 0) {
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
  int pc = f->getPC();
  for (int i = 0; i < 3; ++i)
    if ((ef[i] || ee) && intersectsEE(ed[i]->getE(), e, pc))
      return true;
  return false;
}

bool intersectsEE (Edge *e, Edge *f, int pc)
{
  if (!bboxOverlap(e->getBBox(), f->getBBox()))
    return false;
  Point *et = e->getT()->getP(), *eh = e->getH()->getP(), *ft = f->getT()->getP(),
    *fh = f->getH()->getP();
  int tp1 = TripleProductR(et, ft, fh), tp2 = TripleProductR(eh, ft, fh);
  if (tp1*tp2 > -1)
    return false;
  if (tp1 == 0)
    return onEdge(et, ft, fh, true);
  if (tp2 == 0)
    return onEdge(eh, ft, fh, true);
  int tp3 = TripleProductR(ft, et, eh), tp4 = TripleProductR(fh, et, eh);
  return 
    tp3 == 0 && onEdge(ft, et, eh, true) || tp4 == 0 && onEdge(fh, et, eh, true) ||
    tp3*tp4 == -1;
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
    return collapse(a, d, e->getNext()->getNext()->getE(), fq, octree, vpmap, bad);
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
    double dmax, dr = displacement(vpmap, dmax)/vpmap.size();
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
    double dmax, dprev = displacement(vpmap, dmax), dmin = dprev, bound = 1.0;
    int rcount = 0, dcount = 0;
    while (bound > 1e-5) {
      ++n;
      simplify2s(a, d, fslp, vpmap, false, bound);
      double dbest = displacement(vpmap, dmax);
      while (separation(fslp) < 0.99999*d) {
	++n;
	simplify2s(a, d, fslp, vpmap);
      }
      double dcurr = displacement(vpmap, dmax);
      if (dcurr > 0.999*dmin) {
	++dcount;
	if (dcount == 2)
	  break;
      }
      else
	dcount = 0;
      dmin = min(dmin, dcurr);
      double ratio = dprev == dbest ? 0.0 : fabs((dprev - dcurr)/(dprev - dbest));
      if (ratio < 0.25)
	bound *= 0.5;
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
  return DistancePP(p, q, d) == 1 && f->contains(q);
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

void closestPoints (Vertex *v, Edge *e, PV3 &p, PV3 &q)
{
  Point *a = v->getP(), *t = e->getT()->getP(), *h = e->getH()->getP();
  p = a->getP();
  if (onEdge(a, t, h, true))
    q = LinePoint(a, t, h).getP();
  else
    q = (CloserPair(a, t, a, h) == 1 ? t : h)->getP();
}

void closestPoints (Vertex *v, Face *f, PV3 &p, PV3 &q)
{
  Point *a = v->getP();
  p = a->getP();
  if (a->side(f->getP()) != 0) {
    PlanePoint pp(f->getP(), a);
    if (f->containsConvex(&pp, false)) {
      q = pp.getP();
      return;
    }
  }
  HEdge *h = f->getBoundary(0);
  for (int i = 0; i < 3; ++i) {
    PV3 pi, qi;
    closestPoints(v, h->getE(), pi, qi);
    if (i == 0 || closer(pi, qi, p, q))
      q = qi;
    h = h->getNext();
  }
}

void closestPoints (Edge *e, Edge *f, PV3 &p, PV3 &q)
{
  Vertex *vs[] = {e->getT(), e->getH(), f->getT(), f->getH()};
  Point *et = vs[0]->getP(), *eh = vs[1]->getP(), *ft = vs[2]->getP(),
    *fh = vs[3]->getP();
  if (Orientation(et, eh, ft, fh) != 0) {
    LineLinePoint pe(et, eh, ft, fh), pf(ft, fh, et, eh);
    if (onEdge(&pe, et, eh, true) && onEdge(&pf, ft, fh, true)) {
      p = pe.getP();
      q = pf.getP();
      return;
    }
  }

  Edge *es[] = {f, f, e, e};
  closestPoints(vs[0], es[0], p, q);
  for (int i = 1; i < 4; ++i) {
    PV3 pi, qi;
    closestPoints(vs[i], es[i], pi, qi);
    if (closer(pi, qi, p, q)) {
      p = i == 1 ? pi : qi;
      q = i == 1 ? qi : pi;
    }
  }
}

bool closer (const PV3 &a, const PV3 &b, const PV3 &c, const PV3 &d)
{
  PV3 u = b - a, v = d - c;
  return u.dot(u) < v.dot(v);
}

PV3 orthogonal (const PV3 &u)
{
  Parameter ux = u.x.abs(), uy = u.y.abs(), uz = u.z.abs();
  if ((uy - ux).sign() > -1 && (uz - ux).sign() > -1)
    return PV3(Parameter::constant(0.0), - u.z, u.y);
  if ((ux - uy).sign() > -1 && (uz - uy).sign() > -1)
    return PV3(- u.z, Parameter::constant(0.0), u.x);
  return PV3(- u.y, u.x, Parameter::constant(0.0));
}

FeatureSet smallFeatures (Polyhedron *a, double d, VPMap *vpmap)
{
  FOctree *octree = faceOctree(a, d);
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f)
    if (!(*f)->getBoundary().empty()) {
      (*f)->getPC();
      //(*f)->getP()->getN();
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
	(!vpmap || moved(fa[i], vpmap) || moved(fb[i], vpmap)))
      smallCandidates(fa[i], fb[i], td->d, td->fs);
}

bool moved (Face *f, VPMap *vpmap)
{
  Vertices ve;
  f->boundaryVertices(ve);
  for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
    if (vpmap->find(*v) != vpmap->end())
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
    bool flag = vf[i]->getOutgoing(0)->getF() == f;
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
    bool flag = vg[i]->getOutgoing(0)->getF() == g;
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

bool closeVVT (Vertex *v, Vertex *w, double d)
{
  PV3 u = v->getP()->getApprox(1e-6) - w->getP()->getApprox(1e-6);
  Parameter k = d*d - u.dot(u);
  return k.ub() > 0.0;
}

bool closeVET (Vertex *v, Edge *e, double d)
{
  if (v == e->getT() || v == e->getH())
    return false;
  PV3 a = v->getP()->getApprox(1e-6), t = e->getT()->getP()->getApprox(1e-6),
    h = e->getH()->getP()->getApprox(1e-6), u = h - t;
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
  PV3 p = v->getP()->getApprox(1e-6), n = f->getP()->getApprox(1e-6).n;
  Parameter k1 = n.dot(p) + f->getP()->getApprox(1e-6).k, k2 = d*d*n.dot(n) - k1*k1;
  if (k2.sign(false) == -1)
    return false;
  PV3 pts[3];
  for (int i = 0; i < 3; ++i)
    pts[i] = ve[i]->getP()->getApprox(1e-6);
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
  PV3 aa = a->getApprox(1e-6), u = p->getApprox(1e-6) - aa,
    bb = b->getApprox(1e-6), v = q->getApprox(1e-6) - bb, w = bb - aa;
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

void simplify2s (Polyhedron *a, double d, FeatureSet &fs, VPMap &vpmap,
		 bool vobj, double vbound)
{
  FeatureSet fslp = featuresLP(fs);
  VIMap vimap;
  double *v = solveLP(fslp, d, vpmap, vobj, vbound, vimap);
  moveVertices(a, vpmap, vimap, v);
  delete [] v;
  fs = smallFeatures(a, 3.465*d, fs, &vpmap);
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

FeatureSet featuresLP (const FeatureSet &fs)
{
  FeatureSet fslp;
  for (FeatureSet::const_iterator f = fs.begin(); f != fs.end(); ++f)
    if (f->type == VV)
      featuresLPVV(f->v, f->w, fslp);
    else if (f->type == VE)
      featuresLPVE(f->v, f->e, fslp);
    else
      fslp.insert(*f);
  return fslp;
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

int InnerProduct4::sign ()
{
  return (a->getP() - b->getP()).dot(c->getP() - d->getP()).sign();
}

double separation (const FeatureSet &fs)
{
  double d = 1e10;
  for (FeatureSet::const_iterator f = fs.begin(); f != fs.end(); ++f)
    d = min(d, Separator(*f, 0.0).getApprox(1e-6).pq.ub());
  return d;
}

double * solveLP (const FeatureSet &fs, double d, const VPMap &vpmap,
		  bool vobj, double vbound, VIMap &vimap)
{
  Expander2 expander;
  for (FeatureSet::const_iterator f = fs.begin(); f != fs.end(); ++f) {
    Vertices ve[2];
    if (f->v) {
      ve[0].push_back(f->v);
      if (f->w)
	ve[1].push_back(f->w);
      else if (f->e) {
	ve[1].push_back(f->e->getT());
	ve[1].push_back(f->e->getH());
      }
      else
	f->fa->boundaryVertices(ve[1]);
    }
    else if (f->e) {
      ve[0].push_back(f->e->getT());
      ve[0].push_back(f->e->getH());
      if (f->f) {
	ve[1].push_back(f->f->getT());
	ve[1].push_back(f->f->getH());
      }
      else
	f->fa->boundaryVertices(ve[1]);
    }
    Separator s(*f, d);
    if (!check(&s, ve[0], ve[1]))
      cerr << "bad feature" << endl;
    
    
    SeparatorData data = s.getApprox();
    Expander2::Point u(data.u.x.mid(), data.u.y.mid(), data.u.z.mid());
    Expander2::Point v(data.v.x.mid(), data.v.y.mid(), data.v.z.mid());
    Expander2::Point w(data.w.x.mid(), data.w.y.mid(), data.w.z.mid());

    Expander2::Pair pair(u, v, w, data.pq.mid()/d);
    for (int ab = 0; ab < 2; ++ab)
      for (Vertices::iterator x = ve[ab].begin(); x != ve[ab].end(); ++x) {
	int xi = getVertex(*x, vimap);
	double vwr[3];
	s.getVWR((*x)->getP(), ab == 0, vwr);
	pair.addConstraint(ab, Expander2::Constraint(xi, vwr[0], vwr[1], vwr[2]));
      }
    expander.addPair(pair);
  }
  for (VIMap::iterator i = vimap.begin(); i != vimap.end(); ++i) {
    double dv[3];
    getDV(i->first, d, vpmap, dv);
    expander.addDisplacement(i->second, Expander2::Point(dv));
  }
  expander.expandV(1.0, vobj, vbound);
  double *x = new double [3*vimap.size()];
  for (VIMap::iterator i = vimap.begin(); i != vimap.end(); ++i) {
    int k = i->second;
    Expander2::Point p = expander.motion(k);
    for (int j = 0; j < 3; ++j)
      x[3*k+j] = d*p.x[j];
  }
  return x;
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
  PV3 u = (v->getP()->getApprox() - i->second->getApprox())/d;
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
      insert(v, vpmap);
      PTR<Point> vp = v->getP(),
	dp = new InputPoint(PV3::constant(di[0], di[1], di[2])),
	q = new SumPoint(vp, dp);
      a->moveVertex(v, q);
    }
  }
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

int euler (Polyhedron *a)
{
  int n = 0;
  for (Vertices::iterator v = a->vertices.begin(); v != a->vertices.end(); ++v)
    if ((*v)->EdgesN() > 0)
      ++n;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f)
    if (!(*f)->getBoundary().empty())
      ++n;
  for (Edges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    if ((*e)->HEdgesN() > 0)
      --n;
  return n;
}

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
