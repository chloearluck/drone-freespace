#include "hull.h"

bool HullEdge::horizonEdge (HullVertex *v) const 
{
  return !face->conflict(v) && twin->face->conflict(v);
}

bool HullFace::visible (HullVertex *v) const
{
  HullEdge *f = boundary->twin->next, *g = f->twin->next;
  int s = Orientation(boundary->tail->p, f->tail->p, g->tail->p, v->p);
  return s == 1 || s == 0 && contains(v);
}

bool HullFace::contains (HullVertex *v) const
{
  Point *a = v->p;
  Points b;
  b.push_back(boundary->tail->p);
  b.push_back(boundary->head()->p);
  b.push_back(boundary->twin->next->head()->p);
  TrianglePlane tp(b[0], b[1], b[2]);
  int c = projectionCoordinate(&tp);
  for (int i = 0; i < 3; ++i)
    if (LeftTurn(a, b[i], b[(i+1)%3], c) < 1)
      return false;
  return true;
}

bool HullFace::conflict (HullVertex *v) const
{
  return find(conflicts.begin(), conflicts.end(), v) != conflicts.end();
}

Hull::~Hull ()
{
  for (HullVertices::iterator v = vertices.begin(); v != vertices.end(); ++v)
    delete *v;
  for (HullEdges::iterator e = edges.begin(); e != edges.end(); ++e)
    delete *e;
  for (HullFaces::iterator f = faces.begin(); f != faces.end(); ++f)
    delete *f;
}

HullVertex * Hull::addVertex (Point *p)
{
  HullVertex *v = new HullVertex(p);
  vertices.push_back(v);
  return v;
}

HullEdge * Hull::addEdge (HullVertex *t, HullVertex *h)
{
  HullEdge *e = new HullEdge(t, 0, 0), *et = new HullEdge(h, e, 0);
  e->twin = et;
  edges.push_back(e);
  edges.push_back(et);
  return e;
}

HullFace * Hull::addFace (HullEdge *e)
{
  HullEdge *f = e->twin->next, *g = f->twin->next;
  HullFace *face = new HullFace(e);
  faces.push_back(face);
  e->face = face;
  f->face = face;
  g->face = face;  
  return face;
}

void Hull::formHull ()
{
  initHull();
  for (int i = 4; i < vertices.size(); ++i)
    expandHull(i);
}

void Hull::initHull ()
{
  HullEdge *e[12];
  e[0] = addEdge(vertices[0], vertices[1]);
  e[3] = e[0]->twin;
  e[1] = addEdge(vertices[0], vertices[2]);
  e[6] = e[1]->twin;
  e[2] = addEdge(vertices[0], vertices[3]);
  e[9] = e[2]->twin;
  e[4] = addEdge(vertices[1], vertices[2]);
  e[7] = e[4]->twin;
  e[5] = addEdge(vertices[1], vertices[3]);
  e[10] = e[5]->twin;
  e[8] = addEdge(vertices[2], vertices[3]);
  e[11] = e[8]->twin;
  for (int i = 0; i < 4; ++i) {
    vertices[i]->edge = e[3*i];
    if (Orientation(e[3*i]->head()->p, e[3*i+1]->head()->p, e[3*i+2]->head()->p, e[3*i]->tail->p) == -1)
      for (int j = 0; j < 3; ++j)
        e[3*i+j]->next = e[3*i+(j+1)%3];
    else
      for (int j = 0; j < 3; ++j)
	      e[3*i+(j+1)%3]->next = e[3*i+j];
  }
  for (int i = 0; i < 12; ++i)
    if (!e[i]->face)
      addFace(e[i]);
  for (int i = 4; i < vertices.size(); ++i)
    for (int j = 0; j < 4; ++j)
      if (faces[j]->visible(vertices[i])) {
	faces[j]->conflicts.push_back(vertices[i]);
	vertices[i]->conflicts.push_back(faces[j]);
      }
}

void Hull::expandHull (int i)
{
  HullVertex *v = vertices[i];
  if (v->conflicts.empty())
    return;
  HullEdges h, ev;
  findHorizon(v, h);
  addCone(v, h);
  removeVisible(i);
}

void Hull::findHorizon (HullVertex *v, HullEdges &h) const
{
  HullEdge *e = findHorizonEdge(v), *e0 = e;
  do {
    h.push_back(e);
    e = e->twin;
    while (!e->horizonEdge(v))
      e = e->next;
  }
  while (e != e0);
}

HullEdge * Hull::findHorizonEdge (HullVertex *v) const
{
  for (HullFaces::iterator f = v->conflicts.begin();
       f != v->conflicts.end(); ++f) {
    HullEdge *e = (*f)->boundary;
    for (int i = 0; i < 3; ++i, e = e->twin->next)
      if (e->twin->horizonEdge(v))
	return e->twin;
  }
  return 0;
}

void Hull::addCone (HullVertex *v, HullEdges &h)
{
  HullEdge *ep = *(h.end()-1), *evp = 0;
  for (HullEdges::iterator e = h.begin(); e != h.end(); ++e) {
    HullEdge *ev = addEdge(v, (*e)->tail);
    if (evp)
      evp->next = ev;
    else
      v->edge = ev;
    ev->twin->next = ep->twin;
    (*e)->next = ev->twin;
    (*e)->tail->edge = *e;
    evp = ev;
    ep = *e;
  }
  evp->next = v->edge;
  for (HullEdges::iterator e = h.begin(); e != h.end(); ++e)
    addConeFace(v, *e);
}

void Hull::addConeFace (HullVertex *v, HullEdge *e)
{
  const HullVertices &ec = e->face->conflicts, &etc = e->twin->face->conflicts;
  HullFace *f = addFace(e->twin);
  for (HullVertices::const_iterator w = ec.begin(); w != ec.end(); ++w)
    if (*w != v && f->visible(*w)) {
      f->conflicts.push_back(*w);
      (*w)->conflicts.push_back(f);
    }
  for (HullVertices::const_iterator w = etc.begin(); w != etc.end(); ++w)
    if (*w != v &&
	find(f->conflicts.begin(), f->conflicts.end(), *w) == f->conflicts.end()
	&& f->visible(*w)) {
      f->conflicts.push_back(*w);
      (*w)->conflicts.push_back(f);
    }
}

void Hull::removeVisible (int i) const
{
  const HullFaces &c = vertices[i]->conflicts;
  for (int j = i + 1; j < vertices.size(); ++j) {
    HullVertex *v = vertices[j];
    for (HullFaces::const_iterator f = c.begin(); f != c.end(); ++f) {
      HullFaces::iterator g = remove(v->conflicts.begin(), v->conflicts.end(), *f);
      v->conflicts.erase(g, v->conflicts.end());
    }
  }
  for (HullFaces::const_iterator f = c.begin(); f != c.end(); ++f)
    (*f)->boundary = 0;
}

Polyhedron * Hull::polyhedron () const
{
  Polyhedron *a = new Polyhedron;
  map<HullVertex *, Vertex *> vmap;
  for (HullVertices::const_iterator v = vertices.begin(); v != vertices.end(); ++v)
    vmap.insert(pair<HullVertex *, Vertex *>(*v, a->getVertex((*v)->p)));
  for (HullFaces::const_iterator f = faces.begin(); f != faces.end(); ++f) {
    HullEdge *e = (*f)->boundary;
    if (e) {
      HullVertex *u = e->tail, *v = e->head(), *w = e->twin->next->head();
      a->addTriangle(vmap.find(u)->second, vmap.find(v)->second,
		     vmap.find(w)->second);
    }
  }
  return a;
}

Polyhedron * convexHull (Polyhedron *a)
{
  Points pts;
  for (Vertices::iterator v = a->vertices.begin(); v != a->vertices.end(); ++v)
    pts.push_back((*v)->getP());
  return convexHull(pts);
}

Polyhedron * convexHull (const Points &pts)
{
  int n = pts.size(), *p = new int [n], i = 2;
  randomPermutation(n, p);
  while (pts[p[0]]->onLine(pts[p[1]], pts[p[i]]))
    ++i;
  if (i > 2) {
    int t = p[2];
    p[2] = p[i];
    p[i] = t;
  }
  i = 3;
  while (Orientation(pts[p[0]], pts[p[1]], pts[p[2]], pts[p[i]]) == 0)
    ++i;
  if (i > 3) {
    int t = p[3];
    p[3] = p[i];
    p[i] = t;
  }
  Hull hull;
  for (int i = 0; i < n; ++i) {
    hull.addVertex(pts[p[i]]);
  }
  delete [] p;
  hull.formHull();
  return hull.polyhedron();
}

void randomPermutation (int n, int *p)
{
  for (int i = 0; i < n; ++i)
    p[i] = i;
  for (int i = 0; i < n - 1; ++i) {
    int j = randomInteger(i, n - 1), ti = p[i];
    p[i] = p[j];
    p[j] = ti;
  }
}    

int randomInteger (int lb, int ub)
{
  double r = random()/(double) RAND_MAX;
  int d = (int) floor(r*(ub - lb + 1));
  return lb + d;
}

/*int main (int argc, char *argv[])
{
  if (argc < 2)
    return 0;
  ifstream astr(argv[1]);
  if (!astr.good())
    return 0;
  Parameter::enable();
  Polyhedron *a = readPolyhedronVTK(astr),
    *b = convexHull(a);
  delete b;
  delete a;
  Parameter::disable();
}*/

// debug

void pp (Point *p);

void pv (HullVertex *v)
{
  pp(v->p);
}

void pe (HullEdge *e)
{
  cerr << "(";
  pv(e->tail);
  pv(e->head());
  cerr << ")" << endl;
}

void pf (HullFace *f)
{
  cerr << "(";
  pv(f->boundary->tail);
  pv(f->boundary->twin->next->tail);
  pv(f->boundary->twin->next->twin->next->tail);
  cerr << ")" << endl;
}
