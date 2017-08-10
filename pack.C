#include "pack.h"

double pack3 (Polyhedron *a, Polyhedron *b, Polyhedron *c, double minsep, PTR<Point> t[3])
{
  t[0] = t[1] = t[2] = 0;
  double d[3] = {maxBBoxDimension(a), maxBBoxDimension(b), maxBBoxDimension(c)};
  sort(d, d + 3);
  double lb = d[2], ub = d[1] + d[2] + 0.1, tt = 0.0;
  int k = 0;
  cerr << setprecision(16) << "initial interval: [" << lb << ", " << ub << "]"
       << endl;
  while (ub - lb > minsep) {
    ++k;
    double m = 0.5*(lb + ub), t0 = getTime();
    bool flag = pack3(a, b, c, t, m);
    double t1 = getTime() - t0;
    cerr << k << (flag ? ". yes " : ". no ") << m << "; time: " << t1 << endl;
    tt += t1;
    if (flag)
      ub = m;
    else
      lb = m;
  }
  cerr << "minimum size: " << ub << "; time: " << tt << endl;
  return ub;
}

double maxBBoxDimension (Polyhedron *a)
{
  double *bbox = a->bbox;
  return max(bbox[1] - bbox[0], max(bbox[3] - bbox[2], bbox[5] - bbox[4]));
}

bool pack3 (Polyhedron *a, Polyhedron *b, Polyhedron *c, PTR<Point> t[3], double bs)
{
  double u01[6], u02[6], u03[6], v012[6], v013[6], v023[6],
    bb[6] = {0.0, bs, 0.0, bs, 0.0, bs};
  freeSpace(a, bb, u01);
  freeSpace(b, bb, u02);
  freeSpace(c, bb, u03);
  freeSpace(u01, u02, v012);
  freeSpace(u01, u03, v013);
  freeSpace(u02, u03, v023);
  Polyhedron *mb = b->negative(), *mc = c->negative(),
    *u12 = minkowskiSumFull(a, mb),
    *u13 = minkowskiSumFull(a, mc),
    *u23 = minkowskiSumFull(b, mc),
    *b012 = box(v012), *v12 = b012->boolean(u12, Complement),
    *b013 = box(v013), *v13 = b013->boolean(u13, Complement),
    *b023 = box(v023), *v23 = b023->boolean(u23, Complement),
    *s123 = minkowskiSumFull(v12, v23),
    *w13 = v13->boolean(s123, Intersection);
  bool res = w13->faces.size() > 0;
  if (res) {
    PTR<Point> t13 = interiorPoint(w13);
    double u03t[6], v01[6], v01u02[6], u02t[6], w01[6];
    translateBox(u03, t13, u03t);
    intersectBoxes(u01, u03t, v01);
    freeSpace(v01, u02, v01u02);
    Polyhedron *bv01u02 = box(v01u02),
      *v23mt = v23->negativeTranslate(t13),
      *temp = v23mt->boolean(bv01u02, Intersection),
      *w12 = v12->boolean(temp, Intersection);
    if (w12->faces.empty())
      res = false;
    PTR<Point> t12 = interiorPoint(w12);
    translateBox(u02, t12, u02t);
    intersectBoxes(v01, u02t, w01);
    t[0] = midpointBox(w01);
    t[1] = new SumPoint(t[0], t12);
    t[2] = new SumPoint(t[0], t13);
    delete bv01u02;
    delete v23mt;
    delete temp;
    delete w12;
  }
  delete mb; delete mc; delete u12;  delete u13; delete u23;
  delete b012; delete v12; delete b013; delete v13; delete b023; delete v23;
  delete s123; delete w13;
  return res;
}

void freeSpace (Polyhedron *a, double *box, double *res)
{
  double *bboxa = a->bbox;
  for (int i = 0; i < 6; ++i)
    res[i] = box[i] - bboxa[i];
}

// a moves with respect to b
void freeSpace (double *a, double *b, double *res)
{
  for (int i = 0; i < 3; ++i) {
    res[2*i] = b[2*i] - a[2*i+1];
    res[2*i+1] = b[2*i+1] - a[2*i];
  }
}

void translateBox (double *a, PTR<Point> t, double *res)
{
  PV3 p = t->getP();
  double tx = - p.x.mid(), ty = - p.y.mid(), tz = - p.z.mid();
  res[0] = a[0] + tx;
  res[1] = a[1] + tx;
  res[2] = a[2] + ty;
  res[3] = a[3] + ty;
  res[4] = a[4] + tz;
  res[5] = a[5] + tz;
}

void intersectBoxes (double *a, double *b, double *res)
{
  for (int i = 0; i < 3; ++i) {
    res[2*i] = max(a[2*i], b[2*i]);
    res[2*i+1] = min(a[2*i+1], b[2*i+1]);
  }
}

PTR<Point> midpointBox (double *a)
{
  return new InputPoint(0.5*(a[0] + a[1]), 0.5*(a[2] + a[3]), 0.5*(a[4] + a[5]));
}

PTR<Point> interiorPoint (Polyhedron *a)
{
  Face *f = largestFace(a);
  PTR<Point> p = f->centroid();
  HFaceNormal n(f->getHFace(1));
  PTR<Point> qmin = 0;
  for (Faces::iterator g = a->faces.begin(); g != a->faces.end(); ++g)
    if (*g != f) {
      PTR<Point> q = (*g)->rayIntersection(p, &n);
      if (q && (!qmin || CloserPair(p, q, p, qmin) == 1))
	qmin = q;
    }
  return new MidPoint(p, qmin);
}

Face * largestFace (Polyhedron *a)
{
  double amax = 0.0;
  Face *fmax = 0;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f) {
    double af = area(*f);
    if (amax < af) {
      amax = af;
      fmax = *f;
    }
  }
  return fmax;
}

double area (Face *f)
{
  int c = f->getPC();
  HEdge *e = f->getBoundary(0);
  PV3 u = e->getU(), v = e->getNext()->getU();
  return 0.5*cross(u, v, c).mid();
}

void pack3output (Polyhedron *a, Polyhedron *b, Polyhedron *c,
		  PTR<Point> t[3], ostream &ostr)
{
  Polyhedron *at = a->translate(t[0]), *bt = b->translate(t[1]),
    *ct = c->translate(t[2]);
  Faces fa(at->faces);
  fa.insert(fa.end(), bt->faces.begin(), bt->faces.end());
  fa.insert(fa.end(), ct->faces.begin(), ct->faces.end());
  writePolyhedronVTK(fa, ostr);
}

int main (int argc, char *argv[])
{
  if (argc < 4)
    return 0;
  ifstream astr(argv[1]), bstr(argv[2]), cstr(argv[3]);
  if (!(astr.good() && bstr.good() && cstr.good()))
    return 0;
  inputPerturbed = argc == 4 ? true : *argv[4] != 't';
  Parameter::enable();
  Polyhedron *a = readPolyhedronVTK(astr, inputPerturbed),
    *b = readPolyhedronVTK(bstr, inputPerturbed),
    *c = readPolyhedronVTK(cstr, inputPerturbed);
  double minsep = 1e-6;
  PTR<Point> t[3];
  pack3(a, b, c, minsep, t);
  if (t[0])
    pack3output(a, b, c, t, cout);
  delete a;
  delete b;
  delete c;
}
