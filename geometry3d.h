#ifndef GEOMETRY3D
#define GEOMETRY3D

#include "polyhedron.h"

class PlaneSide : public Primitive {
  Plane *plane;
  Point *point;

  Parameter calculate () {
    PV3 v = plane->get().n, p = point->get();
    Parameter d = plane->get().k;
    return v.dot(p) - d;
  }
 public:
  PlaneSide (Plane *plane, Point *point) : plane(plane), point(point) {}
};

class DiffLength : public Primitive {
  Point *i, *j;

  Parameter calculate () {
    return i->get().dot(i->get()) - j->get().dot(j->get());
  }
 public:
  DiffLength (Point *i, Point *j) : i(i), j(j) {}
};

class IsTangent : public Primitive {
  Point *v0, *v1, *v2;

  int sign () {
    PV3 a = v0->get(), b = v1->get(), c = v2->get(),
      n = (b - a).cross(c - a), nsplit(n.y, - n.x, Parameter::constant(0.0));
    Parameter k = n.dot(a);
    int psa = nsplit.dot(a).sign(), psb = nsplit.dot(b).sign(),
      psc = nsplit.dot(c).sign();
    return psa == psb && psb == psc ? 0 : 1;
  }
 public:
  IsTangent (Point *v0, Point *v1, Point *v2) : v0(v0), v1(v1), v2(v2) {}
};

class Orient3D : public Primitive {
  Point *a, *b, *c, *d;

  Parameter calculate () {
    PV3 dd = d->get(), da = a->get() - dd, db = b->get() - dd,
      dc = c->get() - dd;
    return da.cross(db).dot(dc);
  }
 public:
  Orient3D (Point *a, Point *b, Point *c, Point *d) : a(a), b(b), c(c), d(d) {}
};

class InSphere : public Primitive {
  Point *a, *b, *c, *d, *e;

  Parameter calculate () {
    PV3 ee = e->get(), ea = a->get() - ee, eb = b->get() - ee,
      ec = c->get() - ee, ed = d->get() - ee;
    return ea.cross(eb).dot(ec)*ed.dot(ed) -
      eb.cross(ec).dot(ed)*ea.dot(ea) +
      ec.cross(ed).dot(ea)*eb.dot(eb) -
      ed.cross(ea).dot(eb)*ec.dot(ec);
  }
 public:
  InSphere (Point *a, Point *b, Point *c, Point *d, Point *e)
    : a(a), b(b), c(c), d(d), e(e) {}
};

class RotationPoint : public Point {
  PTR<Point> point, sca;

  PV3 calculate () {
    Parameter sint = sca->get().x, cost = sca->get().y;
    PV3 p = point->get();
    return PV3(cost*p.x - sint*p.y, sint*p.x + cost*p.y, p.z);
  }
  
 public:
  RotationPoint (Point *point, Point *sca) : point(point), sca(sca) {} 
};

class TangentIntersectionPoint : public Point {
  PTR<Point> point, sca;

  PV3 calculate () {
    Parameter alpha = sca->get().z;
    PV3 p = point->get();
    return PV3(p.x - alpha*p.y, p.y + alpha*p.x, p.z);
  }

 public:
  TangentIntersectionPoint (Point *point, Point *sca) : point(point), sca(sca) {}
};

class ZIntercectPoint : public Point {
  PTR<Point> a, b, c;

  PV3 calculate () {
    Parameter t = (c->get().z - a->get().z)/(b->get().z - a->get().z);
    return a->get() + t*(b->get() - a->get());
  }

 public:
 ZIntercectPoint (Point *a, Point *b, Point *c) : a(a), b(b), c(c) {}
};

class SplitPlane : public Plane {
  PTR<Point> v0, v1, v2;
  
  PlaneData calculate () {
    PV3 a = v0->get(), b = v1->get(), c = v2->get(),
      n_tri = (b - a).cross(c - a); 
    Parameter k_tri = n_tri.dot(a); 
    PV3 n(n_tri.y, - n_tri.x, Parameter::constant(0.0));
    return PlaneData(n, Parameter::constant(0.0));
  }

 public:
  SplitPlane (Point *v0, Point *v1, Point *v2) : v0(v0), v1(v1), v2(v2) {}
};

class PointNormalPlane : public Plane {
  PTR<Object<PV3> > point, normal;
  
  PlaneData calculate () {
    PV3 p = point->get(), n = normal->get();
    Parameter d = n.dot(p);
  return PlaneData(n, d);
  }

 public:
  PointNormalPlane (Point *point, Point *normal) : point(point), normal(normal) {}
};

class IntersectionPoint : public Point {
  PTR<Point> tail, head;
  PTR<Plane> plane;

  PV3 calculate () {
    PV3 t = tail->get(), h = head->get(), v = plane->get().n;
    Parameter d = plane->get().k, s = (d - t.dot(v))/(h - t).dot(v);
    return t + (h - t)*s;
  }
  
 public:
  IntersectionPoint (Point *tail, Point *head, Plane *plane)
    : tail(tail), head(head), plane(plane) {}
};

class FaceIntersectionPoint : public Point {
  PTR<Point> head, tail, a, b, c;
  HFace * hface;
  
  PV3 calculate () {
    PV3 t = tail->get(), v = head->get() - t,
      n = (c->get() - b->get()).cross(a->get() - b->get());
    Parameter kt = - n.dot(a->get()), k = - (n.dot(t) + kt)/n.dot(v);
    return t + k*v;
  }

 public:
  FaceIntersectionPoint (Point *tail, Point *head, HFace *hface)
    : tail(tail), head(head), hface(hface)
  { a = hface->getF()->getP()->getA();
    b = hface->getF()->getP()->getB();
    c = hface->getF()->getP()->getC();
  }
  
  HFace * getHFace() { return hface; }
};

#define vol(a,b,c,d) (a-d).cross(b-d).dot(c-d)

class FaceNearestPoint : public Point {
  PTR<Point> point, pa, pb, pc;

  PV3 calculate () {
    PV3 p = point->get(), a = pa->get(), b = pb->get(),
      c = pc->get(), n = (a - b).cross(c - b);
    Parameter d = n.dot(b), s = (d - p.dot(n))/n.dot(n);
    PV3 q = p + n*s;
    bool sideab = vol(q, a, b, b + n).sign() == vol(c, a, b, b + n).sign();
    bool sideac = vol(q, a, c, c + n).sign() == vol(b, a, c, c + n).sign();
    bool sidebc = vol(q, b, c, c + n).sign() == vol(a, b, c, c + n).sign();
    if (sideac && sideab && sidebc) 
      return q;
    if (sideab && sideac) {
      Parameter t = (q-c).dot(b-c);
      if (t < 0.0)
	return c;
      if (t > 1.0)
	return b;
      return c + t*(b-c);
    }
    if (sideab && sidebc) {
      Parameter t = (q-c).dot(a-c);
      if (t < 0.0)
	return c;
      if (t > 1.0)
	return a;
      return c + t*(a-c);
    }
    if (sideac && sidebc) {
      Parameter t = (q-a).dot(b-a);
      if (t < 0.0)
	return a;
      if (t > 1.0)
	return b;
      return a + t*(b-a);    
    }
    if (!sideab && !sideac)
      return a;
    if (!sideab && !sidebc)
      return b;
    if (!sideac && !sidebc)
      return c;
    assert(false);
    return p;
  }
 
public: 
  FaceNearestPoint (Point *point, HFace *hf) : point(point) {
    pa = hf->getF()->getP()->getA();
    pb = hf->getF()->getP()->getB();
    pc = hf->getF()->getP()->getC();
  }
};

#endif
