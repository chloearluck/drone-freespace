#include "geometry3d.h"
#include <iostream>

int PlaneSide::sign () {
  PV3 v = plane->getN();
  Parameter d = plane->getK();
  PV3 p = point->getP();
  return (v.dot(p) - d).sign();
}

// Comment out to test how ACP handles inaccuracy.
#define ACCURATE_ORIENT3D

int Orient3D::sign () {
#ifdef ACCURATE_ORIENT3D
  PV3 d = this->d->getP();
  PV3 da = this->a->getP() - d;
  PV3 db = this->b->getP() - d;
  PV3 dc = this->c->getP() - d;
  return da.cross(db).dot(dc).sign();
#else
  PV3 d = this->d->getP();
  PV3 a = this->a->getP();
  PV3 b = this->b->getP();
  PV3 c = this->c->getP();

  Parameter o = (a.cross(b).dot(c) -
		 d.cross(b).dot(c) +
		 d.cross(a).dot(c) -
		 d.cross(a).dot(b));

  return o.sign();
#endif
}

int InSphere::sign () {
  PV3 e = this->e->getP();
  PV3 ea = a->getP() - e;
  PV3 eb = b->getP() - e;
  PV3 ec = c->getP() - e;
  PV3 ed = d->getP() - e;

  return (ea.cross(eb).dot(ec) * ed.dot(ed) -
          eb.cross(ec).dot(ed) * ea.dot(ea) +
          ec.cross(ed).dot(ea) * eb.dot(eb) -
          ed.cross(ea).dot(eb) * ec.dot(ec)).sign();
}

PlaneData TrianglePlane2::calculate () {
  PV3 a = this->a->get();
  PV3 b = this->b->get();
  PV3 c = this->c->get();
  PV3 v = (b - a).cross(c - a);
  Parameter d = v.dot(a);
  return(PlaneData(v, d));
}


PlaneData PointNormalPlane::calculate () {
  PV3 p = this->point->get();
  PV3 n = this->normal->get();

  Parameter d = n.dot(p);
  return(PlaneData(n, d));
} 

PV3 RotationPoint::calculate () {
  Parameter sint = sin_cos_alpha->getP().getX();
  Parameter cost = sin_cos_alpha->getP().getY();
  PV3 p = point->getP();
  return PV3(cost * p.getX() - sint * p.getY(),
             sint * p.getX() + cost * p.getY(), 
             p.getZ() );
}

PV3 TangentIntersectionPoint::calculate () {
  Parameter alpha = sin_cos_alpha->getP().getZ();
  PV3 p = point->getP();
  return PV3(p.getX()-alpha*p.getY(), 
             p.getY()+alpha*p.getX(), 
             p.getZ() );
}

PV3 AddPoint::calculate () {
  return point1->getP() + point2->getP();
}

PV3 ACrossBPoint::calculate () {
  return point1->getP().cross(point2->getP());
}

PV3 NormalToPlane::calculate () {
  return plane->getN();
}

PV3 PlanePoint::calculate ()  {
  PV3 v = plane->getN();
  Parameter d = plane->getK();
  PV3 p = point->getP();
  return(p - v * ((v.dot(p) - d) / v.dot(v)));
}

PV3 IntersectionPoint::calculate () {
  PV3 t = tail->getP();
  PV3 h = head->getP();
  PV3 v = plane->getN();
  Parameter d = plane->getK();
  // (t + s (h - t)) * v = d
  // t * v + s (h - t) * v = d
  // s = (d - t * v) / ((h - t) * v)

  Parameter s = (d - t.dot(v)) / (h - t).dot(v);
  return(t + (h - t) * s);
}

