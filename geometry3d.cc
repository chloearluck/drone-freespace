#include "geometry3d.h"
#include <iostream>

int PlaneSide::sign () {
  PV3 v = plane->getN();
  Parameter d = plane->getK();
  PV3 p = point->getP();
  return (v.dot(p) - d).sign();
}

int IsTangent::sign () {
  PV3 a  = v0->getP();
  PV3 b  = v1->getP();
  PV3 c  = v2->getP();
  PV3 n = (b - a).cross(c - a); 
  Parameter k = n.dot(a);
  PV3 p = k / (n.dot(n)) * n;

  PV3 n0 = (a-b).cross(n); Parameter k0 = n0.dot(b);
  PV3 n1 = (b-c).cross(n); Parameter k1 = n1.dot(c);
  PV3 n2 = (c-a).cross(n); Parameter k2 = n2.dot(a);

  if ( ((n0.dot(p) - k0).sign() == (n0.dot(c) - k0).sign()) && 
       ((n1.dot(p) - k1).sign() == (n1.dot(a) - k0).sign()) && 
       ((n2.dot(p) - k2).sign() == (n2.dot(b) - k0).sign()) )
    return 1;  
  return 0;
}


int DiffLength::sign() {
  return (i->getP().dot(i->getP()) - j->getP().dot(i->getP())).sign();
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

int SegmentIntersect::sign() {
  PV3 a = pa->getP();
  PV3 b = pb->getP();
  PV3 c = pc->getP();
  PV3 d = pd->getP();

  PV3 n  = (b-a).cross(d-c);
  if ((n.dot(a) - n.dot(c)).sign() != 0) //don't lie in the same plane
    return 0; //don't intersect

  PV3 nab = (b-a).cross(n);
  PV3 ncd = (d-c).cross(n);

  if ( (nab.dot(d-a).sign() != nab.dot(c-a).sign()) && 
       (ncd.dot(a-c).sign() != ncd.dot(b-c).sign()) ) 
    return 1;
  return 0;
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

PlaneData SplitPlane::calculate () {
  PV3 a  = v0->getP();
  PV3 b  = v1->getP();
  PV3 c  = v2->getP();
  PV3 n_tri = (b - a).cross(c - a); 
  Parameter k_tri = n_tri.dot(a); 
  PV3 p = k_tri / (n_tri.dot(n_tri)) * n_tri;

  PV3 n = PV3(n_tri.getY(), -n_tri.getX(), Parameter::constant(0));
  Parameter k = n.dot(p);
  return(PlaneData(n,k));
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

PV3 ZIntercectPoint::calculate () {
  Parameter t = (c->getP().getZ() - a->getP().getZ()) / (b->getP().getZ() - a->getP().getZ());
  return a->getP() + t * (b->getP() - a->getP());
}

PV3 FacePoint::calculate () {
  HFace * hface = cell->getBoundary(0)->getHFaces()[0];
  HEdge * hedge = hface->getF()->getBoundary(0);
  
  PV3 h = hedge->head()->getP()->getP();
  PV3 t = hedge->tail()->getP()->getP();
  PV3 n = hface->getN();

  PV3 mid = t + (h-t)/2;
  PV3 inward_normal = n.cross(h-t);

  return mid + unit * (inward_normal); 
}

PV3 CellInternalPoint::calculate () {
  HFace * hface = cell->getBoundary(0)->getHFaces()[0];
  PV3 fp = facePoint->getP();
  PV3 n = hface->getN();
  return fp - unit * n;
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

