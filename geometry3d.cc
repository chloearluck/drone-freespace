#include "geometry3d.h"

int PlaneSide::sign ()
{
  PV3 v = plane->get().n, p = point->get();
  Parameter d = plane->get().k;
  return (v.dot(p) - d).sign();
}

int DiffLength::sign()
{
  return (i->get().dot(i->get()) - j->get().dot(j->get())).sign();
}

int IsTangent::sign ()
{
  PV3 a = v0->get(), b = v1->get(), c = v2->get(),
    n = (b - a).cross(c - a), nsplit(n.y, - n.x, Parameter::constant(0.0));
  Parameter k = n.dot(a);
  int psa = nsplit.dot(a).sign(), psb = nsplit.dot(b).sign(),
    psc = nsplit.dot(c).sign();
  return psa == psb && psb == psc ? 0 : 1;
}

int Orient3D::sign ()
{
  PV3 dd = d->get(), da = a->get() - dd, db = b->get() - dd,
    dc = c->get() - dd;
  return da.cross(db).dot(dc).sign();
}

int SegmentIntersect::sign()
{
  PV3 a = pa->get(), b = pb->get(), c = pc->get(), d = pd->get(),
    n = (b - a).cross(d - c);
  if (n.dot(a - c).sign() != 0)
    return 0;
  PV3 nab = (b - a).cross(n), ncd = (d - c).cross(n);
  if (nab.dot(d - a).sign() != nab.dot(c - a).sign() && 
      ncd.dot(a - c).sign() != ncd.dot(b - c).sign()) 
    return 1;
  return 0;
}

int InSphere::sign ()
{
  PV3 ee = e->get(), ea = a->get() - ee, eb = b->get() - ee,
    ec = c->get() - ee, ed = d->get() - ee;
  return (ea.cross(eb).dot(ec)*ed.dot(ed) -
          eb.cross(ec).dot(ed)*ea.dot(ea) +
          ec.cross(ed).dot(ea)*eb.dot(eb) -
          ed.cross(ea).dot(eb)*ec.dot(ec)).sign();
}
