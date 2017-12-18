// acp.h provides the acp::Parameter (numerical) class.
#include "acp.h"
// pv.h provides basic 2D and 3D vector of Parameter called acp::PV2 and
// acp::PV3 with the usual vector operations.
#include "pv.h"
// object.h provides acp::Object, acp::Primitive and macros to create
// primitives.
#include "object.h"
// Actually, object.h includes pv.h which includes acp.h, so it is
// sufficient to include object.h
#include "polyhedron.h"

using namespace acp;


// ABSTRACT GEOMETRIC OBJECTS

// class PPlane {
// public:
//   PV3 v;
//   Parameter d;

//   PPlane () {}
//   PPlane (PV3 v, Parameter d) : v(v), d(d) {}

//   // Number of Parameter instances inside this class.
//   int size () { return 4; }

//   // LValue access to the i'th Parameter.
//   Parameter &operator[] (int i) { return i < 3 ? v[i] : d; }
// };

// // Abstract plane class.
// class Plane : public Object<PPlane> {
// public:
//   PV3 getV () { return get().v; }
//   Parameter getD () { return get().d; }
// }; //Plane and PPlane
// Abstract Point class.
// class Point : public Object<PV3> {
// public:
//   PV3 getP () { return get(); }
// };
 //Point

// PRIMITIVES

// Here is our first primitive class.
// Example invocation:  "if (PlaneSide(plane, point) > 0)".
// Sorry this is a clunky macro, requiring commas between types and variables.
// int PlaneSide::side() must be defined in geometry3d.cc.
Primitive2(PlaneSide, Plane*, plane, Point*, point);
Primitive2(DiffLength, Point*, i, Point*, j);
Primitive3(IsTangent, Point*, v0, Point*, v1, Point*, v2);
// Orientation of tetrahedron.
Primitive4(Orient3D, Point*, a, Point*, b, Point*, c, Point*, d);
Primitive4(SegmentIntersect, Point*, pa, Point*, pb, Point*, pc, Point*, pd);
// Insphere test.
Primitive5(InSphere, Point*, a, Point*, b, Point*, c, Point*, d, Point*, e);


// CONCRETE GEOMETRIC OBJECTS

//Points defined by rotating point theta radians around the origin
class RotationPoint : public Point {
  PTR<Point> point, sin_cos_alpha;
  PV3 calculate ();

public:
  RotationPoint(PTR<Point> point, PTR<Point> sin_cos_alpha) : point(point), sin_cos_alpha(sin_cos_alpha) {}
};

class ScalePoint : public Point {
  PTR<Point> p;
  PTR<Object<Parameter> > unit;
  PV3 calculate () { return unit->get() * p->getP(); }
 public:
  ScalePoint(PTR<Point> p, PTR<Object<Parameter> > unit) : p(p), unit(unit) {}
};

//point defined by the intersection of point1 and point2's tangent lines in the xy plane
class TangentIntersectionPoint : public Point {
  PTR<Point> point, sin_cos_alpha;
  PV3 calculate ();

public:
  TangentIntersectionPoint(PTR<Point> point, PTR<Point> sin_cos_alpha) : point(point), sin_cos_alpha(sin_cos_alpha) {}
};

//point along ab which has the same z component as c
class ZIntercectPoint : public Point {
  PTR<Point> a;
  PTR<Point> b;
  PTR<Point> c;
  PV3 calculate ();

public:
  ZIntercectPoint(PTR<Point> a, PTR<Point> b, PTR<Point> c) : a(a), b(b), c(c) {}
};

class FacePoint : public Point {
  Cell * cell;
  double unit;
  PV3 calculate ();

public:
  FacePoint(Cell * cell, double unit) : cell(cell), unit(unit) {}
};

class CellInternalPoint : public Point {
  Cell * cell;
  PTR<Point> facePoint;
  double unit;
  PV3 calculate ();

public:
  CellInternalPoint(Cell * cell, PTR<Point> facePoint, double unit) : cell(cell), facePoint(facePoint), unit(unit) {}
};

class SplitPlane : public Plane {
  PTR<Point> v0, v1, v2;
  PlaneData calculate ();

public:
  SplitPlane(PTR<Point> v0, PTR<Point> v1, PTR<Point> v2) : v0(v0), v1(v1), v2(v2) {}
};

class PointNormalPlane : public Plane {
  ObjPTR<PV3> point, normal;
  PlaneData calculate ();

public:
  PointNormalPlane(PTR<Point> point, PTR<Point> normal) : point(point), normal(normal) {}
};

class IntersectionPoint : public Point {
  PV3 calculate ();

protected:
  PTR<Point> tail;
  PTR<Point> head;
  PTR<Plane> plane;
  
public:
  IntersectionPoint (PTR<Point> tail, PTR<Point> head, PTR<Plane> plane)
    : tail(tail), head(head), plane(plane) {}
};

class FaceIntersectionPoint : public Point {
  PV3 calculate ();

protected:
  PTR<Point> head, tail;
  PTR<Point> a,b,c;
  HFace * hface;

public:
  FaceIntersectionPoint (PTR<Point> tail, PTR<Point> head, HFace * hface) : tail(tail), head(head), hface(hface) { a = hface->getF()->getP()->getA(); b = hface->getF()->getP()->getB(); c = hface->getF()->getP()->getC(); }
  HFace * getHFace() { return hface; }
};

//assume face is triangular
class FaceNearestPoint : public Point {
  PV3 calculate ();

protected:
  PTR<Point> point;
  PTR<Point> pa,pb,pc;

public: 
  FaceNearestPoint (PTR<Point> point, HFace * hf) : point(point) { pa = hf->getF()->getP()->getA(); pb = hf->getF()->getP()->getB(); pc = hf->getF()->getP()->getC(); }
};