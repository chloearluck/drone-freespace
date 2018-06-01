#include "hull.h"
#include "geometry3d.h"
#include "mink.h"
#include "io.h"
#include <cstring>

class InputParameter : public Object<Parameter> {
public:
  InputParameter (double x) { set(Parameter::input(x)); }
};

class SinCosAlpha : public Point {
  Object<Parameter> *tan_theta;
  PV3 calculate () {
    Parameter t = tan_theta->get();
    Parameter sint = 2*t/(1+t*t);
    Parameter cost = (1-t*t)/(1+t*t);
    Parameter alpha = (1-cost)/sint;
    return PV3(sint, cost, alpha);
  }
public:
  SinCosAlpha (Object<Parameter> *t) : tan_theta(t) {}
};

Polyhedron * loadPoly(const char * filename);
void savePoly(Polyhedron * p, const char * filename);

class FreeSpace {
public:
  Polyhedron * robot, * obstacle;
  std::vector<Polyhedron*> blockspaces;
  int numRotations;
  FreeSpace(Polyhedron * robot, Polyhedron * obstacle, PTR<Object<Parameter> > tan_half_angle, int numRotations, bool inner_approximation = true);
};

