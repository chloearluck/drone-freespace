#include "hull.h"
#include "geometry3d.h"
#include "mink.h"
#include "io.h"
#include <cstring>

class InputParameter : public Object<Parameter> {
public:
  InputParameter (double x) { set(Parameter::input(x)); }
};

Polyhedron * loadPoly(const char * filename);
void savePoly(Polyhedron * p, const char * filename);

class FreeSpace {
public:
  Polyhedron * robot, * obstacle;
  std::vector<Polyhedron*> blockspaces;
  int numRotations;
  FreeSpace(Polyhedron * robot, Polyhedron * obstacle, PTR<Object<Parameter> > tan_half_angle, int numRotations);
};

