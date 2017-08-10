#include "cspace.h"
#include "io.h"
#include "simplify.h"

int main (int argc, char *argv[])
{
  if (argc < 3)
    return 0;
  ifstream astr(argv[1]), bstr(argv[2]);
  if (!(astr.good() && bstr.good()))
    return 0;
  Parameter::enable();
  Polyhedron *a = readPolyhedronVTK(astr), *b = readPolyhedronVTK(bstr);
  double t0 = getTime();
  Cspace *c = cspace(a, b);
  c->describe();
  t0 = getTime() - t0;
  cerr << "cspace time: " << t0 << endl;
  double dd = 0.01, ds = 1e-6;
  Parameter::usePrecisionException = false;
  t0 = getTime();
  Polyhedron *p = discretize(c, dd);
  t0 = getTime() - t0;
  cerr << "discretize time: " << t0 << endl;
  p->formCells();
  p->describe();
  simplify(p, ds);
  p->formCells();
  p->describe();
  Parameter::usePrecisionException = true;
  delete p;
  delete c;
  delete b;
  delete a;
  Parameter::disable();
}

/*
int main (int argc, char *argv[])
{
  Parameter::enable();
  bool perturb = true;
  double d = 1e-1, x1 = 1, x2 = 21, x3 = 22, y1 = 1, y2 = 10 - d, y3 = 11 + d,
    y4 = 20, y5 = 21, tx = 5, ty = 0, tz = 2;
  // room1: tx = 2, room2: tx = 4, room3: tx = 5
  InputPoint tb = PV3::constant(tx, ty, tz), ta = PV3::constant(-11, -11, 0);
  double bb1[] = {10, 12, 10.5, 11.5, 0, 1.51}, bb2[] = {10.5, 11.5, 10, 12, 1.5, 3};
  Polyhedron *a1 = box(bb1, perturb), *a2 = box(bb2, perturb),
    *a3 = a1->boolean(a2, Union), *a = a3->translate(&ta);
  Polyhedron *b1 = room(x1, x2, x3, y1, y2, y3, y4, y5, 0, 1, perturb),
    *b2 = b1->translate(&tb), *b = b1->boolean(b2, Union);
  delete a1;
  delete a2;
  delete a3;
  delete a;
  delete b1;
  delete b2;
  delete b;
}
*/
