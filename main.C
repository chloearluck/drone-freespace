#include "mink.h"
#include "simplify.h"

int main (int argc, char *argv[])
{
  if (argc < 3)
    return 0;
  ifstream astr(argv[1]), bstr(argv[2]);
  if (!(astr.good() && bstr.good()))
    return 0;
  inputPerturbed = argc == 3 ? true : *argv[3] != 't';
  Parameter::enable();
  Polyhedron *a = readPolyhedronVTK(astr, inputPerturbed),
    *b = readPolyhedronVTK(bstr, inputPerturbed);
  double t0 = getTime();
  Polyhedron *c = minkowskiSum(a, b);
  t0 = getTime() - t0;
  c->formCells();
  c->describe();
  cerr << "Minkowski time " << t0 << endl;
  double s = 1e-6;
  simplify(c, s);
  c->formCells();
  c->describe();
  delete c;
  delete b;
  delete a;
  Parameter::disable();
}
