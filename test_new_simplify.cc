#include "polyhedron.h"
#include "io.h"
#include <cstring>
#include "simplify.h"
#include "mink.h"

Polyhedron * loadPoly(const char * filename) {
  Polyhedron * poly;
  ifstream infile (filename);
  if (infile.is_open()) {
    poly = readPolyhedronVTK (infile);
    infile.close();
  } else {
    cout<<"could not read from file"<<endl;
    return NULL;
  }

  return poly;
}

void savePoly(Polyhedron * p, const char * filename) {
  ofstream out;
  out.open(filename);
  if (out.is_open()) {
    writePolyhedronVTK (p->faces, out);
    out.close();
  } else {
    cout<<"could not write to file"<<endl;
  }
}

int main (int argc, char *argv[]) {
  Parameter::enable();

  Polyhedron * ball = loadPoly("sphere.vtk");
  Polyhedron * unit_ball = ball->scale(0.05);
  Polyhedron * poly = loadPoly("sum00-out.vtk");
  
  cout<<"without simplify: "<<endl;
  poly->formCells();
  cout<<"poly has "<<poly->cells.size()<<" cells"<<endl;
  Polyhedron * sum = minkowskiSumFull(poly, unit_ball);
  sum->formCells();
  cout<<"sum has "<<sum->cells.size()<<" cells"<<endl;
  savePoly(sum, "nosimplify.vtk");

  cout<<"with simplify"<<endl;
  simplify(poly, 1e-6);

  poly->formCells();
  cout<<"poly has "<<poly->cells.size()<<" cells"<<endl;
  sum = minkowskiSumFull(poly, unit_ball);
  sum->formCells();
  cout<<"sum has "<<sum->cells.size()<<" cells"<<endl;
  savePoly(sum, "simplify.vtk");  
}