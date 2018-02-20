#include "polyhedron.h"
#include "io.h"
#include <cstring>
#include "simplify.h"
#include "hull.h"

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

Primitive2(DiffZ, PTR<Point>, i, PTR<Point>, j);
int DiffZ::sign() { return (i->getP().getZ() - j->getP().getZ()).sign(); }
struct CompareZ {
  bool operator()(PTR<Point> i, PTR<Point> j) {
    return (DiffZ(i, j) < 0); 
  }
};

int main (int argc, char *argv[]) {
  Parameter::enable();

  PTR<Point> top1 = new InputPoint(-1,-1, 1); //top left
  PTR<Point> top2 = new InputPoint(-1, 1, 1); //top right
  PTR<Point> top3 = new InputPoint( 1, 1, 1); //bottom right
  PTR<Point> top4 = new InputPoint( 1,-1, 1); //bottom left

  PTR<Point> bottom1 = new InputPoint(-1,-1,-1); //top left
  PTR<Point> bottom2 = new InputPoint(-1, 1,-1); //top right
  PTR<Point> bottom3 = new InputPoint( 1, 1,-1); //bottom right
  PTR<Point> bottom4 = new InputPoint( 1,-1,-1); //bottom left

  //find the midpoint of the edges you want to split
  PTR<Point> mid1 = new MidPoint(top1, top4);
  PTR<Point> mid2 = new MidPoint(top1, top2);
  //add a (0,0,.000001) to them
  PTR<Point> p1 = new SumPoint(mid1, new InputPoint(0,0,0.0000001));
  PTR<Point> p2 = new SumPoint(mid2, new InputPoint(0,0,0.0000001));

  Points plist;
  plist.push_back(top1);
  plist.push_back(top3);
  plist.push_back(top4);
  plist.push_back(p1);
  plist.push_back(bottom1);
  plist.push_back(bottom2);
  plist.push_back(bottom3);
  plist.push_back(bottom4);
  Polyhedron * hull1 = convexHull(plist, true);
  assert(hull1 != NULL);

  plist.clear();
  plist.push_back(top1);
  plist.push_back(top2);
  plist.push_back(top3);
  plist.push_back(p2);
  plist.push_back(bottom1);
  plist.push_back(bottom2);
  plist.push_back(bottom3);
  plist.push_back(bottom4);
  Polyhedron * hull2 = convexHull(plist, true);
  assert(hull2 != NULL);

  Polyhedron * cube = hull1->boolean(hull2, Union);

  savePoly(cube, "split-faces-cube.vtk");

}