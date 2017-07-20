//minkowski sum of rect-drone and room is incorrect

#include "hull.h"
#include "geometry3d.h"
#include "mink.h"
#include <cstring>

Polyhedron * loadPoly(char * str) {
  Polyhedron * poly;
  ifstream infile (str);
  if (infile.is_open()) {
    poly = readPolyhedronVTK (infile);
    infile.close();
  } else {
    printf("could not read from file\n");
    return NULL;
  }

  return poly;
}


void savePoly(Polyhedron * p, char * filename) {
  int n = strlen(filename);
  char str[n+9];
  strncpy(str, filename, n);
  strncpy(str+n, "-out.vtk", 9);

  ofstream out;
  out.open(str);
  if (out.is_open()) {
    writePolyhedronVTK (p->faces, out);
    out.close();
  } else {
    cout<<"could not write to file"<<endl;
  }
}

void saveShell(HFaces hf,  char * filename) {
  int n = strlen(filename);
  char str[n+9];
  strncpy(str, filename, n);
  strncpy(str+n, "-out.vtk", 9);

  ofstream out;
  out.open(str);
  if (out.is_open()) {
    writePolyhedronVTK (hf, out);
    out.close();
  } else {
    cout<<"could not write to file"<<endl;
  }
}


void saveWithShells(Polyhedron * poly, char * filename) {
  savePoly(poly, filename);
  poly->formCells();
  for (int i=1; i<poly->cells.size(); i++) {
    Cell * c = poly->cells[i];
    for (int j=0; j<c->nShells(); j++) {
      Shell * s = c->getShell(j);
      char str[50];
      sprintf(str, "%s-%d-%d", filename, i, j);
      saveShell(s->getHFaces(), str);
    }
  }
}

void test_cells(Polyhedron * p, char * s) {
  p->formCells();
  cout<<s<<" has "<<p->cells.size()<<" cells"<<endl;
  for (int i=1; i<p->cells.size(); i++) {
    cout<<"cell "<<i<<endl;
    Cell * c = p->cells[i];
    bool valid = !c->getShell(0)->getHFaces()[0]->pos();
    if (valid) {
      for (int j=0; j<c->nShells(); j++) 
        cout<<j<<": "<<( c->getShell(j)->outer() ? "outer" : "inner")<<endl;
      PTR<Point> q = pointInCell(p, i);
      cout<<"p: "<<q->getP().getX().mid()<<" "<<q->getP().getY().mid()<<" "<<q->getP().getZ().mid()<<endl;
    }
  }
  saveWithShells(p, s);
}

int main (int argc, char *argv[]) {
  if (argc == 2) { 
    unsigned seed = atoi(argv[1]);
    srandom(seed);
  }

  // double b1[6] = { 1, 2, 1, 2, 1, 2};
  // Polyhedron * box1 = box(b1);

  // double b2[6] = { -1, -2, -1, -2, -1, -2};
  // Polyhedron * box2 = box(b2);

  // double b3[6] = { -3, 3, -3, 3, -3, 3};
  // Polyhedron * box3 = box(b3);

  // Polyhedron * disjoint = box1->boolean(box2, Union);
  // Polyhedron * nested = complement(box3, box1);

  // test_cells(box3, "box3");
  // test_cells(nested, "nested");
  // test_cells(disjoint, "disjoint");

  Polyhedron * robot = loadPoly("rect-drone.vtk");
  Polyhedron * room = loadPoly("room.vtk");
  Polyhedron * mSum = minkowskiSum(robot, room);
  savePoly(mSum, "mSum");
  
  return 0;
}