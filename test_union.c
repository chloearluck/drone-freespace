#include "hull.h"
#include "geometry3d.h"



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

void savePoly(Polyhedron * p, char * str) {
  ofstream out;
  out.open(str);
  if (out.is_open()) {
    writePolyhedronVTK (p->faces, out);
    out.close();
  } else {
    printf("could not write to file\n");
  }
}

// //take the union of 2 polyhedrons given at vtk files
// //first 2 args are name of vtk files, 3rd is random seed
// int main (int argc, char *argv[]) {
//   if (argc < 3) {
//     printf("not enough arguments\n");
//     return 1;
//   }

//   char * filename1 =  argv[1];
//   char * filename2 =  argv[2];

//   if (argc == 3) { 
//     unsigned seed = atoi(argv[3]);
//     srandom(seed);
//   }

//   Polyhedron * poly1 = loadPoly(filename1);
//   if (poly1 == NULL)
//     return 1;
//   Polyhedron * poly2 = loadPoly(filename2);
//   if (poly2 == NULL)
//     return 1;

//   Polyhedron * out = poly1->boolean(poly2, Union);

//   savePoly(out, "out.vtk");
// }


Polyhedron * union_all(std::vector<Polyhedron*> pList, int start, int end) {
  printf("union_all %d %d\n", start, end);
  if (start == end) 
    return pList[start];
  if ((start+1) == end) 
    return pList[start]->boolean(pList[end], Union);
  int mid = (start+end)/2;
  Polyhedron * p1 = union_all(pList, start, mid);
  Polyhedron * p2 = union_all(pList, mid+1, end);
  printf("taking union of (%d,%d) and (%d,%d)\n", start, mid, mid+1, end);
  return p1->boolean(p2, Union);
}

Polyhedron * union_all(std::vector<Polyhedron*> pList) {
  return union_all(pList, 0, pList.size()-1);
}

double randomInRange(double low, double high) {
  return (random() * (high-low)) + low;
}

int main (int argc, char *argv[]) {
  if (argc == 2) { 
    unsigned seed = atoi(argv[1]);
    srandom(seed);
  }

  bool generateSimple = true;

  if (generateSimple) {
    Points pList;
    for (int i=0; i<4; i++) {
      pList.push_back(new InputPoint(randomInRange(-1,1), randomInRange(-1,1), randomInRange(-1,1)));
    }

    Polyhedron * simple = convexHull(pList);
    savePoly(simple, "simple.vtk");
  } else {
    std::vector<Polyhedron *> pList;
    for (int i=0; i<168; i++) {
      char s[30];
      sprintf(s, "output/%d-out.vtk", i);
      Polyhedron * poly = loadPoly(s);
      pList.push_back(poly);
    }

    Polyhedron * p = union_all(pList);

    savePoly(p, "out.vtk");
  }
}