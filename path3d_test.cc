#include "path3d.h"
#include "io.h"

Polyhedron * loadPoly(const char * filename) {
  int n = strlen(filename);
  char str[n+5];
  strncpy(str, filename, n);
  strncpy(str+n, ".vtk", 5);

  Polyhedron * poly;
  ifstream infile (str);
  if (infile.is_open()) {
    poly = readPolyhedronVTK (infile);
    infile.close();
  } else {
    cout<<"could not read from file"<<endl;
    return NULL;
  }

  return poly;
}

void savePathVTK(Points path, const char * filename) {
  ofstream ostr;
  ostr.open(filename);
  if (ostr.is_open()) { 
    ostr << setprecision(20) << "# vtk DataFile Version 3.0" << endl
         << "vtk output" << endl << "ASCII" << endl
         << "DATASET POLYDATA" << endl 
         << "POINTS " << path.size() << " double" << endl;
    for (int i=0; i<path.size(); i++)
      ostr << path[i]->getApprox().getX().mid() << " " << path[i]->getApprox().getY().mid() << " " << path[i]->getApprox().getZ().mid() << endl;
    ostr<<endl<<"LINES 1 "<<path.size()+1<<endl<<path.size()<<" ";
    for (int i=0; i<path.size(); i++)
      ostr<<i<<" ";
    ostr.close();
  } else {
    cout<<"could not write to file"<<endl;
  }
}

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout<<"not enough arguments"<<endl;
    return 1;
  }

  char * filename =  argv[1];

  if (argc == 3) {
    unsigned seed = atoi(argv[2]);
    srandom(seed);
  }

  Polyhedron * blockspace = loadPoly(filename);

  PTR<Point> start = new InputPoint(5, 5, 2);
  PTR<Point> end = new InputPoint(15, 15, 2);
  
  Points path;

  findPath(blockspace, start, end, path);
  savePathVTK(path, "path.vtk");
}
