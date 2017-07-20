#include "freespace.h"
#include <queue>

Polyhedron * loadPoly(char * filename) {
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

void bfs(FreeSpace * fs, FreeSpace::Node * start, FreeSpace::Node * end) {
  if (start == NULL || end == NULL)
    return;

  queue<FreeSpace::Node *> q;
  start->discovered = true;
  q.push(start);
  while(!q.empty()) {
    FreeSpace::Node * current = q.front();   q.pop();
    if (current == end) {
      cout<<"found path"<<endl;
      return;
    }
    for (int i=0; i<current->edges.size(); i++) {
      FreeSpace::Edge * e = current->edges[i];
      FreeSpace::Node * n = ((e->a == current)? e->b : e->a);
      if (!n->discovered) {
        n->discovered = true;
        n->parent = current;
        q.push(n);
      }
    }
  }
  cout<<"no path found"<<endl;
}

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout<<"not enough arguments"<<endl;
    return 1;
  }

  char * filename =  argv[1];

  // if (argc == 3) { 
    // unsigned seed = atoi(argv[2]);
    srandom(10);
  // }

  Polyhedron * poly = loadPoly(filename);
  if (poly == NULL)
    return 1;

  Polyhedron * obstacle;
  if (argc == 3) 
    obstacle = loadPoly(argv[2]); 
  else {
    double bounds[6] = { 5, 7, 5, 7, 5, 7};
    obstacle = box(bounds);
  }
  
  double theta = M_PI / 10;
  double bb_bounds[6] = {-20, 20, -20, 20, -20, 20};

  FreeSpace * fs  = new FreeSpace(poly, obstacle, theta, bb_bounds);

  cout<<"Finding path"<<endl;
  PTR<Point> startp = new InputPoint(5, 5, 2);
  PTR<Point> endp = new InputPoint(15, 15, 2);
  int starti = fs->cspaces[0]->containingCell(startp);
  int endi = fs->cspaces[0]->containingCell(endp);
  cout<<"need to get from cell "<<starti<<" to cell "<<endi<<" of cspace[0]"<<endl;
  FreeSpace::Node * start =  fs->findNode(0, starti);
  if (start == NULL)
    cout<<"start is null"<<endl;
  FreeSpace::Node * end =  fs->findNode(0, endi);
  if (end == NULL)
    cout<<"end is null"<<endl;

  cout<<"bfs:"<<endl;
  bfs(fs, start, end);

  FreeSpace::Node * current = end;
  while (current != NULL && current->parent != NULL) {
    cout<<"cspaces["<<current->cspace_index<<"] cell "<<current->cell_index<<endl;
    current = current->parent;
  }
}