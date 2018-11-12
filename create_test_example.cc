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

void corkscrewPath() {
  double dtheta = M_PI / 20.0; //angle change at each step
  double dh = 0.65; //height change at each step
  double radius = 3.0; //radius of the corkscrew 
  double radiusx = 1.3; //x thickness of the corkscrew thread
  double radiusy = 0.85; //y thickness of the corkscrew thread
  int nsteps = 20;

  Polyhedron * corkscrew = new Polyhedron(true);
  Vertex * verts[nsteps][4];

  double theta=0.0, height = -dh * (nsteps/2);
  for (int i=0; i<nsteps; i++) {
    Parameter::disable();
    double x = radius * sin(theta);
    double y = radius * cos(theta);
    Parameter::enable();

    verts[i][0] = corkscrew->getVertex(new Point(x - radiusx, y - radiusy, height));
    verts[i][1] = corkscrew->getVertex(new Point(x - radiusx, y + radiusy, height));
    verts[i][2] = corkscrew->getVertex(new Point(x + radiusx, y + radiusy, height));
    verts[i][3] = corkscrew->getVertex(new Point(x + radiusx, y - radiusy, height));

    theta += dtheta;
    height += dh;
  }

  corkscrew->addTriangle(verts[0][0], verts[0][1], verts[0][2]); //top
  corkscrew->addTriangle(verts[0][0], verts[0][2], verts[0][3]);
  for (int i=1; i<nsteps; i++) {
    for (int j=0; j<4; j++) {
      int k= (j+1)%4;
      corkscrew->addTriangle(verts[i-1][k], verts[i-1][j], verts[i][j]);
      corkscrew->addTriangle(verts[i-1][k], verts[i  ][j], verts[i][k]);
    }
  }
  corkscrew->addTriangle(verts[nsteps-1][0], verts[nsteps-1][2], verts[nsteps-1][1]); //bottom
  corkscrew->addTriangle(verts[nsteps-1][0], verts[nsteps-1][3], verts[nsteps-1][2]);

  double d[6] = {-2*radius, 2*radius, -2*radius, 2*radius, -dh * (nsteps/2)-1.0, -dh * (nsteps/2)+1.0};
  Polyhedron * room = box(d);
  Polyhedron * room2 = room->translate(new Point(0.0, 0.0,  dh*nsteps));
  double d2[6] = {-2.1*radius, 2.1*radius, -2.1*radius, 2.1*radius, -dh * (nsteps/2)-1.2, dh * (nsteps/2)+1.2};
  Polyhedron * outer = box(d2);

  Polyhedron * obstacle = outer->boolean(corkscrew, Complement); 
  obstacle = obstacle->boolean(room, Complement);
  obstacle = obstacle->boolean(room2, Complement);
  savePoly(obstacle, "corkscrew.vtk");
}

int main (int argc, char *argv[]) {
  Parameter::enable();
  corkscrewPath();
}