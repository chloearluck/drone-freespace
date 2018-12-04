#include "polyhedron.h"
#include "io.h"
#include <cstring>
#include "simplify.h"
#include "hull.h"
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

Polyhedron * ell(double lx, double ly, double r, double h) {
  int flipX = 1, flipY = 1, flipZ = 1;
  if (lx < 1) {
    lx = -lx; flipX = -1;
  }
  if (ly < 1) {
    ly = -ly; flipY = -1;
  }
  flipZ = (flipX == flipY)? 1 : -1;

  Polyhedron * poly = new Polyhedron(true);
  Vertex * ta = poly->getVertex(new Point(flipX * (lx+r), flipY * (   r), flipZ *  h/2  ));
  Vertex * tb = poly->getVertex(new Point(flipX * (lx+r), flipY * (  -r), flipZ *  h/2  ));
  Vertex * tc = poly->getVertex(new Point(flipX * (  -r), flipY * (  -r), flipZ *  h/2  ));
  Vertex * td = poly->getVertex(new Point(flipX * (   r), flipY * (   r), flipZ *  h/2  ));
  Vertex * te = poly->getVertex(new Point(flipX * (  -r), flipY * (ly+r), flipZ *  h/2  ));
  Vertex * tf = poly->getVertex(new Point(flipX * (   r), flipY * (ly+r), flipZ *  h/2  ));
  Vertex * ba = poly->getVertex(new Point(flipX * (lx+r), flipY * (   r), flipZ * -h/2  ));
  Vertex * bb = poly->getVertex(new Point(flipX * (lx+r), flipY * (  -r), flipZ * -h/2  ));
  Vertex * bc = poly->getVertex(new Point(flipX * (  -r), flipY * (  -r), flipZ * -h/2  ));
  Vertex * bd = poly->getVertex(new Point(flipX * (   r), flipY * (   r), flipZ * -h/2  ));
  Vertex * be = poly->getVertex(new Point(flipX * (  -r), flipY * (ly+r), flipZ * -h/2  ));
  Vertex * bf = poly->getVertex(new Point(flipX * (   r), flipY * (ly+r), flipZ * -h/2  ));
  poly->addTriangle(ta, td, tb);
  poly->addTriangle(td, tc, tb);
  poly->addTriangle(td, te, tc);
  poly->addTriangle(td, tf, te);
  poly->addTriangle(ba, bb, bd);
  poly->addTriangle(bd, bb, bc);
  poly->addTriangle(bd, bc, be);
  poly->addTriangle(bd, be, bf);
  poly->addTriangle(tc, bc, bb);
  poly->addTriangle(tc, bb, tb);
  poly->addTriangle(tb, bb, ba);
  poly->addTriangle(tb, ba, ta);
  poly->addTriangle(ta, ba, bd);
  poly->addTriangle(ta, bd, td);
  poly->addTriangle(td, bd, bf);
  poly->addTriangle(td, bf, tf);
  poly->addTriangle(tf, bf, be);
  poly->addTriangle(tf, be, te);
  poly->addTriangle(te, be, bc);
  poly->addTriangle(te, bc, tc);
  return poly;
}

void ellPath() {
  double l = 0.5, r = 0.05, h = 0.4;
  Polyhedron * robot = ell(l, l, r, h);
  savePoly(robot, "ell.vtk");

  double room_r = 1.25;
  double hall_h = 2*l;
  double hall_w = 0.7*l + 2.1*r;
  double theta = M_PI / 20;

  cout << "fattest point = " << sqrt(2.0)*(l+r) << endl;
  cout << "thinnest point = " << (l+3*r)/sqrt(2.0) << endl;
  cout << "parallel side = " << l+r << endl;
  cout << "hall width = " << hall_w << endl;

  double d[6] = {-room_r, room_r, -room_r, room_r, -room_r, room_r};
  Polyhedron * room = box(d);

  Polyhedron * hall1 = ell(-3.0, -3.0, hall_w/2, hall_h/2);
  Polyhedron * hall2 = ell(-3.0,  3.0, hall_w/2, hall_h/2);

  Polyhedron * room1 = room->translate(new Point(0.0, 0.0, 0.0));
  Polyhedron * room2 = room->translate(new Point(3.0, -3.0, 0.0));
  Polyhedron * room3 = room->translate(new Point(6.0, 0.0, 0.0));
  hall1 = hall1->translate(new Point(3.0, 0.0, 0.5));
  hall2 = hall2->translate(new Point(6.0, -3.0, -0.5));

  Polyhedron * polys[5] = {room1, room2, room3, hall1, hall2};
  Polyhedron * inner = multiUnion(polys, 5);
  savePoly(inner, "ellroom-inner.vtk");

  double d2[6] = {-1.5 , 8, -5, 2, -1.5, 1.5};
  Polyhedron * outer = box(d2);
  outer = outer->boolean(inner, Complement);
  simplify(outer, 1e-6);
  savePoly(outer, "ellroom.vtk");
}

void corkscrewPath() {
  double dtheta = M_PI / 15.0; //angle change at each step
  double dh = 0.65; //height change at each step
  double radius = 3.0; //radius of the corkscrew 
  double radiusx = 1.2; //x thickness of the corkscrew thread
  double radiusy = 1.7; //y thickness of the corkscrew thread
  int nsteps = 20;

  double theta0 = atan(radiusy / radiusx);
  double theta1 = -theta0;
  double theta2 = M_PI + theta0;
  double theta3 = M_PI - theta0;
  double r = sqrt(radiusx*radiusx + radiusy*radiusy);

  Polyhedron * corkscrew = new Polyhedron(true);
  Vertex * verts[nsteps][4];

  double theta=0.0, height = -dh * (nsteps/2);
  for (int i=0; i<nsteps; i++) {
    double x = radius * sin(theta);
    double y = radius * cos(theta);

    verts[i][0] = corkscrew->getVertex(new Point(x + r*cos(theta0), y + r*sin(theta0), height));
    verts[i][1] = corkscrew->getVertex(new Point(x + r*cos(theta1), y + r*sin(theta1), height));
    verts[i][2] = corkscrew->getVertex(new Point(x + r*cos(theta2), y + r*sin(theta2), height));
    verts[i][3] = corkscrew->getVertex(new Point(x + r*cos(theta3), y + r*sin(theta3), height));

    theta += dtheta;
    height += dh;
    theta0 -= dtheta;
    theta1 -= dtheta;
    theta2 -= dtheta;
    theta3 -= dtheta;
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

  Polyhedron * polys[3] = {room, room2, corkscrew};
  Polyhedron * all = multiUnion(polys, 3);
  savePoly(all, "corkscrew-inner.vtk");

  Polyhedron * obstacle = outer->boolean(all, Complement); 
  savePoly(obstacle, "corkscrew.vtk");
}

int main (int argc, char *argv[]) {
  Parameter::enable();
  ellPath();
}