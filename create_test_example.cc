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

  Polyhedron * poly = new Polyhedron();
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

  Polyhedron * corkscrew = new Polyhedron();
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

Polyhedron * get_wire(double x1, double y1, double z1, double x2, double y2, double z2, double wire_r) {
  Polyhedron * p = new Polyhedron();
  
  if (z2 > z1) {
    double x = x1, y = y1, z = z1;
    x1 = x2; y1 = y2; z1 = z2;
    x2 = x;  y2 = y;  z2 = z;
  }

  Vertex * t0 = p->getVertex(new Point(x1-wire_r, y1-wire_r, z1));
  Vertex * t1 = p->getVertex(new Point(x1+wire_r, y1-wire_r, z1));
  Vertex * t2 = p->getVertex(new Point(x1+wire_r, y1+wire_r, z1));
  Vertex * t3 = p->getVertex(new Point(x1-wire_r, y1+wire_r, z1));
  Vertex * b0 = p->getVertex(new Point(x2-wire_r, y2-wire_r, z2));
  Vertex * b1 = p->getVertex(new Point(x2+wire_r, y2-wire_r, z2));
  Vertex * b2 = p->getVertex(new Point(x2+wire_r, y2+wire_r, z2));
  Vertex * b3 = p->getVertex(new Point(x2-wire_r, y2+wire_r, z2));

  p->addTriangle(b0, t1, t0);
  p->addTriangle(b0, b1, t1);
  p->addTriangle(b1, t2, t1);
  p->addTriangle(b1, b2, t2);
  p->addTriangle(b2, t3, t2);
  p->addTriangle(b2, b3, t3);
  p->addTriangle(b3, t0, t3);
  p->addTriangle(b3, b0, t0);

  p->addTriangle(t0, t1, t2);
  p->addTriangle(t0, t2, t3);
  p->addTriangle(b0, b2, b1);
  p->addTriangle(b0, b3, b2);

  return p;
}

Polyhedron * tetrahedron(double A) {
  Polyhedron * p = new Polyhedron();
  Parameter::disable(); double root3 = sqrt(3); double root6 = sqrt(6); Parameter::enable();

  // Vertex * a = p->getVertex(new Point(root3*A/3, 0, -root3*A/3));
  // Vertex * b = p->getVertex(new Point(-root3*A/6, A/2, -root3*A/3));
  // Vertex * c = p->getVertex(new Point(-root3*A/6, -A/2, -root3*A/3));
  // Vertex * d = p->getVertex(new Point(0, 0, root3*A/3));

  // p->addTriangle(b, a, d);
  // p->addTriangle(a, c, d);
  // p->addTriangle(c, b, d);
  // p->addTriangle(a, b, c);

  Vertex * a = p->getVertex(new Point(root3*A/3, -root6*A/9, 0));
  Vertex * b = p->getVertex(new Point(-root3*A/6, -root6*A/9, A/2));
  Vertex * c = p->getVertex(new Point(-root3*A/6, -root6*A/9, -A/2));
  Vertex * d = p->getVertex(new Point(0, root6*A/4.5, 0));

  p->addTriangle(a,d,b);
  p->addTriangle(a,c,d);
  p->addTriangle(a,b,c);
  p->addTriangle(d,c,b);

  return p;
}

Polyhedron * get_spine(double x, double y, double z, double wire_r) {
  Polyhedron * poly = new Polyhedron();
  Vertex * p = poly->getVertex(new Point(x, y, z));
  Vertex * a = poly->getVertex(new Point(0, 0, 0));
  Vertex *b,*c;
  if (z > 0) {
    b = poly->getVertex(new Point(wire_r, 0, 0));
    c = poly->getVertex(new Point(0, wire_r, 0));
  } else {
    b = poly->getVertex(new Point(0, wire_r, 0));
    c = poly->getVertex(new Point(wire_r, 0, 0));
  }

  poly->addTriangle(a,b,p);
  poly->addTriangle(b,c,p);
  poly->addTriangle(c,a,p);
  poly->addTriangle(a,c,b);
  return poly;
}

void spiny() {
  double wire_r = 1e-3;
  double box_r = 2;
  int n_wires = 10;
  double box_length = 6;
  double box_thickness = 0.05;

  Polyhedron * wires = new Polyhedron();
  // for (int i=0; i<n_wires; i++) {
  //   double x[2], y[2], z[2];
  //   for (int j=0; j<2; j++) {
  //     x[j] = -box_r + 2*box_r*random()/double(RAND_MAX);
  //     if (random()/double(RAND_MAX) < 0.5) {
  //       //floor or ceiling point 
  //       y[j] = -box_r + 2*box_r*random()/double(RAND_MAX);
  //       z[j] = (random()/double(RAND_MAX) < 0.5? -box_r : box_r);
  //     } else {
  //       //wall point
  //       y[j] = (random()/double(RAND_MAX) < 0.5? -box_r : box_r);
  //       z[j] = -box_r + 2*box_r*random()/double(RAND_MAX);
  //     }
  //   }

  //   if (y[0] == y[1] || z[0] == z[1]) {
  //     i--; continue;
  //   }
  //   cout << "get_wire("<<x[0]<<", "<<y[0]<<", "<<z[0]<<", "<<x[1]<<", "<<y[1]<<", "<<z[1]<<", wire_r);"<<endl;
  //   Polyhedron * wire = get_wire(x[0], y[0], z[0], x[1], y[1], z[1], wire_r);
  //   wires = wires->boolean(wire, Union);
  // }
  // cout << "-----------" << endl;

  // saved wires
  Polyhedron * wire;
  wire = get_wire(0.245521, -0.427632, -2, -0.859834, 0.254218, 2, wire_r); wires = wires->boolean(wire, Union); delete wire;
  wire = get_wire(-0.129744, -2, -1.59059, 0.45452, 2, 1.61083, wire_r); wires = wires->boolean(wire, Union); delete wire;
  wire = get_wire(0.240393, 2, 0.183991, -0.713537, -2, 0.732517, wire_r); wires = wires->boolean(wire, Union); delete wire;
  wire = get_wire(-1.84819, -0.0252306, -2, 0.721769, 2, 1.77568, wire_r); wires = wires->boolean(wire, Union); delete wire;
  wire = get_wire(1.82971, 0.383239, 2, 0.171406, 2, 0.665213, wire_r); wires = wires->boolean(wire, Union); delete wire;
  wire = get_wire(-0.557526, 2, 0.83767, 0.906497, 0.224287, -2, wire_r); wires = wires->boolean(wire, Union); delete wire;
  wire = get_wire(-0.597933, -2, 1.91441, 1.34802, 1.36169, -2, wire_r); wires = wires->boolean(wire, Union); delete wire;
  wire = get_wire(-1.55025, 2, -1.36723, -0.749517, 1.18442, -2, wire_r); wires = wires->boolean(wire, Union); delete wire;
  wire = get_wire(-0.851031, 2, -1.14023, 0.806354, -1.55504, -2, wire_r); wires = wires->boolean(wire, Union); delete wire;
  wire = get_wire(1.4719, -1.89553, -2, -1.11516, 2, 0.506067, wire_r); wires = wires->boolean(wire, Union); delete wire;

  savePoly(wires, "wires.vtk");

  double d1[6] = {-box_length - box_thickness, box_length + box_thickness, -box_r - box_thickness, box_r + box_thickness, -box_r - box_thickness, box_r + box_thickness};
  double d2[6] = {-box_length + box_thickness, box_length - box_thickness, -box_r + box_thickness, box_r - box_thickness, -box_r + box_thickness, box_r - box_thickness};
  Polyhedron * outer = box(d1);
  Polyhedron * inner = box(d2);
  Polyhedron * cage = outer->boolean(inner, Complement);
  savePoly(cage, "box.vtk");
  Polyhedron * obstacle = cage->boolean(wires, Union);
  savePoly(obstacle, "wirecage.vtk");

  //construction spiny robot
  //question: box or frustum, is frustum better for outerApprox?  box for now
  double robot_r = 0.2;
  double spine_length = 1.5;
  int n_spines = 3;
  double d3[6] = {-robot_r, robot_r, -robot_r, robot_r, -robot_r, robot_r};
  Polyhedron * robot = box(d3);
  robot = tetrahedron(robot_r);

  // for (int i=0; i<n_spines; i++) {
  //   double x = -1 + 2*random()/double(RAND_MAX);
  //   double y = -1 + 2*random()/double(RAND_MAX);
  //   double z = -1 + 2*random()/double(RAND_MAX);
  //   Parameter::disable();
  //   double l = sqrt(x*x + y*y + z*z);
  //   Parameter::enable();
  //   x = x/l; y = y/l; z = z/l;
  //   cout << "get_wire(0.0, 0.0, 0.0, spine_length*"<<x<<", spine_length*"<<y<<", spine_length*"<<z<<", wire_r);"<<endl;
  //   Polyhedron * spine = get_wire(0.0, 0.0, 0.0, x, y, z, wire_r);
  //   robot = robot->boolean(spine, Union);
  // }

  // saved spines
  Polyhedron * spine;
  spine = get_spine(spine_length*0.185313, spine_length*-0.510148, spine_length*0.839886, wire_r); robot = robot->boolean(spine, Union); delete spine;
  spine = get_spine(spine_length*0.0344939, spine_length*0.030168, spine_length*-0.998949, wire_r); robot = robot->boolean(spine, Union); delete spine;
  spine = get_spine(spine_length*0.772395, spine_length*0.28215, spine_length*0.569034, wire_r); robot = robot->boolean(spine, Union); delete spine;
  savePoly(robot, "spiny.vtk");
}

int main (int argc, char *argv[]) {
  Parameter::enable();
  if (argc == 2) { 
    unsigned seed = atoi(argv[1]);
    srandom(seed);
  }

  spiny();
}