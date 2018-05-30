#include "polyhedron.h"
#include "io.h"
#include <cstring>
#include "simplify.h"
#include "hull.h"
#include "mink.h"

const double pole_width = 1e-5;
const double width = 6.5, height = 3.5, depth = 3.0; //lattice dimensions
const double delta_h = 1.5;

const double r = 0.8; //radius of robot
const double a = 0.1; //width of each of the robots arms

const bool CONSTRUCT_ROOM = true; 
const bool CONSTRUCT_ROBOT = true; 
const bool SAVE_ROOM_COMPONENTS = true; //save walls separately form lattice/ramp for visualization

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

Polyhedron * get_pole(double x, double y, double z1, double z2) {//return a pole centered at x,y spanning height [z1,z2]
  double x1 = x - pole_width/2;
  double x2 = x + pole_width/2;
  double y1 = y - pole_width/2;
  double y2 = y + pole_width/2;

  double d[] = {x1, x2, y1, y2, z1, z2};
  return box(d, true);
}

int main (int argc, char *argv[]) {
  Parameter::enable();
  if (argc == 2) { 
    unsigned seed = atoi(argv[1]);
    srandom(seed);
  }


  if (CONSTRUCT_ROOM) {
    //walls
    double ww = 0.1; //wall width
    double d1[] = {-width/2, width/2+2*ww, -depth/2, depth/2, -height/2, height/2};
    double d2[] = {-width/2-ww, width/2+ww, -depth/2-ww, depth/2+ww, -height/2-ww, height/2+ww};
    Polyhedron * outer = box(d2, true);
    Polyhedron * inner = box(d1, true);
    Polyhedron * walls = outer->boolean(inner, Complement);
    //top ramp
    PTR<Point> p0 = new InputPoint(-width/2-ww/2, -depth/2-ww/2, height/2+ww/2);
    PTR<Point> p1 = new InputPoint(-width/2-ww/2,  depth/2+ww/2, height/2+ww/2);
    PTR<Point> p2 = new InputPoint( width/2+ww/2, -depth/2-ww/2, height/2+ww/2);
    PTR<Point> p3 = new InputPoint( width/2+ww/2,  depth/2+ww/2, height/2+ww/2);
    PTR<Point> p4 = new InputPoint(-width/2-ww/2, -depth/2-ww/2, height/2-delta_h);
    PTR<Point> p5 = new InputPoint(-width/2-ww/2,  depth/2+ww/2, height/2-delta_h);
    Points points; 
    points.push_back(p0); points.push_back(p1); points.push_back(p2); 
    points.push_back(p3); points.push_back(p4); points.push_back(p5); 
    Polyhedron * ramp = convexHull(points, true);
    //bottom ramp
    p0 = new InputPoint(-width/2-ww/2, -depth/2-ww/2, -height/2-ww/2);
    p1 = new InputPoint(-width/2-ww/2,  depth/2+ww/2, -height/2-ww/2);
    p2 = new InputPoint( width/2+ww/2, -depth/2-ww/2, -height/2-ww/2);
    p3 = new InputPoint( width/2+ww/2,  depth/2+ww/2, -height/2-ww/2);
    p4 = new InputPoint( width/2+ww/2, -depth/2-ww/2, -height/2+delta_h);
    p5 = new InputPoint( width/2+ww/2,  depth/2+ww/2, -height/2+delta_h);
    points.clear(); 
    points.push_back(p0); points.push_back(p1); points.push_back(p2); 
    points.push_back(p3); points.push_back(p4); points.push_back(p5); 
    Polyhedron * ramp2 = convexHull(points, true);
    Polyhedron * tmp = ramp->boolean(ramp2, Union);
    delete ramp, ramp2;
    ramp = tmp;

    //lattice
    //each lattice pole should be just small enough to not touch either the top or bottom ramp
    Polyhedron * lattice = new Polyhedron();
    for (int x = 0; x<3; x++) {
        //how far am I down the ramp, widthwise
        double d = (x-(-width/2-ww/2)) / (2*(width/2+ww/2));  //on a scale from 0 to 1, how much has the eight changed
        //they both have a delta of delta_h+ww/2
        double delta1 = (delta_h+ww/2) * d;
        double delta2 = (delta_h+ww/2) * (1-d);
        double t = height/2+ww/2 - delta2 -0.01;
        double b = -height/2-ww/2 + delta1 +0.01;

        for (int y = -1; y<=1; y++) {
            Polyhedron * pole = get_pole((double)x, (double)y, b, t);
            Polyhedron * tmp = lattice->boolean(pole, Union);
            delete lattice, pole;
            lattice = tmp;
        }
    }


    Polyhedron * components[] = {lattice, walls, ramp};
    Polyhedron * obstacle = multiUnion(components, 3);
    simplify(obstacle, 1e-6);
    savePoly(obstacle, "lattice-room.vtk");
    if (SAVE_ROOM_COMPONENTS) {
      Polyhedron * ramp_and_lattice = ramp->boolean(lattice, Union);
      simplify(ramp_and_lattice, 1e-6);
      savePoly(ramp_and_lattice, "ramp-and-lattice.vtk");

      simplify(lattice, 1e-6);
      savePoly(lattice, "lattice.vtk");
      simplify(walls, 1e-6);
      savePoly(walls, "walls.vtk");
      simplify(ramp, 1e-6);
      savePoly(ramp, "ramp.vtk");

    }
  }

  if (CONSTRUCT_ROBOT) {
    //manually constructed to limit number of faces
    double h = 0.6;//height for plus, set to a for jack

    PTR<Point>  a0 = new InputPoint(-a,-a, r);
    PTR<Point>  a1 = new InputPoint(-a, a, r);
    PTR<Point>  a2 = new InputPoint( a,-a, r);
    PTR<Point>  a3 = new InputPoint( a, a, r);
    PTR<Point>  b0 = new InputPoint(-a,-a, h);
    PTR<Point>  b1 = new InputPoint(-a, a, h);
    PTR<Point>  b2 = new InputPoint( a,-a, h);
    PTR<Point>  b3 = new InputPoint( a, a, h);
    PTR<Point>  c0 = new InputPoint(-a,-a,-h);
    PTR<Point>  c1 = new InputPoint(-a, a,-h);
    PTR<Point>  c2 = new InputPoint( a,-a,-h);
    PTR<Point>  c3 = new InputPoint( a, a,-h);
    PTR<Point>  d0 = new InputPoint(-a,-a,-r);
    PTR<Point>  d1 = new InputPoint(-a, a,-r);
    PTR<Point>  d2 = new InputPoint( a,-a,-r);
    PTR<Point>  d3 = new InputPoint( a, a,-r);
    PTR<Point> bl0 = new InputPoint(-r,-a, h);
    PTR<Point> bl1 = new InputPoint(-r, a, h);
    PTR<Point> cl0 = new InputPoint(-r,-a,-h);
    PTR<Point> cl1 = new InputPoint(-r, a,-h);
    PTR<Point> br2 = new InputPoint( r,-a, h);
    PTR<Point> br3 = new InputPoint( r, a, h);
    PTR<Point> cr2 = new InputPoint( r,-a,-h);
    PTR<Point> cr3 = new InputPoint( r, a,-h);
    PTR<Point> bb0 = new InputPoint(-a,-r, h);
    PTR<Point> bb2 = new InputPoint( a,-r, h);
    PTR<Point> cb0 = new InputPoint(-a,-r,-h);
    PTR<Point> cb2 = new InputPoint( a,-r,-h);
    PTR<Point> bf1 = new InputPoint(-a, r, h);
    PTR<Point> bf3 = new InputPoint( a, r, h);
    PTR<Point> cf1 = new InputPoint(-a, r,-h);
    PTR<Point> cf3 = new InputPoint( a, r,-h);

    Points points;
    points.push_back(a0); points.push_back(a1); points.push_back(a2); points.push_back(a3);
    points.push_back(b0); points.push_back(b1); points.push_back(b2); points.push_back(b3);
    Polyhedron * top = convexHull(points, true);
    points.clear();
    points.push_back(b0); points.push_back(b1); points.push_back(b2); points.push_back(b3);
    points.push_back(c0); points.push_back(c1); points.push_back(c2); points.push_back(c3);
    Polyhedron * middle = convexHull(points, true);
    points.clear();
    points.push_back(c0); points.push_back(c1); points.push_back(c2); points.push_back(c3);
    points.push_back(d0); points.push_back(d1); points.push_back(d2); points.push_back(d3);
    Polyhedron * bottom = convexHull(points, true);
    points.clear();
    points.push_back(b0); points.push_back(b1); points.push_back(c0); points.push_back(c1);
    points.push_back(bl0); points.push_back(bl1); points.push_back(cl0); points.push_back(cl1);
    Polyhedron * left = convexHull(points, true);
    points.clear();
    points.push_back(b2); points.push_back(b3); points.push_back(c2); points.push_back(c3);
    points.push_back(br2); points.push_back(br3); points.push_back(cr2); points.push_back(cr3);
    Polyhedron * right = convexHull(points, true);
    points.clear();
    points.push_back(b0); points.push_back(b2); points.push_back(c0); points.push_back(c2);
    points.push_back(bb0); points.push_back(bb2); points.push_back(cb0); points.push_back(cb2);
    Polyhedron * back = convexHull(points, true);
    points.clear();
    points.push_back(b1); points.push_back(b3); points.push_back(c1); points.push_back(c3);
    points.push_back(bf1); points.push_back(bf3); points.push_back(cf1); points.push_back(cf3);
    Polyhedron * front = convexHull(points, true);

    Polyhedron * tmp1 = left;
    Polyhedron * tmp2 = tmp1->boolean(middle, Union);
    delete tmp1; tmp1 = tmp2;
    tmp2 = tmp1->boolean(right, Union);
    delete tmp1; tmp1 = tmp2;
    tmp2 = tmp1->boolean(back, Union);
    delete tmp1; tmp1 = tmp2;
    tmp2 = tmp1->boolean(front, Union);
    simplify(tmp2, 1e-6);
    savePoly(tmp2, "plus.vtk");
  }

}