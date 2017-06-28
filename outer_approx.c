#include "hull.h"
#include "geometry3d.h"
#include "mink.h"
#include <cstring>

const double TAN_THETA = tan(M_PI / 36);  //approximately 5 degrees
bool SAVE_CONVEX_HULLS = false;
bool SAVE_TRIANGULATION = false;


class InputParameter : public Object<Parameter> {
public:
  InputParameter (double x) { set(Parameter::input(x)); }
};

class SinCosAlpha : public Point {
  Object<Parameter> *tan_theta;
  PV3 calculate () {
    Parameter t = tan_theta->get();
    Parameter sint = 2*t/(1+t*t);
    Parameter cost = (1-t*t)/(1+t*t);
    Parameter alpha = (1-cost)/sint;
    return PV3(sint, cost, alpha);
  }

public:
  SinCosAlpha (Object<Parameter> *t) : tan_theta(t) {}
};

InputParameter * t = new InputParameter(TAN_THETA);
Point * sin_cos_alpha = new SinCosAlpha(t);



class SimpleTriangle {
public: 
  Point * verts[3];

  SimpleTriangle(Point * a, Point * b, Point * c) {
    verts[0] = a;
    verts[1] = b;
    verts[2] = c;
  }
};

class OuterApproxFace {
  public:
  Point * bottom1, * bottom2, * top1, * top2;
  bool isTrapazoid;

  OuterApproxFace(Point * bottom1, Point * bottom2, Point * top1, Point * top2) {
    this->bottom1 = bottom1;
    this->bottom2 = bottom2;
    this->top1 = top1;
    this->top2 = top2;
    isTrapazoid = (top2 != NULL) && (bottom2 != NULL);
  }
};

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
    printf("could not write to file\n");
  }
}

void save_triangulation(std::vector<OuterApproxFace> tList, char * filename) {
  std::vector<Point *> vertList;
  std::vector<Point *> rotatedVertList;

  int size = 0;
  for (int i=0; i<tList.size(); i++) {
    if(std::find(vertList.begin(), vertList.end(), tList[i].top1) == vertList.end()) {
      vertList.push_back(tList[i].top1);
      rotatedVertList.push_back(new RotationPoint(tList[i].top1, sin_cos_alpha));
    }
    if(tList[i].top2 != NULL && std::find(vertList.begin(), vertList.end(), tList[i].top2) == vertList.end()) {
      vertList.push_back(tList[i].top2);
      rotatedVertList.push_back(new RotationPoint(tList[i].top2, sin_cos_alpha));
    }
    if(std::find(vertList.begin(), vertList.end(), tList[i].bottom1) == vertList.end()) {
      vertList.push_back(tList[i].bottom1);
      rotatedVertList.push_back(new RotationPoint(tList[i].bottom1, sin_cos_alpha));
    }
    if(tList[i].bottom2 != NULL && std::find(vertList.begin(), vertList.end(), tList[i].bottom2) == vertList.end()) {
      vertList.push_back(tList[i].bottom2);
      rotatedVertList.push_back(new RotationPoint(tList[i].bottom2, sin_cos_alpha));
    }
    size = size + 3;
    if (tList[i].top2 != NULL && tList[i].bottom2 != NULL)
      size++;
  }

  int n = strlen(filename);
  char str[n+18];
  char str2[n+13];
  strncpy(str, filename, n);
  strncpy(str2, filename, n);
  strncpy(str+n, "-triangulated.vtk", 18);
  strncpy(str2+n, "-rotated.vtk", 13);
  ofstream ostr;
  ofstream ostr2;
  ostr.open(str);
  ostr2.open(str2);
  if (ostr.is_open()) { 
    ostr << setprecision(20) << "# vtk DataFile Version 3.0" << endl
         << "vtk output" << endl << "ASCII" << endl
         << "DATASET POLYDATA" << endl 
         << "POINTS " << vertList.size() << " double" << endl;
    ostr2 << setprecision(20) << "# vtk DataFile Version 3.0" << endl
         << "vtk output" << endl << "ASCII" << endl
         << "DATASET POLYDATA" << endl 
         << "POINTS " << vertList.size() << " double" << endl;
    for (int i=0; i<vertList.size(); i++)
      ostr << vertList[i]->getP().getX().mid() << " " << vertList[i]->getP().getY().mid() << " " << vertList[i]->getP().getZ().mid() << endl;
    for (int i=0; i<vertList.size(); i++)
      ostr2 << rotatedVertList[i]->getP().getX().mid() << " " << rotatedVertList[i]->getP().getY().mid() << " " << rotatedVertList[i]->getP().getZ().mid() << endl;
 
    ostr << endl << "POLYGONS " << tList.size() << " " << size+tList.size() << endl;
    ostr2 << endl << "POLYGONS " << tList.size() << " " << size+tList.size() << endl;
    for (int i=0; i<tList.size(); i++) {
      int t1,t2,b1,b2;
      t1 = std::find(vertList.begin(), vertList.end(), tList[i].top1) - vertList.begin();
      if (tList[i].top2 != NULL) t2 = std::find(vertList.begin(), vertList.end(), tList[i].top2) - vertList.begin();
      b1 = std::find(vertList.begin(), vertList.end(), tList[i].bottom1) - vertList.begin();
      if (tList[i].bottom2 != NULL) b2 = std::find(vertList.begin(), vertList.end(), tList[i].bottom2) - vertList.begin();

      if (tList[i].top2 == NULL) {
        ostr << "3 " << t1 << " " << b1 << " " << b2 << endl;
        ostr2 << "3 " << t1 << " " << b1 << " " << b2 << endl;
      } else if (tList[i].bottom2 == NULL) {
        ostr << "3 " << t1 << " " << t2 << " " << b1 << endl;
        ostr2 << "3 " << t1 << " " << t2 << " " << b1 << endl;
      } else {
        if (SegmentIntersect(tList[i].top1, tList[i].bottom1, tList[i].top2, tList[i].bottom2)) {
          ostr << "4 " << t1 << " " << t2 << " " << b1 << " " << b2 << endl;
          ostr2 << "4 " << t1 << " " << t2 << " " << b1 << " " << b2 << endl;
        } else {
          ostr << "4 " << t1 << " " << t2 << " " << b2 << " " << b1 << endl;
          ostr2 << "4 " << t1 << " " << t2 << " " << b2 << " " << b1 << endl;
        }
      }
    }
    ostr.close();
    ostr2.close();
  } else {
    printf("could not write to file\n");
  }
}

Primitive2(DiffZ, Point*, i, Point*, j);
int DiffZ::sign() { return (i->getP().getZ() - j->getP().getZ()).sign(); }
struct CompareZ {
  bool operator()(Point * i, Point * j) {
    return (DiffZ(i, j) < 0); 
  }
};

//generate 3 outer approximation points from rotate point p around the origin, add them to pList
void pointOuterApprox(Points &pList, Point * p) {
	pList.push_back(p);
	pList.push_back(new RotationPoint(p, sin_cos_alpha));
	pList.push_back(new TangentIntersectionPoint(p, sin_cos_alpha));
}

// generate the 5 out approximation points from rotating a segment around the origin (assuming the nearest point on the segment is an endpoint)
void segmentOuterApprox(Points &pList, Point * p1, Point * p2) {
	Point * outer;
	Point * inner;

  if (DiffLength(p1,p2) > 0) { 
    outer = p1;
    inner = p2;
  } else {
    outer = p2;
    inner = p1;
  }

  pointOuterApprox(pList, outer);
  pList.push_back(inner);
  pList.push_back(new RotationPoint(inner, sin_cos_alpha));
}

// generate the 8 o 10 convex hull points from rotating a triangle or trapezoid around the origin
Polyhedron * triangleOuterApprox(OuterApproxFace t) { 
	Points pList;
  
  if (t.isTrapazoid) {
    segmentOuterApprox(pList, t.top1, t.top2);
    segmentOuterApprox(pList, t.bottom1, t.bottom2);
  } else if (t.top2 == NULL) {
    pointOuterApprox(pList, t.top1);
    segmentOuterApprox(pList, t.bottom1, t.bottom2);
  } else if (t.bottom2 == NULL) {
    segmentOuterApprox(pList, t.top1, t.top2);
    pointOuterApprox(pList, t.bottom1);
  } else {
    printf("invalid face\n");
  }
  return convexHull(pList);
}

void splitSimple(std::vector<OuterApproxFace> & tList, SimpleTriangle t) {
  Point * verts[3];
	for (int i=0; i<3; i++)
		verts[i] = t.verts[i];
	CompareZ compz;

  std::sort(verts, verts + 3, compz);

	if (DiffZ(verts[0], verts[1]) == 0) { //no splitting needed
    tList.push_back(OuterApproxFace(verts[0], verts[1], verts[2], NULL));
    return;
  }

  if (DiffZ(verts[1], verts[2]) == 0) { //no splitting needed
    tList.push_back(OuterApproxFace(verts[0], NULL, verts[1], verts[2]));
    return;
  }

  Point * newVert = new ZIntercectPoint(verts[0], verts[2], verts[1]);

  tList.push_back(OuterApproxFace(verts[0], NULL, verts[1], newVert));
  tList.push_back(OuterApproxFace(verts[1], newVert, verts[2], NULL));
}

//splits t into sub-triangles that can be rotated, add sub-triangles to tList
void split(std::vector<OuterApproxFace> & tList, SimpleTriangle t) {
  if (IsTangent(t.verts[0], t.verts[1], t.verts[2]) == 0) { 
  	splitSimple(tList, t);
  	return;
  }

  PTR<Plane> split = new SplitPlane(t.verts[0], t.verts[1], t.verts[2]); 

  int validIntersects = 0;
  Point * p0 = new IntersectionPoint(t.verts[0], t.verts[1], split);
  Point * p1 = new IntersectionPoint(t.verts[1], t.verts[2], split);
  Point * p2 = new IntersectionPoint(t.verts[2], t.verts[0], split);
  if (PlaneSide(split, t.verts[0]) != PlaneSide(split, t.verts[1]))
    validIntersects += 1;
  if (PlaneSide(split, t.verts[1]) != PlaneSide(split, t.verts[2]))
    validIntersects += 2;
  if (PlaneSide(split, t.verts[0]) != PlaneSide(split, t.verts[2]))
    validIntersects += 4;
  assert(validIntersects == 3 || validIntersects == 5 || validIntersects == 6);

	Point * intersect1;
	Point * intersect2;
	Point * commonVert;
	Point * otherVert1;
	Point * otherVert2;

	if (validIntersects == 3) {
		intersect1 = p0;
		intersect2 = p1;
		commonVert = t.verts[1];
		otherVert1 = t.verts[0];
		otherVert2 = t.verts[2];
	} else if (validIntersects == 5) {
		intersect1 = p0;
		intersect2 = p2;
		commonVert = t.verts[0];
		otherVert1 = t.verts[1];
		otherVert2 = t.verts[2];
	} else {
		intersect1 = p1;
		intersect2 = p2;
		commonVert = t.verts[2];
		otherVert1 = t.verts[0];
		otherVert2 = t.verts[1];
	} 

	SimpleTriangle simple = SimpleTriangle(intersect1, intersect2, commonVert); 
	splitSimple(tList, simple);

	Point * verts[4];
	verts[0] = intersect1;
	verts[1] = intersect2;
	verts[2] = otherVert1;
	verts[3] = otherVert2;
	CompareZ compz;
	std::sort(verts, verts + 4, compz);

	Point * newVert1;
  Point * newVert2;

  if ((verts[0] == intersect1) || (verts[0] == intersect2) || (verts[3] == intersect1) || (verts[3] == intersect2)) {
    newVert1 = new ZIntercectPoint(verts[0], verts[2], verts[1]);
    newVert2 = new ZIntercectPoint(verts[1], verts[3], verts[2]);
  } else {
    newVert1 = new ZIntercectPoint(verts[0], verts[3], verts[1]);
    newVert2 = new ZIntercectPoint(verts[0], verts[3], verts[2]);
  }

  //if bottom 2 have same z value, there is no bottom triangle 
  if (DiffZ(verts[0], verts[1]) != 0) 
	  tList.push_back(OuterApproxFace(verts[0], NULL, verts[1], newVert1));
 	tList.push_back(OuterApproxFace(verts[1], newVert1, verts[2], newVert2));
  //if top 2 have same z value, there is no top triangle
  if (DiffZ(verts[1], verts[2]) != 0) 
	  tList.push_back(OuterApproxFace(verts[2], newVert2, verts[3], NULL));	
}

Polyhedron * union_all(std::vector<Polyhedron*> pList, int start, int end) {
  if (start == end) 
    return pList[start];
  if ((start+1) == end) 
    return pList[start]->boolean(pList[end], Union);
  int mid = (start+end)/2;
  Polyhedron * p1 = union_all(pList, start, mid);
  Polyhedron * p2 = union_all(pList, mid+1, end);
  printf("taking union of (%d,%d) and (%d,%d)\n", start, mid, mid+1, end);
  Polyhedron * p3 = p1->boolean(p2, Union);
  delete p1;
  delete p2;
  return p3;
}

Polyhedron * union_all(std::vector<Polyhedron*> pList) {
  return union_all(pList, 0, pList.size()-1);
}

int main (int argc, char *argv[]) {
	if (argc < 2) {
    printf("not enough arguments\n");
    return 1;
  }

  char * filename =  argv[1];

  if (argc == 3) { 
    unsigned seed = atoi(argv[2]);
    srandom(seed);
  }

  Polyhedron * poly = loadPoly(filename);
  if (poly == NULL)
    return 1;

  
  Polyhedron * poly2 = poly->triangulate(); 
  delete poly;
  poly = poly2;

  std::vector<SimpleTriangle> tList;
  for (int i=0; i<poly->faces.size(); i++) {
  	HEdges es = poly->faces[i]->getBoundary();
  	Point * p = es[0]->tail()->getP();
  	Point * q = es[0]->getNext()->tail()->getP();
  	Point * r = es[0]->getNext()->getNext()->tail()->getP();

  	tList.push_back(SimpleTriangle(p,q,r));
  }

  std::vector<OuterApproxFace> splitTList;
  //call split on each triangle in tList
  for (int i=0; i< tList.size(); i++) 
  	split(splitTList, tList[i]);

  if (SAVE_TRIANGULATION)
    save_triangulation(splitTList, filename);

  //find the rotated convex hull of each triangle in splitTList
  std::vector<Polyhedron *> polyList;

  for (int i=0; i<splitTList.size(); i++) {
  	polyList.push_back(triangleOuterApprox(splitTList[i]));
  }

  if (SAVE_CONVEX_HULLS) {
    char s[30];
    for (int i=0; i<polyList.size(); i++) {
      sprintf(s, "%d", i);
      savePoly(polyList[i], s);
    }
  }

  Polyhedron * outerApprox = union_all(polyList);
  printf("saving outerApprox\n");
  savePoly(outerApprox, filename);

  return 0;
}
