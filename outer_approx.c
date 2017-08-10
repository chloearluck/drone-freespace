#include "hull.h"
#include "geometry3d.h"
#include <cstring>

// const double THETA = M_PI / 180; //approximately 1 degree  
const double THETA = M_PI / 10;  //temporarily increase degrees for better visibility

/*
 * first argument is the name of the file (no extension) that is the input polygon
 * second argument (optional) is the seed for random perturbation 

 * bugs:
 * precision exception when sorting vertices by z component in splitSimple, 
 *    goes away when I use any random seed value. Ignoring for now.
 *
 * union behaving strangely
 * - union of two non-empty polyhedrons becomes an empty polyhedron
 * - seg fault happening at some unions
 */


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
  char * str = new char[n+4];
  strncpy(str, filename, n);
  strncpy(str+n, ".vtk", 4);

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
  char * str = new char[n+9];
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

struct CompareZ {
  bool operator()(Point * i, Point * j) {
    return ((i->getP().getZ() - j->getP().getZ()).sign() < 0);
  }
};

void printPoint(Point * p) {
	printf("  %f %f %f\n", p->getP().getX().mid(), p->getP().getY().mid(), p->getP().getZ().mid());
}

void printTriangle(OuterApproxFace t) {
  Point * p = ((t.top2 != NULL)? t.top2 : t.bottom2);
  printf("(%f %f %f)(%f %f %f)(%f %f %f)\n",
          t.top1->getP().getX().mid(), t.top1->getP().getY().mid(), t.top1->getP().getZ().mid(), 
          t.bottom1->getP().getX().mid(), t.bottom1->getP().getY().mid(), t.bottom1->getP().getZ().mid(), 
          p->getP().getX().mid(), p->getP().getY().mid(), p->getP().getZ().mid()
        );
}

void printTriangle(SimpleTriangle t) {
  printf("(%f %f %f)(%f %f %f)(%f %f %f)\n",
          t.verts[0]->getP().getX().mid(), t.verts[0]->getP().getY().mid(), t.verts[0]->getP().getY().mid(),
          t.verts[1]->getP().getX().mid(), t.verts[1]->getP().getY().mid(), t.verts[1]->getP().getY().mid(),
          t.verts[2]->getP().getX().mid(), t.verts[2]->getP().getY().mid(), t.verts[2]->getP().getY().mid()
        );
  if ((t.verts[0]->getP().getZ() - t.verts[1]->getP().getZ()).sign() == 0)
    printf("01\n");
  if ((t.verts[0]->getP().getZ() - t.verts[2]->getP().getZ()).sign() == 0)
    printf("02\n");
  if ((t.verts[1]->getP().getZ() - t.verts[2]->getP().getZ()).sign() == 0)
    printf("12\n");
}

void printPolyVerts(Polyhedron * poly) {
  printf("{\n");
  Vertices vlist = poly->vertices;
  for (int i= 0; i<vlist.size(); i++) {
    Point * p = vlist[i]->getP();
    printf("  %f %f %f\n", p->getP().getX().mid(), p->getP().getY().mid(), p->getP().getZ().mid());
  }
  printf("}\n------\n");
}

void printEdge(HEdge * e) {
  Point * t = e->tail()->getP();
  Point * h = e->head()->getP();
  printf("  %f %f %f -> %f %f %f\n", t->getP().getX().mid(), t->getP().getY().mid(), t->getP().getZ().mid(),
                                     h->getP().getX().mid(), h->getP().getY().mid(), h->getP().getZ().mid());
}

void printfFaces(Polyhedron * poly) {
  bool matlab = true;
  for (int i=0; i<poly->faces.size(); i++) {
    HEdges es = poly->faces[i]->getBoundary();
    HEdge * e = es[0];
    if (matlab) printf("T = [\n");
    
    do {
      Point * p = e->tail()->getP();
      printPoint(p);
      e = e->getNext();
    } while (e != es[0]);
    if (matlab) printf("];  Ts = cat(3, Ts, T);\n");
    printf("\n");
  }
}

//generate 3 outer approximation points from rotate point p around the origin, add them to pList
void pointOuterApprox(Points &pList, Point * p) {
	Point * p_theta = new RotationPoint(p, THETA);
	Point * q  = new TangentIntersectionPoint(p, p_theta);

	pList.push_back(p);
	pList.push_back(p_theta);
	pList.push_back(q);
}

// generate the 5 out approximation points from rotating a segment around the origin (assuming the nearest point on the segment is an endpoint)
void segmentOuterApprox(Points &pList, Point * p1, Point * p2) {
	Point * outer;
	Point * inner;

	if (p1->get().dot(p1->get()) > p2->get().dot(p2->get())) { 
    outer = p1;
    inner = p2;
  } else {
    outer = p2;
    inner = p1;
  }

  pointOuterApprox(pList, outer);
  Point * inner_theta = new RotationPoint(inner, THETA);

  pList.push_back(inner);
  pList.push_back(inner_theta);
}

// generate the 8 outer approximation points from rotating a triangle around the origin
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

	if ((verts[0]->getP().getZ() - verts[1]->getP().getZ()).sign() == 0) { //no splitting needed
    tList.push_back(OuterApproxFace(verts[0], verts[1], verts[2], NULL));
    return;
  }

  if ((verts[1]->getP().getZ() - verts[2]->getP().getZ()).sign() == 0) { //no splitting needed
    tList.push_back(OuterApproxFace(verts[0], NULL, verts[1], verts[2]));
    return;
  }

  //get intersection between plane with v=(0,0,1) and d=vert[1]->getP().getZ()
  PTR<Plane> plane = new PointNormalPlane(verts[1], new InputPoint(0,0,1));
  Point * newVert = new IntersectionPoint(verts[0], verts[2], plane);

  tList.push_back(OuterApproxFace(verts[0], NULL, verts[1], newVert));
  tList.push_back(OuterApproxFace(verts[1], newVert, verts[2], NULL));
}

//splits t into sub-triangles that can be rotated, add sub-triangles to tList
void split(std::vector<OuterApproxFace> & tList, SimpleTriangle t) {
  PTR<Plane> plane = new TrianglePlane(t.verts[0], t.verts[1], t.verts[2]);
  Point * n = new NormalToPlane(plane); 
	Point * p = new IntersectionPoint(new InputPoint(0,0,0), n, plane);

	PTR<Plane> edgePlane0 = new TrianglePlane(t.verts[1], t.verts[0], new AddPoint(t.verts[1], n));
	PTR<Plane> edgePlane1 = new TrianglePlane(t.verts[2], t.verts[1], new AddPoint(t.verts[2], n));
	PTR<Plane> edgePlane2 = new TrianglePlane(t.verts[0], t.verts[2], new AddPoint(t.verts[0], n));
  
  if ((PlaneSide(edgePlane0, p) < 0) || (PlaneSide(edgePlane1, p) < 0) || (PlaneSide(edgePlane2, p) < 0)) {
  	splitSimple(tList, t);
  	return;
  }
 
  printf("\ntriangle is tangent\n\n");
  
  Point * z = new InputPoint(0,0,1);
  PTR<Plane> split = new PointNormalPlane(p, new ACrossBPoint(n, z)); 
    //if || n cross z || = 0, then the triangle is in the xy, how should we handle this? 

  int validIntersects = 0;
  Point * p0 = new IntersectionPoint(t.verts[0], t.verts[1], split);
  Point * p1 = new IntersectionPoint(t.verts[1], t.verts[2], split);
  Point * p2 = new IntersectionPoint(t.verts[2], t.verts[0], split);
  if ((PlaneSide(edgePlane1, p0) > 0) && (PlaneSide(edgePlane2, p0) > 0))
		validIntersects += 1;
	if ((PlaneSide(edgePlane0, p1) > 0) && (PlaneSide(edgePlane2, p1) > 0))
		validIntersects += 2;
	if ((PlaneSide(edgePlane0, p2) > 0) && (PlaneSide(edgePlane1, p2) > 0))
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

	PTR<Plane> z1 = new PointNormalPlane(verts[1], z);
  PTR<Plane> z2 = new PointNormalPlane(verts[2], z);
 
  Point * newVert1;
  Point * newVert2;

  Point * tVertSorted[3];
  for (int i=0; i<3; i++)
  	tVertSorted[i] = t.verts[i];
  std::sort(tVertSorted, tVertSorted + 3, compz);

	if ((verts[1] == t.verts[0]) || (verts[1] == t.verts[1]) || (verts[1] == t.verts[2]))
		newVert1 = new IntersectionPoint(intersect1, intersect2, z1);
	else {
		if (verts[1]->getP().getZ() < tVertSorted[1]->getP().getZ())
			newVert1 = new IntersectionPoint(tVertSorted[0], tVertSorted[1], z1);
		else 
			newVert1 = new IntersectionPoint(tVertSorted[1], tVertSorted[2], z1);
	}

	if ((verts[2] == t.verts[0]) || (verts[2] == t.verts[1]) || (verts[2] == t.verts[2]))
		newVert2 = new IntersectionPoint(intersect1, intersect2, z2);
	else {
		if (verts[2]->getP().getZ() < tVertSorted[1]->getP().getZ())
			newVert2 = new IntersectionPoint(tVertSorted[0], tVertSorted[1], z2);
		else 
			newVert2 = new IntersectionPoint(tVertSorted[1], tVertSorted[2], z2);
	}

  //if bottom 2 have same z value, there is no bottom triangle 
  if ((verts[0]->getP().getZ() - verts[1]->getP().getZ()).sign() != 0) //test me!!!
	  tList.push_back(OuterApproxFace(verts[0], NULL, verts[1], newVert1));
 	tList.push_back(OuterApproxFace(verts[1], newVert1, verts[2], newVert2));
  //if top 2 have same z value, there is no top triangle
 	if ((verts[1]->getP().getZ() - verts[2]->getP().getZ()).sign() != 0) //test me!!!
	  tList.push_back(OuterApproxFace(verts[2], newVert2, verts[3], NULL));	
}

Polyhedron * union_all(std::vector<Polyhedron*> pList) {
  Polyhedron * outerApprox = pList[0];
  for (int i=1; i< pList.size(); i++) {
    printf("union %d of %d\n", i, pList.size()-1);
    outerApprox = outerApprox->boolean(pList[i], Union);
  }
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

  Parameter::enable();

  Polyhedron * poly = loadPoly(filename);
  if (poly == NULL)
    return 1;

  poly = poly->triangulate(); 

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

  //find the rotated convex hull of each triangle in splitTList
  std::vector<Polyhedron *> polyList;

  for (int i=0; i<splitTList.size(); i++) {
  	polyList.push_back(triangleOuterApprox(splitTList[i]));
  }

  char s[30];
  for (int i=0; i<polyList.size(); i++) {
    sprintf(s, "%d", i);
    savePoly(polyList[i], s);
  }

  Polyhedron * outerApprox = union_all(polyList);
  printf("saving outerApprox\n");
  savePoly(outerApprox, filename);

  Parameter::disable();

  return 0;
}
