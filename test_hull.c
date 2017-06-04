#include "hull.h"
#include "geometry3d.h"

// const double THETA = M_PI / 180; //approximately 1 degree  
const double THETA = M_PI / 10;  //temporarily increase degrees for better visibility

/*
 * notes: 
 * first argument is the name of the file (no extension, should I change this?) that is the input polygon
 * second argument (optional) is the seed for random perturbation 
 * the program generates the convex hull polygons and saves them as '0-out.vtk', '1-out.vtk', ...
 * getting a precision exception from the when splitting triangles,
 *  likely the same issue I was getting and ignoring in the original version of the code. Debug.
 * -> problem is when sorting vertices by z component in splitSimple, goes away when I use a random seed value. Ignore for now.
 */


class TriangleOrTrapezoid {
  public:
  Point * verts[4];
  bool isTriangle;

  TriangleOrTrapezoid(Point * a, Point * b, Point * c) {
    verts[0] = a;
    verts[1] = b;
    verts[2] = c;
    isTriangle = true;
  }

  TriangleOrTrapezoid(Point * a, Point * b, Point * c, Point * d) {
  	verts[0] = a;
    verts[1] = b;
    verts[2] = c;
    verts[3] = d;
    isTriangle = false;
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
  char * str = new char[n+8];
  strncpy(str, filename, n);
  strncpy(str+n, "-out.vtk", 8);

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
  printf("pointOuterApprox:\n"); // !!! precision exception somewhere in here
	Point * p_theta = new RotationPoint(p, THETA);
	Point * q  = new TangentIntersectionPoint(p, p_theta);

	pList.push_back(p);
	pList.push_back(p_theta);
	pList.push_back(q);
}

// generate the 5 out approximation points from rotating a segment around the origin (assuming the nearest point on the segment is an endpoint)
void segmentOuterApprox(Points &pList, Point * p1, Point * p2) {
  printf("segmentOuterApprox:\n"); // !!! precision exception somewhere in here
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
Polyhedron * triangleOuterApprox(TriangleOrTrapezoid t) { 
	printf("traingleOuterApprox:\n"); // !!! precision exception somewhere in here
  Points pList;
  bool b;
  printf("comp1\n");
  b = ((t.verts[0]->getP().getZ() - t.verts[1]->getP().getZ()).sign() == 0);
  printf("comp2\n");
  b = ((t.verts[0]->getP().getZ() - t.verts[2]->getP().getZ()).sign() == 0);
  printf("done comp\n");

  if ((t.verts[0]->getP().getZ() - t.verts[1]->getP().getZ()).sign() == 0) {
		segmentOuterApprox(pList, t.verts[0], t.verts[1]);
		pointOuterApprox(pList, t.verts[2]);
	} else if ((t.verts[0]->getP().getZ() - t.verts[2]->getP().getZ()).sign() == 0) {
		segmentOuterApprox(pList, t.verts[0], t.verts[2]);
		pointOuterApprox(pList, t.verts[1]);
	} else {
		segmentOuterApprox(pList, t.verts[1], t.verts[2]);
		pointOuterApprox(pList, t.verts[0]);
	}
	return convexHull(pList);
}

void splitSimple(std::vector<TriangleOrTrapezoid> & tList, TriangleOrTrapezoid t) {
  Point * verts[3];
	for (int i=0; i<3; i++)
		verts[i] = t.verts[i];
	CompareZ compz;

  std::sort(verts, verts + 3, compz);

	if (((verts[0]->getP().getZ() - verts[1]->getP().getZ()).sign() == 0) || ((verts[1]->getP().getZ() - verts[2]->getP().getZ()).sign() == 0)) {
    tList.push_back(t); //no splitting is needed
    return;
  }

  //get intersection between plane with v=(0,0,1) and d=vert[1]->getP().getZ()
  PTR<Plane> plane = new PointNormalPlane(verts[1], new InputPoint(0,0,1));
  Point * newVert = new IntersectionPoint(verts[0], verts[2], plane);

  //TO DO: does this maintain the face direction of triangle t? what if the middle vert is on the other side? does it matter?
  tList.push_back(TriangleOrTrapezoid(verts[0], verts[1], newVert));
  tList.push_back(TriangleOrTrapezoid(verts[1], verts[2], newVert));
}

//splits t into sub-triangles that can be rotated, add sub-triangles to tList
void split(std::vector<TriangleOrTrapezoid> & tList, TriangleOrTrapezoid t) {
  PTR<Plane> plane = new TrianglePlane(t.verts[0], t.verts[1], t.verts[2]);
  Point * n = new NormalToPlane(plane); 
	Point * p = new IntersectionPoint(new InputPoint(0,0,0), n, plane);

	PTR<Plane> edgePlane0 = new TrianglePlane(t.verts[1], t.verts[0], new AddPoint(t.verts[1], n));
	PTR<Plane> edgePlane1 = new TrianglePlane(t.verts[2], t.verts[1], new AddPoint(t.verts[2], n));
	PTR<Plane> edgePlane2 = new TrianglePlane(t.verts[0], t.verts[2], new AddPoint(t.verts[0], n));
  
  // bool b = (PlaneSide(edgePlane0, p) < 0);
  // printf("PlaneSide1 %s\n", (b? "true" : "false"));
  // b = (PlaneSide(edgePlane1, p) < 0);
  // printf("PlaneSide2 %s\n", (b? "true" : "false"));
  // b = (PlaneSide(edgePlane2, p) < 0);
  // printf("PlaneSide3 %s\n", (b? "true" : "false"));

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

	TriangleOrTrapezoid simple = TriangleOrTrapezoid(intersect1, intersect2, commonVert); //TODO: check direction of triangle faces
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

  //if bottom 2 have same z value, this is no bottom triangle 
  if ((verts[0]->getP().getZ() - verts[1]->getP().getZ()).sign() != 0) //test me!!!
	  tList.push_back(TriangleOrTrapezoid(verts[0], verts[1], newVert1));
 	tList.push_back(TriangleOrTrapezoid(verts[1], newVert1, verts[2], newVert2));
  //if top 2 have same z value, there is no top triangle
 	if ((verts[1]->getP().getZ() - verts[2]->getP().getZ()).sign() != 0) //test me!!!
	  tList.push_back(TriangleOrTrapezoid(verts[2], verts[3], newVert2));	
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

  poly = poly->triangulate(); 

  std::vector<TriangleOrTrapezoid> tList;
  for (int i=0; i<poly->faces.size(); i++) {
  	HEdges es = poly->faces[i]->getBoundary();
  	Point * p = es[0]->tail()->getP();
  	Point * q = es[0]->getNext()->tail()->getP();
  	Point * r = es[0]->getNext()->getNext()->tail()->getP();

  	tList.push_back(TriangleOrTrapezoid(p,q,r));
  }

  std::vector<TriangleOrTrapezoid> splitTList;
  //call split on each triangle in tList
  for (int i=0; i< tList.size(); i++) 
  	split(splitTList, tList[i]);

  //find the rotated convex hull of each triangle in splitTList
  std::vector<Polyhedron *> polyList;

  for (int i=0; i<splitTList.size(); i++) {
  	polyList.push_back(triangleOuterApprox(splitTList[i]));
  }

  // char s[3];
  // for (int i=0; i<polyList.size(); i++) {
  //   sprintf(s, "%d", i);
  //   savePoly(polyList[i], s);
  // }

  // // repeatedly take the union 
  // for (int i=0; i< polyList.size(); i++) 
  //   outerApprox = outerApprox->boolean(polyList[i], Union);

  return 0;
}
