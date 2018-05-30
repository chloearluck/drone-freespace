#include "freespace.h"
#include "simplify.h"
#include <map>

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

void savePoly(Polyhedron * p, const char * filename) {
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

PTR<Point> sin_cos_alpha;

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

class SimpleTriangle {
public: 
  PTR<Point> verts[3];

  SimpleTriangle(PTR<Point> a, PTR<Point> b, PTR<Point> c) {
    verts[0] = a;
    verts[1] = b;
    verts[2] = c;
  }
};

class OuterApproxVertex {
public:
  PTR<Point> p;
  PTR<Point> r;
  PTR<Point> t;
  OuterApproxVertex(PTR<Point> p) {
    this->p = p; r = 0; t = 0;
  }
};

class OuterApproxFace {
  public:
  OuterApproxVertex *bottom1, *bottom2, *top1, *top2;
  bool isTrapazoid;

  OuterApproxFace(PTR<Point> bottom1, PTR<Point> bottom2, PTR<Point> top1, PTR<Point> top2, std::map<PTR<Point>, OuterApproxVertex*> & pvmap) {
    if (pvmap.find(bottom1) == pvmap.end())
      pvmap[bottom1] = new OuterApproxVertex(bottom1);
    if (bottom2 != NULL && pvmap.find(bottom2) == pvmap.end())
      pvmap[bottom2] = new OuterApproxVertex(bottom2);
    if (pvmap.find(top1) == pvmap.end())
      pvmap[top1] = new OuterApproxVertex(top1);
    if (top2 != NULL && pvmap.find(top2) == pvmap.end())
      pvmap[top2] = new OuterApproxVertex(top2);

    this->bottom1 = pvmap[bottom1];
    this->bottom2 = (bottom2 == NULL? NULL: pvmap[bottom2]);
    this->top1 = pvmap[top1];
    this->top2 = (top2 == NULL? NULL: pvmap[top2]);
    isTrapazoid = (top2 != NULL) && (bottom2 != NULL);
  }
};

Primitive2(DiffZ, PTR<Point>, i, PTR<Point>, j);
int DiffZ::sign() { return (i->getP().getZ() - j->getP().getZ()).sign(); }
struct CompareZ {
  bool operator()(PTR<Point> i, PTR<Point> j) {
    return (DiffZ(i, j) < 0); 
  }
};

//returns a point in the ith cell of polyhedron poly
/*PTR<Point> pointInCell(Polyhedron * poly, int i) {
  Cell * cell =  poly->cells[i];
  Face * face = cell->getShell(0)->getHFaces()[0]->getF();

  PTR<Point> fp;
  double unit = 1;

  do {
    fp = new FacePoint(cell, unit);
    unit /= 2;
  } while (!face->contains(fp));

  PTR<Point> p;
  unit = 1;
  
  do {
    p = new CellInternalPoint(cell, fp, unit);
    unit /= 2;
  } while (!cell->contains(p));
  
  return p;
}*/

//generate 3 outer approximation points from rotate point p around the origin, add them to pList
void pointOuterApprox(Points &pList, OuterApproxVertex * p) {
  pList.push_back(p->p);
  if (p->r == NULL)
    p->r = new RotationPoint(p->p, sin_cos_alpha);
  pList.push_back(p->r);
  if (p->t == NULL)
    p->t = new TangentIntersectionPoint(p->p, sin_cos_alpha);
  pList.push_back(p->t);
}

// generate the 5 out approximation points from rotating a segment around the origin (assuming the nearest point on the segment is an endpoint)
void segmentOuterApprox(Points &pList, OuterApproxVertex * p1, OuterApproxVertex * p2) {
  OuterApproxVertex * outer;
  OuterApproxVertex * inner;

  if (DiffLength(p1->p,p2->p) > 0) { 
    outer = p1;
    inner = p2;
  } else {
    outer = p2;
    inner = p1;
  }

  pointOuterApprox(pList, outer);
  pList.push_back(inner->p);
  if (inner->r == NULL)
    inner->r = new RotationPoint(inner->p, sin_cos_alpha);
  pList.push_back(inner->r);
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
  } 
  return convexHull(pList, true);
}

void splitSimple(std::vector<OuterApproxFace> & tList, SimpleTriangle t, std::map<PTR<Point>, OuterApproxVertex*> & pvmap) {
  PTR<Point> verts[3];
  for (int i=0; i<3; i++)
    verts[i] = t.verts[i];
  CompareZ compz;

  std::sort(verts, verts + 3, compz);

  if (DiffZ(verts[0], verts[1]) == 0) { //no splitting needed
    tList.push_back(OuterApproxFace(verts[0], verts[1], verts[2], NULL, pvmap));
    return;
  }

  if (DiffZ(verts[1], verts[2]) == 0) { //no splitting needed
    tList.push_back(OuterApproxFace(verts[0], NULL, verts[1], verts[2], pvmap));
    return;
  }

  PTR<Point> newVert = new ZIntercectPoint(verts[0], verts[2], verts[1]);

  tList.push_back(OuterApproxFace(verts[0], NULL, verts[1], newVert, pvmap));
  tList.push_back(OuterApproxFace(verts[1], newVert, verts[2], NULL, pvmap));
}

//splits t into sub-triangles that can be rotated, add sub-triangles to tList
void split(std::vector<OuterApproxFace> & tList, SimpleTriangle t, std::map<PTR<Point>, OuterApproxVertex*> & pvmap) {
  if (IsTangent(t.verts[0], t.verts[1], t.verts[2]) == 0) { 
    splitSimple(tList, t, pvmap);
    return;
  }

  PTR<Plane> split = new SplitPlane(t.verts[0], t.verts[1], t.verts[2]); 

  int validIntersects = 0;
  PTR<Point> p0 = new IntersectionPoint(t.verts[0], t.verts[1], split);
  PTR<Point> p1 = new IntersectionPoint(t.verts[1], t.verts[2], split);
  PTR<Point> p2 = new IntersectionPoint(t.verts[2], t.verts[0], split);
  if (PlaneSide(split, t.verts[0]) != PlaneSide(split, t.verts[1]))
    validIntersects += 1;
  if (PlaneSide(split, t.verts[1]) != PlaneSide(split, t.verts[2]))
    validIntersects += 2;
  if (PlaneSide(split, t.verts[0]) != PlaneSide(split, t.verts[2]))
    validIntersects += 4;
  assert(validIntersects == 3 || validIntersects == 5 || validIntersects == 6);

  PTR<Point> intersect1;
  PTR<Point> intersect2;
  PTR<Point> commonVert;
  PTR<Point> otherVert1;
  PTR<Point> otherVert2;

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
  splitSimple(tList, simple, pvmap);

  PTR<Point> verts[4];
  verts[0] = intersect1;
  verts[1] = intersect2;
  verts[2] = otherVert1;
  verts[3] = otherVert2;
  CompareZ compz;
  std::sort(verts, verts + 4, compz);

  PTR<Point> newVert1;
  PTR<Point> newVert2;

  if ((verts[0] == intersect1) || (verts[0] == intersect2) || (verts[3] == intersect1) || (verts[3] == intersect2)) {
    newVert1 = new ZIntercectPoint(verts[0], verts[2], verts[1]);
    newVert2 = new ZIntercectPoint(verts[1], verts[3], verts[2]);
  } else {
    newVert1 = new ZIntercectPoint(verts[0], verts[3], verts[1]);
    newVert2 = new ZIntercectPoint(verts[0], verts[3], verts[2]);
  }

  //if bottom 2 have same z value, there is no bottom triangle 
  if (DiffZ(verts[0], verts[1]) != 0) 
    tList.push_back(OuterApproxFace(verts[0], NULL, verts[1], newVert1, pvmap));
  tList.push_back(OuterApproxFace(verts[1], newVert1, verts[2], newVert2, pvmap));
  //if top 2 have same z value, there is no top triangle
  if (DiffZ(verts[1], verts[2]) != 0) 
    tList.push_back(OuterApproxFace(verts[2], newVert2, verts[3], NULL, pvmap)); 
}

Polyhedron * rotate(Polyhedron * p) {
  Polyhedron *a = new Polyhedron;
  PVMap pvmap;
  for (Vertices::const_iterator v = p->vertices.begin(); v != p->vertices.end(); ++v)
    pvmap.insert(PVPair((*v)->getP(), a->getVertex(new RotationPoint((*v)->getP(), sin_cos_alpha))));
  for (Faces::const_iterator f = p->faces.begin(); f != p->faces.end(); ++f)
    a->getTriangle(*f, pvmap);
  return a;
}

FreeSpace::FreeSpace(Polyhedron * robot, Polyhedron * obstacle, PTR<Object<Parameter> > tan_half_angle, int numRotations) {
  this->robot = robot->triangulate();
  this->obstacle = obstacle;

  sin_cos_alpha = new SinCosAlpha(tan_half_angle);

  std::vector<SimpleTriangle> tList;
  for (int i=0; i<this->robot->faces.size(); i++) {
    HEdges es = this->robot->faces[i]->getBoundary();
    PTR<Point> p = es[0]->tail()->getP();
    PTR<Point> q = es[0]->getNext()->tail()->getP();
    PTR<Point> r = es[0]->getNext()->getNext()->tail()->getP();
    tList.push_back(SimpleTriangle(p,q,r));
  }

  cout<<"loaded triangles"<<endl;

  std::map<PTR<Point>, OuterApproxVertex*> pvmap;
  std::vector<OuterApproxFace> splitTList;
  for (int i=0; i< tList.size(); i++) 
    split(splitTList, tList[i], pvmap);
  tList.clear();

  cout<<"split triangles"<<endl;

  std::vector<Polyhedron *> polyList;
  for (int i=0; i<splitTList.size(); i++) {
    polyList.push_back(triangleOuterApprox(splitTList[i]));
  }
  cout<<"found triangle polyhedrons"<<endl;

  Polyhedron * outerApprox = multiUnion(&polyList[0], polyList.size());
  cout<<"found outerApprox"<<endl;
  Polyhedron * tmp = outerApprox->boolean(robot, Union);
  delete outerApprox;
  outerApprox = tmp;
  simplify(outerApprox, 1e-6, false);
  savePoly(outerApprox, "outerApprox");

  for (int i=0; i<polyList.size(); i++)
    delete polyList[i];
  polyList.clear();

  Polyhedron * reflectedOuterApprox = outerApprox->negative();
  vector<Polyhedron *> allRotations;
  allRotations.push_back(reflectedOuterApprox);
  for (int i=0; i< numRotations; i++)
    allRotations.push_back(rotate(allRotations[i]));
  cout<<"found all rotations"<<endl;

  for (int i=0; i< allRotations.size(); i++) {
    cout<<"minkowskiSum "<<i<<" of "<<allRotations.size()-1<<endl;
    Polyhedron * mSum = minkowskiSumFull(allRotations[i], obstacle);
    simplify(mSum, 1e-6, false);
    // bool selfInt = mSum->intersectsEdges(mSum);
    // cout << "selfInt " << selfInt << endl;
    char s[25];
    sprintf(s, "sum%02d", i);
    savePoly(mSum, s);
    blockspaces.push_back(mSum);
  }

  cout<<"cspaces populated"<<endl;

  for (int i=0; i< numRotations; i++)
    delete allRotations[i];
  allRotations.clear(); 
}
