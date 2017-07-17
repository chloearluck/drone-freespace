#include "freespace.h"

void savePolyTmp(Polyhedron * p, char * filename) {
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

void saveShell(HFaces hf,  char * filename) {
  int n = strlen(filename);
  char str[n+9];
  strncpy(str, filename, n);
  strncpy(str+n, "-out.vtk", 9);

  ofstream out;
  out.open(str);
  if (out.is_open()) {
    writePolyhedronVTK (hf, out);
    out.close();
  } else {
    cout<<"could not write to file"<<endl;
  }
}

void saveWithShells(Polyhedron * poly, char * filename) {
  savePolyTmp(poly, filename);
  poly->formCells();
  for (int i=0; i<poly->cells.size(); i++) {
    Cell * c = poly->cells[i];
    for (int j=0; j<c->nBoundary(); j++) {
      Shell * s = c->getBoundary(j);
      char str[50];
      sprintf(str, "%s-%d-%d", filename, i, j);
      saveShell(s->getHFaces(), str);
    }
  }
}


PTR<Point> sin_cos_alpha;

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

class SimpleTriangle {
public: 
  PTR<Point> verts[3];

  SimpleTriangle(PTR<Point> a, PTR<Point> b, PTR<Point> c) {
    verts[0] = a;
    verts[1] = b;
    verts[2] = c;
  }
};

class OuterApproxFace {
  public:
  PTR<Point> bottom1, bottom2, top1, top2;
  bool isTrapazoid;

  OuterApproxFace(PTR<Point> bottom1, PTR<Point> bottom2, PTR<Point> top1, PTR<Point> top2) {
    this->bottom1 = bottom1;
    this->bottom2 = bottom2;
    this->top1 = top1;
    this->top2 = top2;
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
PTR<Point> pointInCell(Polyhedron * poly, int i) {
  Cell * cell =  poly->cells[i];
  Face * face = cell->getBoundary(0)->getHFaces()[0]->getF();

  PTR<Point> fp;
  double unit = 1;

  while (true) {
    fp = new FacePoint(cell, unit);
    if (face->contains(fp))
      break;
    fp = new FacePoint(cell, -unit);
    if (face->contains(fp))
      break;
    unit= unit/2;
  }

  PTR<Point> p;
  unit = 1;
  
  while (true) {
    p = new CellInternalPoint(cell, fp, unit);
    if (cell->contains(p))
      break;
    p = new CellInternalPoint(cell, fp, -unit);
    if (cell->contains(p))
      break;
    unit = unit/2;
  }
  return p;
}

int containingCell(Polyhedron * poly, PTR<Point> p) {
  poly->formCells(); 
  for (int i=1; i< poly->cells.size(); i++)
    if (poly->cells[i]->contains(p))
      return i;
  return -1;
}

//generate 3 outer approximation points from rotate point p around the origin, add them to pList
void pointOuterApprox(Points &pList, PTR<Point> p) {
  pList.push_back(p);
  pList.push_back(new RotationPoint(p, sin_cos_alpha));
  pList.push_back(new TangentIntersectionPoint(p, sin_cos_alpha));
}

// generate the 5 out approximation points from rotating a segment around the origin (assuming the nearest point on the segment is an endpoint)
void segmentOuterApprox(Points &pList, PTR<Point> p1, PTR<Point> p2) {
  PTR<Point> outer;
  PTR<Point> inner;

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
  } 
  return convexHull(pList);
}

void splitSimple(std::vector<OuterApproxFace> & tList, SimpleTriangle t) {
  PTR<Point> verts[3];
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

  PTR<Point> newVert = new ZIntercectPoint(verts[0], verts[2], verts[1]);

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
  splitSimple(tList, simple);

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
  Polyhedron * p3 = p1->boolean(p2, Union);
  if (start != mid)
    delete p1;
  if ((mid+1) != end)
    delete p2;
  return p3;
}

Polyhedron * union_all(std::vector<Polyhedron*> pList) {
  return union_all(pList, 0, pList.size()-1);
}

Polyhedron * rotate(Polyhedron * p) {
  Polyhedron * poly = p->triangulate(); 
  for (int i = 0; i< poly->vertices.size(); i++) {
    Vertex * v = poly->vertices[i];
    poly->moveVertex(v, new RotationPoint(v->getP(), sin_cos_alpha));
  }
  return poly;
}

FreeSpace::Node * FreeSpace::findOrAddNode(int cspace_index, int cell_index) {
  for (int i=0; i<nodes[cspace_index].size(); i++) {
    if (nodes[cspace_index][i]->cell_index == cell_index)
      return nodes[cspace_index][i];
  }
  FreeSpace::Node * n = new Node(cspaces[cspace_index], cell_index);
  nodes[cspace_index].push_back(n);
  return n;
}

FreeSpace::FreeSpace(Polyhedron * robot, Polyhedron * obstacle, double theta, double * bounding_box) {
  this->robot = robot->triangulate();
  this->obstacle = obstacle;
  this->bb = box(bounding_box);
  numRotations = floor(2*M_PI/theta);

  const double TAN_THETA = tan(theta/2);
  InputParameter * t = new InputParameter(TAN_THETA);
  sin_cos_alpha = new SinCosAlpha(t);

  std::vector<SimpleTriangle> tList;
  for (int i=0; i<this->robot->faces.size(); i++) {
    HEdges es = this->robot->faces[i]->getBoundary();
    PTR<Point> p = es[0]->tail()->getP();
    PTR<Point> q = es[0]->getNext()->tail()->getP();
    PTR<Point> r = es[0]->getNext()->getNext()->tail()->getP();

    tList.push_back(SimpleTriangle(p,q,r));
  }

  std::vector<OuterApproxFace> splitTList;
  for (int i=0; i< tList.size(); i++) 
    split(splitTList, tList[i]);
  tList.clear();

  std::vector<Polyhedron *> polyList;
  for (int i=0; i<splitTList.size(); i++) {
    polyList.push_back(triangleOuterApprox(splitTList[i]));
  }

  Polyhedron * outerApprox = union_all(polyList);
  Polyhedron * tmp = outerApprox->boolean(robot, Union);
  delete outerApprox;
  outerApprox = tmp;

  for (int i=0; i<polyList.size(); i++)
    delete polyList[i];
  polyList.clear();

  Polyhedron * reflectedOuterApprox = outerApprox->negative();
  vector<Polyhedron *> allRotations;
  allRotations.push_back(reflectedOuterApprox);
  for (int i=0; i< numRotations; i++)
    allRotations.push_back(rotate(allRotations[i]));

  std::vector<Polyhedron*> blockspaces;
  for (int i=0; i<allRotations.size(); i++) {
    Polyhedron * mSum = minkowskiSum(allRotations[i], obstacle); 
    Polyhedron * mSum_complement = complement(bb, mSum);
    cspaces.push_back(mSum_complement); 
    blockspaces.push_back(mSum);
  }

  cout<<"cspaces populated"<<endl;

  for (int i=0; i< numRotations; i++)
    delete allRotations[i];
  allRotations.clear();

  for (int i=0; i<numRotations; i++) {
    cspaces[i]->formCells();
    nodes.push_back(std::vector<Node *>());
  }

  cout<<"creating nodes and edges"<<endl;

  // create all edges
  for (int i=0; i<numRotations; i++) {
    int j = (i+1)%numRotations;
    cout<<"i,j = "<<i<<","<<j<<endl;
    Polyhedron * block_union = blockspaces[i]->boolean(blockspaces[j], Union);
    Polyhedron * intersection = complement(bb, block_union);
    intersection->formCells();
    // char s[50];
    // sprintf(s, "%d-%d", i, j);
    // saveWithShells(intersection, s);
    for (int k=1; k<intersection->cells.size(); k++) {
      PTR<Point> p = pointInCell(intersection, k);
      // bool valid = intersection->contains(p); 
      bool valid = !intersection->cells[k]->getBoundary(0)->getHFaces()[0]->pos();
      if (valid) {
        cout<<"  "<<k<<": "<<p->getP().getX().mid()<<" "<<p->getP().getY().mid()<<" "<<p->getP().getZ().mid()<<endl;
        int ci = containingCell(cspaces[i], p);
        int cj = containingCell(cspaces[j], p);
        FreeSpace::Node * ni = findOrAddNode(i, ci);
        FreeSpace::Node * nj = findOrAddNode(j, cj);
        edges.push_back(new Edge(ni, nj, p));
      }
    }
    delete block_union;
    delete intersection;
    cout<<endl;
  }

  for (int i=0; i<numRotations; i++)
    delete blockspaces[i];
  blockspaces.clear();
}