#include "freespace.h"
#include "simplify.h"

const bool USE_SIMPLIFY = false;

void savePolyTmp(Polyhedron * p, const char * filename) {
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
  for (int i=0; i<poly->cells.size(); i++) {
    Cell * c = poly->cells[i];
    for (int j=0; j<c->nShells(); j++) {
      Shell * s = c->getShell(j);
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

Primitive3(SaveTriangulation, std::vector<OuterApproxFace>, tList, char*, filename, PTR<Point>, sin_cos_alpha);
int SaveTriangulation::sign() {
  std::vector< PTR<Point> > vertList;
  std::vector< PTR<Point> > rotatedVertList;

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
    cout<<"could not write to file"<<endl;
  }
  return 0;
}

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

Polyhedron * rotate(Polyhedron * p) {
  Polyhedron * poly = p->triangulate(); 
  for (int i = 0; i< poly->vertices.size(); i++) {
    Vertex * v = poly->vertices[i];
    poly->moveVertex(v, new RotationPoint(v->getP(), sin_cos_alpha));
  }
  return poly;
}

FreeSpace::Node * FreeSpace::findOrAddNode(int space_index, int cell_index) {
  for (int i=0; i<nodes[space_index].size(); i++) {
    if (nodes[space_index][i]->cell_index == cell_index)
      return nodes[space_index][i];
  }
  FreeSpace::Node * n = new Node(space_index, cell_index);
  nodes[space_index].push_back(n);
  return n;
}

FreeSpace::Node * FreeSpace::findNode(int space_index, int cell_index) {
  for (int i=0; i<nodes[space_index].size(); i++) {
    if (nodes[space_index][i]->cell_index == cell_index)
      return nodes[space_index][i];
  }
  return NULL;
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
  this->robot->computeWindingNumbers();
  for (int i=0; i<this->robot->faces.size(); i++) {
    HEdges es = this->robot->faces[i]->getBoundary();
    PTR<Point> p = es[0]->tail()->getP();
    PTR<Point> q = es[0]->getNext()->tail()->getP();
    PTR<Point> r = es[0]->getNext()->getNext()->tail()->getP();

    if (this->robot->faces[i]->getHFace(0)->getS()->getC()->getWN() == 0)
      tList.push_back(SimpleTriangle(p,q,r));
    else
      tList.push_back(SimpleTriangle(r,q,p));
  }

  cout<<"loaded triangles"<<endl;

  std::vector<OuterApproxFace> splitTList;
  for (int i=0; i< tList.size(); i++) 
    split(splitTList, tList[i]);
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
  // savePolyTmp(outerApprox, "outerApprox");

  for (int i=0; i<polyList.size(); i++)
    delete polyList[i];
  polyList.clear();

  Polyhedron * reflectedOuterApprox = outerApprox->negative();
  vector<Polyhedron *> allRotations;
  allRotations.push_back(reflectedOuterApprox);
  for (int i=0; i< numRotations; i++)
    allRotations.push_back(rotate(allRotations[i]));
  cout<<"found all rotations"<<endl;


  const char *mSumName[] = { "sum00",
		       "sum01",
		       "sum02",
		       "sum03",
		       "sum04",
		       "sum05",
		       "sum06",
		       "sum07",
		       "sum08",
		       "sum09",
		       "sum10",
		       "sum11",
		       "sum12",
		       "sum13",
		       "sum14",
		       "sum15",
		       "sum16",
		       "sum17",
		       "sum18",
		       "sum19",
		       "sum20",
		       "sum21",
		       "sum22",
		       "sum23",
		       "sum24",
		       "sum25",
		       "sum26",
		       "sum27",
		       "sum28",
		       "sum29",
		       "sum30",
		       "sum31",
		       "sum32",
		       "sum33",
		       "sum34",
		       "sum35",
		       "sum36",
		       "sum37",
		       "sum38",
		       "sum39",
		       "sum40",
		       "sum41",
		       "sum42",
		       "sum43",
		       "sum44",
		       "sum45",
		       "sum46",
		       "sum47",
		       "sum48",
		       "sum49" };

  for (int i=0; i< allRotations.size(); i++) {
    cout<<"minkowskiSum "<<i<<" of "<<allRotations.size()-1<<endl;
    Polyhedron * mSum = minkowskiSumFull(allRotations[i], obstacle);
    if (USE_SIMPLIFY) {
      simplify(mSum, 1e-6);
      bool selfInt = mSum->intersectsEdges(mSum);
      cout << "selfInt " << selfInt << endl;
      savePolyTmp(mSum, mSumName[i]);
    }
    Polyhedron * mSum_complement = bb->boolean(mSum, Complement);
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
    // savePolyTmp(blockspaces[i], "blockspacesi");
    // savePolyTmp(blockspaces[j], "blockspacesj");
    Polyhedron * block_union = blockspaces[i]->boolean(blockspaces[j], Union);
    if (COMPUTE_USING_BLOCK_SPACE) {
      block_union->computeWindingNumbers();
      // cout << "union has " << block_union->cells.size() << " cells" << endl;
      char s[50];
      sprintf(s, "%d-%d", i, j);
      // saveWithShells(block_union, s);
      // savePolyTmp(block_union, "blockunionij");

      for (int k=0; k<block_union->cells.size(); k++) {
        bool valid = (block_union->cells[k]->getWN() == 0);
        if (valid) {
          PTR<Point> p = pointInCell(block_union, k);
	  PV3 pp = p->getApprox(1e-6);
          // cout<<"  "<<k<<": "<<pp.x.mid()<<" "<<pp.y.mid()<<" "<<pp.z.mid()<<endl;
          int ci = blockspaces[i]->containingCell(p);
          int cj = blockspaces[j]->containingCell(p);
          FreeSpace::Node * ni = findOrAddNode(i, ci);
          FreeSpace::Node * nj = findOrAddNode(j, cj);
          edges.push_back(new Edge(ni, nj, p));
        }
      }

    } else {
      Polyhedron * intersection = bb->boolean(block_union, Complement);
      intersection->computeWindingNumbers();
      char s[50];
      sprintf(s, "%d-%d", i, j);
      // saveWithShells(intersection, s);

      for (int k=0; k<intersection->cells.size(); k++) {
        bool valid = (intersection->cells[k]->getWN() == 1);
        if (valid) {
          PTR<Point> p = pointInCell(intersection, k);
          // cout<<"  "<<k<<": "<<p->getP().getX().mid()<<" "<<p->getP().getY().mid()<<" "<<p->getP().getZ().mid()<<endl;
          int ci = cspaces[i]->containingCell(p);
          int cj = cspaces[j]->containingCell(p);
          FreeSpace::Node * ni = findOrAddNode(i, ci);
          FreeSpace::Node * nj = findOrAddNode(j, cj);
          edges.push_back(new Edge(ni, nj, p));
        }
      }
      delete intersection;
    }
    delete block_union;
  }
}
