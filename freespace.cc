#include "freespace.h"
#include "simplify.h"
#include <map>

Polyhedron * loadPoly2(const char * filename) {
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

void savePoly2(Polyhedron * p, const char * filename) {
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
    this->p = p;
    r = new RotationPoint(p, sin_cos_alpha);
    t = new TangentIntersectionPoint(p, sin_cos_alpha);
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
    setInnerOuter();
  }
  OuterApproxFace(OuterApproxVertex* bottom1, OuterApproxVertex* bottom2, OuterApproxVertex* top1, OuterApproxVertex* top2) {
    this->bottom1 = bottom1;
    this->bottom2 = (bottom1 != bottom2? bottom2 : NULL);
    this->top1 = top1;
    this->top2 = (top1 != top2? top2 : NULL);
    isTrapazoid = (this->top2 != NULL) && (this->bottom2 != NULL);
    setInnerOuter();
  }
  void setInnerOuter() {
    if (top2 != NULL && DiffLength(top1->p, top2->p) > 0) {
      OuterApproxVertex * tmp = top1;
      top1 = top2;
      top2 = tmp;
    }
    if (bottom2 != NULL && DiffLength(bottom1->p, bottom2->p) > 0) {
      OuterApproxVertex * tmp = bottom1;
      bottom1 = bottom2;
      bottom2 = tmp;
    }
  }
};

Primitive2(DiffZ, PTR<Point>, i, PTR<Point>, j);
int DiffZ::sign() { return (i->get().z - j->get().z).sign(); }
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

Primitive4(InOrderVertices, PTR<Point>, p0, PTR<Point>, p1, PTR<Point>, p2, PTR<Point>, p3);
int InOrderVertices::sign() {
  PV3 n1  = (p1->get()-p0->get()).cross(p2->get()-p0->get());
  PV3 n2  = (p2->get()-p1->get()).cross(p3->get()-p1->get());

  if (n1.dot(n2) > 0) {
    return 1;
  }
  return 0;
}

void saveTriangulation(const char * filename, std::vector<OuterApproxFace> & splitTList) {
  ofstream out;
  out.open(filename);


  int ntrapezoids = 0;
  for (std::vector<OuterApproxFace>::iterator it=splitTList.begin(); it!=splitTList.end(); ++it)
    if (it->isTrapazoid)
      ntrapezoids++;

  out << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS "<<3*(splitTList.size()-ntrapezoids)+4*ntrapezoids<<" double"<<endl;

  for (std::vector<OuterApproxFace>::iterator it=splitTList.begin(); it!=splitTList.end(); ++it) {
    //if the polyhedron is a trapezoid, we need to make sure we're printing them in the right order
    bool swap = it->isTrapazoid && (InOrderVertices(it->top1->p, it->top2->p, it->bottom1->p, it->bottom2->p) == 0);

    out << it->top1->p->getApprox().getX().mid() << " "
        << it->top1->p->getApprox().getY().mid() << " "
        << it->top1->p->getApprox().getZ().mid() << endl;
    if (it->top2 != NULL)
      out << it->top2->p->getApprox().getX().mid() << " "
          << it->top2->p->getApprox().getY().mid() << " "
          << it->top2->p->getApprox().getZ().mid() << endl;
    if (!swap)
      out << it->bottom1->p->getApprox().getX().mid() << " "
          << it->bottom1->p->getApprox().getY().mid() << " "
          << it->bottom1->p->getApprox().getZ().mid() << endl;
    if (it->bottom2 != NULL)
      out << it->bottom2->p->getApprox().getX().mid() << " "
          << it->bottom2->p->getApprox().getY().mid() << " "
          << it->bottom2->p->getApprox().getZ().mid() << endl;
    if (swap)
      out << it->bottom1->p->getApprox().getX().mid() << " "
          << it->bottom1->p->getApprox().getY().mid() << " "
          << it->bottom1->p->getApprox().getZ().mid() << endl;
  }

  out << "POLYGONS "<<splitTList.size()<<" "<<4*(splitTList.size()-ntrapezoids)+5*ntrapezoids<<endl;
  int index = 0;
  for (std::vector<OuterApproxFace>::iterator it=splitTList.begin(); it!=splitTList.end(); ++it) {
    out << (it->isTrapazoid? 4 : 3) << " ";
    out << index++ << " " << index++ << " " << index++;
    if (it->isTrapazoid)
      out << " " << index++;
    out << endl;
  }

  out.close();
}

Primitive4(FrontDiagonal, Vertex*,top1,Vertex*,top2,Vertex*,bottom1,Vertex*,bottom2);
int FrontDiagonal::sign() {
  //return -1 if the front diagonal is top2 bottom1
  //return  1 if the front diagonal is top1 bottom2
  //note: top1 top2 and bottom1 bottom2 are in coutnerclockwise order around their respective pentagons
  PV3 t1 = top1->getP()->get();
  PV3 t2 = top2->getP()->get();
  PV2 vt(t2.getX()-t1.getX(), t2.getY()-t1.getY());
  PV3 b1 = bottom1->getP()->get();
  PV3 b2 = bottom2->getP()->get();
  PV2 vb(b2.getX()-b1.getX(), b2.getY()-b1.getY());
  return vb.cross(vt).sign();
}

Primitive3(Area2D, PTR<Point>, a, PTR<Point>, b, PTR<Point>, c);
int Area2D::sign() {
  PV3 v1 = b->get()-a->get();
  PV3 v2 = c->get()-a->get();
  return (v1.getX()*v2.getY() - v1.getY() * v2.getX()).sign();
}

//find the top and bottom pentagon's for trapezoid 2 and connect each of the sides by their convex hulls 
Polyhedron * trapezoidOuterApprox(OuterApproxFace t) {
  bool topIsTriangle = (t.top2 == NULL);
  bool bottomIsTriangle = (t.bottom2 == NULL);

  Polyhedron *a = new Polyhedron;

  Vertex * tp = a->getVertex(t.top1->p);
  Vertex * tq = topIsTriangle? tp : a->getVertex(t.top2->p);
  Vertex * tq_tanint = topIsTriangle? a->getVertex(t.top1->t) : a->getVertex(t.top2->t);
  Vertex * tp_rot = a->getVertex(t.top1->r);
  Vertex * tq_rot = topIsTriangle? tp_rot : a->getVertex(t.top2->r);

  Vertex * bp = a->getVertex(t.bottom1->p);
  Vertex * bq = bottomIsTriangle? bp : a->getVertex(t.bottom2->p);
  Vertex * bq_tanint = bottomIsTriangle? a->getVertex(t.bottom1->t) : a->getVertex(t.bottom2->t);
  Vertex * bp_rot = a->getVertex(t.bottom1->r);
  Vertex * bq_rot = bottomIsTriangle? bp_rot : a->getVertex(t.bottom2->r);

  //top pentagon, triangulated
  if (topIsTriangle) {
    a->addTriangle(tq, tq_tanint, tq_rot);
  } else {
    Vertices polygon; 
    polygon.push_back(tp); polygon.push_back(tq); polygon.push_back(tq_tanint); polygon.push_back(tq_rot); polygon.push_back(tp_rot);
    while (polygon.size() > 3) {
      int n=polygon.size();
      for (int i=0; i<n; i++) {
        bool isEar = Area2D(polygon[i]->getP(), polygon[(i+1)%n]->getP(), polygon[(i+2)%n]->getP()) > 0;
        if (isEar)
          for (int j=3; j<n; j++)
            if (Area2D(polygon[i]->getP(),       polygon[(i+1)%n]->getP(), polygon[(i+j)%n]->getP()) > 0  &&
                Area2D(polygon[(i+1)%n]->getP(), polygon[(i+2)%n]->getP(), polygon[(i+j)%n]->getP()) > 0  &&
                Area2D(polygon[(i+2)%n]->getP(), polygon[i]->getP(),       polygon[(i+j)%n]->getP()) > 0) {
              isEar = false;
              break;
            }
        if (isEar) {
          a->addTriangle(polygon[i], polygon[(i+1)%n], polygon[(i+2)%n]);
          polygon.erase(polygon.begin() + (i+1)%n ); 
          break;     
        }
      }
      assert(polygon.size() == n-1);
    }
    a->addTriangle(polygon[0], polygon[1], polygon[2]);
  }

  //bottom pentagon, triangulated
  if (bottomIsTriangle) {
    a->addTriangle(bq_rot, bq_tanint, bq);
  } else {
    Vertices polygon;
    polygon.push_back(bp); polygon.push_back(bq); polygon.push_back(bq_tanint); polygon.push_back(bq_rot); polygon.push_back(bp_rot);
    while (polygon.size() > 3) {
      int n=polygon.size();
      for (int i=0; i<n; i++) {
        bool isEar = Area2D(polygon[i]->getP(), polygon[(i+1)%n]->getP(), polygon[(i+2)%n]->getP()) > 0;
        if (isEar)
          for (int j=3; j<n; j++)
            if (Area2D(polygon[i]->getP(),       polygon[(i+1)%n]->getP(), polygon[(i+j)%n]->getP()) > 0  &&
                Area2D(polygon[(i+1)%n]->getP(), polygon[(i+2)%n]->getP(), polygon[(i+j)%n]->getP()) > 0  &&
                Area2D(polygon[(i+2)%n]->getP(), polygon[i]->getP(),       polygon[(i+j)%n]->getP()) > 0) {
              isEar = false;
              break;
            }
        if (isEar) {
          a->addTriangle(polygon[i], polygon[(i+2)%n], polygon[(i+1)%n]);
          polygon.erase(polygon.begin() + (i+1)%n ); 
          break;     
        }
      }
      assert(polygon.size() == n-1);
    }
    a->addTriangle(polygon[0], polygon[2], polygon[1]);
  }


  //p q side (is parallel, diagonal doesn't matter)
  if (!bottomIsTriangle)
    a->addTriangle(tp, bp, bq);
  if (!topIsTriangle)
    a->addTriangle(tp, bq, tq);

  //q t side
  if (FrontDiagonal(tq, tq_tanint, bq, bq_tanint) < 0) {
    a->addTriangle(tq, bq, tq_tanint);
    a->addTriangle(bq, bq_tanint, tq_tanint);
  } else {
    a->addTriangle(tq, bq, bq_tanint);
    a->addTriangle(tq, bq_tanint, tq_tanint);
  }

  //t q' side
  if (FrontDiagonal(tq_tanint, tq_rot, bq_tanint, bq_rot) < 0) {
    a->addTriangle(tq_tanint, bq_tanint, tq_rot);
    a->addTriangle(bq_tanint, bq_rot, tq_rot);
  } else {
    a->addTriangle(tq_tanint, bq_tanint, bq_rot);
    a->addTriangle(tq_tanint, bq_rot, tq_rot);
  }

  //q' p' side (parallel, doesn't matter which diagonal)
  if (!topIsTriangle)
    a->addTriangle(tq_rot, bq_rot, tp_rot);
  if (!bottomIsTriangle)
    a->addTriangle(bq_rot, bp_rot, tp_rot);

  //p' p side
  if (FrontDiagonal(tp_rot, tp, bp_rot, bp) < 0) {
    a->addTriangle(tp_rot, bp_rot, tp);
    a->addTriangle(bp_rot, bp, tp);
  } else {
    a->addTriangle(tp_rot, bp_rot, bp);
    a->addTriangle(tp_rot, bp, tp);
  }

  //DEBUG 
  a->computeWindingNumbers();
  if (a->cells.size() != 2) {
    cout << "trapezoids outer approximation has "<<a->cells.size()<<" cells"<<endl;
    std::vector<OuterApproxFace> tmp;
    tmp.push_back(t);
    saveTriangulation("trap.vtk", tmp);
    savePoly2(a, "trapApprox"); 
    assert(false);
  }
  return a;
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
  for (Faces::const_iterator f = p->faces.begin(); f != p->faces.end(); ++f) {
    Vertices ve = (*f)->getBoundary()->loop();
    a->addTriangle(ve[0]->getP(), ve[1]->getP(), ve[2]->getP(), pvmap);
  }
  return a;
}

Primitive3(RequiresVerticalSplit, OuterApproxFace, f, PTR<Object<Parameter> >, sin2theta, PTR<Object<Parameter > >, Rtheta);
int RequiresVerticalSplit::sign() {
  Parameter dh = f.top1->p->get().getZ() - f.bottom1->p->get().getZ();
  if (dh > Rtheta->get()) {
    return 1;
  }
  //p p' side
  PV3 t1 = f.top1->p->get();
  PV3 t2 = f.top1->r->get();
  PV3 b1 = f.bottom1->p->get();
  PV3 b2 = f.bottom1->r->get();
  PV3 vt = t2-t1;
  PV3 vb = b2-b1;
  Parameter vtxvb = vt.getX()*vb.getY() - vt.getY()*vb.getX();
  Parameter sin2alpha = vtxvb*vtxvb / vt.dot(vt) / vb.dot(vb);
  if (sin2alpha > sin2theta->get())
    return 1;

  //q t side
  t1 = (f.top2 != NULL? f.top2->p->get() : f.top1->p->get());
  t2 = (f.top2 != NULL? f.top2->t->get() : f.top1->t->get());
  b1 = (f.bottom2 != NULL? f.bottom2->p->get() : f.bottom1->p->get());
  b2 = (f.bottom2 != NULL? f.bottom2->t->get() : f.bottom1->t->get());
  vt = t2-t1;
  vb = b2-b1;
  vtxvb = vt.getX()*vb.getY() - vt.getY()*vb.getX();
  Parameter vt2 = vt.dot(vt);
  Parameter vb2 = vb.dot(vb);
  sin2alpha = vtxvb*vtxvb / vt2 / vb2;
  if (sin2alpha > sin2theta->get()) 
    return 1;

  //q' t side (QUESTION: identical to q t? hence redundant)
  t1 = (f.top2 != NULL? f.top2->r->get() : f.top1->r->get());
  b1 = (f.bottom2 != NULL? f.bottom2->r->get() : f.bottom1->r->get());
  vt = t2-t1;
  vb = b2-b1;
  vtxvb = vt.getX()*vb.getY() - vt.getY()*vb.getX();
  sin2alpha = vtxvb*vtxvb / vt.dot(vt) / vb.dot(vb);
  if (sin2alpha > sin2theta->get())
    return 1;

  return 0;
}

void splitVertically(OuterApproxFace f, std::vector<OuterApproxFace> & list, PTR<Object<Parameter> > sin2theta, PTR<Object<Parameter> > Rtheta) {
  if (RequiresVerticalSplit(f, sin2theta, Rtheta) == 1) {
    OuterApproxVertex* mid1 = new OuterApproxVertex(new MidPoint(f.top1->p, f.bottom1->p));
    OuterApproxVertex* mid2 = new OuterApproxVertex(new MidPoint((f.top2 != NULL? f.top2->p : f.top1->p), (f.bottom2 != NULL? f.bottom2->p : f.bottom1->p)));

    OuterApproxFace top(mid1, mid2, f.top1, f.top2);
    OuterApproxFace bottom(f.bottom1, f.bottom2, mid1, mid2);
    splitVertically(top, list, sin2theta, Rtheta);
    splitVertically(bottom, list, sin2theta, Rtheta);
  } else
    list.push_back(f);
}

Polyhedron * FreeSpace::generateSweep(std::vector<OuterApproxFace> & trapezoids) {
  std::vector<Polyhedron *> polyList;
  for (int i=0; i<trapezoids.size(); i++)
    polyList.push_back(trapezoidOuterApprox(trapezoids[i]));

  cout<<"found polyhedral outer approximations of "<<polyList.size()<<" trapezoids"<<endl;
  Polyhedron *temp = multiUnion(&polyList[0], polyList.size()),
    *outerApproxShell = coalesce(temp);
  delete temp;
  cout <<"done multiUnion"<<endl;
  Polyhedron * robotApprox = inner_approximation? robot->boolean(outerApproxShell, Union) : robot->boolean(outerApproxShell, Complement);
  for (int i=0; i<polyList.size(); i++)
    delete polyList[i];
  polyList.clear();
  Polyhedron *reflectedRobotApprox = robotApprox->negative();
  delete robotApprox, outerApproxShell;
  return reflectedRobotApprox;
}

void FreeSpace::generateFreeSpaces(std::vector<Polyhedron*> & spaces, Polyhedron * sweep, const char * basename) {
  std::vector<Polyhedron *> allRotations;
  allRotations.push_back(sweep);
  for (int i=0; i< numRotations; i++)
    allRotations.push_back(rotate(allRotations[i]));
  cout<<"found all rotations"<<endl;

  for (int i=0; i< allRotations.size(); i++) {
    cout<<"minkowskiSum "<<i<<" of "<<allRotations.size()-1<<endl;
    Polyhedron * mSum = minkowskiSum(allRotations[i], obstacle);
    // simplify(mSum, 1e-6, false);
    char s[25];
    sprintf(s, "%s%02d", basename, i); 
    savePoly2(mSum, s);
    spaces.push_back(mSum);
  }

  cout<<"cspaces populated"<<endl;

  for (int i=0; i< numRotations; i++)
    delete allRotations[i];
  allRotations.clear(); 
}

FreeSpace::FreeSpace(Polyhedron * robot, Polyhedron * obstacle, PTR<Object<Parameter> > tan_half_angle, int numRotations, bool inner_approximation) {
  this->robot = robot;
  this->obstacle = obstacle;
  this->inner_approximation = inner_approximation;
  this->numRotations = numRotations;

  sin_cos_alpha = new SinCosAlpha(tan_half_angle);

  std::vector<SimpleTriangle> tList;
  for (int i=0; i<this->robot->faces.size(); i++) {
    HEdge *es = this->robot->faces[i]->getBoundary();
    PTR<Point> p = es->tail()->getP();
    PTR<Point> q = es->getNext()->tail()->getP();
    PTR<Point> r = es->getNext()->getNext()->tail()->getP();
    tList.push_back(SimpleTriangle(p,q,r));
  }
  cout<<"loaded triangles"<<endl;

  std::map<PTR<Point>, OuterApproxVertex*> pvmap;
  std::vector<OuterApproxFace> splitTList;
  for (int i=0; i< tList.size(); i++)
    split(splitTList, tList[i], pvmap);
  tList.clear();

  cout<<"split triangles"<<endl;
  saveTriangulation("triangulation.vtk", splitTList);

  //find R * theta and sin2theta
  PTR<Point> furthestVertex = robot->vertices[0]->getP();
  for (int i=1; i<robot->vertices.size(); i++) {
    if (DiffLength(robot->vertices[i]->getP(), furthestVertex) > 0)
      furthestVertex = robot->vertices[i]->getP();
  }
  double r2 = furthestVertex->getApprox().getX().mid() * furthestVertex->getApprox().getX().mid() 
            + furthestVertex->getApprox().getY().mid() * furthestVertex->getApprox().getY().mid() 
            + furthestVertex->getApprox().getZ().mid() * furthestVertex->getApprox().getZ().mid();
  Parameter::disable();
  double rt = sqrt(r2) * 2 * M_PI / numRotations;
  double sin2t = sin(2 * M_PI / numRotations) * sin(2 * M_PI / numRotations);
  Parameter::enable();
  PTR<Object<Parameter> > Rtheta = new InputParameter(rt);
  PTR<Object<Parameter> > sin2theta = new InputParameter(sin2t);


  //split trapezoid horizontally to ensure excess is O(theta^2)
  std::vector<OuterApproxFace> splitTList2;
  for (std::vector<OuterApproxFace>::iterator f = splitTList.begin(); f != splitTList.end(); f++) {
    splitVertically(*f, splitTList2, sin2theta, Rtheta);
  }
  cout <<splitTList2.size() << " trapezoids" << endl;

  saveTriangulation("triangulation2.vtk", splitTList2);

  Polyhedron * closeApprox = generateSweep(splitTList2);
  savePoly2(closeApprox, "closeApprox");
  // Polyhedron * roughApprox = generateSweep(splitTList);
  // savePoly2(roughApprox, "roughApprox");

  time_t start,end;

  time (&start);
  generateFreeSpaces(blockspaces_close, closeApprox, "close");
  time (&end);
  double closetime = difftime (end,start);

  time (&start);
  // generateFreeSpaces(blockspaces_rough, roughApprox, "rough");
  time (&end);
  double roughtime = difftime (end,start);

  cout << "close approximate free spaces generated in " << closetime/60.0 << " minutes" << endl;
  cout << "rough approximate free spaces generated in " << roughtime/60.0 << " minutes" << endl;
}
