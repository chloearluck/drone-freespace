#include "sumfaces.h"

const bool DEBUG = false;

void save(PTR<Feature> feature, const char * filename) {
  ofstream out;
  out.open(filename);
  out << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS "<<feature->p.size()<<" double"<<endl;
  for (int i=0; i<feature->p.size(); i++)
    out << feature->p[i]->getApprox().getX().mid() << " " <<
           feature->p[i]->getApprox().getY().mid() << " " <<
           feature->p[i]->getApprox().getZ().mid() << endl;
  out << "POLYGONS "<<"1 "<<feature->p.size()+1<<endl;
  out << feature->p.size() << " 0 1 2" << (feature->p.size() == 4? " 3" : "") << endl;
  out.close();
}

void save(Edge * e, const char * filename) {
  ofstream out;
  out.open(filename);
  out << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS 2 double"<<endl;
  out << e->getH()->getP()->getApprox().getX().mid() << " "
      << e->getH()->getP()->getApprox().getY().mid() << " "
      << e->getH()->getP()->getApprox().getZ().mid() << endl;
  out << e->getT()->getP()->getApprox().getX().mid() << " "
      << e->getT()->getP()->getApprox().getY().mid() << " "
      << e->getT()->getP()->getApprox().getZ().mid() << endl;
  out << "LINES 1 3"<< endl;
  out << "2 0 1" << endl;
  out.close();
}

void save(Face * f, const char * filename) {
  ofstream out;
  out.open(filename);
  out << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS 3 double"<<endl;
  Vertices vs;
  f->boundaryVertices(vs);
  assert(vs.size() == 3);
  for (int i=0; i<3; i++)
    out << vs[i]->getP()->getApprox().getX().mid() << " "
        << vs[i]->getP()->getApprox().getY().mid() << " "
        << vs[i]->getP()->getApprox().getZ().mid() << endl;
  out << "POLYGONS 1 4"<<endl;
  out << "3 0 1 2" << endl;
  out.close();
}

void debug(PTR<Feature> feature, Edge * e, Face * f, int angle_range, int debug_msg) {
  if (random()/((double)RAND_MAX) < 0.05 && debug_msg != 0) {
    cout<<"angle range: "<<angle_range<<endl;
    save(feature, "sumface.vtk");
    if (e != NULL)
      save(e, "block.vtk");
    if (f != NULL)
      save(f, "block.vtk");

    if (debug_msg == 0)
      cout << "sumface vertex outside block space" << endl;
    else if (debug_msg == 1)
      cout << "blockspace edge pieces sumface" << endl;
    else if (debug_msg == 2)
      cout << "sumface edge pierce a blockspace face" << endl;
    else if (debug_msg == 3)
      cout << "blockspace edge pieces sumface (sumface not a triangle)" << endl;
    cout<< endl;
  }
}

Primitive5(FaceEdgeIntersect, PTR<Point>, tail, PTR<Point>, head, PTR<Point>, pa, PTR<Point>, pb, PTR<Point>, pc);
int FaceEdgeIntersect::sign() {
  PV3 a = pa->getP();
  PV3 b = pb->getP();
  PV3 c = pc->getP();
  PV3 t = tail->getP();
  PV3 v = head->getP()-t;

  PV3 n = (b-a).cross(c-a);
  Parameter k = a.dot(n);

  Parameter s = (k-t.dot(n))/v.dot(n);

  if (s < 0 || s > 1)
    return 0;

  PV3 q = t + s*v;

  bool sideab = vol(q, a, b, b+n).sign() == vol(c, a, b, b+n).sign();//q is on c's side of ab
  bool sideac = vol(q, a, c, c+n).sign() == vol(b, a, c, c+n).sign();
  bool sidebc = vol(q, b, c, c+n).sign() == vol(a, b, c, c+n).sign();

  if (sideab && sideac && sidebc)
    return 1;
  else
    return 0;
}

bool * candidatePairs(PTR<Feature> fobs, PTR<Feature> frob, PTR<Point> sin_cos_alpha, PTR<Point> sin_cos_alpha_sample, std::vector<Polyhedron*> & blockspaces, int n_samples) {
  bool * candidate = new bool[blockspaces.size()];

  for (int i=0; i<blockspaces.size(); i++) {
    candidate[i] = false;
    frob = frob->rotate(sin_cos_alpha);
    PTR<Feature> g = frob;

    for (int j=0; j< n_samples; j++) {
      if (j>0)
        g = g->rotate(sin_cos_alpha_sample);

      PTR<Feature> sumface = g->sum(fobs);

      if (!blockspaces[i]->contains(sumface->p[0]) ||
          !blockspaces[i]->contains(sumface->p[1]) ||
          !blockspaces[i]->contains(sumface->p[2]) ||
          (sumface->p.size() == 4 && !blockspaces[i]->contains(sumface->p[3]))) {
        candidate[i] = true;
        debug(sumface, NULL, NULL, i, 0);
        break;
      }

      //for every edge of the block space, does that edge pierce sumface
      PTR<Plane> sumfaceplane = new TrianglePlane(sumface->p[0], sumface->p[1], sumface->p[2]);
      for (Edges::iterator edge = blockspaces[i]->edges.begin(); edge != blockspaces[i]->edges.end(); edge++) {
        if (FaceEdgeIntersect((*edge)->getT()->getP(), (*edge)->getH()->getP(),
                              sumface->p[0], sumface->p[1], sumface->p[2]) == 1) {
          candidate[i] = true;
          debug(sumface, (*edge), NULL, i, 1);
          break;
        }

        if (sumface->p.size() == 4 &&
            FaceEdgeIntersect((*edge)->getT()->getP(), (*edge)->getH()->getP(),
                              sumface->p[0], sumface->p[2], sumface->p[3]) == 1) {
          candidate[i] = true;
          debug(sumface, (*edge), NULL, i, 3);
          break;
        }
      }

      if (candidate[i])
        break;

      //for ever face of the block space, does any edge of sumface pierce it
      for (Faces::iterator face=blockspaces[i]->faces.begin(); face != blockspaces[i]->faces.end(); face++) {
        Vertices vs;
        (*face)->boundaryVertices(vs);
        assert(vs.size() == 3);
        PTR<Plane> plane = new TrianglePlane(vs[0]->getP(), vs[1]->getP(), vs[2]->getP());

        for (int j=0; j<sumface->p.size(); j++) {
          PTR<Point> head = sumface->p[j];
          PTR<Point> tail = sumface->p[(j+1)%sumface->p.size()];

          if (FaceEdgeIntersect(tail, head, vs[0]->getP(), vs[1]->getP(), vs[2]->getP()) == 1) {
            candidate[i] = true;
            debug(sumface, NULL, (*face), i, 2);
            break;
          }
        }
      }

      if (candidate[i])
        break;
    }
  }
  return candidate;
}
