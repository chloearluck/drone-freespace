#include "io.h"

Polyhedron * readPolyhedron (istream &istr, bool tflag, bool check)
{
  Polyhedron *a = new Polyhedron;
  readVertices(istr, a);
  readAttributes(istr, a);
  FaceRecords frs;
  readFaceRecords(istr, frs, a);
  formFaces(frs, a, tflag, check);
  return a;
}

void readVertices (istream &istr, Polyhedron *a)
{
  Parameter::disable();
  skipComments(istr);
  int nv;
  istr >> nv;
  for (int i = 0; i < nv; ++i) {
    skipComments(istr);
    ID id;
    double x, y, z;
    istr >> id >> x >> y >> z;
    a->getVertex(x, y, z, false);
  }
  Parameter::enable();
}

void readAttributes (istream &istr, Polyhedron *a)
{
  skipComments(istr);
  int n;
  istr >> n;
  for (int i = 0; i < n; ++i) {
    ID id;
    string name, value;
    istr >> id >> name >> value;
    a->attributes.push_back(Attribute(id, name, value));
  }
}

void skipComments (istream &istr)
{
  char s[10000];
  ws(istr);
  while (!istr.eof() && (istr.peek() == '#' || istr.peek() == '\n'))
    istr.getline(s, 1000);
}

void readFaceRecords (istream &istr, FaceRecords &frs, Polyhedron *a)
{
  skipComments(istr);
  int nf, nh;
  istr >> nf;
  for (int i = 0; i < nf; ++i)
    frs.push_back(readFaceRecord(istr));
  skipComments(istr);
  istr >> nh;
  for (int i = 0; i < nh; ++i) {
    ID hid;
    int fid, cid, na;
    istr >> hid >> fid >> na >> cid;
    FaceRecord &fr = frs[abs(fid)-1];
    for (int j = 0; j < na; ++j) {
      ID aid;
      istr >> aid;
      if (aid == 0 || aid > a->attributes.size())
	continue;
      Attribute &att = a->attributes[aid-1];
      if (att.name == "offset") {
	if (fid > 0)
	  fr.o1 = true;
	else
	  fr.o2 = true;
      }
    }
  }
}

FaceRecord readFaceRecord (istream &istr)
{
  skipComments(istr);
  ID id;
  int nb;
  istr >> id >> nb;
  FaceRecord fs;
  for (int i = 0; i < nb; ++i) {
    skipComments(istr);
    int nv, k;
    istr >> nv;
    IVector iv;
    for (int j = 0; j < nv; ++j) {
      istr >> k;
      iv.push_back(k - 1);
    }
    fs.b.push_back(iv);
  }
  return fs;
}

void formFaces (const FaceRecords &frs, Polyhedron *a, bool tflag, bool check)
{
  for (FaceRecords::const_iterator i = frs.begin(); i != frs.end(); ++i) {
    VVertices reg;
    faceVertices(i->b, a, reg);
    if (tflag || check && !facePlane(reg, a)) {
      Triangles tr;
      int c = projectionCoordinate(*reg[0]);
      triangulate(reg, c, tr);
      for (Triangles::iterator t = tr.begin(); t != tr.end(); ++t)
	a->addTriangle(t->a, t->b, t->c);
    }
    else
      a->addFace(reg);
    deleteRegion(reg);
  }
}

void faceVertices (const IVectors &iv, Polyhedron *a, VVertices &reg)
{
  for (IVectors::const_iterator i = iv.begin(); i != iv.end(); ++i) {
    Vertices *ve = new Vertices;
    reg.push_back(ve);
    for (IVector::const_iterator j = i->begin(); j != i->end(); ++j)
      ve->push_back(a->vertices[*j]);
  }
}

bool facePlane (const VVertices &reg, Polyhedron *a)
{
  Vertices ve;
  for (VVertices::const_iterator r = reg.begin(); r != reg.end(); ++r)
    ve.insert(ve.end(), (*r)->begin(), (*r)->end());
  int i = 0, n = ve.size();
  Plane *p = 0;
  while (i + 2 < n) {
    Vertex *v1 = ve[i], *v2 = ve[i+1], *v3 = 0;
    for (int j = i + 2; j < n; ++j)
      if (!ve[j]->getP()->onLine(v1->getP(), v2->getP())) {
	v3 = ve[j];
	i = j + 1;
	break;
      }
    if (!v3)
      break;
    Plane *q = new TrianglePlane(v1->getP(), v2->getP(), v3->getP());
    if (!p)
      p = q;
    else if (SamePlane(p, q) == 1)
      delete q;
    else {
      delete p;
      delete q;
      return false;
    }
  }
  delete p;
  return true;
}

int projectionCoordinate (const Vertices &ve)
{
  Vertex *v1 = ve[0], *v2 = ve[1], *v3 = 0;
  for (int i = 2; !v3 && i < ve.size(); ++i) {
    Vertex *vi = ve[i];
    if (!vi->getP()->onLine(v1->getP(), v2->getP()))
      v3 = vi;
  }
  TrianglePlane p(v1->getP(), v2->getP(), v3->getP());
  int c = ProjectionCoordinate(&p);
  return outerLoop(ve, c) ? c : - c;
}

void writePolyhedron (Polyhedron *a, ostream &ostr)
{
  VIMap vimap;
  ostr << setprecision(20) << a->vertices.size() << endl;
  for (unsigned int i = 0u; i < a->vertices.size(); ++i) {
    Vertex *v = a->vertices[i];
    PV3 p = v->getP()->getP();
    unsigned int j = i + 1u;
    ostr << j << " " << p.x.mid() << " " << p.y.mid() << " " << p.z.mid() << endl;
    vimap.insert(VIPair(v, j));
  }
  ostr << endl << a->attributes.size() << endl;
  for (Attributes::iterator i = a->attributes.begin(); i != a->attributes.end(); ++i)
    ostr << i->id << " " << i->name << " " << i->value << endl;
  ostr << endl << a->faces.size() << endl;
  for (unsigned int i = 0; i < a->faces.size(); ++i) {
    Face *f = a->faces[i];
    const HEdges &b = f->getBoundary();
    ostr << i + 1u << " " << b.size() << endl;
    for (HEdges::const_iterator e = b.begin(); e != b.end(); ++e) {
      Vertices ve;
      (*e)->loop(ve);
      ostr << ve.size() << endl;
      for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
	ostr << vimap.find(*v)->second << " ";
      ostr << endl;
    }
    ostr << endl;
  }
  writeHFaces(a, ostr);
}

void writeHFaces (Polyhedron *a, ostream &ostr)
{
  if (a->cells.empty())
    a->formCells();
  map<Cell *, int> cimap;
  for (int i = 0; i < a->cells.size(); ++i)
    cimap.insert(pair<Cell *, int>(a->cells[i], i));
  int n = a->faces.size();
  ostr << endl << 2*n << endl;
  for (int i = 0; i < a->faces.size(); ++i)
    for (int j = 0; j < 2; ++j) {
      HFace *f = a->faces[i]->getHFace(j);
      int k = 2*i+j+1, fid = j == 0 ? i + 1u : - (i + 1u),
	cid = cimap.find(f->twin()->getS()->getC())->second; // Intel convention
      ostr << k << " " << fid << " " << cid << " ";
    }
}

Polyhedron * readPolyhedronVTK (istream &istr, bool perturb)
{
  Parameter::disable();
  Polyhedron *a = new Polyhedron;
  skipComments(istr);
  string dummy;
  do {
    istr >> dummy;
  }
  while (!(dummy == "POINTS" || dummy == "points"));
  int nv;
  istr >> nv >> dummy;
  for (int i = 0; i < nv; ++i) {
    double x, y, z;
    istr >> x >> y >> z;
    a->getVertex(x, y, z, perturb);
  }
  Parameter::enable();
  int nt;
  istr >> dummy >> nt >> dummy;
  for (int i = 0; i < nt; ++i) {
    int k, u, v, w;
    istr >> k >> u >> v >> w;
    assert(k == 3);
    Vertex *uv = a->vertices[u], *vv = a->vertices[v], *wv = a->vertices[w];
    a->addTriangle(uv, vv, wv);
  }
  return a;
}

void writePolyhedronVTK (const Faces &fa, ostream &ostr)
{
  PV3s pts;
  IVector data;
  VIMap vimap;
  for (Faces::const_iterator f = fa.begin(); f != fa.end(); ++f) {
    const HEdges &fb = (*f)->getBoundary();
    if (!fb.empty()) {
      HEdge *e = fb[0], *e0 = e;
      Vertices ve;
      do {
  ve.push_back(e->tail());
  e = e->getNext();
      }
      while (e != e0);
      data.push_back(ve.size());
      for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
  data.push_back(getPoint(vimap, pts, *v));
    }
  }
  outputVTK(pts, data, true, ostr);
}

void writePolyhedronOBJ (const Faces &fa, ostream &ostr)
{
  PV3s pts;
  IVector data;
  VIMap vimap;
  for (Faces::const_iterator f = fa.begin(); f != fa.end(); ++f) {
    const HEdges &fb = (*f)->getBoundary();
    if (!fb.empty()) {
      HEdge *e = fb[0], *e0 = e;
      Vertices ve;
      do {
  ve.push_back(e->tail());
  e = e->getNext();
      }
      while (e != e0);
      data.push_back(ve.size());
      for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
  data.push_back(getPoint(vimap, pts, *v));
    }
  }
  outputOBJ(pts, data, true, ostr);
}

int getPoint (VIMap &vimap, PV3s &pts, Vertex *v)
{
  VIMap::iterator iter = vimap.find(v);
  if (iter != vimap.end())
    return iter->second;
  int k = pts.size();
  pts.push_back(v->getP()->getApprox());
  vimap.insert(VIPair(v, k));
  return k;
}

void outputVTK (const PV3s &pts, const IVector &data, bool pflag, 
		ostream &ostr)
{
  int nv = pts.size();
  ostr << setprecision(20) << "# vtk DataFile Version 3.0" << endl
       << "vtk output" << endl << "ASCII" << endl
       << "DATASET POLYDATA" << endl 
       << "POINTS " << nv << " double" << endl;
  for (PV3s::const_iterator p = pts.begin(); p != pts.end(); ++p) 
    ostr << p->getX().mid() << " " << p->getY().mid() << " " 
	 << p->getZ().mid() << endl;
  int np = 0, i = 0;
  while (i < data.size()) {
    ++np;
    i += data[i] + 1;
  }
  ostr << endl << (pflag ? "POLYGONS " : "LINES ") << np << " " << 
    data.size() << endl;
  i = 0;
  while (i < data.size()) {
    ostr << data[i] << " ";
    for (int j = 0; j < data[i]; ++j)
      ostr << data[i+j+1] << " ";
    ostr << endl;
    i += data[i] + 1;
  }
}

void outputOBJ (const PV3s &pts, const IVector &data, bool pflag, 
    ostream &ostr)
{
  int nv = pts.size();
  ostr <<setprecision(20) << "g"<<endl;
  for (PV3s::const_iterator p = pts.begin(); p != pts.end(); ++p) 
    ostr << "v "<< p->getX().mid() << " " << p->getY().mid() << " " << p->getZ().mid() << endl;
  
  int i = 0;
  while (i < data.size()) {
    ostr << "f ";
    for (int j = 0; j < data[i]; ++j)
      ostr << data[i+j+1]+1 << " ";
    ostr << endl;
    i += data[i] + 1;
  }
}

void writePolyhedronVTK (const HFaces &fa, ostream &ostr)
{
  PV3s pts;
  IVector data;
  VIMap vimap;
  for (HFaces::const_iterator f = fa.begin(); f != fa.end(); ++f) {
    Face *ff = (*f)->getF();
    const HEdges &fb = ff->getBoundary();
    if (!fb.empty()) {
      HEdge *e = fb[0], *e0 = e;
      Vertices ve;
      do {
	ve.push_back(e->tail());
	e = e->getNext();
      }
      while (e != e0);
      if (!(*f)->pos())
	reverse(ve.begin(), ve.end());
      data.push_back(ve.size());
      for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
	data.push_back(getPoint(vimap, pts, *v));
    }
  }
  outputVTK(pts, data, true, ostr);
}

Polyhedron * readPolyhedronSTL (istream &istr)
{
  Parameter::disable();
  Polyhedron *a = new Polyhedron;
  char cd[1000];
  istr.getline(cd, 1000);
  while (istr.peek() != EOF) {
    istr.getline(cd, 1000);
    istr.getline(cd, 1000);
    Vertex *v[3];
    for (int i = 0; i < 3; ++i) {
      string dummy;
      double x, y, z;
      istr >> dummy >> x >> y >> z;
      v[i] = a->getVertex(x, y, z, false);
    }
    istr.getline(cd, 1000);
    istr.getline(cd, 1000);
    istr.getline(cd, 1000);
    if (!v[0]->getP()->onLine(v[1]->getP(), v[2]->getP()))
      a->addTriangle(v[0], v[1], v[2]);
  }
  Parameter::disable();
  return a;
}

void ptriangles (const Triangles &tr, ostream &ostr)
{
  PV3s pts;
  IVector data;
  VIMap vimap;
  for (Triangles::const_iterator t = tr.begin(); t != tr.end(); ++t) {
    data.push_back(3);
    data.push_back(getPoint(vimap, pts, t->a));
    data.push_back(getPoint(vimap, pts, t->b));
    data.push_back(getPoint(vimap, pts, t->c));
  }
  outputVTK(pts, data, true, ostr);
}

string outi (int i)
{
  return "out" + static_cast<ostringstream*>( &(ostringstream() << i) )->str()
    + ".vtk";
}

void plines (const vector<PV3s> &lines, int i)
{
  PV3s pts;
  IVector data;
  int k = 0;
  for (vector<PV3s>::const_iterator l = lines.begin(); l != lines.end(); ++l) {
    data.push_back(l->size());
    for (PV3s::const_iterator p = l->begin(); p != l->end(); ++p, ++k) {
      pts.push_back(*p);
      data.push_back(k);
    }
  }
  ofstream ostr(outi(i).c_str());
  outputVTK(pts, data, false, ostr);
}

// debug

void ptriangles (const Triangles &tr, int i)
{
  ofstream ostr(outi(i).c_str());
  ptriangles(tr, ostr);
}

void pfaces (const Faces &faces, int i)
{
  ofstream ostr(outi(i).c_str());
  writePolyhedronVTK(faces, ostr);
}

void pfaces (Polyhedron *a, int i)
{
  pfaces(a->faces, i);
}

void pface (Face *f, int i)
{
  Faces fa;
  fa.push_back(f);
  pfaces(fa, i);
}

void pfaces (const FaceSet &fs, int i)
{
  Faces fa(fs.begin(), fs.end());
  pfaces(fa, i);
}

void pfaces (const HFaces &hfaces, int i)
{
  ofstream ostr(outi(i).c_str());
  writePolyhedronVTK(hfaces, ostr);
}

void pfaces (ID pid, const Faces &faces, int i)
{
  Faces fa;
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f)
    if ((*f)->getP()->getid() == pid)
      fa.push_back(*f);
  pfaces(fa, i);
}

void pcells (Polyhedron *a)
{
  if (a->cells.empty())
    a->formCells();
  for (int i = 1; i < a->cells.size(); ++i)
    pfaces(a->cells[i]->getShell(0)->getHFaces(), i);
}

void pfaces (const Faces &faces)
{
  writePolyhedronVTK(faces, cout);
}

void pfaces (const HFaces &hfaces)
{
  Faces faces;
  for (HFaces::const_iterator f = hfaces.begin(); f != hfaces.end(); ++f)
    faces.push_back((*f)->getF());
  writePolyhedronVTK(faces, cout);
}

void pface (Face *f)
{
  Faces faces;
  faces.push_back(f);
  writePolyhedronVTK(faces, cout);
}

void pedges (const Edges &edges, ostream &ostr)
{
  VIMap vidmap;
  PV3s pts;
  IVector data;
  for (Edges::const_iterator e = edges.begin(); e != edges.end(); ++e) {
    data.push_back(2);
    data.push_back(getPoint(vidmap, pts, (*e)->getT()));
    data.push_back(getPoint(vidmap, pts, (*e)->getH()));
  }
  outputVTK(pts, data, false, ostr);
}

void pedges (const Edges &edges)
{
  pedges(edges, cout);
}

void pedges (const Edges &edges, int i)
{
  ofstream ostr(outi(i).c_str());
  pedges(edges, ostr);
}

void pedges (const HEdges &hedges)
{
  Edges edges;
  for (HEdges::const_iterator e = hedges.begin(); e != hedges.end(); ++e)
    edges.push_back((*e)->getE());
  pedges(edges);
}

void pedge (Edge *e)
{
  Edges ed;
  ed.push_back(e);
  pedges(ed);
}

void pedge (HEdge *e)
{
  pedge(e->getE());
}

void pvertices (const Vertices &ve)
{
  VIMap vidmap;
  PV3s pts;
  IVector data;
  data.push_back(ve.size());
  for (Vertices::const_iterator v = ve.begin(); v != ve.end(); ++v)
    data.push_back(getPoint(vidmap, pts, *v));
  outputVTK(pts, data, false, cout);
}
