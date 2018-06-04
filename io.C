#include "io.h"

Polyhedron * readPolyhedron (istream &istr, bool perturbed)
{
  Polyhedron *a = new Polyhedron(perturbed);
  int nv, nf;
  istr >> nv >> nf;
  Parameter::disable();
  for (int i = 0; i < nv; ++i) {
    double x, y, z;
    istr >> x >> y >> z;
    a->getVertex(x, y, z);
  }
  Parameter::enable();
  for (int i = 0; i < nf; ++i) {
    int u, v, w;
    istr >> u >> v >> w;
    Vertex *uv = a->vertices[u], *vv = a->vertices[v], *wv = a->vertices[w];
    a->addTriangle(uv, vv, wv);
  }
  istr >> ws;
  if (istr.peek() != EOF)
    readCells(istr, a);
  return a;
}

void readCells (istream &istr, Polyhedron *a)
{
  int nc;
  istr >> nc;
  for (int i = 0; i < nc; ++i)
    a->cells.push_back(readCell(istr, a));
}

Cell * readCell (istream &istr, Polyhedron *a)
{
  int ns;
  istr >> ns;
  Shells sh;
  for (int i = 0; i < ns; ++i)
    sh.push_back(readShell(istr, a));
  int i0 = a->cells.empty() ? 0 : 1;
  Cell *c = new Cell(i0 == 0 ? 0 : sh[0]);
  for (int i = i0; i < ns; ++i)
    c->addInner(sh[i]);
  return c;
}

Shell * readShell (istream &istr, Polyhedron *a)
{
  int nf;
  istr >> nf;
  HFaces hf;
  for (int i = 0; i < nf; ++i) {
    int f, h = 0;
    istr >> f;
    if (f < 0) {
      f = - f;
      h = 1;
    }    
    hf.push_back(a->faces[f-1]->getHFace(h));
  }
  return new Shell(hf);
}

void writePolyhedron (Polyhedron *a, ostream &ostr)
{
  VIMap vimap;
  ostr << setprecision(17);
  ostr << a->vertices.size() << " " << a->faces.size() << endl << endl;
  for (int i = 0; i < a->vertices.size(); ++i) {
    Vertex *v = a->vertices[i];
    vimap.insert(VIPair(v, i));
    PV3 p = v->getP()->getApprox();
    ostr << p.x.mid() << " " << p.y.mid() << " " << p.z.mid() << endl;
  }
  ostr << endl;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f) {
    Vertices ve;
    (*f)->boundaryVertices(ve);
    ostr << vimap.find(ve[0])->second << " " << vimap.find(ve[1])->second
	 << " " << vimap.find(ve[2])->second << endl;
  }
  ostr << endl;
  if (!a->cells.empty())
    writeCells(a, ostr);  
}

void writeCells (Polyhedron *a, ostream &ostr)
{
  FIMap fimap;
  for (int i = 0; i < a->faces.size(); ++i)
    fimap.insert(FIPair(a->faces[i], i + 1));
  int nc = a->cells.size();
  ostr << nc << endl << endl;
  for (Cells::iterator c = a->cells.begin(); c != a->cells.end(); ++c) {
    int ns = (*c)->nShells();
    ostr << ns << endl;
    for (int j = 0; j < ns; ++j)
      writeShell((*c)->getShell(j), fimap, ostr);
    ostr << endl;
  }
}

void writeShell (Shell *s, FIMap &fimap, ostream &ostr)
{
  const HFaces &hf = s->getHFaces();
  ostr << hf.size() << endl;
  for (HFaces::const_iterator h = hf.begin(); h != hf.end(); ++h) {
    int s = (*h)->pos() ? 1 : -1, f = fimap.find((*h)->getF())->second;
    ostr << s*f << " ";
  }
  ostr << endl;
}

Polyhedron * readPolyhedronVTK (istream &istr, bool perturbed)
{
  Polyhedron *a = new Polyhedron(perturbed);
  skipComments(istr);
  string dummy;
  do {
    istr >> dummy;
  }
  while (!(dummy == "POINTS" || dummy == "points"));
  int nv;
  istr >> nv >> dummy;
  Parameter::disable();
  for (int i = 0; i < nv; ++i) {
    double x, y, z;
    istr >> x >> y >> z;
    a->getVertex(x, y, z);
  }
  Parameter::enable();
  int nf;
  istr >> dummy >> nf >> dummy;
  for (int i = 0; i < nf; ++i) {
    int k, u, v, w;
    istr >> k >> u >> v >> w;
    assert(k == 3);
    Vertex *uv = a->vertices[u], *vv = a->vertices[v], *wv = a->vertices[w];
    a->addTriangle(uv, vv, wv);
  }
  return a;
}

void skipComments (istream &istr)
{
  char s[10000];
  ws(istr);
  while (!istr.eof() && (istr.peek() == '#' || istr.peek() == '\n'))
    istr.getline(s, 1000);
}

void writePolyhedronVTK (const Polyhedron *a, ostream &ostr)
{
  writePolyhedronVTK(a->faces, ostr);
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
  outputVTK(pts, data, 1, ostr);
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

void outputVTK (const PV3s &pts, const IVector &data, int ptype,
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
  ostr << endl;
  if (ptype == 1) ostr << "POLYGONS";
  else if (ptype == 2) ostr << "LINES";
  else ostr << "VERTICES";
  ostr << " " << np << " " << data.size() << endl;
  i = 0;
  while (i < data.size()) {
    ostr << data[i] << " ";
    for (int j = 0; j < data[i]; ++j)
      ostr << data[i+j+1] << " ";
    ostr << endl;
    i += data[i] + 1;
  }
}

// debug

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
  outputVTK(pts, data, 1, ostr);
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
  outputVTK(pts, data, 2, ostr);
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
  outputVTK(pts, data, 1, ostr);
}

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
  outputVTK(pts, data, 2, ostr);
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

void pvertices (const Vertices &ve, int i)
{
  VIMap vidmap;
  PV3s pts;
  IVector data;
  data.push_back(ve.size());
  for (Vertices::const_iterator v = ve.begin(); v != ve.end(); ++v)
    data.push_back(getPoint(vidmap, pts, *v));
  ofstream ostr(outi(i).c_str());
  outputVTK(pts, data, 3, ostr);
}
