/*
  ACP (Adaptive Controlled Precision) Library
  for robust computational geometry

  Copyright (c) 2013-07-15
  Victor Joseph Milenkovic
  University of Miami
  vjm@miami.edu
  Elisha Sacks
  Purdue University
  eps@cs.purdue.edu

  This file contains code described in

Robust Complete Path Planning in the Plane
Victor Milenkovic, Elisha Sacks, and Steven Trac
Proceedings of the Workshop on the Algorithmic Foundations of Robotics (WAFR)
pages 37-52, 2012

   This code is under development.
   It is free for use for academic and research purposes.
   Permission from the authors is required for commercial use.
*/

#include "poly.h"
using namespace std;
using namespace acp;

namespace acp {

bool Root2::polishing = false;

typedef Object<PPoly> Poly;

bool zerop (const Parameter &p)
{
  return p.sign() == 0;
}

int PolySolver::Degree::sign () {
  const PPoly &p = poly->get();
  int d = p.deg();
  while (d >= 0 && zerop(p.a[d]))
    d--;
  return d;
}
/*
int Poly::deg ()
{
  PPoly<Parameter> p = getApprox(1.0);
  while (!p.zero() && p.lc().sign(false) == 0)
    p.a.pop_back();
  return p.deg();
}
*/

Roots PolySolver::getRoots ()
{
  static int count;
  count++;
  int d = deg();
  if (d < 1)
    return Roots();
  if (d == 1) {
    Roots res;
    res.push_back(new LinearRoot(poly));
    return res;
  }
  if (d == 2)
    return quadraticRoots();
  double b = CauchyBound(poly).getApprox(1.0).ub();
  // DEBUG
  return getRoots(new Object<Parameter>(Parameter::constant(-b * 1.123456789)), 
                  new Object<Parameter>(Parameter::constant(b * 1.2143645987)));
}

Roots PolySolver::getRoots (Object<Parameter>* l, Object<Parameter>* u)
{
  int d = deg();
  if (d < 1)
    return Roots();
  if (d == 1)
    return linearRoots(l, u);
  if (d == 2)
    return quadraticRoots(l, u);
  return descartesRoots(l, u);
}

Roots PolySolver::linearRoots (Object<Parameter>* l, Object<Parameter>* u)
{
  Roots res;
  RootPTR r = new LinearRoot(poly);
  if (Order(l, r) == 1 && Order(r, u) == 1)  
    res.push_back(r);
  return res;
}

Roots PolySolver::quadraticRoots ()
{
  Roots res;
  if (QuadraticRoots(poly) == 1) {
    RootPTR r1 = new QuadraticRoot(poly, false), r2 = new QuadraticRoot(poly, true);
    if (Order(r1, r2) == 1) {
      res.push_back(r1);
      res.push_back(r2);
    }
    else {
      res.push_back(r2);
      res.push_back(r1);
    }
  }
  return res;
}

Roots PolySolver::quadraticRoots (Object<Parameter>* l, Object<Parameter>* u)
{
  Roots res;
  if (QuadraticRoots(poly) == 1) {
    RootPTR r1 = new QuadraticRoot(poly, false), r2 = new QuadraticRoot(poly, true);
    if (Order(l, r1) == 1 && Order(r1, u) == 1)
      res.push_back(r1);
    if (Order(l, r2) == 1 && Order(r2, u) == 1)
      res.push_back(r2);
    if (res.size() == 2 && Order(res[0], res[1]) == -1)
      reverse(res.begin(), res.end());
  }
  return res;
}

Roots PolySolver::descartesRoots (Object<Parameter>* l, Object<Parameter>* u)
{
  Roots res;
  vector<DInt> st;
  st.push_back(DInt(l, u, 4));
  int k = 0;
  while (!st.empty()) {
    ++k;
    DInt i = *st.rbegin();
    st.pop_back();
    int v = Descartes(poly, i.l, i.u);
    if (v == 1)
      res.push_back(new PolyRoot(poly, i.l, i.u));
    else if (v > 1 // && !descartes1(st, i, v, true) && !descartes1(st, i, v, false)
	     && !descartes2(st, i, v)) {
      PTR<Object<Parameter> > m = new SubKnot(i.l, i.u, 1, 2);
      unsigned int n = max(4u, (unsigned int) sqrt(i.n));
      st.push_back(DInt(i.l, m, n));
      st.push_back(DInt(m, i.u, n));
    }
  }
  //cerr << "descartes iterations = " << k << endl;
  reverse(res.begin(), res.end());

  if (res.size() == 0)
    return res;

#ifdef USE_STURM
  Object<Parameter> *l_ub = new Object<Parameter>(Parameter::constant(l->get().ub()), true);
  Object<Parameter> *u_lb = new Object<Parameter>(Parameter::constant(u->get().lb()), true);
  Object<Parameter> *prev = l_ub;
  Object<Parameter> *next;

  for (int i = 0; i < res.size(); i++) {
    if (i < res.size() - 1)
      next = new Object<Parameter>(Parameter::constant((res[i]->get().ub() + res[i+1]->get().lb()) * 0.5), true);
    else
      next = u_lb;
    dynamic_cast<PolyRoot*>((Root*)res[i])->updateKnots(prev, next);
    prev = next;
  }
#else
  for (int i = 0; i < res.size(); i++) {
    dynamic_cast<PolyRoot*>((Root*)res[i])->updateKnots();
  }
#endif

  return res;
}

void PolyRoot::updateKnots () {
  PTR<Object<Parameter> > l2 = 
    new Object<Parameter>(Parameter::constant(l->get().mid()), true);
  PTR<Object<Parameter> > u2 = 
    new Object<Parameter>(Parameter::constant(u->get().mid()), true);
  updateKnots(l2, u2);

#ifdef BLEEN
  double ld = l->get().mid();
  double ud = u->get().mid();
  double rd = getApprox().mid();

  ld = ld - (ld - rd) / 7;
  ud = ud - (ud - rd) / 7;

  PTR<Object<Parameter> > lo = 
    new Object<Parameter>(Parameter::constant(ld), true);
  PTR<Object<Parameter> > uo = 
    new Object<Parameter>(Parameter::constant(ud), true);

  int des = 0;
  assert(::BaseObject::throwSignException() == false);
  ::BaseObject::throwSignException() = true;
  try {
    des = Descartes(p, lo, uo);
  }
  catch (SignException se) {}
  ::BaseObject::throwSignException() = false;
  assert(des == 1);

  bool lStep = true, uStep = true;
  while (lStep || uStep) {
    if (lStep) {
      double ld2 = ld + (ld - rd) / 2;
      PTR<Object<Parameter> > lo2 = 
        new Object<Parameter>(Parameter::constant(ld2), true);
      des = 0;
      ::BaseObject::throwSignException() = true;
      try {
        des = Descartes(p, lo2, uo);
      }
      catch (SignException se) {}
      ::BaseObject::throwSignException() = false;
      if (des == 1) {
        ld = ld2;
        lo = lo2;
        cout << "lStep " << lo->get().mid() << " " << uo->get().mid() << endl;
        if (ld < -1000)
          lStep = false;
      }
      else
        lStep = false;
    }
    if (uStep) {
      double ud2 = ud + (ud - rd) / 2;
      PTR<Object<Parameter> > uo2 = 
        new Object<Parameter>(Parameter::constant(ud2), true);
      des = 0;
      ::BaseObject::throwSignException() = true;
      try {
        des = Descartes(p, lo, uo2);
      }
      catch (SignException se) {}
      ::BaseObject::throwSignException() = false;
      if (des == 1) {
        ud = ud2;
        uo = uo2;
        cout << "uStep " << lo->get().mid() << " " << uo->get().mid() << endl;
        if (ud2 > 1000)
          uStep = false;
      }
      else
        uStep = false;
    }
  }
  cout << "updateKnots" << endl;
  updateKnots(lo, uo);
#endif
}  

bool PolySolver::descartes1 (vector<DInt> &st, const DInt &i, int v, bool lflag)
{
  PTR<Object<Parameter> > m = new SubKnot(i.l, i.u, lflag ? 1 : i.n - 1, i.n),
    nl = lflag ? i.l : m, nu = lflag ? m : i.u;
  if (Descartes(poly, nl, nu) != v)
    return false;
  unsigned long int ln = i.n, nn = min(ln*ln, 4294967295ul);
  st.push_back(DInt(nl, nu, nn));
  return true;
}

bool PolySolver::descartes2 (vector<DInt> &st, const DInt &i, int v)
{
  int k = descartes2int(i, v);
  if (k == 0)
    return false;
  PTR<Object<Parameter> > nl = new SubKnot(i.l, i.u, k - 2, 4*i.n),
    nu = new SubKnot(i.l, i.u, k + 2, 4*i.n);
  if (Descartes(poly, nl, nu) != v)
    return false;
  unsigned long int ln = i.n, nn = min(ln*ln, 4294967295ul);
  st.push_back(DInt(nl, nu, nn));
  return true;
}

int PolySolver::descartes2int (const DInt &i, int v)
{
  static int count;
  PPoly p = getApprox(1.0);
  Parameter il = i.l->getApprox(1.0);
  Parameter iu = i.u->getApprox(1.0);
  Parameter x = 0.5*(il + iu);
  Parameter dp = p.der(x);
  if (dp.sign(false) == 0)
    return 0;
  Parameter nx = x - v*p.value(x)/dp, w = (nx - il)/(iu - il);
  int k = min(max(4*i.n*w.mid(), 2.0), 4*i.n - 2.0);
  return k;
}

Parameter QuadraticRoots::calculate ()
{
  const PPoly &q = p->get();
  Parameter d = q[1]*q[1] - 4.0*q[0]*q[2];
  return d;
}

int Descartes::sign ()
{ 
  PPoly q = p->get().moebius(l->get(), u->get());
  int v = 0u, d = q.deg();
  int s = q.a[d].sign();
  for (int i = d - 1u; i >= 0; --i)
    if (q.a[i].sign() == - s) {
      ++v;
      s = - s;
    }
  return v;
}

int Sturm::sign () { 
#ifdef BLEEN
  static int count;
  if (++count == 0)
    cout << "this is it" << endl;
#endif
  vector<PPoly> ss;
  ss.push_back(p->get());
  ss.push_back(ss.back().der());
  while (ss.back().deg() > 0) {
    PPoly q;
    ss.push_back(ss[ss.size()-2].rem(ss[ss.size()-1], q));
    vector<Parameter> &a = ss.back().a;
    for (int i = 0; i < a.size(); i++)
      a[i] = -a[i];
  }
  assert(ss.back().deg() == 0);

  Parameter a = l->get();
  int na = 0;
  int sa = ss[0].value(a).sign();
  for (int i = 1; i < ss.size(); i++)
    if (ss[i].value(a).sign() == -sa) {
      na++;
      sa = -sa;
    }

  Parameter b = u->get();
  int nb = 0;
  int sb = ss[0].value(b).sign();
  for (int i = 1; i < ss.size(); i++)
    if (ss[i].value(b).sign() == -sb) {
      nb++;
      sb = -sb;
    }

  assert(na >= nb);
  return na - nb;
}

PV2 Root2::polish(PV2 p) {
 
  /*
  vector<PV2> b;
  b.push_back(p);
  vector< ObjPTR<PPoly2> > func;
  func.push_back(new Object<PPoly2>(ff));
  func.push_back(new Object<PPoly2>(gg));
  */
  //debugDraw(b, func);

  polishing = true;

  //printf("------------------------------------------------------------------------------------------\n\n");

  RootBoundary rb = {.I = p, .v = vector< vector<RootPTR> >()};

  //run Newton's until it fails then try to subdivide
  do {
    newton(rb);
  } while(subdivide(rb));

  polishing = false;

  //printf("I [%.16lf, %.16lf] [%.16lf, %.16lf]\n", rb.I.x.lb(), rb.I.x.ub(), rb.I.y.lb(), rb.I.y.ub());

  //printf("------------------------------------------------------------------------------------------\n\n");

  return rb.I;
}

//sign() == -1 if a < b
Primitive2(LessThan, Object<Parameter>* , a, Object<Parameter>* , b);
int LessThan::sign() {
  return (a->get() - b->get()).sign();
}

vector<int> Root2::intersections(vector<vector<RootPTR> > &rbv) {
 
  vector< vector<RootPTR> > fv;
  vector< vector<RootPTR> > gv;

  for(int i = 0; i < 4; i++) fv.push_back(rbv[i]);
  for(int i = 4; i < 8; i++) gv.push_back(rbv[i]);

  vector<int> alt;

  for(int k = 0; k < 4; k++) {

    vector<RootPTR> frr = fv[k];
    vector<RootPTR> grr = gv[k];

    for(int i = 0, j = 0; i < frr.size() || j < grr.size(); ) {

      RootPTR fr = (i < frr.size() ? frr[i] : nullptr);
      RootPTR gr = (j < grr.size() ? grr[j] : nullptr);

      //k < 2 is bottom and right which are ascending roots, k >= 2 are top and left which are descending
      if(fr != nullptr && (gr == nullptr || (k < 2 ? LessThan(fr, gr) == -1 : LessThan(fr, gr) == 1))) {
        alt.push_back(0);
        i++;
      } else {
        alt.push_back(1);
        j++;
      }

    }

  }

  return alt;

}

int Root2::parity(const vector<int> &alt) {
 
  //find the first f intersection
  int k = 0;
  int j = 0;
  while(j < alt.size()) {
    if(alt[j] == 0) break;
    j++;
  }
  
  //if f doesn't intersect this cell, keep going
  if(j == alt.size()) {
    return 0;
  }

  int count = 0;
  k = (j+1)%alt.size();

  vector<int> par;

  //compute x1, x2, x3, x4, ...
  //count the number of g intersections inbetween pairs of f intersections
  while(k != j) {
    if(alt[k] == 1) {
      count++;
    } else {
      par.push_back(count);
      count = 0;
    }

    k = (k+1)%alt.size();
  }

  //don't forget the last pair
  par.push_back(count);

  //parity of intersections == parity of sum of odd indexed xj's
  int parity = 0;
  for(j = 1; j < par.size(); j += 2) {
    parity += par[j];
  }

  return parity % 2;

}

class SubXPoly : public Object<PPoly> {
  PTR<Poly2> f;
  PTR<Object<Parameter> > x;
  PPoly calculate () { return f->get().subX(x->get()); }
public:
  SubXPoly (Poly2* f, Object<Parameter>* x) : f(f), x(x) {}
};

class SubYPoly : public Object<PPoly> {
  PTR<Poly2> f;
  PTR<Object<Parameter> > y;
  PPoly calculate () { return f->get().subY(y->get()); }
public:
  SubYPoly (Poly2* f, Object<Parameter>* y) : f(f), y(y) {}
};

void Root2::setRoots(RootBoundary &rb) {
  PV2 b = rb.I;

  rb.v.clear();

  PTR<Object<Parameter> > xl = new Object<Parameter>(b.x.lbP());
  PTR<Object<Parameter> > xu = new Object<Parameter>(b.x.ubP());
  PTR<Object<Parameter> > yl = new Object<Parameter>(b.y.lbP());
  PTR<Object<Parameter> > yu = new Object<Parameter>(b.y.ubP());
  
  PTR<Object<PPoly> > bot_f   = new SubYPoly(f, yl);
  PTR<Object<PPoly> > right_f = new SubXPoly(f, xu);
  PTR<Object<PPoly> > top_f   = new SubYPoly(f, yu);
  PTR<Object<PPoly> > left_f  = new SubXPoly(f, xl);

  rb.v.push_back(PolySolver(bot_f).getRoots(xl, xu));
  rb.v.push_back(PolySolver(right_f).getRoots(yl, yu));
  rb.v.push_back(PolySolver(top_f).getRoots(xl, xu));
  rb.v.push_back(PolySolver(left_f).getRoots(yl, yu));

  PTR<Object<PPoly> > bot_g   = new SubYPoly(g, yl);
  PTR<Object<PPoly> > right_g = new SubXPoly(g, xu);
  PTR<Object<PPoly> > top_g   = new SubYPoly(g, yu);
  PTR<Object<PPoly> > left_g  = new SubXPoly(g, xl);

  rb.v.push_back(PolySolver(bot_g).getRoots(xl, xu));
  rb.v.push_back(PolySolver(right_g).getRoots(yl, yu));
  rb.v.push_back(PolySolver(top_g).getRoots(xl, xu));
  rb.v.push_back(PolySolver(left_g).getRoots(yl, yu));

  //top and left need to be reversed for parity check
  reverse(rb.v[2].begin(), rb.v[2].end());
  reverse(rb.v[3].begin(), rb.v[3].end());
  reverse(rb.v[6].begin(), rb.v[6].end());
  reverse(rb.v[7].begin(), rb.v[7].end());

}

Root2::RootBoundary Root2::splitHoriz(RootBoundary &rb, Parameter c) {

  PTR<Object<Parameter> > s = new Object<Parameter>(c);
  
  vector< vector<RootPTR> > bot(8);
  vector< vector<RootPTR> > top(8);

  //bot side for bot f and g roots are same
  bot[0] = rb.v[0];
  bot[4] = rb.v[4];

  //top side for top f and g roots are the same 
  top[2] = rb.v[2];
  top[6] = rb.v[6];

  
  //Separate roots around the split point, give left to left and right to right
  //for f [0, 4) and g [4, 8)
  for(int i = 0; i < 8; i += 4) {

    int li = 0;
    while(li < rb.v[i+1].size() && LessThan(rb.v[i+1][li], s) == -1) { li++; }
    //[0, li) are < s [li, size) are >= s


    int ui = 0;
    while(ui < rb.v[i+3].size() && LessThan(rb.v[i+3][ui], s) == 1) { ui++; }
    //[0, ui) are > s [ui, size) are <= s  

    bot[i+1] = vector<RootPTR>(rb.v[i+1].begin(), rb.v[i+1].begin()+li);
    top[i+1] = vector<RootPTR>(rb.v[i+1].begin()+li, rb.v[i+1].end());

    top[i+3] = vector<RootPTR>(rb.v[i+3].begin(), rb.v[i+3].begin()+ui);
    bot[i+3] = vector<RootPTR>(rb.v[i+3].begin()+ui, rb.v[i+3].end());

  }

  PTR<Object<Parameter> > xl = new Object<Parameter>(rb.I.x.lbP());
  PTR<Object<Parameter> > xu = new Object<Parameter>(rb.I.x.ubP());
 
  PTR<Poly> mid_line_f = new SubYPoly(f, s);
  PTR<Poly> mid_line_g = new SubYPoly(g, s);

  vector<RootPTR> mid_roots_f_asc = PolySolver(mid_line_f).getRoots(xl, xu);
  vector<RootPTR> mid_roots_g_asc = PolySolver(mid_line_g).getRoots(xl, xu);

  vector<RootPTR> mid_roots_f_desc(mid_roots_f_asc.begin(), mid_roots_f_asc.end());
  vector<RootPTR> mid_roots_g_desc(mid_roots_g_asc.begin(), mid_roots_g_asc.end());

  reverse(mid_roots_f_desc.begin(), mid_roots_f_desc.end());
  reverse(mid_roots_g_desc.begin(), mid_roots_g_desc.end());

  top[0] = mid_roots_f_asc;
  top[4] = mid_roots_g_asc;

  bot[2] = mid_roots_f_desc;
  bot[6] = mid_roots_g_desc;

  PV2 I = rb.I;

  Parameter yl = I.y.lbP();
  Parameter yu = I.y.ubP();

  rb = {.I = PV2(I.x, Parameter(c, yu)), .v = top};

  return {.I = PV2(I.x, Parameter(yl, c)), .v = bot};

}

Root2::RootBoundary Root2::splitVert(RootBoundary &rb, Parameter c) {

  PTR<Object<Parameter> > s = new Object<Parameter>(c);
  
  vector< vector<RootPTR> > left(8);
  vector< vector<RootPTR> > right(8);

  //left side for left f and g roots are same
  left[3] = rb.v[3];
  left[7] = rb.v[7];

  //right side for right f and g roots are the same 
  right[1] = rb.v[1];
  right[5] = rb.v[5];

  
  //Separate roots around the split point, give left to left and right to right
  //for f [0, 4) and g [4, 8)
  for(int i = 0; i < 8; i += 4) {
    int li = 0;
    while(li < rb.v[i+0].size() && LessThan(rb.v[i+0][li], s) == -1) { li++; }
    //[0, li) are < s [li, size) are >= s


    int ui = 0;
    while(ui < rb.v[i+2].size() && LessThan(rb.v[i+2][ui], s) == 1) { ui++; }
    //[0, ui) are > s [ui, size) are <= s  

    left[i+0] = vector<RootPTR>(rb.v[i+0].begin(), rb.v[i+0].begin()+li);
    right[i+0] = vector<RootPTR>(rb.v[i+0].begin()+li, rb.v[i+0].end());

    right[i+2] = vector<RootPTR>(rb.v[i+2].begin(), rb.v[i+2].begin()+ui);
    left[i+2] = vector<RootPTR>(rb.v[i+2].begin()+ui, rb.v[i+2].end());
  }

  PTR<Object<Parameter> > yl = new Object<Parameter>(rb.I.y.lbP());
  PTR<Object<Parameter> > yu = new Object<Parameter>(rb.I.y.ubP());
 
  PTR<Poly> mid_line_f = new SubXPoly(f, s);
  PTR<Poly> mid_line_g = new SubXPoly(g, s);

  vector<RootPTR> mid_roots_f_asc = PolySolver(mid_line_f).getRoots(yl, yu);
  vector<RootPTR> mid_roots_g_asc = PolySolver(mid_line_g).getRoots(yl, yu);

  vector<RootPTR> mid_roots_f_desc(mid_roots_f_asc.begin(), mid_roots_f_asc.end());
  vector<RootPTR> mid_roots_g_desc(mid_roots_g_asc.begin(), mid_roots_g_asc.end());

  reverse(mid_roots_f_desc.begin(), mid_roots_f_desc.end());
  reverse(mid_roots_g_desc.begin(), mid_roots_g_desc.end());

  left[1] = mid_roots_f_asc;
  left[5] = mid_roots_g_asc;

  right[3] = mid_roots_f_desc;
  right[7] = mid_roots_g_desc;

  PV2 I = rb.I;

  Parameter xl = I.x.lbP();
  Parameter xu = I.x.ubP();

  rb = {.I = PV2(Parameter(c, xu), I.y), .v = right};

  return {.I = PV2(Parameter(xl, c), I.y), .v = left};;;;

}

void foo () {
  bool hse = BaseObject::throwSignException();
}

//subdivide the cell once along major axis
bool Root2::subdivide(RootBoundary &rb) {

  PV2 I = rb.I;

  //printf("SUBDIVIDE\n\n");

  bool hse = ::BaseObject::throwSignException();
  ::BaseObject::throwSignException() = true;
  unsigned int hp = ::BaseObject::getPrecision();

  int even_count = 0;
  int odd_count  = 0;
  int odd_index = -1;

  vector<RootBoundary> b;

  try {

    if(rb.v.size() == 0) {
      setRoots(rb);
    }

    if(I.x.intervalWidth() > I.y.intervalWidth()) {

      Parameter xl = I.x.lbP();
      Parameter xu = I.x.ubP();
      Parameter xl75 = 0.75*xl + 0.25*xu;
      Parameter xl25 = 0.25*xl + 0.75*xu;
   
      b.push_back(splitVert(rb, xl75));
      b.push_back(splitVert(rb, xl25));
      b.push_back(rb);

    } else {
      Parameter yl = I.y.lbP();
      Parameter yu = I.y.ubP();
      Parameter yl75 = 0.75*yl + 0.25*yu;
      Parameter yl25 = 0.25*yl + 0.75*yu;

      b.push_back(splitHoriz(rb, yl75));
      b.push_back(splitHoriz(rb, yl25));
      b.push_back(rb);

    }

    for(int it = 0; it < b.size(); it++) {

      //vector<PV2> list;
      //list.push_back(bb);

      //vector< ObjPTR<PPoly2> > fs;
      //fs.push_back(new Object<PPoly2>(f));
      //fs.push_back(new Object<PPoly2>(g));

      //debugDraw(list, fs);

      if(parity(intersections(b[it].v))) {
        odd_count++;
        odd_index = it;
      } else {
        even_count++;
      }

    }

  } catch(SignException e) {
    ::BaseObject::setPrecision(hp);
    ::BaseObject::throwSignException() = hse;
    return false;
  }
  
  ::BaseObject::throwSignException() = hse;

  assert(odd_count >= 1);
  //if(odd_count < 1)
  //  return false;

  if(odd_count >= 1) {
    //I = b[odd_index];
    rb = b[odd_index];
    return true;
  } else {
    return false;
  }

}

void Root2::newton(RootBoundary &rb) {

  PV2 I = rb.I;

  PPoly2 ff = f->get();
  PPoly2 gg = g->get();
  PPoly2 fx = ff.derX();
  PPoly2 fy = ff.derY();
  PPoly2 gx = gg.derX();
  PPoly2 gy = gg.derY();

  int count = 0;

  do {

    PV2 y(I.x.midP(), I.y.midP());
    
    PV2 foy(ff.value(&y.x), gg.value(&y.x));

    PV2 fr(fx.value(&I.x), fy.value(&I.x));
    PV2 gr(gx.value(&I.x), gy.value(&I.x));

    PV2 x;


    //let this be an exception in the end code
    if(!solve(fr, gr, foy, x)) {
      //printf("NEWTON determinant 0, need to subdivide\n\n");
      break;
    }

    PV2 newI = y - x;

    //printf("Newton Result [%20.16g, %20.16g] [%20.16g, %20.16g]\n\n", newI.x.lb(), newI.x.ub(), newI.y.lb(), newI.y.ub());


    //this should just be an assertion, not exception in the end code
    assert(newI.x.intersects(I.x) && newI.y.intersects(I.y));

    if(!newI.x.intersects(I.x) || !newI.y.intersects(I.y)) {
      //printf("NEWTON result disjoint\n\n");
      break;
    }

    /*
    vector<PV2> b;
    b.push_back(I);
    b.push_back(newI);
    vector< ObjPTR<PPoly2> > fs;
    fs.push_back(new Object<PPoly2>(f));
    fs.push_back(new Object<PPoly2>(g));
    */
    //debugDraw(b, fs);
 
    PV2 inter(newI.x.intersect(I.x), newI.y.intersect(I.y));

    //this is a termination condition
    //if the new interval, (result of intersecting the result and old interval, is not a proper subset of the old interval
    //if(!inter.x.subset(I.x) && !inter.y.subset(I.y))

    //this is less efficient
    //if(inter.x.intervalWidth() >= I.x.intervalWidth() && inter.y.intervalWidth() >= I.y.intervalWidth()) {
    if(!inter.x.subset(I.x) && !inter.y.subset(I.y)) {
      //printf("NEWTON no progress\n\n");
      break;
    }
 /*
    int pI = parity(intersections(I));
    int pIn = parity(intersections(inter));

    if(pI == 0 || pIn == 0) {
      printf("here\n");
    }
*/
    I = inter;
  
    rb = {.I = I, .v = vector< vector<RootPTR> >()};

    //printf("      I       [%20.16g, %20.16g] [%20.16g, %20.16g]\n\n", I.x.lb(), I.x.ub(), I.y.lb(), I.y.ub());
    
    //printf("NEWTON STEP\n\n");

  } while(true);
  //} while(++count < 100);

  //return I;

}

bool Root2::solve(PV2 fr, PV2 gr, PV2 b, PV2 &sol) {

  //pivot on the smallest element!
  //compare two things: sing(false) for comparison

  //if det contains zero, need to subdivide
  Parameter d = fr.x*gr.y - fr.y*gr.x;
  if(!d.sign(false)) return false;

  if(!fr.x.sign(false) || !gr.x.sign(false)) return false;
 
  //pivot first row 
  if( (fr.x.abs() - gr.x.abs()).sign(false) >= 0) {  

    Parameter c = gr.x / fr.x;

    gr.x = gr.x - c*fr.x;
    gr.y = gr.y - c*fr.y;
    b.y  = b.y - c*b.x;

    if(!gr.y.sign(false)) return false;
    Parameter y = b.y / gr.y;

    if(!fr.x.sign(false)) return false;
    Parameter x = (b.x - fr.y*y) / fr.x;

    sol = PV2(x, y); 

  } 
  //pivot on the second row
  else {

    Parameter c = fr.x / gr.x;

    fr.x = fr.x - c*gr.x;
    fr.y = fr.y - c*gr.y;
    b.x = b.x - c*b.y;

    if(!fr.y.sign(false)) return false;
    Parameter y = b.x / fr.y;

    if(!gr.x.sign(false)) return false;
    Parameter x = (b.y - gr.y*y) / gr.x;

    sol = PV2(x, y);

  }

  return true;
}

}

// debug

void pp (Object<Parameter> *p)
{
  cerr << p->getApprox().mid() << endl;
}

void pp (Object<PPoly> *p)
{
  p->getApprox(1.0).print();
}

