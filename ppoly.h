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

#ifndef PPOLY
#define PPOLY

#include "pv.h"
#include <vector>
#include <map>
#include <unordered_map>
using namespace std;
using namespace acp;

namespace acp {

class PPoly {
 public:
  PPoly () {}
  PPoly (double c) { add(Parameter::constant(c), 0); }
  PPoly (const vector<Parameter> &a) : a(a) {}

  int size () const { return a.size(); }
  const Parameter &operator[] (int i) const { return a[i]; }
  Parameter &operator[] (int i) { return a[i]; }
  Parameter lc () const { return *a.rbegin(); }
  int deg () const { return a.size() - 1; }
  bool zero () const { return a.empty(); }

  void add (const Parameter &c, int e) {
    int d = deg();
    if (e <= d)
      a[e] = a[e] + c;
    else {
      while (d + 1 < e) {
	a.push_back(Parameter::constant(0.0));
	++d;
      }
      a.push_back(c);
    }
  }

  PPoly operator+ (const PPoly &p) const {
    int d = deg(), e = p.deg();
    PPoly q;
    if (d <= e) {
      for (int i = 0; i <= d; ++i)
	q.add(a[i] + p.a[i], i);
      for (int i = d + 1; i <= e; ++i)
	q.add(p.a[i], i);
    }
    else {
      for (int i = 0; i <= e; ++i)
	q.add(a[i] + p.a[i], i);
      for (int i = e + 1; i <= d; ++i)
	q.add(a[i], i);
    }
    return q;
  }

  PPoly operator- (const PPoly &p) const { return (*this) + (- p); }

  PPoly operator- () const {
    PPoly q;
    for (int i = 0; i <= deg(); ++i)
      q.add(- a[i], i);
    return q;
  }

  PPoly operator* (const PPoly &p) const {
    int d = deg(), e = p.deg();
    PPoly q;
    for (int i = 0; i <= d; ++i)
      for (int j = 0; j <= e; ++j)
	q.add(a[i]*p.a[j], i + j);
    return q;
  }
  
  PPoly monic () const {
    PPoly p;
    int d = deg();
    Parameter ad = a[d];
    for (unsigned int i = 0u; i < d; ++i)
      p.a.push_back(a[i]/ad);
    p.a.push_back(Parameter::constant(1.0));
    return p;
  }

  PPoly quotient (const PPoly &divisor, PPoly &remainder) {
    if (divisor.deg() > deg()) {
      remainder = *this;
      return PPoly();
    }
    PPoly f = *this;
    PPoly g = divisor.monic();
    int gd = gd;
    PPoly q;
    q.a.resize(deg() - gd + 1);
    for (int i = q.deg(); i >= 0; i--) {
      q.a[i] = f.a[gd + i];
      for (int j = 0; j < gd; j++)
        f.a[i + j] = f.a[i + j] - q.a[i] * g.a[j];
    }
    f.a.resize(gd - 1);
    while (f.deg() >= 0 && f.a[f.deg()].sign() == 0)
      f.a.pop_back();
    remainder = f;
    return q;
  }

  PPoly rem (const PPoly &p, PPoly &q) const {
    Parameter lp = p.lc();
    int dp = p.deg();
    PPoly r(*this);
    while (dp <= r.deg()) {
      Parameter k = r.lc()/lp;
      int dk = r.deg() - dp;
      q.add(k, dk);
      r.a.pop_back();
      for (int i = 0; i < dp; ++i)
	r.add(- k*p.a[i], i + dk);
      r.checkDegree();
    }
    return r;
  }

  void checkDegree () {
    while (!zero() && lc().sign() == 0)
      a.pop_back();
  }

  PPoly gcd (const PPoly &p, PPoly &s, PPoly &t) const
  {
    s = PPoly(1);
    t = PPoly(0);
    PPoly r(*this), nr(p), ns, nt(1);
    r.checkDegree();
    nr.checkDegree();
    while (!nr.zero()) {
      PPoly q, nnr = r.rem(nr, q), nns = s - q*ns, nnt = t - q*nt;
      r = nr;
      nr = nnr;
      s = ns;
      ns = nns;
      t = nt;
      nt = nnt;
    }
    return r;
  }

  Parameter value (const Parameter &x) const
  {
    int d = deg();
    Parameter y = a[d];
    for (int i = d - 1; i >= 0; --i)
      y = x*y + a[i];
    return y;
  }

  Parameter der (const Parameter &x) const
  {
    int d = deg();
    Parameter y = d*a[d];
    for (int i = d - 1; i > 0; --i)
      y = x*y + i*a[i];
    return y;
  }

  PPoly der () const
  {
    PPoly g;
    for (int i = 1; i <= deg(); ++i)
      g.a.push_back(i*a[i]);
    return g;
  }
  
  PPoly neg () const
  {
    PPoly p;
    for (unsigned int i = 0u; i <= deg(); ++i)
      p.a.push_back(i%2 == 0 ? a[i] : - a[i]);
    return p;
  }

  PPoly dual () const
  {
    PPoly p;
    int d = deg();
    for (unsigned int i = 0u; i <= d; ++i)
      p.a.push_back(a[d-i]);
    return p;
  }

  Parameter shrinkBracket (Parameter x, int sub) const {
    bool bflag = false;
    while (true) {
      bool flag = true;
      Parameter xm = x.midP(), fm = value(xm), df = der(x);
      if (df.sign(false) == 0)
	flag = false;
      else {
	Parameter nx = (xm - fm/df).intersect(x);
	if (nx.subset(x))
	  x = nx;
	else
	  flag = false;
      }
      if (!flag) {
	int sm = fm.sign(false);
	if (sm == sub) {
	  Parameter nx = x.interval(xm);
	  if (!nx.subset(x))
	    break;
	  x = nx;
	}
	else if (sm == - sub) {
	  Parameter nx = xm.interval(x);
	  if (!nx.subset(x))
	    break;
	  x = nx;
	}
	else if (!bflag) {
	  bflag = true;
	  x = shrinkBracketB(x, sub);
	}
	else
	  break;
      }
    }
    return x;
  }

  Parameter shrinkBracketB (const Parameter &x, int sub) const {
    Parameter xlb = x.lbP(), xm = x.midP(), xub = x.ubP();
    while ((xlb - xm).sign(false) < 0) {
      Parameter nx = xlb.interval(xm).midP();
      if ((nx - xlb).sign(false) == 0 || (nx - xm).sign(false) == 0)
	break;
      int sm = value(nx).sign(false);
      if (sm == 0)
	xm = nx;
      else if (sm == sub)
	xub = nx;
      else
	xlb = nx;
    }
    xm = x.midP();
    while ((xm - xub).sign(false) < 0) {
      Parameter nx = xm.interval(xub).midP();
      if ((nx - xm).sign(false) == 0 || (nx - xub).sign(false) == 0)
	break;
      int sm = value(nx).sign(false);
      if (sm == 0)
	xm = nx;
      else if (sm == -sub)
	xlb = nx;
      else
	xub = nx;
    }
    Parameter nx = xlb.interval(xub);
    return nx.subset(x) ? nx : x;
  }

  PPoly moebius (const Parameter &l, const Parameter &u) const
  {
    return shift(l).dual().shift(1.0/(u - l));
  }

  PPoly shift (const Parameter &s) const
  {
    PPoly p = *this;
    if (s.sign() == 0)
      return p;
    int d = deg();
    for (int i = 0; i < d; ++i)
      for (int j = d - 1; j >= i; --j)
	p.a[j] = p.a[j] + s * p.a[j+1];
    return p;
  }

  PPoly shiftX (const Parameter &s) const
  {
    PPoly p = *this;
    if (s.sign() == 0)
      return p;
    p.scale(s);
    p.shift1();
    p.scale(1.0/s);
    return p;
  }

  void scale (const Parameter &s)
  {
    Parameter ss = s;
    for (int i = 1; i <= deg(); ++i) {
      a[i] = ss*a[i];
      ss = ss*s;
    }
  }

  void shift1 ()
  {
    int d = deg();
    for (int i = 0; i < d; ++i)
      for (int j = d - 1; j >= i; --j)
	a[j] = a[j] + a[j+1];
  }

  void print () const {
    int d = deg();
    cerr << "( ";
    for (int i = 0; i <= d; ++i)
      cerr << a[i].mid() << " ";
    cerr << ")" << endl;
  }

  vector<Parameter> a;
};

typedef std::pair<Parameter, Parameter> Int;
typedef std::vector<Int> Ints;
class Mobius;

class PPoly1 : public PPoly {
public:

  static bool print_all;
  int d;

  PPoly1 () : d(-1) {}
  PPoly1 (int d, const Parameter *a = 0) : d(d) {
    if (a == 0)
      add(Parameter::constant(0), d);
    else
      for (int i = 0; i <= d; i++)
        add(a[i], i);
  }
  PPoly1 (const Parameter &c, int d) : d(d) { add(c, d); }
  PPoly1 (const PPoly &ppoly) : PPoly(ppoly), d(ppoly.deg()) {}

  bool increased () const { return false; /* a[0].increased(); */ }
  bool decreased () const { return true; /* a[0].decreased(); */ }
  Parameter value (const Parameter &x) const;
  Parameter der (const Parameter &x) const;
  PPoly1 der () const;

  static PPoly1 XminusR (const Parameter &r) {
    Parameter a[2] = { -r, Parameter::constant(1) };
    return PPoly1(1, a);
  }

  std::vector<Parameter> roots (const Parameter &lb, const Parameter &ub) const
    { return rootsR(lb, ub); }
  
  Parameter polish (const Parameter &x) const;
  Parameter shrinkBracket (Parameter x, int sub) const;
  Parameter shrinkBracketB (const Parameter &x, int sub) const;
  std::vector<Parameter> rootsR (const Parameter &lb, const Parameter &ub) const;
  std::vector<Parameter> roots1 (const Parameter &lb, const Parameter &ub) const;
  std::vector<Parameter> roots2 (const Parameter &lb, const Parameter &ub) const;

  PPoly1 operator* (const PPoly1 &g) const;
  PPoly1 operator* (const Parameter &p) const;
  PPoly1 operator* (const double p) const;
  PPoly1 plusCtimes (double c, const PPoly1 &b) const;
  PPoly1 plusCtimes (const Parameter &c, const PPoly1 &b) const;
  PPoly1 operator+ (const PPoly1 &b) const { return plusCtimes(1.0, b); }
  PPoly1 operator- (const PPoly1 &b) const { return plusCtimes(-1.0, b); }

  PPoly1 karatsuba (const PPoly1 &g) const;

  PPoly1 normalized() const;

  PPoly1 quotient (const PPoly1 &b, PPoly1 &remainder);
  PPoly1 gcd (const PPoly1 &b);
  int degree () { return deg(); }

  void print();

  int size () const { return a.size(); }
  Parameter &operator[] (int i) { return a[i]; }

  PPoly1 monic () const;
  PPoly1 normal () const;
  PPoly1 shift (double s) const;
  PPoly1 shift (const Parameter &s) const;
  PPoly1 neg () const;
  PPoly1 dual () const;
  PPoly1 removeZeroRoots () const;
  unsigned int descartes () const;
  Parameter ubk () const;

  std::vector<Parameter> rootsD (const Parameter &lb, const Parameter &ub) const;
  Ints rootInts (const Parameter &lb, const Parameter &ub) const;
  void vas (PPoly1 p, Mobius m, Ints &ints) const;
  Int vasInt (const PPoly1 &p, const Mobius &m) const;
};

class PPoly2 {
 public:
  PPoly2 () : d(-2), degx(-2), degy(-2) {}
  PPoly2 (int nt) : a(nt), m(2*nt), d(-2), degx(-2), degy(-2) {}
  PPoly2 (int nt, Parameter *a, int *m) : a(nt), m(2*nt), d(-2), degx(-2), degy(-2) {
    int k = nt*2;
    for (int i = 0; i < nt; ++i)
      this->a[i] = a[i];
    for (int i = 0; i < k; ++i)
      this->m[i] = m[i];
    setDegree();
    degx = degX();
    degy = degY();
  }

  void setDegree () {
    d = 0;
    for (int i = 0; i < a.size(); ++i) {
      int di = 0;
      for (int j = 0; j < 2; ++j)
        di += m[2*i+j];
      if (d < di)
        d = di;
    }
    degx = degX();
    degy = degY();
  }

  int degree () const { return d; }
  bool increased () const { assert(0); return false; /* a[0].increased(); */ }
  bool decreased () const { assert(0); return true; /* a[0].decreased(); */ }
  Parameter value (Parameter *x) const {
    Parameter xp[degx+1];
    Parameter yp[degy+1];
    
    Parameter p = x[0];
    for(int i = 1; i <= degx; i++) {
      xp[i] = p;
      p = p*x[0];
    }
    
    p = x[1];
    for(int i = 1; i <= degy; i++) {
      yp[i] = p;
      p = p*x[1];
    }
    
    Parameter y;
    auto mp = m.begin();
    for(int i = 0; i < a.size(); ++i) {
      Parameter z = a[i];
      if(*mp) z = z*xp[*mp];
      mp++;
      if(*mp) z = z*yp[*mp];
      mp++;
      if(i == 0)
        y = z;
      else
        y = y + z;
    }
    
    return y;
  }
  
  Parameter value (const std::initializer_list<Parameter>& x) const {
    Parameter xy[2] = {*(x.begin()), *(x.begin()+1)};
    return value(xy);
  }

  Parameter value (const PV2 p) {
    Parameter xy[2] = {p.x, p.y};
    return value(xy);
  }

  static PPoly2 one() {
    int nt = 1;
    Parameter a[1];
    a[0] = Parameter::constant(1);
    int m[] = {0, 0};
    return PPoly2(nt, a, m);
  }

  static PPoly2 zero() {
    int nt = 1;
    Parameter a[1];
    a[0] = Parameter::constant(0);
    int m[] = {0, 0};
    return PPoly2(nt, a, m);
  }

  PPoly1 subX(Parameter x) const  {
    vector<Parameter> coef(degy + 1);
    
    for(int i = 0; i < coef.size(); i++) {
      coef[i] = Parameter::constant(0);
    }
    
    vector<Parameter> powx;
    powx.push_back(Parameter::constant(1));
    powx.push_back(x);
    for(int i = 2; i <= degx; i++) {
      powx.push_back(powx[powx.size()-1] * x);
    }
    
    for(int i = 0; i < a.size(); i++) {

      int dx = m[2*i];
      int dy = m[2*i+1];
      coef[dy] = coef[dy] + (a[i] * powx[dx]);
      
    }
    
    return PPoly1(degx, &coef[0]);
  }

  PPoly1 subY(Parameter y) const {
    vector<Parameter> coef(degx + 1);
    
    for(int i = 0; i < coef.size(); i++) {
      coef[i] = Parameter::constant(0);
    }
    
    vector<Parameter> powy;
    powy.push_back(Parameter::constant(1));
    powy.push_back(y);
    for(int i = 2; i <= degy; i++) {
      powy.push_back(powy[powy.size()-1] * y);
    }
    
    for(int i = 0; i < a.size(); i++) {
      
      int dx = m[2*i];
      int dy = m[2*i+1];
      coef[dx] = coef[dx] + (a[i] * powy[dy]);
      
    }

    return PPoly1(degy, &coef[0]);
  }

  PPoly2 derX () const  {
    vector<Parameter> coef;
    vector<int> pow;
    
    for(int i = 0; i < a.size(); i++) {
      if(m[2*i] < 1) continue;
      coef.push_back(a[i] * m[2*i]);
      pow.push_back(m[2*i] - 1);
      pow.push_back(m[2*i + 1]);
    }
    
    return PPoly2(coef.size(), &coef[0], &pow[0]);
  }

  PPoly2 derY () const  {
    vector<Parameter> coef;
    vector<int> pow;
    
    for(int i = 0; i < a.size(); i++) {
      if(m[2*i + 1] < 1) continue;
      coef.push_back(a[i] * m[2*i + 1]);
      pow.push_back(m[2*i]);
      pow.push_back(m[2*i + 1] - 1);
    }
    
    return PPoly2(coef.size(), &coef[0], &pow[0]);
    
  }
  
  int degX() const  {
    int d = 0;
    
    for(int i = 0; i < a.size(); i++) {
      int x = m[2*i];
      d = std::max(d, x);
    }
    
    return d;
  }

  int degY() const {
    int d = 0;
    
    for(int i = 0; i < a.size(); i++) {
      int x = m[2*i + 1];
      d = std::max(d, x);
    }
    
    return d;
  }

  void print_input() const {
    for(int i = 0; i < a.size(); i++) {
      printf("%.16f %d %d%s", a[i].mid(), m[2*i], m[2*i+1], (i == a.size()-1 ? "\n" : " " ));
    }
  }

  void print() const {
    for(int i = 0; i < a.size(); i++) {
      printf("%.16f*x^(%d)*y^(%d) %s", a[i].mid(), m[2*i], m[2*i+1], (i == a.size()-1 ? "\n" : "+ "));
    }
  }

  PV2 polish (const PV2 &r) const {
    assert(0);
    return r;
  }

  PPoly2 operator* (const PPoly2 &b) const {
    int cnt = 0;
    // Parameter ca[nt*b.nt];
    vector<Parameter> ca(a.size()*b.a.size());
    int cm[2*a.size()*b.a.size()];
    int exps[2];
    for (int i = 0; i < a.size(); i++) {
      auto mp = m.begin() + i * 2;
      for (int j = 0; j < b.a.size(); j++) {
        Parameter coef = a[i] * b.a[j];
        auto bmp = b.m.begin() + j * 2;
        for (int k = 0; k < 2; k++)
          exps[k] = mp[k] + bmp[k];
        int l;
        for (l = 0; l < cnt; l++) {
          int *cmp = cm + l * 2;
          int h;
          for (h = 0; h < 2; h++)
            if (cmp[h] != exps[h])
              break;
          if (h == 2) {
            ca[l] = ca[l] + coef;
            break;
          }
        }
        if (l == cnt) {
          int *cmp = cm + l * 2;
          for (int h = 0; h < 2; h++)
            cmp[h] = exps[h];
          ca[l] = coef;
          cnt++;
        }
      }
    }
    // return PPoly2(cnt, ca, cm);
    return PPoly2(cnt, &ca[0], cm);
  }

  PPoly2 operator* (const Parameter &p) const {
    vector<Parameter> ca(a.size());
    int cm[2*a.size()];
    
    for(int i = 0; i < a.size(); i++) {
      ca[i] = a[i] * p;
      cm[2*i] = m[2*i];
      cm[2*i + 1] = m[2*i + 1];
    }
    
    return PPoly2(a.size(), &ca[0], cm);
  }

  PPoly2 operator* (const double p) const {
    vector<Parameter> ca(a.size());
    int cm[2*a.size()];
    
    for(int i = 0; i < a.size(); i++) {
      ca[i] = a[i] * p;
      cm[2*i] = m[2*i];
      cm[2*i + 1] = m[2*i + 1];
    }
    
    return PPoly2(a.size(), &ca[0], cm);
  }

  PPoly2 plusCtimes (double c, const PPoly2 &b) const  {
    int pnt = a.size();
    vector<Parameter> pa(a.size() + b.a.size());
    int pm[2 * (a.size() + b.a.size())];
    
    for (int i = 0; i < a.size(); i++)
      pa[i] = a[i];
    int ntnv = a.size() * 2;
    for (int i = 0; i < ntnv; i++)
      pm[i] = m[i];
    
    auto bmp = b.m.begin();
    for (int i = 0; i < b.a.size(); i++) {
      int j;
      for (j = 0; j < a.size(); j++) {
        auto mp = m.begin() + j * 2;
        int k;
        for (k = 0; k < 2; k++)
          if (mp[k] != bmp[k])
            break;
        if (k == 2) {
          pa[j] = pa[j] + c * b.a[i];
          break;
        }
      }
      if (j == a.size()) {
        pa[pnt] = c * b.a[i];
        int *pmp = pm + pnt * 2;
        for (int k = 0; k < 2; k++)
          pmp[k] = bmp[k];
        pnt++;
      }
      bmp += 2;
    }
    
    return PPoly2(pnt, &pa[0], pm);
  }

  PPoly2 operator+ (const PPoly2 &b) const { return plusCtimes(1.0, b); }
  PPoly2 operator- (const PPoly2 &b) const { return plusCtimes(-1.0, b); }

  int size () const { return a.size(); }
  const Parameter &operator[] (int i) const { return a[i]; }
  Parameter &operator[] (int i) { return a[i]; }

  int d;
  vector<int> m;
  vector<Parameter> a;
  int degx, degy;
};

class Index3 {
public:
  int i[3];
  Index3 (int i0, int i1, int i2) { i[0] = i0; i[1] = i1; i[2] = i2; }

  bool operator< (const Index3 &b) const {
    const Index3 &a = *this;
    for (int j = 3; --j >= 0;)
      if (a.i[j] != b.i[j])
        return a.i[j] < b.i[j];
    return false;
  }
    
  class Compare {
  public:
    bool operator() (const Index3 &a, const Index3 &b) {
      for (int j = 2; --j >= 0;)
        if (a.i[j] != b.i[j])
          return a.i[j] < b.i[j];
      return false;
    }
  };

  class Hasher {
  public:
    size_t operator() (const Index3 &a) const {
      return a.i[0] * 779230947 + a.i[1] * 247091631 + a.i[2] * 1194289623;
    }
  };
  
  class Equals {
  public:
    bool operator() (const Index3 &a, const Index3 &b) const {
      return a.i[0] == b.i[0] && a.i[1] == b.i[1] && a.i[2] == b.i[2];
    }
  };
};
 
class PPoly3 {
public:
  PPoly3 () {}

  void add (int i, int j, int k, const Parameter& c) { 
    auto it = ia.find(Index3(i, j, k));
    if (it == ia.end()) {
      ia[Index3(i, j, k)] = a.size();
      a.push_back(c);
    }
    else
      a[it->second] = a[it->second] + c;
  }

  Parameter value (const Parameter& x, const Parameter& y, const Parameter& z) const {
    const Parameter *xyz[3] = { &x, &y, &z };
    
    vector<Parameter> powers[3];
    for (int i = 0; i < 3; i++)
      powers[i].push_back(Parameter::constant(1));
    
    Parameter sum = Parameter::constant(0);
    
    for (auto pair = ia.begin(); pair != ia.end(); pair++) {
      Index3 ind = pair->first;
      Parameter term = a[pair->second];
      
      for (int i = 0; i < 3; i++) {
        while (powers[i].size() <= ind.i[i])
          powers[i].push_back(*xyz[i] * powers[i].back());
        term = term * powers[i][ind.i[i]];
      }
      
      sum = sum + term;
    }
    
    return sum;
  }

  // xyz = 0, 1, 2:  substituting for yz, xz, xy.
  PPoly substitute2
    (int ixyz, const Parameter& x, const Parameter& y, const Parameter& z) const {
    const Parameter *xyz[3] = { &x, &y, &z };
    
    vector<Parameter> powers[3];
    for (int i = 0; i < 3; i++)
      powers[i].push_back(Parameter::constant(1));
    
    PPoly ppoly;
    
    for (auto pair = ia.begin(); pair != ia.end(); pair++) {
      Index3 ind = pair->first;
      Parameter term = a[pair->second];
      
      for (int i = 0; i < 3; i++) 
	if (i != ixyz) {
	  while (powers[i].size() <= ind.i[i])
	    powers[i].push_back(*xyz[i] * powers[i].back());
	  term = term * powers[i][ind.i[i]];
	}
      
      ppoly.add(term, ind.i[ixyz]);
    }
    
    return ppoly;
  }
  
  // unordered_map<Index3, int, Index3::Hasher, Index3::Equals> ia;
  map<Index3, int> ia;
  vector<Parameter> a;
  int size () const { return a.size(); }
  const Parameter &operator[] (int i) const { return a[i]; }
  Parameter &operator[] (int i) { return a[i]; }
};

}

#endif
