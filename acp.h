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

#ifndef ACP
#define ACP

// #define NO_MODE

#include <gmp.h>
#include <mpfr.h>
#include <assert.h>
#include <iostream>
#include <fenv.h>
#include <exception>
#include <math.h>
#include <algorithm>
#include <float.h>
#include <vector>

using namespace std;

namespace acp {

double randomNumber (double rmin, double rmax);

inline double nextD (double x) {
  static double p = 1.0/(1 << 26)/(1 << 26);
  double x2 = x + p*fabs(x);
  if (x2 == x)
    return nextafter(x, 1.0);
  double x3 = 0.5*(x + x2);
  return x3 > x ? x3 : x2;
}
 
inline double prevD (double x) {
  static double p = 1.0/(1 << 26)/(1 << 26);
  double x2 = x - p*fabs(x);
  if (x2 == x)
    return nextafter(x, -1.0);
 double x3 = 0.5*(x + x2);
 return x3 < x ? x3 : x2;
}

class SignException : public std::exception {
 public:
  virtual const char* what() const throw() {
    return "Not enough precision";
  }
};

extern SignException signException;
 
enum RoundMode { RoundUp=1, RoundDown=-1, RoundNearest=0 };

class MValue {
 public:
  MValue (unsigned int ip) : p(ip) { mpfr_init2(m, p); }
  MValue (double x, unsigned int ip) : p(ip) { 
    mpfr_init2(m, p);
    mpfr_set_d(m, x, GMP_RNDN);
  }
  MValue (const MValue &v, unsigned int ip, mpfr_rnd_t round = GMP_RNDN) : p(ip) { 
    mpfr_init2(m, ip); 
    mpfr_set(m, v.m, round); 
  }
  MValue (const MValue &v) : p(v.p) { 
    mpfr_init2(m, v.p); 
    mpfr_set(m, v.m, GMP_RNDN); 
  }
  ~MValue () { mpfr_clear(m); }
  void operator= (const MValue &v) { mpfr_set(m, v.m, GMP_RNDN); }
  double value () const { return mpfr_get_d(m, GMP_RNDN); }
  MValue plus (const MValue &b, mpfr_rnd_t round) const {
    MValue res(max(p, b.p));
    mpfr_add(res.m, m, b.m, round);
    return res;
  }
  MValue plus (double b, mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_add_d(res.m, m, b, round);
    return res;
  }
  MValue minus () const {
    MValue res(p);
    mpfr_neg(res.m, m, GMP_RNDN);
    return res;
  }
  MValue minus (const MValue &b, mpfr_rnd_t round) const {
    MValue res(max(p, b.p));
    mpfr_sub(res.m, m, b.m, round);
    return res;
  }
  MValue times (const MValue &b, mpfr_rnd_t round) const {
    MValue res(max(p, b.p));
    mpfr_mul(res.m, m, b.m, round);
    return res;
  }
  MValue times (double b, mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_mul_d(res.m, m, b, round);
    return res;
  }
  MValue divide (const MValue &b, mpfr_rnd_t round) const {
    MValue res(max(p, b.p));
    mpfr_div(res.m, m, b.m, round);
    return res;
  }
  int sign () const { return mpfr_sgn(m); }
  bool operator< (const MValue &b)  const { return mpfr_less_p(m, b.m); }
  MValue sqrt (mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_sqrt(res.m, m, round);
    return res;
  }
  MValue root (unsigned long int k, mpfr_rnd_t round) const {
    MValue res(p);
    //mpfr_rootn_ui(res.m, m, k, round);
    mpfr_root(res.m, m, k, round);
    return res;
  }

  mpfr_t m;
  unsigned int p;
};

typedef unsigned int uint;

typedef unsigned long ulong;

unsigned int inverse (unsigned int a, unsigned int n);

#define ulong(a) ((unsigned long)(a))
 
class ModP {
  static const int eShift;
  static const int eMax;
  static const int eMin;

  vector<unsigned int> pow2v;
  unsigned int *pow2;
  unsigned int p;

public:
  ModP (unsigned int p) : p(p), pow2v(eMax - eMin + 1) {
    pow2 = &pow2v[0] - eMin;
    pow2[0] = 1;
    unsigned long x = 1;
    for (int e = 1; e <= eMax; e++) {
      x = (2 * x) % p;
      pow2[e] = x;
    }
    x = 1;
    unsigned int i2 = (p + 1) / 2;
    for (int e = -1; e >= eMin; e--) {
      x = (i2 * x) % p;
      pow2[e] = x;
    }
  }

  unsigned int mod (double x) {
    if (x == (long) x)
      return x >= 0 ? ((long) x) % p : p + ((long) x) % p;
    int e;
    double m = frexp(x, &e) * (1l << eShift);
    assert(m == (long) m);
    e -= eShift;
    long mp = ((long) m) % p;
    if (mp < 0)
      mp += p;
    assert(eMin <= e && e <= eMax);
    return ((unsigned long) mp * pow2[e]) % p;
  }
};

class ModInt {
 public:
  static unsigned int ps[7];
  static int pIndex;
  static unsigned int p1, p2;
  static ModP modP1, modP2;

  static void changePrime (unsigned int p);

  ModInt () {}
  
  ModInt (double r, unsigned int ip) : p(ip) {
    if (p == p1)
      a = modP1.mod(r);
    else if (p == p2)
      a = modP2.mod(r);
    else
      assert(0);
  }

  ModInt operator- () const {
    long l = (-long(a))%p;
    return ModInt(l < 0 ? l + p : l, p);
  }
  
  ModInt operator+ (const ModInt &x) const {

    return ModInt((long(a) + long(x.a))%p, p);
  }

  ModInt operator+ (double x) const {
    return *this + ModInt(x, p);
  }

  ModInt operator- (const ModInt &x) const {
    long l = (long(a) - long(x.a))%p;
    return ModInt(l < 0 ? l + p : l, p);
  }

  ModInt operator- (double x) const {
    return *this - ModInt(x, p);
  }

  ModInt operator* (const ModInt &x) const {
    return ModInt((ulong(a)*ulong(x.a))%p, p);
  }

  ModInt operator* (double x) const {
    return *this*ModInt(x, p);
  }
  
  ModInt operator/ (const ModInt &x) const {
    if (x.a == 0)
      throw p;
    return ModInt((ulong(a)*ulong(inverse(x.a, p)))%p, p);
  }

  unsigned int a, p;
};

extern pthread_key_t mikey;
extern pthread_key_t cpkey;

class Parameter;

class MInt {
  static unsigned int & curPrecision () {
    void *v = pthread_getspecific(cpkey);
    return v ? *(unsigned int *) v : curPrecisions[0];
  }

public:
  static unsigned int curPrecisions[128];

  static unsigned int getPrecision () {
    return curPrecision();
  }

  static void setPrecision (unsigned int p) {
    curPrecision() = p;
    if (p == 53u)
      clearMInts();
  }

  static vector<MInt*> mints[128];

  static vector<MInt*>& getMInts () {
    void *v = pthread_getspecific(mikey);
    return v ? *(vector<MInt*>*) v : mints[0];
  }

  const ModInt m1, m2;

  MInt ()
    : m1(ModInt(random(), ModInt::p1)), m2(ModInt(random(), ModInt::p2)) {
    getMInts().push_back(this);
  }
  MInt (const ModInt &m1, const ModInt &m2) : m1(m1), m2(m2) {
    getMInts().push_back(this);
  }
  virtual ~MInt () {}

  static void clearMInts () {
    vector<MInt*> &mints = getMInts();
    for (int i = 0; i < mints.size(); i++)
      delete mints[i];
    mints.clear();
  }

  bool zeroMod () { return m1.a == 0u && m2.a == 0u; }

  static MInt* make (unsigned int precision, const Parameter &p);
  static MInt* make (unsigned int precision, const MInt *m);
  MInt *clone () { 
    MInt *m = make(getPrecision(), this);
    getMInts().pop_back();
    return m;
  }
  virtual unsigned int precision () const = 0;
  virtual double intervalWidth () const = 0;
  virtual double lb () const = 0;
  virtual double ub () const = 0;
  virtual MInt* lbP () const = 0;
  virtual MInt* ubP () const = 0;
  virtual MInt* plus (const MInt* b) const = 0;
  // virtual MInt* plus (double b) const = 0;
  virtual MInt* minus (const MInt* b) const = 0;
  virtual MInt* minus () const = 0;
  virtual MInt* times (const MInt* b) const = 0;
  // virtual MInt* times (double b) const = 0;
  virtual MInt* divide (const MInt* b) const = 0;
  virtual int sign (bool fail = true) const = 0;
  virtual MInt* mid () const = 0;
  virtual bool subset (const MInt* b) const = 0;
  virtual MInt* interval (const MInt* b) const = 0;
  virtual MInt* innerInterval (const MInt* b) const = 0;
  virtual bool intersects (const MInt* b) const = 0;
  virtual MInt* intersect (const MInt* b) const = 0;
  virtual MInt* sqrt () const = 0;
  virtual MInt* root (unsigned long int k) const = 0;
  Parameter par () const;
};

class Parameter {
  friend class MInt;
  friend class PInt;
  friend class EInt;
  template<class P> friend class Object;
  friend class Primitive;
  static bool penabled, enabled;
  static double no_optimize (volatile double x) { return x; }
  
  static const double sentinel;

  static Parameter sqrt (double);
  static Parameter root (double a, unsigned int n);
  static Parameter pow (const Parameter &x, unsigned long int n);

  double l;
  union {
    double r;
    MInt *m;
  } u;

  Parameter (double l, double u) : l(l) { this->u.r = u; }

  static Parameter interval (double x, double y) { 
    Parameter p(x, y);
    if (MInt::getPrecision() > 53u)
      p.increasePrecision();
    return p;
  }
  
  Parameter (double x) : l(x) { u.r = x; }

public:
  Parameter (const Parameter &il, const Parameter &iu) {
#ifdef USE_ASSERT
    assert(il.low() == iu.low());
#endif
    if (il.low()) {
      l = il.l;
      u.r = iu.u.r;
    }
    else {
      l = sentinel;
      u.m = il.u.m->interval(iu.u.m);
    }
  }

private:
  Parameter (MInt *m) : l(sentinel) { u.m = m; }

  void increasePrecision () {
#ifdef USE_ASSERT
    assert(low());
#endif
    u.m = MInt::make(MInt::getPrecision(), *this);
    l = sentinel;
  }

  void decreasePrecision () {
#ifdef USE_ASSERT
    assert(high());
#endif
    l = u.m->lb();
    u.r = u.m->ub();
  }

 public:
  static double delta; // default is 2^{-26}

  static bool isEnabled () { return enabled; }

  static void enable () {
    extern pthread_key_t mikey;
    extern pthread_key_t cpkey;
    extern pthread_key_t hsekey;
    if (!penabled) {
      penabled = true;
      pthread_key_create(&mikey, NULL);
      pthread_key_create(&cpkey, NULL);
      pthread_key_create(&hsekey, NULL);
    }
#ifndef NO_MODE
    if (!enabled) { enabled = true; fesetround(FE_UPWARD); }
#endif
  }
  
  static void disable () {
#ifndef NO_MODE
    if (enabled) { enabled = false; fesetround(FE_TONEAREST); }
#endif
  }

  static void exit () {
    extern pthread_key_t mikey;
    extern pthread_key_t cpkey;
    extern pthread_key_t hsekey;
    pthread_key_delete(mikey);
    pthread_key_delete(cpkey);
    pthread_key_delete(hsekey);
    disable();
  }
  
  int size () const { return 1;}
  const Parameter &operator[] (int i) const {return *this; }
  Parameter &operator[] (int i) {return *this; }

  Parameter () : l(0.0) { u.r = -1.0; }
  bool uninitialized () const { return l == 0.0 && u.r == -1.0; }
  bool initialized () const { return l <= u.r || high(); }
  
  static Parameter constant (double x) {
    Parameter p(x);
    if (MInt::getPrecision() > 53u)
      p.increasePrecision();
    return p;
  }
  
  static Parameter input (double x) { 
    double y = x + delta*(1.0 + fabs(x))*randomNumber(-1.0, 1.0);
    return Parameter::constant(y);
  }
  
  bool low () const { return l != sentinel; }
  bool high () const { return l == sentinel; }

  bool hasTheRightPrimes () {
    // assert(u.m->m1.p == ModInt::p1 || u.m->m2.p == ModInt::p2);
    return u.m->m1.p == ModInt::p1 && u.m->m2.p == ModInt::p2;
  }

  unsigned int precision () const { 
    return uninitialized() ? 0 : low() ? 53u : u.m->precision(); 
  }
  
  bool zero () { return lb() == 0 && ub() == 0; }

  bool zeroMod () {
#ifdef USE_ASSERT
    assert(high());
#endif
    return u.m->zeroMod();
  }

  int sign (bool fail = true) const {
    if (high())
      return u.m->sign(fail);
    if (l > 0.0)
      return 1;
    if (u.r < 0.0)
      return -1;
    if (!fail || (l == 0.0 && u.r == 0.0))
      return 0;
    throw signException;
  }

  double mid () const { return 0.5 * (lb() + ub()); }
  double lb () const { return low() ? l : u.m->lb(); }
  double ub () const { return low() ? u.r : u.m->ub(); }
  double intervalWidth () const { return low() ? u.r - l : u.m->intervalWidth(); }
  Parameter lbP () const {
    return low() ? Parameter(l, l) : Parameter(u.m->lbP());
  }
  Parameter ubP () const {
    return low() ? Parameter(u.r, u.r) : Parameter(u.m->ubP());
  }
  Parameter midP () const { 
    return low() ? Parameter(0.5*(l + u.r)) : Parameter(u.m->mid()); 
  }

  bool subset (const Parameter &b) const {
#ifdef USE_ASSERT
    if (low() && b.low())
#else
    if (low())
#endif
      return !(l < b.l) && !(b.u.r < u.r) && (b.l < l || u.r < b.u.r);
#ifdef USE_ASSERT      
    assert(high() && b.high());
#endif
    return u.m->subset(b.u.m);
  }

  Parameter interval (const Parameter &b) const { 
#ifdef USE_ASSERT
    if (low() && b.low()) {
#else
    if (low()) {
#endif
      assert(l <= b.u.r);
      return Parameter(l, b.u.r);
    }
#ifdef USE_ASSERT      
    assert(high() && b.high());
#endif
    return Parameter(u.m->interval(b.u.m));
  }

  Parameter innerInterval (const Parameter &b) const { 
#ifdef USE_ASSERT
    if (low() && b.low()) {
#else
    if (low()) {
#endif
      assert(u.r <= b.l);
      return Parameter(u.r, b.l);
    }
#ifdef USE_ASSERT      
    assert(high() && b.high());
#endif
    return Parameter(u.m->innerInterval(b.u.m));
  }

  bool intersects (const Parameter &b) const { 
#ifdef USE_ASSERT
    if (low() && b.low())
#else
    if (low())
#endif
      return !(u.r < b.l || b.u.r < l); 
#ifdef USE_ASSERT      
    assert(high() && b.high());
#endif
    return u.m->intersects(b.u.m);
  }
  
  Parameter intersect (const Parameter &b) const {
#ifdef USE_ASSERT
    if (low() && b.low()) {
#else
    if (low()) {
#endif
      assert(!(u.r < b.l || b.u.r < l));
      double il = l < b.l ? b.l : l;
      double iu = u.r < b.u.r ? u.r : b.u.r;
      return Parameter(il, iu);
    }
#ifdef USE_ASSERT      
    assert(high() && b.high());
#endif
    return Parameter(u.m->intersect(b.u.m));
  }

  Parameter operator+ (const Parameter &b) const {
#ifdef USE_ASSERT
    if (low() && b.low())
#else
    if (low())
#endif
#ifndef NO_MODE
      return Parameter(- no_optimize(no_optimize(- l) - b.l), u.r + b.u.r);
#else
      return Parameter(prevD(l + b.l), nextD(u.r + b.u.r));
#endif    
#ifdef USE_ASSERT      
    assert(high() && b.high());
#endif
    return Parameter(u.m->plus(b.u.m));
  }

  Parameter & operator+= (const Parameter &p) {
    return *this = *this + p;
  }

  Parameter operator+ (double b) const {
    return *this + Parameter::constant(b);
  }

  Parameter operator- (const Parameter &b) const {
#ifdef USE_ASSERT
    if (low() && b.low())
#else
    if (low())
#endif
#ifndef NO_MODE
      return Parameter(- no_optimize(b.u.r - l), u.r - b.l);
#else
      return Parameter(prevD(l - b.u.r), nextD(u.r - b.l));
#endif    
#ifdef USE_ASSERT      
    assert(high() && b.high());
#endif
    return Parameter(u.m->minus(b.u.m));
  }
  
  Parameter operator- (double b) const {
    return *this - Parameter::constant(b);
  }

  Parameter operator- () const {
    if (low())
      return Parameter(- u.r, - l);
    return Parameter(u.m->minus());
  }
  
  Parameter operator* (const Parameter &b) const {
#ifndef NO_MODE
#ifdef USE_ASSERT
    if (low() && b.low()) {
#else
    if (low()) {
#endif
      Parameter s = u.r < 0.0 ? - *this : *this, t = u.r < 0.0 ? - b : b;
      if (s.l > 0.0) {
	volatile double k = t.l > 0.0 ? - s.l : - s.u.r;
	return Parameter(- no_optimize(k*t.l), t.u.r > 0.0 ? s.u.r*t.u.r : s.l*t.u.r);
      }
      if (t.l > 0.0) {
	volatile double k = - s.l;
	return Parameter(- no_optimize(k*t.u.r), s.u.r*t.u.r);
      }
      if (t.u.r < 0.0) {
	volatile double k = - s.u.r;
	return Parameter(- no_optimize(k*t.l), s.l*t.l);
      }
      volatile double k1 = - s.l, k2 = - s.u.r;
      double cl1 = no_optimize(k1*t.u.r), cl2 = no_optimize(k2*t.l), 
	cu1 = s.l*t.l, cu2 = s.u.r*t.u.r,
	cl = cl1 < cl2 ? - cl2 : - cl1, 
        cu = cu1 < cu2 ? cu2 : cu1;
      return Parameter(cl, cu);
    }
#else
#ifdef USE_ASSERT
    if (low() && b.low()) {
#else
    if (low()) {
#endif
      Parameter s = u.r < 0.0 ? - *this : *this, t = u.r < 0.0 ? - b : b;
      if (s.l >= 0.0)
	if (t.l >= 0.0)
	  return Parameter(prevD(s.l*t.l), nextD(s.u.r*t.u.r));
	else if (t.u.r <= 0.0)
	  return Parameter(prevD(s.u.r*t.l), nextD(s.l*t.u.r));
	else
	  return Parameter(prevD(s.u.r*t.l), nextD(s.u.r*t.u.r));
      if (t.l >= 0.0)
	return Parameter(prevD(s.l*t.u.r), nextD(s.u.r*t.u.r));
      if (t.u.r <= 0.0)
	return Parameter(prevD(s.u.r*t.l), nextD(s.l*t.l));
      double k1 = s.l*t.u.r, k2 = s.u.r*t.l, nl = k1 < k2 ? k1 : k2,
	k3 = s.l*t.l, k4 = s.u.r*t.u.r, nu = k3 < k4 ? k4 : k3;
      return Parameter(prevD(nl), nextD(nu));
    }
#endif
#ifdef USE_ASSERT      
    assert(high() && b.high());
#endif
    return Parameter(u.m->times(b.u.m));
  }
  
  Parameter operator* (double b) const {
    return *this * Parameter::constant(b);
  }

  Parameter operator/ (const Parameter &b) const {
    int bs = b.sign();
    assert(bs != 0);
#ifndef NO_MODE
#ifdef USE_ASSERT
    if (low() && b.low()) {
#else
    if (low()) {
#endif
      if (bs == 1) {
	if (l >= 0.0)
	  return Parameter(- no_optimize(no_optimize(- l)/b.u.r), u.r/b.l);
	if (u.r <= 0.0)
	  return Parameter(- no_optimize(no_optimize(- l)/b.l), u.r/b.u.r);
	return Parameter(- no_optimize(no_optimize(- l)/b.l), u.r/b.l);
      }
      if (l >= 0.0)
	return Parameter(- no_optimize(no_optimize(- u.r)/b.u.r), l/b.l);
      if (u.r <= 0.0)
	return Parameter(- no_optimize(no_optimize(- u.r)/b.l), l/b.u.r);
      return Parameter(- no_optimize(no_optimize(- u.r)/b.u.r), l/b.u.r);
    }
#else
#ifdef USE_ASSERT
    if (low() && b.low()) {
#else
    if (low()) {
#endif
      if (bs == 1)
	if (l >= 0.0)
	  return Parameter(prevD(l/b.u.r), nextD(u.r/b.l));
	else if (u.r <= 0.0)
	  return Parameter(prevD(l/b.l), nextD(u.r/b.u.r));
	else
	  return Parameter(prevD(l/b.l), nextD(u.r/b.l));
      if (l >= 0.0)
	return Parameter(prevD(u.r/b.u.r), nextD(l/b.l));
      if (u.r <= 0.0)
	return Parameter(prevD(u.r/b.l), nextD(l/b.u.r));
      return Parameter(prevD(u.r/b.u.r), nextD(l/b.u.r));
    }
#endif
#ifdef USE_ASSERT      
    assert(high() && b.high());
#endif
    return Parameter(u.m->divide(b.u.m));
  }

  Parameter operator/ (double b) const {
    return *this / Parameter::constant(b);
  }
    
  bool operator< (const Parameter &b) const { return (b - *this).sign() == 1; }
  bool operator< (double b) const { return (*this - b).sign() == -1; }
  bool operator> (const Parameter &b) const { return (b - *this).sign() == -1; }
  bool operator> (double b) const { return (*this - b).sign() == 1; }

  static Parameter max (const Parameter &a, const Parameter &b) { return a < b ? b : a; }

  Parameter abs () const { 
    return sign() == 1 ? *this : -*this; 
  }

  Parameter sqrt () const;

  Parameter root (unsigned long int n) const;

  Parameter pow (long int n) const {
    if (n == 0)
      return Parameter(1.0);
    if (n < 0)
      return Parameter(1.0)/Parameter::pow(*this, -n);
    return Parameter::pow(*this, n);
  }
};

inline Parameter MInt::par () const { return Parameter(lb(), ub()); }

inline Parameter operator+ (double a, const Parameter &b)
{
  return b + a;
}

inline Parameter operator- (double a, const Parameter &b)
{
  return (- b) + a;
}

inline Parameter operator* (double a, const Parameter &b)
{
  return b*a;
}

inline Parameter operator/ (double a, const Parameter &b)
{
  return Parameter::constant(a) / b;
}

inline bool operator< (double a, const Parameter &b)
{
  return (b - a).sign() == 1;
}

inline bool operator> (double a, const Parameter &b)
{
  return (b - a).sign() == -1;
}

class PInt : public MInt {
  Parameter p;
  friend class Parameter;
 public:
  PInt (const Parameter &p) { this->p = p; }
  
  PInt (const Parameter &p, const ModInt &m1, const ModInt &m2) 
    : MInt(m1, m2) { this->p = p; }

  unsigned int precision () const { return 106u; }

  double intervalWidth () const { return p.intervalWidth(); }
  
  double lb () const { return p.lb(); }
  
  double ub () const { return p.ub(); }
  
  MInt* lbP () const { return new PInt(p.lbP()); }

  MInt* ubP () const { return new PInt(p.ubP()); }

  MInt* plus (const MInt* that) const {
    return new PInt(p + that->par(), m1 + that->m1, m2 + that->m2);
  }
  
  MInt* plus (double b) const { return new PInt(p + Parameter(b), m1 + b, m2 + b); }
  
  MInt* minus (const MInt* that) const {
    return new PInt(p - that->par(), m1 - that->m1, m2 - that->m2);
  }
  
  MInt* minus () const { return new PInt(- p, - m1, - m2); }
  
  MInt* times (const MInt* that) const {
    return new PInt(p * that->par(), m1*that->m1, m2*that->m2);
  }
  
  MInt* times (double b) const { return new PInt(p * Parameter(b), m1*b, m2*b); }

  MInt* divide (const MInt* that) const {
    return new PInt(p / that->par(), m1/that->m1, m2/that->m2); 
  }

  int sign (bool fail = true) const {
    if (p.u.r < 0.0)
      return -1;
    if (p.l > 0.0)
      return 1;
    if (!fail || (p.l == 0.0 && p.u.r == 0.0) || (m1.a == 0u && m2.a == 0u))
      return 0;
    throw signException;
  }

  MInt* mid () const { return new PInt(p.midP()); }

  bool subset (const MInt* that) const {
    return p.subset(that->par());
  }

  MInt* interval (const MInt* that) const {
    return new PInt(p.interval(that->par()));
  }

  MInt* innerInterval (const MInt* that) const {
    return new PInt(p.innerInterval(that->par()));
  }

  bool intersects (const MInt* that) const {
    return p.intersects(that->par());
  }

  MInt* intersect (const MInt* that) const {
    return new PInt(p.intersect(that->par()));
  }

  MInt* sqrt () const { return new PInt(p.sqrt()); }

  MInt* root (unsigned long int k) const { assert(0); return new PInt(p.root(k)); }

};

 class EInt : public MInt {
   int prec;
 public:
   EInt (const MValue &l, const MValue &u, int prec) : lm(l), um(u), prec(prec) {}
   EInt (const MValue &l, const MValue &u, const ModInt &m1, const ModInt &m2, int prec)
     : MInt(m1, m2), lm(l), um(u), prec(prec) {}
  ~EInt () {}

  PInt *pInt () const { return new PInt(Parameter(lb(), ub()), m1, m2); }
 
  unsigned int precision () const { return lm.p; }

  double intervalWidth () const { return um.minus(lm, GMP_RNDN).value(); }

  double lb () const { return mpfr_get_d(lm.m, GMP_RNDD); }

  double ub () const { return mpfr_get_d(um.m, GMP_RNDU); }

  MInt* lbP () const { return new EInt(lm, lm, prec); }

  MInt* ubP () const { return new EInt(um, um, prec); }

  MInt* plus (const MInt* that) const {
    const EInt *b = dynamic_cast<const EInt*>(that);
    if (!b) return pInt()->plus(that);
    return new EInt(lm.plus(b->lm, GMP_RNDD), um.plus(b->um, GMP_RNDU),
		    m1 + b->m1, m2 + b->m2, min(prec, b->prec));
  }

  MInt* plus (double b) const {
    return new EInt(lm.plus(b, GMP_RNDD), um.plus(b, GMP_RNDU),
		    m1 + b, m2 + b, prec);
  }

  MInt* minus (const MInt* that) const {
    const EInt *b = dynamic_cast<const EInt*>(that);
    if (!b) return pInt()->minus(that);
    return new EInt(lm.minus(b->um, GMP_RNDD), um.minus(b->lm, GMP_RNDU),
		    m1 - b->m1, m2 - b->m2, min(prec, b->prec));
  }

  MInt* minus () const { 
    return new EInt(um.minus(), lm.minus(), - m1, - m2, prec);
  }

  MInt* times (const MInt* that) const;

  MInt* times (double b) const {
    return b > 0.0 ? new EInt(lm.times(b, GMP_RNDD), um.times(b, GMP_RNDU),
			      m1*b, m2*b, prec)
      : new EInt(um.times(b, GMP_RNDD), lm.times(b, GMP_RNDU),
		 m1*b, m2*b, prec);
  }

  MInt* divide (const MInt* b) const;

  int sign (bool fail = true) const {
    if (m1.a == 0u && m2.a == 0u)
      return 0;
    int su = um.sign();
    if (su == -1)
      return -1;
    int sl = lm.sign();
    if (sl == 1)
      return 1;
    if (!fail || (sl == 0 && su == 0))
      return 0;
    throw signException;
  }

  MInt* mid () const {
    MValue m = lm.plus(um, GMP_RNDN).times(0.5, GMP_RNDN);
    return new EInt(m, m, prec);
  }

  bool subset (const MInt* that) const {
    const EInt *b = dynamic_cast<const EInt*>(that);
    if (!b) return pInt()->subset(that);
    return !(lm < b->lm) && !(b->um < um) && (b->lm < lm || um < b->um);
  }

  MInt* interval (const MInt* that) const {
    const EInt *b = dynamic_cast<const EInt*>(that);
    if (!b) return pInt()->interval(that);
    return new EInt(lm, b->um, prec);
  }

  MInt* innerInterval (const MInt* that) const {
    const EInt *b = dynamic_cast<const EInt*>(that);
    if (!b) return pInt()->innerInterval(that);
    return new EInt(um, b->lm, prec);
  }

  bool intersects (const MInt* that) const {
    const EInt *b = dynamic_cast<const EInt*>(that);
    if (!b) return pInt()->intersects(that);
    return !(um < b->lm || b->um < lm);
  }

  MInt* intersect (const MInt* that) const {
    const EInt *b = dynamic_cast<const EInt*>(that);
    if (!b) return pInt()->intersect(that);
    assert(!(um < b->lm || b->um < lm));
    MValue l = lm < b->lm ? b->lm : lm;
    MValue u = um < b->um ? um : b->um;
    return new EInt(l, u, prec);
  }

  MInt* sqrt () const {
    return new EInt(lm.sqrt(GMP_RNDD), um.sqrt(GMP_RNDU), prec);
  }

  MInt* root (unsigned long int k) const {
    return new EInt(lm.root(k, GMP_RNDD), um.root(k, GMP_RNDU), prec);
  }

  MValue lm, um;
};

inline MInt* MInt::make (unsigned int precision, const Parameter &p) {
  double l = p.lb(), u = p.ub();
  ModInt m1 = l == u ? ModInt(l, ModInt::p1) : ModInt(random(), ModInt::p1);
  ModInt m2 = l == u ? ModInt(l, ModInt::p2) : ModInt(random(), ModInt::p2);
  if (precision == 106u)
    return new PInt(p, m1, m2);
  else
    return new EInt(MValue(l, precision), MValue(u, precision), m1, m2, precision);
}

inline MInt* MInt::make (unsigned int precision, const MInt *m) {
  if (m->precision() == 106u)
    return new EInt(MValue(m->lb(), precision), MValue(m->ub(), precision), 
                    m->m1, m->m2, precision);
  else {
    const EInt* e = dynamic_cast<const EInt*>(m);
    return new EInt(MValue(e->lm, precision, GMP_RNDD),
		    MValue(e->um, precision, GMP_RNDU),
                    m->m1, m->m2, precision);
  }
}

template <class T>
T pow (T x, long int n)
{
  if (n == 0)
    return T(1);
  if (n < 0)
    return T(1)/pow(x, - n);
  return pow1(x, n);
 }

template <class T> 
T pow1 (T x, unsigned long int n)
{
  if (n == 1)
    return x;
  T y = pow(x, n/2);
  if (n %2 == 0)
    return y*y;
  else
    return x*y*y;
}
}
#endif
