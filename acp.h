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

//#define NO_MODE
//#define USE_ASSERT

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

unsigned int random32bitPrime ();

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
  const bool algebraic;
  SignException (bool algebraic=false) : algebraic(algebraic) {}
  virtual const char* what() const throw() {
    return "Not enough precision";
  }
};

class MixedModException : public std::exception {
public:
  MixedModException () {}
  virtual const char* what() const throw() {
    return "Zero modulo some but not all primes";
  }
};

class MixedHomotopyException : public std::exception {
public:
  MixedHomotopyException () {}
  virtual const char* what() const throw() {
    return "Some but not all homotopies are zero";
  }
};

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
  bool operator== (const MValue &b)  const { return mpfr_equal_p(m, b.m) != 0; }
  bool operator!= (const MValue &b)  const { return mpfr_equal_p(m, b.m) == 0; }
  bool operator< (const MValue &b)  const { return mpfr_less_p(m, b.m); }
  MValue sqrt (mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_sqrt(res.m, m, round);
    return res;
  }
  MValue root (unsigned long int k, mpfr_rnd_t round) const {
    MValue res(p);
    // mpfr_rootn_ui(res.m, m, k, round);
    mpfr_root(res.m, m, k, round);
    return res;
  }

  mpfr_t m;
  unsigned int p;
};

//typedef unsigned int uint;

//typedef unsigned long ulong;

unsigned int inverse (unsigned int a, unsigned int n);

#define ulong(a) ((unsigned long)(a))
 
class Mod {
 public:
  Mod (unsigned int a=0, unsigned int p=0) : a(a), p(p) {}
  
  Mod operator- () const {
    long l = (-long(a))%p;
    return Mod(l < 0 ? l + p : l, p);
  }
  
  Mod operator+ (const Mod &x) const {
    return Mod((long(a) + long(x.a))%p, p);
  }

  Mod operator- (const Mod &x) const {
    long l = (long(a) - long(x.a))%p;
    return Mod(l < 0 ? l + p : l, p);
  }

  Mod operator* (const Mod &x) const {
    return Mod((ulong(a)*ulong(x.a))%p, p);
  }

  Mod operator/ (const Mod &x) const {
    if (x.a == 0)
      throw p;
    return Mod((ulong(a)*ulong(inverse(x.a, p)))%p, p);
  }

  unsigned int a, p;
};

class Modder {
  static const int eShift;
  static const int eMax;
  static const int eMin;

  vector<unsigned int> pow2v;
  unsigned int *pow2;
  unsigned int p;

public:
  Modder () : pow2v(0), pow2(0), p(0) {}
  Modder (unsigned int p) : p(p), pow2v(eMax - eMin + 1) {
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

  Mod mod (double x) {
    if (x == (long) x)
      return Mod((x >= 0 ? ((long) x) % p : p + ((long) x) % p), p);
    int e;
    double m = frexp(x, &e) * (1l << eShift);
    assert(m == (long) m);
    e -= eShift;
    long mp = ((long) m) % p;
    if (mp < 0)
      mp += p;
    assert(eMin <= e && e <= eMax);
    return Mod(((unsigned long) mp * pow2[e]) % p, p);
  }
};

class Parameter;
class EInt;

#define NPrimes 7
#define NMods 2

/*
class Mods {
 public:
  static unsigned int primes[NPrimes];
  static int primeIndex;
  static unsigned int ps[NMods];
  static Modder modder[NMods];

  static void changePrime (unsigned int p);

  Mod mod[NMods];
  const bool algebraic;

  Mods (int start, int up) : algebraic(false) {
    for (int i = 0; i < NMods; i++) {
      assert(ps[i] == 0);
      ps[i] = primes[i];
      modder[i] = Modder(ps[i]);
    }
  }

  Mods (bool algebraic) : algebraic(algebraic) {}

  Mods (double r) : algebraic(false) {
    for (int i = 0; i < NMods; i++)
      mod[i] = modder[i].mod(r);
  }

  Mods (const Parameter &p);
  Mods (const MValue &l, const MValue &u);
  
  bool hasTheRightPrimes () const {
    for (int i = 0; i < NMods; i++)
      if (mod[i].p != ps[i])
        return false;
    return true;
  }

  bool mixed () const {
    for (int i = 1; i < NMods; i++)
      if ((mod[i].a == 0) != (mod[0].a == 0))
        return true;
    return false;
  }

  bool zero () const {
    for (int i = 1; i < NMods; i++)
      if (mod[i].a != 0)
        return false;
    return true;
  }

  Mods operator- () const {
    Mods m(algebraic);
    for (int i = 0; i < NMods; i++)
      m.mod[i] = -mod[i];
    return m;
  }
  
  Mods operator+ (const Mods &x) const {
    Mods m(algebraic || x.algebraic);
    for (int i = 0; i < NMods; i++)
      m.mod[i] = mod[i] + x.mod[i];
    return m;
  }

  Mods operator+ (double x) const {
    return *this + Mods(x);
  }

  Mods operator- (const Mods &x) const {
    Mods m(algebraic || x.algebraic);
    for (int i = 0; i < NMods; i++)
      m.mod[i] = mod[i] - x.mod[i];
    return m;
  }

  Mods operator- (double x) const {
    return *this - Mods(x);
  }

  Mods operator* (const Mods &x) const {
    Mods m(algebraic || x.algebraic);
    for (int i = 0; i < NMods; i++)
      m.mod[i] = mod[i] * x.mod[i];
    return m;
  }

  Mods operator* (double x) const {
    return *this * Mods(x);
  }

  Mods operator/ (const Mods &x) const {
    Mods m(algebraic || x.algebraic);
    for (int i = 0; i < NMods; i++)
      m.mod[i] = mod[i] / x.mod[i];
    return m;
  }

  Mods operator/ (double x) const {
    return *this / Mods(x);
  }
};
*/

class Mods {
 public:
  static unsigned int primes[NPrimes];
  static int primeIndex;
  static unsigned int ps[NMods];
  static Modder modder[NMods];
  static unsigned long nMods;
  static unsigned long nMixed;
  static double maxBits;

  static void changePrime (unsigned int p);

  Mod mod[NMods];
  const bool algebraic;
  bool constant;
  double bitsA, bitsB;

  void noBits () {
    bitsA = 0;
    bitsB = 0;
  }

  Mods (int start, int it, int up) : algebraic(false), constant(false) {
    for (int i = 0; i < NMods; i++) {
      assert(ps[i] == 0);
      // ps[i] = primes[i]; // PRIMES
      ps[i] = random32bitPrime();
      // cout << "p[" << i << "] = " << ps[i] << endl;
      modder[i] = Modder(ps[i]);
    }
  }

  Mods (bool algebraic, bool constant, double bitsA, double bitsB) : algebraic(algebraic), constant(constant), bitsA(bitsA), bitsB(bitsB) {
    if (bitsA > maxBits)
      maxBits = bitsA;
  }

  Mods (double r) : algebraic(false), constant(true), bitsA(53), bitsB(0) {
    for (int i = 0; i < NMods; i++)
      mod[i] = modder[i].mod(r);
  }

  Mods (const Parameter &p);
  Mods (const MValue &l, const MValue &u);
  
  bool hasTheRightPrimes () const {
    for (int i = 0; i < NMods; i++)
      if (mod[i].p != ps[i])
        return false;
    return true;
  }

  bool mixed () const {
    for (int i = 1; i < NMods; i++)
      if ((mod[i].a == 0) != (mod[0].a == 0))
        return true;
    return false;
  }

  void checkMixed () {
    nMods++;
    if (mixed()) {
      nMixed++;
      cout << (mod[1].a==0) << " nMods " << nMods << " nMixed " << nMixed
	   << " bits " << bitsA << endl;
    }
  }

  bool zero () const {
    for (int i = 0; i < NMods; i++)
      if (mod[i].a != 0)
        return false;
    return true;
  }

  Mods operator- () const {
    Mods m(algebraic, constant, bitsA, bitsB);
    for (int i = 0; i < NMods; i++)
      m.mod[i] = -mod[i];
    return m;
  }
  
  Mods operator+ (const Mods &x) const {
    Mods m(algebraic || x.algebraic, constant && x.constant,
	   std::max(bitsA + x.bitsB, bitsB + x.bitsA), bitsB + x.bitsB);
    for (int i = 0; i < NMods; i++)
      m.mod[i] = mod[i] + x.mod[i];
    if (!(mixed() && x.mixed()))
      m.checkMixed();
    else
      nMods++;
    return m;
  }

  Mods operator+ (double x) const {
    Mods ret = *this + Mods(x);
    ret.bitsA = bitsA;
    ret.bitsB = bitsB;
    return ret;
  }

  Mods operator- (const Mods &x) const {
    Mods m(algebraic || x.algebraic, constant && x.constant,
	   std::max(bitsA + x.bitsB, bitsB + x.bitsA), bitsB + x.bitsB);
    for (int i = 0; i < NMods; i++)
      m.mod[i] = mod[i] - x.mod[i];
    if (!(mixed() && x.mixed()))
      m.checkMixed();
    else
      nMods++;
    return m;
  }

  Mods operator- (double x) const {
    Mods ret = *this - Mods(x);
    ret.bitsA = bitsA;
    ret.bitsB = bitsB;
    return ret;
  }

  Mods operator* (const Mods &x) const {
    Mods m(algebraic || x.algebraic,  constant && x.constant,
	   bitsA + x.bitsA, bitsB + x.bitsB);
    for (int i = 0; i < NMods; i++)
      m.mod[i] = mod[i] * x.mod[i];
    // nMods++;
    return m;
  }

  Mods operator* (double x) const {
    Mods ret = *this * Mods(x);
    ret.bitsA = bitsA;
    ret.bitsB = bitsB;
    return ret;
  }

  Mods operator/ (const Mods &x) const {
    Mods m(algebraic || x.algebraic,  constant && x.constant,
	   bitsA + x.bitsB, bitsB + x.bitsA);
    for (int i = 0; i < NMods; i++)
      m.mod[i] = mod[i] / x.mod[i];
    // nMods++;
    return m;
  }

  Mods operator/ (double x) const {
    Mods ret = *this / Mods(x);
    ret.bitsA = bitsA;
    ret.bitsB = bitsB;
    return ret;
  }
};



extern pthread_key_t idkey;
extern pthread_key_t mikey;
extern pthread_key_t cpkey;
extern pthread_key_t hsekey;

class MInt {
  static unsigned int & curPrecision () {
    void *v = pthread_getspecific(cpkey);
    return v ? *(unsigned int *) v : curPrecisions[0];
  }

public:
  void noBits () { mods.noBits(); } 

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

  static unsigned int threadIds[128];

  static unsigned int threadId () {
    void *i = pthread_getspecific(idkey);
    return i ? *(unsigned int *) i : 0;
  }

  static unsigned int algebraicId;

  // const Mods mods;
  Mods mods; // vjm debug

  /*
  MInt () : m1(Mod(random(), Mod::p1)), m2(Mod(random(), Mod::p2)),
            algebraic(true) {
    getMInts().push_back(this);
  }
  */
  MInt (const Mods &mods) : mods(mods) {
    getMInts().push_back(this);
  }

  virtual ~MInt () {}

  static void clearMInts () {
    vector<MInt*> &mints = getMInts();
    for (int i = 0; i < mints.size(); i++)
      delete mints[i];
    mints.clear();
  }

  bool zeroMod () { assert(!mods.mixed()); return mods.zero(); }

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

public:  // DEBUG
static double algT, algR, sumBits;
static int algI, algP;
static unsigned long nSign;
static unsigned int nAZ, nANZ;

private:

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
  void noBits () { u.m->noBits(); }

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
    if (!penabled) {
      penabled = true;
      pthread_key_create(&idkey, NULL);
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
    pthread_key_delete(idkey);
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
    if (MInt::getPrecision() > 53u) {
      p.increasePrecision();
      p.u.m->mods.constant = true; // debug
    }
    return p;
  }
  
  static Parameter input (double x) { 
    double y = x + delta*(1.0 + fabs(x))*randomNumber(-1.0, 1.0);
    return Parameter::constant(y);
  }
  
  bool low () const { return l != sentinel; }
  bool high () const { return l == sentinel; }

  bool hasTheRightPrimes () {
    return u.m->mods.hasTheRightPrimes();
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
    nSign++;
    if (high())
      return u.m->sign(fail);
    if (l > 0.0)
      return 1;
    if (u.r < 0.0)
      return -1;
    if (!fail || (l == 0.0 && u.r == 0.0))
      return 0;
    throw SignException();
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
  PInt (const Parameter &p, const Mods &mods) : MInt(mods), p(p) {}

  unsigned int precision () const { return 106u; }

  double intervalWidth () const { return p.intervalWidth(); }
  
  double lb () const { return p.lb(); }
  
  double ub () const { return p.ub(); }
  
  MInt* lbP () const { return new PInt(p.lbP(), Mods(p.lbP())); }

  MInt* ubP () const { return new PInt(p.ubP(), Mods(p.ubP())); }

  MInt* mid () const { return new PInt(p.midP(), Mods(p.midP())); }

  MInt* plus (const MInt* that) const {
    return new PInt(p + that->par(), mods + that->mods);
  }
  
  MInt* plus (double b) const { 
    return new PInt(p + Parameter(b), mods + b);
  }
  
  MInt* minus (const MInt* that) const {
    return new PInt(p - that->par(), mods - that->mods);
  }
  
  MInt* minus () const { return new PInt(-p, -mods); }
  
  MInt* times (const MInt* that) const {
    return new PInt(p * that->par(), mods *that->mods);
  }
  
  MInt* times (double b) const { 
    return new PInt(p * Parameter(b), mods * b);
  }

  MInt* divide (const MInt* that) const {
    return new PInt(p / that->par(), mods / that->mods);
  }

  int sign (bool fail = true) const {
    if (p.l > 0.0) {
      if (fail && Parameter::algT > 0 && MInt::threadId() == MInt::algebraicId) {
        double rat = p.l == p.u.r ? Parameter::algR : p.l/(p.u.r-p.l);
        if (rat < Parameter::algR) {
          // cerr << "rat1 " << rat << endl;
          Parameter::algR = rat;
        }
      }
      return 1;
    }
    if (p.u.r < 0.0) {
      if (fail && Parameter::algT > 0 && MInt::threadId() == MInt::algebraicId) {
        double rat = p.l == p.u.r ? Parameter::algR : p.u.r/(p.l-p.u.r);
        if (rat < Parameter::algR) {
          // cerr << "rat2 " << rat << endl;
          Parameter::algR = rat;
        }
      }
      return -1;
    }
    if (!fail || (p.l == 0 && p.u.r == 0))
      return 0;
    if (mods.mixed())
      throw MixedModException();
    if (mods.zero()) {
      Parameter::nAZ++;
      Parameter::sumBits += mods.bitsA;
      return 0;
    }
    Parameter::nANZ++;
    throw SignException(mods.algebraic);
  }

  bool subset (const MInt* that) const {
    return p.subset(that->par());
  }

  MInt* interval (const MInt* that) const {
    Parameter q = p.interval(that->par());
    return new PInt(q, Mods(q));
  }

  MInt* innerInterval (const MInt* that) const {
    Parameter q = p.innerInterval(that->par());
    return new PInt(q, Mods(q));
  }

  bool intersects (const MInt* that) const {
    return p.intersects(that->par());
  }

  MInt* intersect (const MInt* that) const {
    Parameter q = p.intersect(that->par());
    return new PInt(q, Mods(q));
  }

  MInt* sqrt () const {
    Parameter q = p.sqrt();
    return new PInt(q, Mods(q));
  }

  MInt* root (unsigned long int k) const { 
    assert(0);
    Parameter q = p.root(k);
    return new PInt(q, Mods(q));
  }
};

class EInt : public MInt {
  int prec;
public:
  EInt (const MValue &l, const MValue &u, int prec) 
    : MInt(Mods(l, u)), lm(l), um(u), prec(prec) {}
  EInt (const MValue &l, const MValue &u, const Mods &mods, int prec)
    : MInt(mods), lm(l), um(u), prec(prec) {}
  ~EInt () {}

  PInt *pInt () const { return new PInt(Parameter(lb(), ub()), mods); }
 
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
		    mods + b->mods, min(prec, b->prec));
  }

  MInt* plus (double b) const {
    return new EInt(lm.plus(b, GMP_RNDD), um.plus(b, GMP_RNDU),
		    mods + b, prec);
  }

  MInt* minus (const MInt* that) const {
    const EInt *b = dynamic_cast<const EInt*>(that);
    if (!b) return pInt()->minus(that);
    return new EInt(lm.minus(b->um, GMP_RNDD), um.minus(b->lm, GMP_RNDU),
		    mods - b->mods, min(prec, b->prec));
  }

  MInt* minus () const { 
    return new EInt(um.minus(), lm.minus(), -mods, prec);
  }

  MInt* times (const MInt* that) const;

  MInt* times (double b) const {
    return b > 0.0 ? new EInt(lm.times(b, GMP_RNDD), um.times(b, GMP_RNDU),
			      mods*b, prec)
      : new EInt(um.times(b, GMP_RNDD), lm.times(b, GMP_RNDU),
		 mods*b, prec);
  }

  MInt* divide (const MInt* b) const;
   
  int sign (bool fail = true) const {
    int su = um.sign();
    if (su == -1) {
      if (mods.mixed() || mods.zero()) {
	cerr << "nonzero sign with zero mod" << endl;
	//eps removed this: exit(0);
      }
      if (fail && Parameter::algT > 0 && MInt::threadId() == MInt::algebraicId) {
        double rat = mpfr_get_d(um.divide(lm.minus(um, GMP_RNDD), GMP_RNDU).m, GMP_RNDU);
        if (rat < Parameter::algR) {
          // cerr << "rat3 " << rat << endl;
          Parameter::algR = rat;
        }
      }
      return -1;
    }
    int sl = lm.sign();
    if (sl == 1) {
      if (mods.mixed() || mods.zero()) {
	cerr << "nonzero sign with zero mod" << endl;
	//exit(0);
      }
      if (fail && Parameter::algT > 0 && MInt::threadId() == MInt::algebraicId) {
        double rat = mpfr_get_d(lm.divide(um.minus(lm, GMP_RNDU), GMP_RNDU).m, GMP_RNDU);
        if (rat < Parameter::algR) {
          // cerr << "rat4 " << rat << endl;
          Parameter::algR = rat;
        }
      }
      return 1;
    }
    if (!fail || (sl == 0 && su == 0))
      return 0;
    if (mods.mixed())
      throw MixedModException();
    if (mods.zero())
      return 0;
    throw SignException(mods.algebraic);
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
  if (precision == 106u)
    return new PInt(p, Mods(p));
  else
    return new EInt(MValue(l, precision), MValue(u, precision), Mods(p), precision);
}

inline MInt* MInt::make (unsigned int precision, const MInt *m) {
  if (m->precision() == 106u)
    return new EInt(MValue(m->lb(), precision), MValue(m->ub(), precision), 
                    m->mods, precision);
  else {
    const EInt* e = dynamic_cast<const EInt*>(m);
    return new EInt(MValue(e->lm, precision, GMP_RNDD),
		    MValue(e->um, precision, GMP_RNDU),
                    m->mods, precision);
  }
}

inline Mods::Mods (const Parameter &p)
  : algebraic(p.lb()<p.ub()), constant(false), bitsA(53), bitsB(0) {
  //inline Mods::Mods (const Parameter &p) : algebraic(p.lb()<p.ub()) {
  double l = p.lb(), u = p.ub();
  if (l == u)
    for (int i = 0; i < NMods; i++)
      mod[i] = modder[i].mod(l);
  else
    for (int i = 0; i < NMods; i++)
      mod[i] = Mod(random() % ps[i], ps[i]);
}

inline Mods::Mods (const MValue &l, const MValue &u)
  : algebraic(l != u), constant(false), bitsA(53), bitsB(0) {
  //inline Mods::Mods (const MValue &l, const MValue &u) : algebraic(l != u) {
  if (l != u)
    for (int i = 0; i < NMods; i++)
      mod[i] = Mod(random() % ps[i], ps[i]);
  else
    for (int i = 0; i < NMods; i++) {
      MValue m = l; 
      unsigned int ret = 0;
      while (m.sign() != 0) {
        double x = m.value();
        m = m.minus(MValue(x, m.p), GMP_RNDN);
        ret = (ret + modder[i].mod(x).a) % ps[i];
      }
      mod[i] = Mod(ret, ps[i]);
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
