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

#include "acp.h"
using namespace acp;

namespace acp {

double randomNumber (double rmin, double rmax)
{
  return rmin + (rmax - rmin)*random()/double(RAND_MAX);
}

const double Parameter::sentinel = 1e300;
double Parameter::delta = ::pow(2.0, -27);
bool Parameter::enabled = false;
bool Parameter::penabled = false;
unsigned int MInt::curPrecisions[128] = { 53u };
unsigned int MInt::threadIds[128];
vector<MInt*> MInt::mints[128];
unsigned int MInt::algebraicId = 0u;

double Parameter::algT = 1e-13;
double Parameter::algR;
double Parameter::sumBits = 0.0;
int Parameter::algI = -1;
int Parameter::algP;
unsigned long Parameter::nSign = 0;
unsigned int Parameter::nAZ = 0;
unsigned int Parameter::nANZ = 0;
Parameter Parameter::sqrt () const
{
  if (high())
    return Parameter(u.m->sqrt());
  int s = sign();
  assert(s > -1);
  return s == 0 ? Parameter(0, 0) : 
    Parameter(Parameter::sqrt(l).l, Parameter::sqrt(u.r).u.r);
}

Parameter Parameter::sqrt (double x)
{
  double s = ::sqrt(x);
  Parameter xx(x);
  Parameter lr = xx/Parameter(s);
  if (s < lr.l)
    return Parameter(s, lr.u.r);
  if (s > lr.u.r)
    return Parameter(lr.l, s);
  return lr;
}

Parameter Parameter::pow (const Parameter &x, unsigned long int n) {
  if (n == 1)
    return x;
  Parameter y = pow(x, n / 2);
  if (n %2 == 0)
    return y * y;
  else
    return x * y * y;
}

Parameter Parameter::root (unsigned long int n) const
{
  int s = sign();
  assert(s > -1);
  return s == 0 ? Parameter(0.0, 0.0) :
    Parameter(root(l, n).l, root(u.r, n).u.r);
}
  
Parameter Parameter::root (double a, unsigned int n)
{
  double x = a < 1 ? 1 : a;
  double x_old;
  do {
    x_old = x;
    Parameter xx(x);
    x = (xx - (xx - a / pow(xx, n-1)) / n).u.r;
  } while (x < x_old);
  Parameter xx(x);
  Parameter l = a/pow(xx, n-1);
  return Parameter(l.l, x);
}

unsigned int inverse (unsigned int a, unsigned int n)
{
  long t = 0, nt = 1, r = n, nr = a;
  while (nr != 0) {
    long q = r/nr, pt = t, pr = r;
    t = nt;
    nt = pt - q*nt;
    r = nr;
    nr = pr - q*nr;
  }
  return t < 0 ? t + n : t;
}

const int Modder::eShift = 53;
const int Modder::eMax = 1024 - eShift;
const int Modder::eMin = -1073 - eShift;

  //unsigned int Mods::primes[NPrimes] = { 4294967291u, 4294967279u, 4228232747u, 4197064799u, 3691300979u, 3116510503u, 4278467023u };
  //  unsigned int Mods::primes[NPrimes] = { 4228232747u, 4197064799u, 3691300979u, 3116510503u, 4278467023u, 4294967291u, 4294967279u };
  unsigned int Mods::primes[NPrimes] = { 4228232747u, 4197064799u, 3691300979u, 4228232747u, 4197064799u, 3691300979u, 3116510503u };
int Mods::primeIndex = NMods;
unsigned int Mods::ps[NMods];
Modder Mods::modder[NMods];
unsigned long Mods::nMods;
unsigned long Mods::nMixed;
double Mods::maxBits = 0;

  Mods initializeMods(0, 0, 0);

void Mods::changePrime (unsigned int p) {
  cout << "changePrime" << endl;
  for (int i = 0; i < NMods; i++)
    if (p == ps[i]) {
      // ps[i] = primes[primeIndex]; // PRIMES
      ps[i] = random32bitPrime();
      cout << "changeprime " << ps[i] << endl;
      modder[i] = Modder(ps[i]);
      primeIndex = (primeIndex + 1) % NPrimes;
      return;
    }
  assert(0);
}

MInt* EInt::times (const MInt* that) const
{
  const EInt *b = dynamic_cast<const EInt*>(that);
  if (!b) return pInt()->times(that);
  bool pflag = um.sign() == -1;
  MValue sl = pflag ? um.minus() : lm, su = pflag ? lm.minus() : um,
    tl = pflag ? b->um.minus() : b->lm, tu = pflag ? b->lm.minus() : b->um;
  if (sl.sign() == 1) {
    MValue &l1 = tl.sign() == 1 ? sl : su, &u1 = tu.sign() == 1 ? su : sl;
    return new EInt(l1.times(tl, GMP_RNDD), u1.times(tu, GMP_RNDU),
		    mods*b->mods, min(prec, b->prec));
  }
  if (tl.sign() == 1)
    return new EInt(sl.times(tu, GMP_RNDD), su.times(tu, GMP_RNDU),
		    mods*b->mods, min(prec, b->prec));
  if (tu.sign() == -1)
    return new EInt(su.times(tl, GMP_RNDD), sl.times(tl, GMP_RNDU),
		    mods*b->mods, min(prec, b->prec));
  MValue cl1 = sl.times(tu, GMP_RNDD), cl2 = su.times(tl, GMP_RNDD),
    cu1 = sl.times(tl, GMP_RNDU), cu2 = su.times(tu, GMP_RNDU);
  return new EInt(cl1 < cl2 ? cl1 : cl2, cu1 < cu2 ? cu2 : cu1,
		  mods*b->mods, min(prec, b->prec));
}

MInt* EInt::divide (const MInt* that) const
{
  const EInt *b = dynamic_cast<const EInt*>(that);
  if (!b) return pInt()->divide(that);
  int as = sign(false), bs = b->sign();
  if (bs == 1)
    switch (as) {
    case 1:
      return new EInt(lm.divide(b->um, GMP_RNDD), um.divide(b->lm, GMP_RNDU),
		      mods/b->mods, min(prec, b->prec));
    case 0:
      return new EInt(lm.divide(b->lm, GMP_RNDD), um.divide(b->lm, GMP_RNDU),
		      mods/b->mods, min(prec, b->prec));
    case -1:
      return new EInt(lm.divide(b->lm, GMP_RNDD), um.divide(b->um, GMP_RNDU),
		      mods/b->mods, min(prec, b->prec));
    }
  switch (as) {
  case 1:
    return new EInt(um.divide(b->um, GMP_RNDD), lm.divide(b->lm, GMP_RNDU),
		    mods/b->mods, min(prec, b->prec));
  case 0:
    return new EInt(um.divide(b->um, GMP_RNDD), lm.divide(b->um, GMP_RNDU),
		    mods/b->mods, min(prec, b->prec));
  case -1:
    return new EInt(um.divide(b->lm, GMP_RNDD), lm.divide(b->um, GMP_RNDU),
		    mods/b->mods, min(prec, b->prec));
  }
  return 0;
}

/* 
 * calculates (a * b) % c taking into account that a * b might overflow 
 */
unsigned long mulmod(unsigned long a, unsigned long b, unsigned long mod)
{
  unsigned long x = 0,y = a % mod;
  while (b > 0)
    {
      if (b % 2 == 1)
        {    
	  x = (x + y) % mod;
        }
      y = (y * 2) % mod;
      b /= 2;
    }
  return x % mod;
}
/* 
 * modular exponentiation
 */
unsigned long modulo(unsigned long base, unsigned long exponent, unsigned long mod)
{
  unsigned long x = 1;
  unsigned long y = base;
  while (exponent > 0)
    {
      if (exponent % 2 == 1)
	x = (x * y) % mod;
      y = (y * y) % mod;
      exponent = exponent / 2;
    }
  return x % mod;
}
 
/*
 * Miller-Rabin Primality test, iteration signifies the accuracy
 */
int Miller(unsigned long p,int iteration)
{
 
  int i;
  unsigned long s;
  if (p < 2)
    {
      return 0;
    }
  if (p != 2 && p % 2==0)
    {
      return 0;
    }
  s = p - 1;
  while (s % 2 == 0)
    {
      s /= 2;
    }
  for (i = 0; i < iteration; i++)
    {
      unsigned long a = random() % (p - 1) + 1, temp = s;
      unsigned long mod = modulo(a, temp, p);
      while (temp != p - 1 && mod != 1 && mod != p - 1)
        {
	  mod = mulmod(mod, mod, p);
	  temp *= 2;
        }
      if (mod != p - 1 && temp % 2 == 0)
        {
	  return 0;
        }
    }
  return 1;
}


bool isPrime (unsigned int p) {
  return Miller(p, 40);
}

unsigned int random32bitPrime () {
  while (true) {
    unsigned int p = random();
    p |= (1 << 31);
    if (!isPrime(p))
      continue;
    return p;
  }
}


}


