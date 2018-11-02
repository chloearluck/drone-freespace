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
SignException signException;
unsigned int MInt::curPrecisions[128] = { 53u };
vector<MInt*> MInt::mints[128];

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

const int ModP::eShift = 53;
const int ModP::eMax = 1024 - eShift;
const int ModP::eMin = -1073 - eShift;

unsigned int ModInt::ps[7] = { 4294967291u, 4294967279u, 4228232747u, 4197064799u, 3691300979u, 3116510503u, 4278467023u };
int ModInt::pIndex = 2;
unsigned int ModInt::p1 = ps[0];
unsigned int ModInt::p2 = ps[1];
ModP ModInt::modP1(ModInt::ps[0]);
ModP ModInt::modP2(ModInt::ps[1]);

void ModInt::changePrime (unsigned int p) {
  if (p == p1) {
    p1 = ps[pIndex];
    modP1 = ModP(p);
  }
  else if (p == p2) {
    p2 = ps[pIndex];
    modP2 = ModP(p);
  }
  else
    assert(0);
  pIndex = (pIndex + 1) % 7;
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
		    m1*b->m1, m2*b->m2, min(prec, b->prec));
  }
  if (tl.sign() == 1)
    return new EInt(sl.times(tu, GMP_RNDD), su.times(tu, GMP_RNDU),
		    m1*b->m1, m2*b->m2, min(prec, b->prec));
  if (tu.sign() == -1)
    return new EInt(su.times(tl, GMP_RNDD), sl.times(tl, GMP_RNDU),
		    m1*b->m1, m2*b->m2, min(prec, b->prec));
  MValue cl1 = sl.times(tu, GMP_RNDD), cl2 = su.times(tl, GMP_RNDD),
    cu1 = sl.times(tl, GMP_RNDU), cu2 = su.times(tu, GMP_RNDU);
  return new EInt(cl1 < cl2 ? cl1 : cl2, cu1 < cu2 ? cu2 : cu1,
		  m1*b->m1, m2*b->m2, min(prec, b->prec));
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
		      m1/b->m1, m2/b->m2, min(prec, b->prec));
    case 0:
      return new EInt(lm.divide(b->lm, GMP_RNDD), um.divide(b->lm, GMP_RNDU),
		      m1/b->m1, m2/b->m2, min(prec, b->prec));
    case -1:
      return new EInt(lm.divide(b->lm, GMP_RNDD), um.divide(b->um, GMP_RNDU),
		      m1/b->m1, m2/b->m2, min(prec, b->prec));
    }
  switch (as) {
  case 1:
    return new EInt(um.divide(b->um, GMP_RNDD), lm.divide(b->lm, GMP_RNDU),
		    m1/b->m1, m2/b->m2, min(prec, b->prec));
  case 0:
    return new EInt(um.divide(b->um, GMP_RNDD), lm.divide(b->um, GMP_RNDU),
		    m1/b->m1, m2/b->m2, min(prec, b->prec));
  case -1:
    return new EInt(um.divide(b->lm, GMP_RNDD), lm.divide(b->um, GMP_RNDU),
		    m1/b->m1, m2/b->m2, min(prec, b->prec));
  }
  return 0;
}

}


