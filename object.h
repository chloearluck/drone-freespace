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

#ifndef OBJECT
#define OBJECT

#include "pv.h"
#include <map>
#include <limits.h>
#include<typeinfo>
#include <sys/time.h>

double getTime (); 

namespace acp {
  
class RefCnt {
  template<class T> friend class PTR;
  int refCnt;
  
  void incRef () {
    pthread_mutex_lock(&mutex);
    refCnt++;
    pthread_mutex_unlock(&mutex);
  }

  void decRef () { 
    pthread_mutex_lock(&mutex);
    bool flag = --refCnt == 0;
    pthread_mutex_unlock(&mutex);
    if (flag)
      delete this;
  }
protected:
  pthread_mutex_t mutex;
public:
  RefCnt () : refCnt(0) { mutex = PTHREAD_MUTEX_INITIALIZER; }
  virtual ~RefCnt () { assert(refCnt == 0); }
};

template<class T>
class PTR {
  T *t;
public:
  PTR () : t(0) {}
  PTR (T *t) : t(t) { incRef(); }
  PTR (const PTR &p) : t(p.t) { incRef(); }
  const PTR &operator= (const PTR &p) { 
    p.incRef(); decRef(); t = p.t; return *this; 
  }
  ~PTR () { decRef(); }
  void incRef () const { if (t != 0) t->incRef(); }
  void decRef () const { if (t != 0) t->decRef(); }
  operator T* () const { return t; }
  T *operator-> () const { return t; }
};

extern pthread_key_t idkey;
extern pthread_key_t mikey;
extern pthread_key_t cpkey;
extern pthread_key_t hsekey;

class BaseObject : public RefCnt {
protected:
  static unsigned int deltaPrecision;
  static unsigned int maxPrecision;
  static bool throwSignExceptions[128];

  static unsigned int getPrecision () { return MInt::getPrecision(); }
  static void setPrecision (unsigned int p) { return MInt::setPrecision(p); }

public:
  static bool usePrecisionException;

  static unsigned int nAlg, nAmb, nAmbNZ, nHomo;
  static double minHomo, maxWidth, tHom;

  static bool & throwSignException () {
    void *v = pthread_getspecific(hsekey);
    return v ? *(bool *) v : throwSignExceptions[0];
  }

  BaseObject () {}
  virtual ~BaseObject () {}

  static void addThread (unsigned int i) {
    assert(i < 128);
    pthread_setspecific(idkey, (void *) (MInt::threadIds + i));
    MInt::threadIds[i] = i;
    pthread_setspecific(mikey, (void *) (MInt::mints + i));
    pthread_setspecific(cpkey, (void *) (MInt::curPrecisions + i));
    pthread_setspecific(hsekey, (void *) (throwSignExceptions + i));
    MInt::curPrecisions[i] = 53u;
  }   
};

class PrecisionException : public std::exception {
public:
  virtual const char* what() const throw() {
    return "Maximum precision exceeded";
  }
};

extern PrecisionException precisionException;

template<class P>
bool checkAccuracy (const P &p, double acc)
{
  if (acc == 1.0) return true;  
  for (int i = 0; i < p.size(); i++) {
    double l = p[i].lb(), u = p[i].ub();  
    if ((l < 0.0 && u >= 0.0) || (l <= 0.0 && u > 0.0))
      return false;
    if ((l > 0.0 && u - l > l*acc) || (u < 0.0 && u - l > - u*acc)) {
      double ln = nextafter(l, 1.0 + u);
      if (ln < u) {
	ln = nextafter(ln, 1.0 + u);
	if (ln < u)
	  return false;
      }
    }
  }
  return true;
}

template<class P>
class Object :  public BaseObject {
public:
  template<class T> friend class PTR;
  template<class T> friend class ObjPTR;

  class PP {
  public:
    P p[2];
    double t;
    PP () : t(0) {}
    ~PP () {
      for (int i = 0; i < 2; i++)
	for (int j = 0; j < p[i].size(); j++)
	  if (p[i][j].high())
	    delete p[i][j].u.m;
    }
  };

private:
  friend class Parameter;
  
  P p;

  PP *pp;

  virtual P calculate () {
    if (!input())
      cerr << "Object<P> missing calculate method" << endl;
    assert(0);
  }

public:
  Object () : pp(0) {}
  
  Object (const P& q, bool constant=false) : p(q), pp(0) {
    for (int i = 0; i < p.size(); i++) {
      if (p[i].high())
        p[i].decreasePrecision();
      if (p[i].lb() != p[i].ub()) {
        std::cerr << "Input object has non-trivial interval." << std::endl;
        assert(0);
      }
    }
    
    if (constant)
      pp = (PP*)1;
  }
  
  ~Object () { 
    for (int i = 0; i < p.size(); i++)
      if (p[i].high())
        delete p[i].u.m;
    if (pp != 0 && pp != (PP*)1)
      delete[] pp;
  }
  
  int precision () const { return p.size() == 0 ? 0 : p[0].precision(); }
  bool uninitialized () const {
    return p.size() == 0 ? 0 : p[0].uninitialized();
  }
  
  bool input () {
    if (p.size() == 0)
      return false;
    for (int i = 0; i < p.size(); i++)
      if (p[i].high() || p[i].lb() != p[i].ub())
        return false;
    return true;
  }
  
  P getCurrentP () const {
    P q = p;
    if (MInt::getPrecision() > 53u) {
      if (precision() == 53u)
        for (int i = 0; i < q.size(); i++)
          q[i].increasePrecision();
      else {
        unsigned int precObject = precision();
        for (int i = 0; i < q.size(); i++)
          q[i].u.m = MInt::make(precObject, q[i].u.m);
      }
    }
    
    return q;
  }
  
  P get () {
    try {
      pthread_mutex_lock(&mutex);
      P q = get1();
      pthread_mutex_unlock(&mutex);
      return q;
    }
    catch (SignException se) {
      pthread_mutex_unlock(&mutex);
      throw se;
    }
    catch (unsigned int p) {
      pthread_mutex_unlock(&mutex);
      throw p;
    }
  }
  
private:
  P get1 () {
    if (Parameter::algI >= 0 && MInt::threadId() == MInt::algebraicId && pp != (PP*)1) {
      double algT = Parameter::algT;
      int algI = Parameter::algI;
      int algP = Parameter::algP;
      if (input()) {
        if (pp == 0 || pp[0].p[0][0].precision() < MInt::getPrecision()) {
	  if (pp != 0)
	    delete[] pp;
	  pp = new PP[2];
	  for (int ipp = 0; ipp < 2; ipp++) {
	    pp[ipp].p[0] = p;
	    for (int i = 0; i < pp[ipp].p[0].size(); i++)
	      pp[ipp].p[0][i] = Parameter::constant(randomNumber(-1, 1));
	    if (pp[ipp].p[0][0].high())
	      for (int i = 0; i < pp[ipp].p[0].size(); ++i)
		pp[ipp].p[0][i].u.m = pp[ipp].p[0][i].u.m->clone();
	    pp[ipp].t = 0;
	  }
        }

        Parameter::algI = -1;
        P q = get1();
        Parameter::algI = algI;
        if (Parameter::algP == 0) {
          for (int i = 0; i < p.size(); i++)
            if (pp[algI].p[0][i].sign() > 0)
              q[i] = q[i].interval(q[i] + pp[algI].p[0][i] * algT);
            else
              q[i] = (q[i] + pp[algI].p[0][i] * algT).interval(q[i]);
        }
        else {
          for (int i = 0; i < p.size(); i++)
            q[i] = q[i] + pp[algI].p[0][i] * algT;
        }
        return q;
      }
      else { // not input
        if (pp != 0 &&
	    ((algP == 0 && (pp[algI].t == -algT || pp[algI].t == algT)) ||
	     (algP == 1 && pp[algI].t == algT)))
          return pp[algI].p[algP];
        if (pp == 0)
          pp = new PP[2];
        pp[algI].p[algP] = calculate();
	if (pp[algI].p[algP][0].high())
	  for (int i = 0; i < pp[algI].p[algP].size(); ++i) {
	    pp[algI].p[algP][i].u.m = pp[algI].p[algP][i].u.m->clone();
	    if (algI == 1 && algP == 1) {
	      Parameter d0 = pp[0].p[algP][i] - p[i];
	      Parameter d1 = pp[1].p[algP][i] - p[i];
	      int s0 = d0.sign();
	      int s1 = d1.sign();
	      assert((s0 == 0) == (s1 == 0));
	      if (s0 != 0) {
		++nHomo;
		double h = s0 > 0 ? d0.lb() : -d0.ub();
		if (h < minHomo)
		  minHomo = h;
		h = s1 > 0 ? d1.lb() : -d1.ub();
		if (h < minHomo)
		  minHomo = h;
	      }
	      else
		assert(p[i].sign(false) == 0 || p[i].u.m->mods.constant);
	    }
	  }

        pp[algI].t = algP == 0 ? -algT : algT;

        return pp[algI].p[algP];
      }
    }

    unsigned int precObject = precision();
    unsigned int precNeeded = MInt::getPrecision();
    if (precObject >= precNeeded && 
        (precObject == 53u || p[0].hasTheRightPrimes())) {
      if (precObject == 53u)
        return p;
      if (precNeeded > 53u) {
	P q = p;
	for (int i = 0; i < q.size(); i++)
	  q[i].u.m = MInt::make(precObject, q[i].u.m);
	return q;
      }
      P q = p;
      for (int i = 0; i < q.size(); i++)
        q[i].decreasePrecision();
      return q;
    }

    if (input()) {
      P q = p;
      for (int i = 0; i < p.size(); i++)
        q[i].increasePrecision();
      if (pp == (PP*)1)
	for (int i = 0; i < p.size(); i++)
	  q[i].u.m->mods.constant = true;
      return q;
    }

    P q = calculate();
    if (precObject > 53u)
      for (int i = 0; i < p.size(); i++)
        delete p[i].u.m;
    if (precNeeded == 53u)
      return p = q;
    for (int i = 0; i < q.size(); i++)
      if (q[i].sign(false) == 0 && q[i].zeroMod())
	q[i] = Parameter::constant(0.0);
    p = q;
    for (int i = 0; i < p.size(); i++)
      p[i].u.m = p[i].u.m->clone();
    return q;
  }

public:
  int sign (int i); // sign of i'th coordinate

  P getApprox (double acc = 1e-17) {
    assert(MInt::getPrecision() == 53u);
    if (acc > 1)
      acc = 1;
    if (acc < 1) {
      getApprox(1); // make sure p is set so p.size() is correct
      bool *z = new bool [p.size()];
      for (int i = 0; i < p.size(); i++)
	z[i] = sign(i) == 0;
      pthread_mutex_lock(&mutex);
      unsigned int prec = precision();
      for (int i = 0; i < p.size(); ++i)
	if (z[i])
	  if (prec > 53u) {
	    p[i].l = Parameter::sentinel;
	    delete p[i].u.m;
	    p[i].u.m = MInt::make(prec, Parameter::constant(0));
	    MInt::getMInts().pop_back();
	  }
	  else
	    p[i].l = p[i].u.r = 0.0;
      pthread_mutex_unlock(&mutex);
      delete [] z;
    }

    try {
      P q = get();
      if (checkAccuracy(q, acc))
        return q;
    }
    catch (SignException se) {}
    MInt::setPrecision(212u);
    while (MInt::getPrecision() <= maxPrecision) {
      try {
        P q = get();
        if (checkAccuracy(q, acc)) {	
          MInt::setPrecision(53u);
	  return get();
        }
      }
      catch (SignException se) {}
      MInt::setPrecision(MInt::getPrecision() + deltaPrecision);
    }
    if (usePrecisionException || precision() == 0)
      throw precisionException;
    setPrecision(53u);
    return get();
  }
};

template<class P>
class ObjPTR : public PTR< Object<P> > {
public:
  ObjPTR () {}
  ObjPTR (Object<P> *o) : PTR< Object<P> >(o) {}
  ObjPTR (const ObjPTR &p) : PTR< Object<P> >(p) {}
  const ObjPTR &operator= (const ObjPTR &p) {
    PTR< Object<P> >::operator=(p);
    return p;
  }
};

class Primitive : public BaseObject {
  virtual Parameter calculate () {
    cerr << "Primitive::calculate not overridden" << endl;
    assert(0);
    return Parameter::constant(0);
  }

  static pthread_mutex_t algM;

  virtual int sign () { return calculate().sign(); }
public:
  static double tMod;

  operator int () {
    if (MInt::getPrecision() > 53u || throwSignException())
      return sign();
    try {
      return sign();
    }
    catch (SignException se) {}
    setPrecision(106u);
    while (MInt::getPrecision() <= maxPrecision) {
      double t0;
      try {
	t0 = getTime();
        int s = sign();
	if (MInt::getPrecision() == 106u)
	  tMod += getTime() - t0;
        MInt::setPrecision(53u);
        return s;
      }
      catch (unsigned int p) {
	if (MInt::getPrecision() == 106u)
	  tMod += getTime() - t0;
        Mods::changePrime(p);
      }
      catch (SignException se) {
	if (MInt::getPrecision() == 106u)
	  tMod += getTime() - t0;
        if (MInt::getPrecision() == 106u && se.algebraic) {
          pthread_mutex_lock(&algM);
	  MInt::algebraicId = MInt::threadId();
	  double t0 = getTime();
          bool a1 = algebraicIdentity(0);
          bool a2 = algebraicIdentity(1);
	  tHom += getTime() - t0;
	  MInt::algebraicId = 0u;
          pthread_mutex_unlock(&algM);
          if (a1 != a2)
            throw MixedHomotopyException();
          if (a1) {
	    MInt::setPrecision(53u);
            return 0;
          }
        }
        if (getPrecision() == 106u)
          MInt::setPrecision(212u);
        else
          MInt::setPrecision(MInt::getPrecision() + deltaPrecision);
      }
    }
    if (usePrecisionException)
      throw precisionException;
    MInt::setPrecision(53u);
    return 0;
  }

  bool algebraicIdentity (int algI) {
    ++nAlg;
    //return false;

    bool fail106 = true;

    if (!fail106) {
      Parameter::algT = 1e-14;
      Parameter::algR = 1e9;
      throwSignException() = true;
      try {
        calculate();
      } catch (SignException se) {
        //cerr << "signAlgebraic:  algT " << Parameter::algT << " failed." << endl;
        Parameter::algT = 0;
        Parameter::algR = 0;
        fail106 = true;
      }
    }

    if (!fail106) {
      //cerr << "algR " << Parameter::algR << endl;
      Parameter::algT *= Parameter::algR / 10;  // larger but conservative perturbation
      //cerr << "algT " << Parameter::algT << endl;
      while (true) // hopefully only once
        try {
          Parameter::algR = 1;
          calculate();
          break;
        } catch (SignException se) {
        //cerr << "algT shrunk to " << Parameter::algT << endl;
        Parameter::algT /= 10;
      }
      Parameter::algR = -1; // signal final get
      Parameter p;
      try {
        p = calculate();
      } catch (SignException se) {
        cerr << "signAlgebraic:  impossible" << endl;
        assert(0);
      }
      int sf = p.sign(false);
      Parameter::algT = 0;
      Parameter::algR = 0;
      if (sf) {
	throwSignException() = false;
	MInt::setPrecision(106u);
	return false;
      }
    }
    
    throwSignException() = false;
    MInt::setPrecision(212u);
    Parameter p = calculate();
    int sf = p.sign(false);
    if (sf != 0) {
      MInt::setPrecision(106u);
      return false;
    }
    ++nAmb;

#ifdef SET_ALGT
    Parameter::algT = 1e-62;
    // Parameter::algT = 1e-48; // vjm debug
    Parameter::algR = 1e56;
    // Parameter::algR = 1; // vjm debug
    throwSignException() = true;
    try {
      calculate();
    } catch (SignException se) {
      cerr << "signAlgebraic:  algT " << Parameter::algT << " failed." << endl;
      assert(false);
    }
    
    //cerr << "algR " << Parameter::algR << endl;
    Parameter::algT *= Parameter::algR / 10;  // larger but conservative perturbation
    //cerr << "algT " << Parameter::algT << endl;
    
    static double algTmin = 1e100;
    if (Parameter::algT < algTmin) {
      cout << "algTmin " << Parameter::algT << endl;
      algTmin = Parameter::algT;
    }
#endif

    // eps debug
    // Parameter::algT = 1.001 * 1e-62;
    Parameter::algI = algI;
    while (true) // hopefully only once
      try {
        calculate();
        break;
      } catch (SignException se) {
      //cerr << "algT shrunk to " << Parameter::algT << endl;
      Parameter::algT /= 10;
    }

    // eps debug 
    /*
    sf = p.sign(false);
    Parameter::algT = 0;
    Parameter::algR = 0;
    throwSignException() = false;
    MInt::setPrecision(106u);
    return sf == 0;
    */

    Parameter::algP = 1; // signal final get
    try {
      p = calculate();
    } catch (SignException se) {
      cerr << "signAlgebraic:  impossible" << endl;
      assert(0);
    }
    Parameter::algI = -1;
    Parameter::algP = 0;
    //cout << "interval " << p[0].intervalWidth() << " "
    //   << p[0].lb() << " " << p[0].ub() << endl;
    throwSignException() = false;
    sf = p.sign(false);
    if (sf == 0) {
      double w = p.intervalWidth();
      if (w > maxWidth)
	maxWidth = w;
    }
    MInt::setPrecision(106u);
    if (sf != 0)
      ++nAmbNZ;
    return sf == 0;
  }
};

template<class P>
class IthCoordinate : public Primitive {
  Object<P> *o;
  int i;
  Parameter calculate () { return o->get()[i]; }
public:
  IthCoordinate (Object<P> *o, int i) : o(o), i(i) {}
};

template<class P>
int Object<P>::sign (int i) {
  return IthCoordinate<P>(this, i);
}

void report ();

}

#define DeclareSign				\
  private:					\
    int sign ()

#define Primitive1(P, t1, v1)		\
  class P : public acp::Primitive {	\
    int sign();				\
    t1 v1;				\
  public:				\
    P (t1 v1) : v1(v1) {}		\
  };

#define Primitive2(P, t1, v1, t2, v2)		\
  class P : public acp::Primitive {		\
    int sign();					\
    t1 v1; t2 v2;				\
  public:					\
    P (t1 v1, t2 v2) : v1(v1), v2(v2) {}	\
  };

#define Primitive3(P, t1, v1, t2, v2, t3, v3)		\
  class P : public acp::Primitive {			\
    int sign();						\
    t1 v1; t2 v2; t3 v3;				\
  public:						\
    P (t1 v1, t2 v2, t3 v3) : v1(v1), v2(v2), v3(v3) {}	\
  };

#define Primitive4(P, t1, v1, t2, v2, t3, v3, t4, v4)			\
  class P : public acp::Primitive {					\
    int sign();								\
    t1 v1; t2 v2; t3 v3; t4 v4;						\
  public:								\
    P (t1 v1, t2 v2, t3 v3, t4 v4) : v1(v1), v2(v2), v3(v3), v4(v4) {}	\
  };

#define Primitive5(P, t1, v1, t2, v2, t3, v3, t4, v4, t5, v5)	\
  class P : public acp::Primitive {				\
    int sign();							\
    t1 v1; t2 v2; t3 v3; t4 v4; t5 v5;				\
  public:							\
    P (t1 v1, t2 v2, t3 v3, t4 v4, t5 v5)			\
      : v1(v1), v2(v2), v3(v3), v4(v4), v5(v5) {}		\
  };

#define Primitive1c(P, t1, v1)		\
  class P : public acp::Primitive {	\
    Parameter calculate();				\
    t1 v1;				\
  public:				\
    P (t1 v1) : v1(v1) {}		\
  };

#define Primitive2c(P, t1, v1, t2, v2)		\
  class P : public acp::Primitive {		\
    Parameter calculate();					\
    t1 v1; t2 v2;				\
  public:					\
    P (t1 v1, t2 v2) : v1(v1), v2(v2) {}	\
  };

#define Primitive3c(P, t1, v1, t2, v2, t3, v3)		\
  class P : public acp::Primitive {			\
    Parameter calculate();						\
    t1 v1; t2 v2; t3 v3;				\
  public:						\
    P (t1 v1, t2 v2, t3 v3) : v1(v1), v2(v2), v3(v3) {}	\
  };

#define Primitive4c(P, t1, v1, t2, v2, t3, v3, t4, v4)			\
  class P : public acp::Primitive {					\
    Parameter calculate();								\
    t1 v1; t2 v2; t3 v3; t4 v4;						\
  public:								\
    P (t1 v1, t2 v2, t3 v3, t4 v4) : v1(v1), v2(v2), v3(v3), v4(v4) {}	\
  };

#define Primitive5c(P, t1, v1, t2, v2, t3, v3, t4, v4, t5, v5)	\
  class P : public acp::Primitive {				\
    Parameter calculate();							\
    t1 v1; t2 v2; t3 v3; t4 v4; t5 v5;				\
  public:							\
    P (t1 v1, t2 v2, t3 v3, t4 v4, t5 v5)			\
      : v1(v1), v2(v2), v3(v3), v4(v4), v5(v5) {}		\
  };

#endif
