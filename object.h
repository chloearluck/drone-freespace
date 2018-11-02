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

  static bool & throwSignException () {
    void *v = pthread_getspecific(hsekey);
    return v ? *(bool *) v : throwSignExceptions[0];
  }

  BaseObject () {}
  virtual ~BaseObject () {}
  
  static void addThread (unsigned int i) {
    assert(i < 128);
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
class Object :  public BaseObject {
  friend class Parameter;
  
  P p;

  virtual P calculate () {
    if (!input())
      cerr << "Object<P> missing calculate method" << endl;
    assert(0);
  }

  bool checkAccuracy (double acc) {
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

public:
  Object () {}
  
  Object (const P& q) : p(q) {
    for (int i = 0; i < p.size(); i++) {
      if (p[i].high())
        p[i].decreasePrecision();
      if (p[i].lb() != p[i].ub()) {
        std::cerr << "Input object has non-trivial interval." << std::endl;
        assert(0);
      }
    }
  }

  ~Object () { 
    for (int i = 0; i < p.size(); i++)
      if (p[i].high())
        delete p[i].u.m;
  }
  
  int precision () const { return p.size() == 0 ? 0 : p[0].precision(); }
  bool uninitialized () const { return p.size() == 0 ? 0 : p[0].uninitialized(); }

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
    if (precision() == 53u && MInt::getPrecision() > 53u)
      for (int i = 0; i < q.size(); i++)
        q[i].increasePrecision();
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
    if (MInt::getPrecision() == 53u && p.size() > 0 && p[0].low() && !p[0].uninitialized())
      return p;        

    unsigned int precObject = precision();
    unsigned int precNeeded = MInt::getPrecision();

    if (precObject >= precNeeded && 
        (precObject == 53u || p[0].hasTheRightPrimes())) {
      if (precObject == 53u || precNeeded > 53u)
        return p;
      P q = p;
      for (int i = 0; i < q.size(); i++)
        q[i].decreasePrecision();
      return q;
    }

    if (input()) {
      P q = p;
      for (int i = 0; i < p.size(); i++)
        q[i].increasePrecision();
      return q;
    }

    if (precNeeded == 53u)
      return p = calculate();

    P q = calculate();
    for (int i = 0; i < q.size(); i++)
      if (q[i].zeroMod())
        q[i] = Parameter::constant(0.0);

    if (precObject > 53u)
      for (int i = 0; i < p.size(); i++)
        delete p[i].u.m;

    p = q;

    for (int i = 0; i < p.size(); i++)
      p[i].u.m = p[i].u.m->clone();

    return p;
  }

public:
  P getApprox (double acc = 1e-17) {
    try {
      get();
      assert(precision() > 0); // DEBUG
      if (checkAccuracy(acc))
        return get();
    }
    catch (SignException se) {}
    MInt::setPrecision(212u);
    while (MInt::getPrecision() <= maxPrecision) {
      try {
        get();
        assert(precision() > 0); // DEBUG
        if (checkAccuracy(acc)) {	
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
    assert(precision() > 0); // DEBUG
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
  virtual int sign () = 0;
public:
  operator int () {
    if (MInt::getPrecision() > 53u || throwSignException())
      return sign();
    try {
      return sign();
    }
    catch (SignException se) {}
    setPrecision(106u);
    while (MInt::getPrecision() <= maxPrecision) {
      try {
        int s = sign();
        MInt::setPrecision(53u);
        return s;
      }
      catch (unsigned int p) {
        ModInt::changePrime(p);
      }
      catch (SignException se) {
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
};

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

#endif
