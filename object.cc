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

#include "object.h"

double getTime ()
{
  timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1e-6*tv.tv_usec;
}

namespace acp {
  unsigned int getIndent;
  pthread_mutex_t Primitive::algM = PTHREAD_MUTEX_INITIALIZER;
  pthread_key_t idkey;
  pthread_key_t mikey;
  pthread_key_t cpkey;
  pthread_key_t hsekey;
  PrecisionException precisionException;
  bool BaseObject::usePrecisionException = true;
  unsigned int BaseObject::deltaPrecision = 53u;
  unsigned int BaseObject::maxPrecision = 424u;
  bool BaseObject::throwSignExceptions[128];
  unsigned int BaseObject::nAlg = 0;
  unsigned int BaseObject::nAmb = 0;
  unsigned int BaseObject::nAmbNZ = 0;
  unsigned int BaseObject::nHomo = 0;
  double BaseObject::minHomo = 1e100;
  double BaseObject::maxWidth = 0.0;
  double BaseObject::tHom = 0.0;
  double Primitive::tMod = 0.0;

void report ()
{
  cerr << "nSign = " << Parameter::nSign << "; nAZ = " << Parameter::nAZ
       << "; nANZ = " << Parameter::nANZ << "; nMods = " << Mods::nMods
       << "; tMod = " << Primitive::tMod << endl << "maxBits = " << Mods::maxBits
       << "; avBits = " << Parameter::sumBits/max(1u, Parameter::nAZ) << endl;
  if (BaseObject::nAlg) 
    cerr << "nAlg = " << BaseObject::nAlg << "; nAmb = " << BaseObject::nAmb
	 << "; nAmbNZ = " << BaseObject::nAmbNZ << "; nHomo " << BaseObject::nHomo
	 << "; minHomo = " << BaseObject::minHomo << endl << "maxWidth = "
	 << BaseObject::maxWidth << "; tHom = " << BaseObject::tHom << endl;
}

}
