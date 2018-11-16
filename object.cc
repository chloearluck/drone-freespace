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

namespace acp {
  pthread_mutex_t Primitive::algM = PTHREAD_MUTEX_INITIALIZER;
  pthread_key_t idkey;
  pthread_key_t mikey;
  pthread_key_t cpkey;
  pthread_key_t hsekey;
  PrecisionException precisionException;
  bool BaseObject::usePrecisionException = true;
  unsigned int BaseObject::deltaPrecision = 53u;
  unsigned int BaseObject::maxPrecision = 424u;
  unsigned int BaseObject::algebraicId = 0u;
  bool BaseObject::throwSignExceptions[128];
}

