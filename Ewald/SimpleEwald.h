#ifndef SIMPLE_EWALD_H
#define SIMPLE_EWALD_H

#include "EwaldBase.h"

class SimpleEwald : Ewald
{
public:
  void BreakUp (double cutoff);
  double Vshort (double r);
  double Vlong_k (double k);
};

#endif
