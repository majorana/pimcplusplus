#ifndef SPLINE_POT_H
#define SPLINE_POT_H

#include "PotentialBase.h"
#include "../Splines/CubicSpline.h"

class SplinePot : public Potential
{
protected:
  // This is an optionally set potential that kicks in outside the
  // maximum value of the grid point.
  Potential *Vouter;
public:
  CubicSpline Spline;
  double V      (double r);
  double dVdr   (double r);
  double d2Vdr2 (double r);
  void Write(IOSectionClass &out);
  void Read (IOSectionClass &in);
  SplinePot() : Vouter(NULL)
  { 
    /* No nothing else for now */
  }
};

#endif
