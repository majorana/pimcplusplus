#ifndef SPLINE_POT_H
#define SPLINE_POT_H

#include "PotentialBase.h"
#include "../Splines/CubicSpline.h"

class SplinePot : public Potential
{
public:
  CubicSpline Spline;
  double V      (double r);
  double dVdr   (double r);
  double d2Vdr2 (double r);
  void Write(IOSectionClass &out);
  void Read (IOSectionClass &in);
};

#endif
