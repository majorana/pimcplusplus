#ifndef SPLINE_POT_H
#define SPLINE_POT_H

#include "PotentialBase.h"
#include "../Splines/CubicSpline.h"

/// This class stores a tabulated potential and interpolates the data
/// with a cubic spline.  In case r is outside the tabulated grid, it
/// optionally calls Vouter.
class SplinePot : public Potential
{
protected:

public:
  /// This stores 
  CubicSpline Spline;
  /// This is an optionally set potential that kicks in outside the
  /// maximum value of the grid point.
  Potential *Vouter;

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
