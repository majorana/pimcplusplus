#ifndef GENERAL_POT_H
#define GENERAL_POT_H

#include "PotentialBase.h"
#include "../Splines/CubicSpline.h"

class GeneralPot : public Potential
{
protected:
  Grid *PotGrid;
  CubicSpline PotSpline;
  double Z;
public:
  double V      (double r);
  double dVdr   (double r);
  double d2Vdr2 (double r);

  void Write (IOSectionClass &out);
  void Read (IOSectionClass &in);
  GeneralPot();
  ~GeneralPot();
};

#endif
