#ifndef SCREENED_POT_H
#define SCREENED_POT_H

#include "PotentialBase.h"
#include "../Splines/CubicSpline.h"

class ScreenedPot : public Potential
{
public:
  double Charge;
  CubicSpline HXC;
  Potential *BarePot;

  bool IsPH();
  double CoreRadius();
  double A      (double r);
  double B      (double r);
  double dAdr   (double r);
  double d2Adr2 (double r);

  double V      (double r);
  double dVdr   (double r);
  double d2Vdr2 (double r);
  void   Write(IOSectionClass &out);
  void   Read (IOSectionClass &in);
};



#endif
