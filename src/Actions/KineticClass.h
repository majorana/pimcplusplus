#ifndef KINETIC_CLASS_H
#define KINETIC_CLASS_H

#include "ActionBase.h"

/// The KineticClass calculates the kinetic part of the action.  This
/// is the "spring term", of the form
/// \f$ K \equiv \left(4\pi\lambda\tau\right)^{-\frac{ND}{2}} 
/// \exp\left[-\frac{(R-R')^2}{4\lambda\tau}\right] \f$
class KineticClass : public ActionBaseClass
{
  int NumImages;
public:
  void Read (IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  inline void SetNumImages (int num) { NumImages = num; }
  KineticClass (PathDataClass &pathData);
};


#endif
