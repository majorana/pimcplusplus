#ifndef TIP5PWATER_CLASS_H
#define TIP5PWATER_CLASS_H

#include "ActionBase.h"

/// The KineticClass calculates the kinetic part of the action.  This
/// is the "spring term", of the form
/// \f$ K \equiv \left(4\pi\lambda\tau\right)^{-\frac{ND}{2}} 
/// \exp\left[-\frac{(R-R')^2}{4\lambda\tau}\right] \f$
class TIP5PWaterClass : public ActionBaseClass
{
public:
  void Read (IOSectionClass &in);
  double Action (int slice1, int slice2, 
		 const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  double OOSeparation (int slice,int ptcl1,int ptcl2);
  TIP5PWaterClass (PathDataClass &pathData);
};


#endif
