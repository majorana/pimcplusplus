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
  double RotationalKinetic(int startSlice, int endSlice, const Array<int,1> &activeParticles,int level);
  double RotationalEnergy(int startSlice, int endSlice, int level);
  double GetAngles(dVec disp);
  dVec COMVelocity (int slice1,int slice2,int ptcl);
  TIP5PWaterClass (PathDataClass &pathData);
};

const double O_H_moment_arm = 0.9572;
const double lambda_p = 0.047848;

#endif
