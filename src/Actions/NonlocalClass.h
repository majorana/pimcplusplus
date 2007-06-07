#ifndef NONLOCAL_ACTION_H
#define NONLOCAL_ACTION_H

#include "ActionBase.h"
#include <Common/PH/NLPP.h>

class FixedPhaseClass;

class NonlocalClass : public ActionBaseClass
{
protected:
  NLPPClass *NLPP;
  FixedPhaseClass *FixedPhase;
  int IonSpecies, UpSpecies, DownSpecies;
  // Points for angular quadrature
  Array<Vec3,1> QuadPoints, ScaledPoints;
  // Weights for angular quadrature
  Array<double,1> QuadWeights, Legendre, DeltaV;
  Array<complex<double>,1> WFratios;
  void SetQuadratureRule (int nrule);
  void CheckQuadratureRule(int lexact);
  void NearestIon (int slice, int ptcl, Vec3 &ionpos, Vec3 &disp, double &dist);
  void ScaleQuadPoints (Vec3 ionpos, double dist);
  Array<int,1> Electrons;
public:
  double SingleAction (int slice1, int slice2,
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  void GradAction (int slice1, int slice2, const Array<int,1> &ptcls,
		   int level, Array<dVec,1> &gradVec);
  string GetName();
  void Setup (FixedPhaseClass *fixedPhase);
  void Read  (IOSectionClass &in);
  

  NonlocalClass (PathDataClass &pathData);
};



#endif
