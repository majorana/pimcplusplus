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
  Array<CubicSpline,1> DeltaV;
  // Points for angular quadrature
  vector<Vec3> QuadPoints;
  // Weights for angular quadrature
  vector<double> QuadWeights;
  void SetQuadratureRule (int nrule);
  void CheckQuadratureRule(int lexact);
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
