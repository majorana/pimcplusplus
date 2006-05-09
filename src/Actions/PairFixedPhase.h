#ifndef PAIR_FIXED_PHASE_CLASS_H
#define PAIR_FIXED_PHASE_CLASS_H

#include "NodalActionClass.h"
#include <Common/PlaneWavePHDFT/PlaneWavesMPI.h>
#include <Common/Splines/ComplexMultiTricubicSpline.h>
#include "../MirroredClass.h"

class PairFixedPhaseClass  : public ActionBaseClass
{
public:
  Array<complex<double>,2> DetMatrix;
  void BuildRefDet();
  int Species1Num;
  int Species2Num;
  int NumberAssymetry;
  PathDataClass &PathData;
  PathClass &Path;
  double SingleAction (int slice1, int slice2,
		 const Array<int,1> &activeParticles, 
		 int level);
  PairFixedPhaseClass (PathDataClass &pathData);
  double d_dBeta(int slice1, int slice2, int level);
  complex<double> Det (int slice,dVec kVec);
  Array<complex<double>,1> DegenerateRefSliceDeterminates;
  complex<double>  CalcDegenerateDet(int slice);

  void Read (IOSectionClass &in);
};

#endif
