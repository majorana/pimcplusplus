#ifndef GROUND_STATE_NODAL_ACTION_CLASS_H
#define GROUND_STATE_NODAL_ACTION_CLASS_H

#include "NodalActionClass.h"
#include "../Common/PlaneWavePHDFT/PlaneWaves.h"
#include "../Common/Splines/MultiTricubicSpline3.h"

class GroundStateNodalActionClass : NodalActionClass
{
private:
  SystemClass *System;
  int IonSpeciesNum, UpSpeciesNum, DownSpeciesNum;
  double kCut;
  MultiTricubicSpline BandSplines;
  Potential &PH;

  Array<double,1> Workspace;
  Array<double,2> Matrix, Cofactors;
  Array<double,1> UpDists, DownDists;

  Array<Vec3,1> Gradient, Temp, Rions;
  int NumUp, NumDown, NumIons, NumBands;
  // This stores the real space grid dimensions
  LinearGrid xGrid, yGrid, zGrid;
  double Det         (int slice, int speciesNum);
  double GradientDet (int slice, int speciesNum);
  bool IonsHaveMoved();
  void UpdateBands();
  double SimpleDistance (int slice, int species);
  double LineSearchDistance (int slice, int species);
public:
  double Action (int slice1, int slice2,
		 const Array<int,1> &activeParticles, int level);
  
  double d_dBeta(int slice1, int slice2, int level);
  
  void Read (IOSectionClass &in);
  GroundStateNodalActionClass (PathDataClass &pathData,
			       Potential &ph) :
    NodalActionClass (pathData), PH(ph)
  {
    
  }
};


#endif
