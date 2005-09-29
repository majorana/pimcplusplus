#ifndef GROUND_STATE_NODAL_ACTION_CLASS_H
#define GROUND_STATE_NODAL_ACTION_CLASS_H

#include "NodalActionClass.h"
#include <Common/PlaneWavePHDFT/PlaneWaves.h>
#include <Common/Splines/MultiTricubicSpline.h>
#include "../MirroredClass.h"

class GroundStateClass 
{
private:
  SystemClass *System;
  PathDataClass &PathData;
  PathClass &Path;
  double kCut;
  MultiTricubicSpline BandSplines;
  Potential *PH;

  Array<double,1> Workspace;
  Array<double,2> Matrix, Cofactors;
  Array<Vec3,2>   GradMat;
  Mirrored1DClass<double> UpDists, DownDists;

  Array<Vec3,1> Gradient, Temp, Rions;
  int NumUp, NumDown, NumIons, NumBands;
  // This stores the real space grid dimensions
  LinearGrid xGrid, yGrid, zGrid;
  double GradientDet   (int slice, int speciesNum);
  double GradientDetFD (int slice, int speciesNum);
  bool IonsHaveMoved();
  void UpdateBands();
  double SimpleDistance (int slice, int species);
  double LineSearchDistance (int slice, int species);
public:
  int IonSpeciesNum, UpSpeciesNum, DownSpeciesNum;
  double Action (int slice1, int slice2,
		 const Array<int,1> &activeParticles, 
		 int level, int speciesNum);
  
  double d_dBeta(int slice1, int slice2, int level, int speciesNum);

  bool IsPositive (int slice, int speciesNum);
  double Det (int slice, int speciesNum);
  Array<double,2> GetMatrix (int slice, int speciesNum);
  void Read (IOSectionClass &in);
  void ShiftData (int slices2Shift, int speciesNum);
  void AcceptCopy (int slice1, int slice2);
  void RejectCopy (int slice1, int slice2);
  void Init(int speciesNum);
  GroundStateClass (PathDataClass &pathData);

};


// This is a wrapper for the above class.
class GroundStateNodalActionClass : public NodalActionClass
{
private:
  GroundStateClass &GroundState;
  int SpeciesNum;
public:
  double Action (int slice1, int slice2,
		 const Array<int,1> &activeParticles, int level);
  
  double d_dBeta(int slice1, int slice2, int level);

  bool IsPositive (int slice);
  double Det (int slice);
  //  Array<double,2> GetMatrix (int slice);
  void ShiftData (int slices2Shift);  
  void AcceptCopy (int slice1, int slice2);
  void RejectCopy (int slice1, int slice2);
  void Init();
  bool IsGroundState();
  NodeType Type();
  void WriteInfo(IOSectionClass &out);

  GroundStateNodalActionClass (PathDataClass &pathData, GroundStateClass &GS,
			       int speciesNum) :
    GroundState(GS), NodalActionClass (pathData), SpeciesNum(speciesNum)
  {
    
  }

};


#endif
