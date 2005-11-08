#ifndef FIXED_PHASE_ACTION_CLASS_H
#define FIXED_PHASE_ACTION_CLASS_H

#include "NodalActionClass.h"
#include <Common/PlaneWavePHDFT/PlaneWaves.h>
#include <Common/Splines/ComplexMultiTricubicSpline.h>
#include "../MirroredClass.h"

class FixedPhaseActionClass;

class FixedPhaseClass
{
private:
  SystemClass *System;
  PathDataClass &PathData;
  PathClass &Path;
  double kCut;
  ComplexMultiTricubicSpline BandSplines;
  Potential *PH;
  Vec3 kVec;

  Array<complex<double>,1> Workspace;
  Array<complex<double>,2> Matrix, Cofactors;
  Array<cVec3,2>   GradMat;
  Mirrored1DClass<double> UpGrad2, DownGrad2;

  Array<cVec3,1> Gradient;
  Array<Vec3,1> Rions;
  int NumUp, NumDown, NumIons, NumBands;
  // This stores the real space grid dimensions
  LinearGrid xGrid, yGrid, zGrid;
  complex<double> GradientDet   (int slice, int speciesNum);
  complex<double> GradientDetFD (int slice, int speciesNum);
  bool IonsHaveMoved();
  void UpdateBands();
  /// This updates GUp, gUp, and vUp, or GDown, gDown, and vDown,
  /// depending on the species.
  friend class FixedPhaseActionClass;
public:
  double CalcGrad2 (int slice, int species);
  void   CalcGrad2 (int slice, int species, Array<double,1> &grad2);
  int IonSpeciesNum, UpSpeciesNum, DownSpeciesNum;
  double Action (int slice1, int slice2,
		 const Array<int,1> &activeParticles, 
		 int level, int speciesNum);
  
  double d_dBeta(int slice1, int slice2, int level, int speciesNum);

  bool IsPositive (int slice, int speciesNum);
  complex<double> Det (int slice, int speciesNum);
  void Read (IOSectionClass &in);
  void ShiftData (int slices2Shift, int speciesNum);
  void AcceptCopy (int slice1, int slice2);
  void RejectCopy (int slice1, int slice2);
  void Init(int speciesNum);
  void Setk(Vec3 k);
  inline Vec3 Getk() { return kVec; }
  FixedPhaseClass (PathDataClass &pathData);
};


// This is a wrapper for the above class.
class FixedPhaseActionClass : public NodalActionClass
{
private:
  FixedPhaseClass &FixedPhaseA, &FixedPhaseB;
  int SpeciesNum;
public:
  double CalcGrad2 (int slice);
  void   CalcGrad2 (int slice, Array<double,1> &grad2);

  double SingleAction (int slice1, int slice2,
		       const Array<int,1> &activeParticles, int level);
  
  double d_dBeta(int slice1, int slice2, int level);

  bool IsPositive (int slice);
  complex<double> Det (int slice);
  //  Array<double,2> GetMatrix (int slice);
  /// Shifts internal time-slice dependent data
  void ShiftData (int slices2Shift);  
  void AcceptCopy (int slice1, int slice2);
  void RejectCopy (int slice1, int slice2);

  /// Initializes internals.
  void Init();

  /// Returns true
  bool IsGroundState();

  /// Returns FIXED_PHASE
  NodeType Type();

  /// Sets the k-vector or "twist angle"
  void Setk (Vec3 k);
  inline Vec3 Getk() { return FixedPhaseA.Getk(); }
  void WriteInfo (IOSectionClass &out);

  FixedPhaseActionClass (PathDataClass &pathData, 
			 FixedPhaseClass &FPA, FixedPhaseClass &FPB,
			 int speciesNum) :
    FixedPhaseA(FPA), FixedPhaseB(FPB), NodalActionClass (pathData), SpeciesNum(speciesNum)
  {
    
  }

};


#endif
