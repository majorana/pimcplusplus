/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef FIXED_PHASE_ACTION_CLASS_H
#define FIXED_PHASE_ACTION_CLASS_H

#include "NodalActionClass.h"
#include <Common/PlaneWavePHDFT/PlaneWavesMPI.h>
#include <Common/Splines/ComplexMultiTricubicSpline.h>
#include "../MirroredClass.h"

class FixedPhaseActionClass;

typedef enum {UPDATE_NONE, UPDATE_ALL, UPDATE_ACTIVE} UpdateType;


class FixedPhaseClass
{
private:
  MirroredRefClass<MPISystemClass> System;
  PathDataClass &PathData;
  PathClass &Path;
  double kCut;
  MirroredClass<ComplexMultiTricubicSpline> BandSplines;
  Potential *V_elec_ion, *V_ion_ion;
  Vec3 kVec;

  Array<complex<double>,1> Workspace;
  Array<complex<double>,2> Matrix, Cofactors;
  Array<cVec3,2>   GradMat;
  Mirrored1DClass<double> UpGrad2, DownGrad2;
  Mirrored1DClass<double> UpAction, DownAction;
  Mirrored2DClass<Vec3> UpGrad, DownGrad;
  Mirrored1DClass<double> UpPhase, DownPhase;

  Mirrored3DClass<complex<double> > UpMatrixCache, DownMatrixCache;
  Mirrored3DClass<cVec3> UpGradMatCache, DownGradMatCache;

  Array<cVec3,1> Gradient;
  Array<complex<double>,1> OrbitalValues;
  Mirrored1DClass<Vec3> Rions;
  int NumUp, NumDown, NumIons, NumBands, NumFilled;
  Array<int,1> UpParticles, DownParticles;
  // This stores the real space grid dimensions
  LinearGrid xGrid, yGrid, zGrid;
//   complex<double> GradientDet   (int slice, int speciesNum);
  complex<double> GradientDetFD (int slice, int speciesNum);
  complex<double> GradientDet   (int slice, int speciesNum,
				 const Array<int,1> &activeParticles,
				 UpdateType update=UPDATE_ALL);
  bool IonsHaveMoved();
  /// Recomputes the electron orbitals for the new ion positions.
  void UpdateBands();
  /// Fills all the cached slater matrices with the current electron
  /// positions. 
  void UpdateCache();
  /// This updates GUp, gUp, and vUp, or GDown, gDown, and vDown,
  /// depending on the species.
  friend class FixedPhaseActionClass;
  bool UseMDExtrap;
  /// Include DFT in local density approximation
  bool UseLDA;
  // Note this action DOES NOT INCLUDE TAU FACTOR
  double CalcAction (Array<Vec3,1> &G1, Array<Vec3,1> &G2,
		     double phase1, double phase2, Array<Vec3,1> &dR);
  bool MonteCarloIons;
public:
  //  double CalcGrad2 (int slice, int species);
  double CalcGrad2 (int slice, int species, 
		    const Array<int,1> &activeParticles,
		    UpdateType update);
  // This returns the phase gradient.
  double PhaseGrad (int slice, int species, 
		    const Array<int,1> &activeParticles,
		    Array<Vec3,1> &gradPhase, UpdateType update);
//   void   CalcGrad2 (int slice, int species, Array<double,1> &grad2,
// 		    const Array<int,1> &activeParticles);
  int IonSpeciesNum, UpSpeciesNum, DownSpeciesNum;
  double Action (int slice1, int slice2,
		 const Array<int,1> &activeParticles, 
		 int level, int speciesNum);
  void CalcWFratios (int slice, int ptcl, const Array<Vec3,1> &pos, 
		     Array<complex<double>,1> &ratios);
  
  double d_dBeta(int slice1, int slice2, int level, int speciesNum);
  void CalcDensity     (Array<double,3> &rho);
  inline const Array<double,3>&  GetDensity()
  { return System().GetDensity(); }
  void CalcBandDensity (Array<double,4> &rho);

  bool IsPositive (int slice, int speciesNum);
  complex<double> Det (int slice, int speciesNum);
  void Read (IOSectionClass &in);
  void ShiftData (int slices2Shift, int speciesNum);
  void MoveJoin  (int oldJoinPos, int newJoinPos, int speciesNum);
  void AcceptCopy (int slice1, int slice2, int speciesNum);
  void RejectCopy (int slice1, int slice2, int speciesNum);
  void Init(int speciesNum);
  void Setk(Vec3 k);
  inline Vec3 Getk() { return kVec; }
  inline int GetNumBands() { return NumBands; }
  void GetBandEnergies(Array<double,1> &energies);
  void GetIonForces (Array<Vec3,1> &F);
  FixedPhaseClass (PathDataClass &pathData);
};


// This is a wrapper for the above class.
class FixedPhaseActionClass : public NodalActionClass
{
private:
  FixedPhaseClass &FixedPhaseA, &FixedPhaseB;
  int SpeciesNum;
public:
  //  double CalcGrad2 (int slice);
  double CalcGrad2 (int slice, const Array<int,1> &activeParticles,
		    UpdateType update);
//   void   CalcGrad2 (int slice, Array<double,1> &grad2,
// 		    const Array<int,1> &activeParticles);

  double SingleAction (int slice1, int slice2,
		       const Array<int,1> &activeParticles, int level);
  
  double d_dBeta(int slice1, int slice2, int level);

  bool IsPositive (int slice);
  complex<double> Det (int slice);
  //  Array<double,2> GetMatrix (int slice);
  /// Shifts internal time-slice dependent data
  void ShiftData (int slices2Shift);
  void MoveJoin (int oldJoinPos, int newJoinPos);
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
  inline int GetNumBands() { return FixedPhaseA.GetNumBands(); }
  void WriteInfo (IOSectionClass &out);
  /// Updates the bands if I'm  the ion species
  void Update();
  void CalcDensity (Array<double,3> &rho);
  const Array<double,3>& GetDensity ();
  void CalcBandDensity (Array<double,4> &rho);
  void GetBandEnergies(Array<double,1> &energies);
  void GetIonForces (Array<Vec3,1> &F);
  string GetName();

  FixedPhaseActionClass (PathDataClass &pathData, 
			 FixedPhaseClass &FPA, FixedPhaseClass &FPB,
			 int speciesNum) :
    FixedPhaseA(FPA), FixedPhaseB(FPB), NodalActionClass (pathData), 
    SpeciesNum(speciesNum)
  {
    
  }

};


#endif
