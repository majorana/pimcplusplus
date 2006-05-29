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

#ifndef LANGEVIN_MOVE_H
#define LANGEVIN_MOVE_H

#include "MoveBase.h"
#include <deque>

typedef enum { VERLET, LANGEVIN } MDIntegratorType;

class LangevinMoveClass : public MoveClass
{
protected:
  MDIntegratorType Integrator;

  Array<dVec,1> V, Vold, R, Rp, OldFShort, OldFLong, 
    FShortSum, FLongSum, FShortTmp, FLongTmp;

  ///////////////////////////
  // Friction calculations //
  ///////////////////////////
  /// These stores a deque of the force evaluations during the
  /// accumlation stage of the present ion move.
  Array<dVec,1> FShort, FLong, FTmp;
  deque<Array<dVec,1> > FDeque;
  /// This stores the current friction coefficient matrix
  Array<double,2> A, CoVar;
  /// L stores the eigenvectors of A.  Ltran is its transpose
  Array<double,2> L, Ltrans;
  /// This stores the eigenvalues of A;
  Array<double,1> Lambda, ExpLambda, Expm1Lambda, RArray, RpArray, TempArray;
  /// FF stores <F_i F_j>.  This is used in computing the covariance
  /// matrix for the forces, which in turn is used to compute friction
  /// coefficents.  
  /// This stores the long-averaged means
  Array<double,1> Fmean, FallSum;
  /// Kappa is the autocorrelation time for the force calculation
  double kappa;
  /// This stores the number of means we have summed.
  int NumFs;
  void CalcFriction();
  void CalcCovariance();
  void CalcCovariance2();


  //////////////////////
  /// Misc functions ///
  //////////////////////

  /// The species we are doing Langevin dynamics on
  int LDSpecies;
  /// This holds the list of particles with which we're doing dynamics
  Array<int,1> Particles;
  /// The time step for the langevin dynamics integration
  double TimeStep;
  /// Number of MC steps to equilibrate before accumulating forces
  int NumEquilSteps;
  /// Number of MC steps to average forces before making an LD step
  int NumAccumSteps;
  /// The number of MC steps made in the current LD step.  Reset to 0
  /// after each LD step is made
  int MCSteps;
  /// The number of LD steps made so far
  int LDSteps;
  /// This is the friction coefficient to use
  double gamma;
  /// The mass of the LD species and its inverse
  double Mass, MassInv;
  /// The LD time
  double Time;
  /// A fudge-factor to increase/decrease the friction to get the
  /// ionic temperature right.
  double FrictionFactor;
  /// The standard deviation of the extra noise constant
  double ExtraNoiseSigma;
  void InitVelocities();
  void AccumForces();
  void VerletStep();
  void LangevinStep();

  /// This optionally writes the charge density to the file
  ObservableVecDouble3 RhoVar;
  /// This bool tells us whether or not to dump Rho
  bool DumpRho;
  /// And this store the actual density
  Array<double,3> Rho;
  ObservableVecDouble2 Rvar, Vvar, VOldVar, FVar, FShortVar, FLongVar,
    CoVarVar;
  ObservableVecDouble1 LambdaVar, BandEnergiesVar;
  Array<double,2> WriteArray;
  ObservableDouble TimeVar;
public:
  void MakeMove();
  void Read(IOSectionClass &in);
  LangevinMoveClass (PathDataClass &pathData, IOSectionClass &outSection) :
    MoveClass (pathData, outSection), MCSteps(0), LDSteps(0),
    Rvar      ("R",                 IOSection, pathData.Path.Communicator),
    Vvar      ("V",                 IOSection, pathData.Path.Communicator),
    VOldVar   ("Vold",              IOSection, pathData.Path.Communicator),
    FVar      ("F",                 IOSection, pathData.Path.Communicator),
    FShortVar ("FShort",            IOSection, pathData.Path.Communicator),
    FLongVar  ("FLong",             IOSection, pathData.Path.Communicator),
    LambdaVar ("Lambda",            IOSection, pathData.Path.Communicator),
    CoVarVar  ("CoVar",             IOSection, pathData.Path.Communicator),
    TimeVar   ("Time",              IOSection, pathData.Path.Communicator),
    BandEnergiesVar("BandEnergies", IOSection, pathData.Path.Communicator),
    RhoVar    ("Rho",               IOSection, pathData.Path.Communicator),
    Time(0.0), Integrator (VERLET), ExtraNoiseSigma(0.0), FrictionFactor(1.0),
    DumpRho(false)
  {
    // do nothing for now
  }
};


#endif
