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
  void InitVelocities();
  void AccumForces();
  void VerletStep();
  void LangevinStep();
  ObservableVecDouble2 Rvar, Vvar, VOldVar, FVar, FShortVar, FLongVar;
  ObservableVecDouble1 LambdaVar;
  Array<double,2> WriteArray;
  ObservableDouble TimeVar;
public:
  void MakeMove();
  void Read(IOSectionClass &in);
  LangevinMoveClass (PathDataClass &pathData, IOSectionClass &outSection) :
    MoveClass (pathData, outSection), MCSteps(0), LDSteps(0),
    Rvar("R", IOSection, pathData.Path.Communicator),
    Vvar("V", IOSection, pathData.Path.Communicator),
    VOldVar("Vold", IOSection, pathData.Path.Communicator),
    FVar("F", IOSection, pathData.Path.Communicator),
    FShortVar("FShort", IOSection, pathData.Path.Communicator),
    FLongVar ("FLong",  IOSection, pathData.Path.Communicator),
    LambdaVar ("Lambda",  IOSection, pathData.Path.Communicator),
    TimeVar("Time", IOSection, pathData.Path.Communicator),
    Time(0.0), Integrator (VERLET)
  {
    // do nothing for now
  }
};


#endif
