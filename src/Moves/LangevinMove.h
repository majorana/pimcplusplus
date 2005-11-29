#ifndef LANGEVIN_MOVE_H
#define LANGEVIN_MOVE_H

#include "MoveBase.h"

class LangevinMoveClass : public MoveClass
{
protected:
  Array<dVec,1> V, R, OldF, Fsum;
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
  void LDStep();
  ObservableVecDouble2 Rvar, Vvar;
  Array<double,2> WriteArray;
  ObservableDouble TimeVar;
public:
  void MakeMove();
  void Read(IOSectionClass &in);
  LangevinMoveClass (PathDataClass &pathData, IOSectionClass &outSection) :
    MoveClass (pathData, outSection), MCSteps(0), LDSteps(0),
    Rvar("R", OutSection, pathData.Path.Communicator),
    Vvar("V", OutSection, pathData.Path.Communicator),
    TimeVar("Time", OutSection, pathData.Path.Communicator),
    Time(0.0)
  {
    // do nothing for now
  }
};


#endif
