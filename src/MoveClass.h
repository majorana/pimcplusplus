#ifndef MOVE_CLASS_H
#define MOVE_CLASS_H


#include "Common.h"
#include "PathDataClass.h"

/// This is the generic parent class for all moves, including "real moves"
/// which actually move particles and "pseudo moves", which just shift around
/// data, but don't move anything physical.
class MoveClass 
{
 public:
  PathDataClass &PathData;
  virtual void makeMove()=0;
  MoveClass(PathDataClass &myPathData) : PathData(myPathData);
};
  
/// This is a specialization of MoveClass which actually physically moves
/// particles.
class ParticleMoveClass : public MoveClass
{
protected:
  /// Stores which species of particles this moves will work on
  Array<int,1> ActiveSpecies;
  /// Total number of particles in the active species
  int TotalParticles;
  /// Scratch Array holding a random subset of particles
  Array<int,1> MyParticleIndices; 
  /// A mapping from integers to particle ids
  Array<ParticleID,1> MyParticles;

 public:
  int NumMoves, NumAccepted;
  double AcceptanceRatio();
  virtual void makeMove()=0;
  int NumParticlesToMove;
  Array<ParticleID,1> ActiveParticles;
  void SetActiveSpecies(Array<int,1> ActSpecies);
  inline void SetNumParticlesToMove(int i)
  {
    NumParticlesToMove = i;
    MyParticleIndices.resize(i);
    ActiveParticles.resize(i);
  }

  void ChooseParticles();
 
};
















#endif
