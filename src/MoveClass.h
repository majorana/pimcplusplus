#ifndef MOVE_CLASS_H
#define MOVE_CLASS_H


#include "Common.h"
#include "PathDataClass.h"


class MoveClass 
{
 public:
  PathDataClass *PathData;
  virtual void makeMove()=0;
};
  

class ParticleMoveClass : public MoveClass

{
 public:
  /// Stores which species of particles this moves will work on
  Array<int,1> ActiveSpecies;
  ///Total number of particles in the active species
  int TotalParticles;
  ///Scratch Array holding a random subset of particles
  Array<int,1> MyParticleIndices; 
  ///A mapping from integers to particle ids
  Array<ParticleID,1> MyParticles;

 public:
  virtual void makeMove()=0;
  int NumParticlesToMove;
  Array<ParticleID,1> ActiveParticles;
  void SetActiveSpecies(Array<int,1> ActSpecies);
  inline void SetNumParticlesToMove(int i)
  {
    NumParticlesToMove = i;
    MyParticleIndices.resize(i);
  }

  void ChooseParticles();
 
};
















#endif
