#ifndef MOVE_CLASS_H
#define MOVE_CLASS_H

#include "Common.h";

///ParticleID=(species,particle number)
typedef TinyVec<int,2> ParticleID;


class MoveClass
{
private:
  /// Stores which species of particles this moves will work on
  Array<int,i> ActiveSpecies;
  ///Total number of particles in the active species
  int TotalParticles;
  ///Scratch Array holding a random subset of particles
  Array<int,1> MyParticleIndices; 
  ///A mapping from integers to particle ids
  Array<ParticleID,1> MyParticles;
protected:
  PathDataClass *PathData;
  int NumParticlesToMove;
  Array<ParticleID,1> ActiveParticles;
  void ChooseParticles ();


public:
  void SetActiveSpecies(Array<int,1> ActSpecies);
  inline void SetNumParticlesToMove(int i)
  {
    NumParticlesToMove = i;
    MyParticleIndices.resize(i);
  }
  virtual void makeMove()=0;
 
};
















#endif
