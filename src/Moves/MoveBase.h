#ifndef MOVE_BASE_H
#define MOVE_BASE_H

#include "../PathDataClass.h"

/// This is the generic parent class for all moves, including "real moves"
/// which actually move particles and "pseudo moves", which just shift around
/// data, but don't move anything physical.
class MoveClass 
{
protected:
  /// The first time you write to an observable you have to do the
  /// write a little differently and you might need to write additional
  /// info like the description, etc.
  bool FirstTime;
  int TimesCalled, DumpFreq;
  ///You can add more IOVar pointers to inhereted classes if necessary
  VarClass *IOVar;
  IOSectionClass OutSection;  
public:
  /// This hold a reference to the Path Data
  PathDataClass &PathData;
  /// Call this in order to make a move.
  virtual void MakeMove();
  ///Moves have a name by which they can be referenced
  string Name;
  ///All moves ought to be able to read
  virtual void Read(IOSectionClass &input)=0;
  virtual double AcceptanceRatio() {return sqrt((double)-1.0);}
  virtual void WriteRatio();
  
  /// MoveClass constructor. Sets reference to the PathData object
  MoveClass(PathDataClass &myPathData, IOSectionClass outSection) : 
    PathData(myPathData), FirstTime(true), OutSection(outSection),
    TimesCalled(0), DumpFreq(1000)
    {Name="";}
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
  ////  /// A mapping from integers to particle ids
  //  Array<ParticleID,1> MyParticles;
 ///Our move class takes the number of particles to move at a
/// time. This is stored here. 
  int NumParticlesToMove;

 public:
  /// Stores the number of moves made and the number accepted
  int NumMoves, NumAccepted;
  /// This returns the Acceptance Ratio.
  inline double AcceptanceRatio() {return (double)(NumAccepted)/(double)NumMoves;}
  /// Call this to make a move
  virtual void MakeMove()=0;

  /// This array contains the int's of particles that you are 
  /// currently moving (i.e. NumParticlesToMove of them
  Array<int,1> ActiveParticles;
  /// When we choose particles we select the  particles (randomly)
  /// from the set of species enumerated in ActiveSpecies 
  void SetActiveSpecies(Array<int,1> ActSpecies);
  /// Function that sets the number of particles to move
  inline void SetNumParticlesToMove(int i)
  {
    NumParticlesToMove = i;
    MyParticleIndices.resize(i);
    ActiveParticles.resize(i);
  }
  /// Function that chooses the particles that you should move and
  /// places them in ActiveParticles; 
  void ChooseParticles();
  void ChooseParticlesOpen();
  inline int RandInt(int x);
  ParticleMoveClass(PathDataClass &myPathData, IOSectionClass outSection) : 
    MoveClass (myPathData, outSection)
  { 
    NumAccepted=0;
    NumMoves=0;
    /* Do nothing for now.*/  
  }
};











#endif
