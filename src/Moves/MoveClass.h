#ifndef MOVE_CLASS_H
#define MOVE_CLASS_H


#include "Common.h"
#include "PathDataClass.h"
#include "RandomPermClass.h"
#include "BisectionClass.h" 

/// This is the generic parent class for all moves, including "real moves"
/// which actually move particles and "pseudo moves", which just shift around
/// data, but don't move anything physical.
class MoveClass 
{
 public:
  /// This hold a reference to the Path Data
  PathDataClass &PathData;
  /// Call this in order to make a move.
  virtual void MakeMove()=0;
  ///Moves have a name by which they can be referenced
  string Name;
  ///All moves ought to be able to read
  virtual void Read(IOSectionClass &input)=0;
  virtual double AcceptanceRatio() {return sqrt((double)-1.0);}
  
  /// MoveClass constructor. Sets reference to the PathData object
  MoveClass(PathDataClass &myPathData) : PathData(myPathData)
    {Name="";}

};
  
///This is the loop move that allows you to loop


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
  ParticleMoveClass(PathDataClass &myPathData) : MoveClass (myPathData)
  { 
    NumAccepted=0;
    NumMoves=0;
    /* Do nothing for now.*/  
  }
};


class PermMove : public ParticleMoveClass
{
 public:
  RandomPermClass myPerm;
  BisectionClass Bisection;
  int NumLevels;
  void MakeMove()
  {
    // First, decide on the chunk of slices we're working on
    int slice1;
    int slice2;
    int numSlices=PathData.NumTimeSlices();
    
    if (!PathData.Path.OpenPaths){
      double xi=PathData.Path.Random.Local();
      slice2=(int)(xi*(double)(numSlices-(1<<NumLevels)))+(1<<NumLevels);
      slice1=slice2-(1<<NumLevels);
    }
    else {
      // First, decide on the chunk of slices we're working on
      do{
	double xi=PathData.Path.Random.Local();
	slice2=(int)(xi*(double)(numSlices-(1<<NumLevels)))+(1<<NumLevels);
	slice1=slice2-(1<<NumLevels);
	//      if (slice2>=PathData.Path.NumTimeSlices()-1){
	//	cerr<<"ERROR!"<<slice2<<endl;
	//	slice2--;
	//      }
      } while ((slice1<=(int)(PathData.Path.OpenLink) && (int)(PathData.Path.OpenLink)<=slice2));
    }
    //  if (slice2>=PathData.Path.NumTimeSlices()-1){
    //    cerr<<"ERROR! ERROR! ERROR!";
    //  }
    //  cerr<<slice1<<" "<<slice2<<" "<<PathData.Path.NumTimeSlices()<<endl;
    PathData.MoveJoin(slice2);
    myPerm.SpeciesNum=0;
    myPerm.Slice1=slice1;
    myPerm.Slice2=slice2;
    myPerm.GetActionPairs();

    for (int counter=0;counter<100;counter++){
      myPerm.GeneratePerm();
      //    firstIndex i;
      //    Array<int,1> currentParticles(PathData.Path.NumParticles());
      //    currentParticles=i;
      myPerm.Apply();
      bool acceptBisect;
      if (myPerm.CurrentParticles.size()!=0){
	acceptBisect = 
	  Bisection.Bisect(slice1, NumLevels, myPerm.CurrentParticles,0);
      }
      else 
	acceptBisect=false;
      if (acceptBisect){
	//	cerr<<"I've accepted!"<<endl;
	PathData.AcceptMove(slice1,slice2,myPerm.CurrentParticles);
      }
      else 
	PathData.RejectMove(slice1, slice2,myPerm.CurrentParticles);
    }
  }    
  void Read(IOSectionClass &in)
  {
    assert(in.ReadVar("name",Name));
    assert(in.ReadVar("NumLevels",NumLevels));
  }
  PermMove(PathDataClass &myPathData) : 
    ParticleMoveClass (myPathData) , myPerm(myPathData),
    Bisection(myPathData)
  {
    /* do nothing for now */
  }


};












#endif
