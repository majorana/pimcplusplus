#ifndef BISECTIONMOVE_CLASS_H
#define BISECTIONMOVE_CLASS_H

#include "Common.h"
#include "MoveClass.h"
#include "ActionClass.h"



class JoinMoveClass : public MoveClass
{
  int JoinLocation;
  void MoveJoin()
  {
    //    PathData.MoveJoin();
  }
  JoinMoveClass (PathDataClass &myPathData ) : MoveClass(myPathData)
  {/* Do nothing for now. */ }

};

/// This is the bisection move class inherited from ParticleMoveClass
/// Explanation of how bisection moves work is in  
/// Path Integrals in the theory of condensed helium
/// by D.M. Ceperley  (Review of Modern Physics 1995) section V.H
class BisectionMoveClass : public ParticleMoveClass
{
 public:
  ///This is the slice in which the bisection move starts.  It ends up
  ///going to StartTimeSlice+2^NumLevels
  int StartTimeSlice; 
  ///Number of levels the bisection move works on 
  int NumLevels;
  void Read(IOSectionClass *moveInput);
  /// Function to actually make a bisection move.
  void MakeMove();
  BisectionMoveClass(PathDataClass &myPathData ) : ParticleMoveClass(myPathData)
  { 
    /* Do nothing for now. */ 
  }
};


///This is a "psuedo-move" that is inherited from MoveClass.
///"Pseudo-moves do not increment the monte-carlo time slice. It
///shifts the data in memory by numTimeSlicesToShift. We do this to
///shift data between processors, etc.
class ShiftMoveClass : public MoveClass
{
 public:
  /// Contains the number of time slices to shift at a time. We have
///this change itself randomly from other objects.
  int numTimeSlicesToShift;
  /// Function to actually make a shift move. 
  void MakeMove();
  ShiftMoveClass (PathDataClass &myPathData) : MoveClass(myPathData)
  { /* Do nothing for now. */ }
};





       
  














#endif
