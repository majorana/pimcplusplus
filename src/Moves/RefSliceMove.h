#ifndef REF_SLICE_MOVE_H
#define REF_SLICE_MOVE_H

#include "MoveBase.h"
#include "../PathDataClass.h"
#include "PermuteStageClass.h"
#include "BisectionStageClass.h"

/// This is the bisection move class inherited from ParticleMoveClass
/// Explanation of how bisection moves work is in  
/// Path Integrals in the theory of condensed helium
/// by D.M. Ceperley  (Review of Modern Physics 1995) section V.H
class RefSliceMoveClass : public MultiStageClass
{
private:
  /// Number of bisection stage levels
  int NumLevels;

  /// Holds the current master processor
  int MasterProc;

  /// The species this particular instance is working on
  int SpeciesNum;

  void WriteRatio()
  {
    //do nothing for now
  }

  /// This function is called if I have the reference slice
  void MakeMoveMaster();

  /// This move is called if I don't have the reference slice
  void MakeMoveSlave();

public:
  /// Read in the parameters this class needs from the input file.
  void Read(IOSectionClass &in);
  

  /// Override base class MakeMove to do a block of moves
  void MakeMove();

  RefSliceMoveClass(PathDataClass &pathData, IOSectionClass &out) : 
    MultiStageClass(pathData, out)
  { 
    // do nothing for now
  }
};



#endif
