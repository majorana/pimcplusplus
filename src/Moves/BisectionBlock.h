#ifndef BISECTION_BLOCK_CLASS_H
#define BISECTION_BLOCK_CLASS_H


#include "MoveBase.h"
#include "../PathDataClass.h"
#include "PermuteStageClass.h"
#include "BisectionStageClass.h"


/// This is the bisection move class inherited from ParticleMoveClass
/// Explanation of how bisection moves work is in  
/// Path Integrals in the theory of condensed helium
/// by D.M. Ceperley  (Review of Modern Physics 1995) section V.H
class BisectionBlockClass : public MultiStageLocalClass
{
private:
  int NumLevels;
  int StepsPerBlock;
  bool IsFermion;
  int SpeciesNum;
  void WriteRatio()
  {
    //do nothing for now
  }

  void ChooseTimeSlices();
public:

  /// Number of levels the bisection move works on 
  void Read(IOSectionClass &in);
  
  /// Override base class MakeMove to do a block of moves
  void MakeMove();

  BisectionBlockClass(PathDataClass &pathData, IOSectionClass &out) : 
    MultiStageLocalClass(pathData, out)
  { 
    // do nothing for now
  }
};


#endif
