#ifndef BISECTIONBLOCK_CLASS_H
#define BISECTIONBLOCK_CLASS_H


#include "MoveBase.h"
#include "MultiStageClass.h"
#include "../PathDataClass.h"
#include "BisectionClass.h"

/// This is the bisection move class inherited from ParticleMoveClass
/// Explanation of how bisection moves work is in  
/// Path Integrals in the theory of condensed helium
/// by D.M. Ceperley  (Review of Modern Physics 1995) section V.H
class BisectionBlockClass : public MultiStageLocalClass
{
 public:
  ///Number of levels the bisection move works on 
  int NumLevels;
  void Read(IOSectionClass &in);

  BisectionBlockClass(PathDataClass &pathData, IOSectionClass &out) : 
    MultiStageLocalClass(pathData, outSection), 
  { 
    ///Defaults to the 0'th time slice but shouldn't matter because it
    ///chooses a different time slice each time.
    StartTimeSlice=0;
  }
};


#endif
