#ifndef COUPLINGMOVE_CLASS_H
#define COUPLINGMOVE_CLASS_H


#include "MoveBase.h"
#include "../PathDataClass.h"
#include "MultiStage.h"




/// This is the bisection move class inherited from ParticleMoveClass
/// Explanation of how bisection moves work is in  
/// Path Integrals in the theory of condensed helium
/// by D.M. Ceperley  (Review of Modern Physics 1995) section V.H
class CouplingMoveClass : public MultiStageClass
{
 private:
  //  EndType Open;
  int SpeciesNum;
 public:
  ///Number of levels the bisection move works on 
  int NumLevels;
  void Read(IOSectionClass &moveInput);
  void WriteRatio();
  void MakeMove();
  CouplingMoveClass(PathDataClass &myPathData, IOSectionClass iosection) :
    MultiStageClass(myPathData,iosection) 
  { 
    /* Do nothing for now. */ 
  }
};


#endif
