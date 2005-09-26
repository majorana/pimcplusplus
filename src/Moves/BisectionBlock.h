#ifndef BISECTION_BLOCK_CLASS_H
#define BISECTION_BLOCK_CLASS_H


#include "MoveBase.h"
#include "../PathDataClass.h"
#include "PermuteStageClass.h"
#include "CoupledPermuteStageClass.h"
#include "BisectionStageClass.h"
#include "../Observables/ObservableBase.h"

/// This is the bisection move class inherited from ParticleMoveClass
/// Explanation of how bisection moves work is in  
/// Path Integrals in the theory of condensed helium
/// by D.M. Ceperley  (Review of Modern Physics 1995) section V.H
class BisectionBlockClass : public MultiStageClass
{
private:
  int StepNum;
  int NumLevels;
  int StepsPerBlock;
  bool HaveRefslice;
  int SpeciesNum;
  void ChooseTimeSlices();
  StageClass* PermuteStage;
  void WriteRatio();
  //  ObservableDouble AcceptanceRatioVar;
public:


  /// Number of levels the bisection move works on 
  void Read(IOSectionClass &in);
  
  /// Override base class MakeMove to do a block of moves
  void MakeMove();

  BisectionBlockClass(PathDataClass &pathData, IOSectionClass &out) : 
    MultiStageClass(pathData, out),StepNum(0)

  { 
    // do nothing for now
  }
};


#endif
