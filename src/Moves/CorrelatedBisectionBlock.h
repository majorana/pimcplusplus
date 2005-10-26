#ifndef CORRELATED_BISECTION_BLOCK_CLASS_H
#define CORRELATED_BISECTION_BLOCK_CLASS_H


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
class CorrelatedBisectionBlockClass : public MultiStageClass
{

private:
  int StepNum;
  int NumLevels;
  int StepsPerBlock;
  bool HaveRefslice;
  int SpeciesNum;
  void ChooseTimeSlices();
  PermuteStageClass* PermuteStage;
  Array<BisectionStageClass*,1> BisectionStages;
  list<ActionBaseClass*> LocalActions;
  list<ActionBaseClass*> CommonActions;
  
  FILE *EAout, *EBout, *SAout, *SBout;
  ObservableDouble wAEAvar, wBEBvar, wASAvar, wBSBvar, wAvar, wBvar;
  //  ObservableDouble AcceptanceRatioVar;
public:
  /// Number of levels the bisection move works on 
  void Read(IOSectionClass &in);

  /// Override base class MakeMove to do a block of moves
  void MakeMove();

  CorrelatedBisectionBlockClass(PathDataClass &pathData, IOSectionClass &out) : 
    MultiStageClass(pathData, out),StepNum(0), 
    wAEAvar("wAEA", OutSection, pathData.Path.Communicator),
    wBEBvar("wBEB", OutSection, pathData.Path.Communicator),
    wASAvar("wASA", OutSection, pathData.Path.Communicator),
    wBSBvar("wBSB", OutSection, pathData.Path.Communicator),
    wAvar  ("wA",   OutSection, pathData.Path.Communicator),
    wBvar  ("wB",   OutSection, pathData.Path.Communicator)
  { 
    EAout = fopen ("EA.dat", "w");
    SAout = fopen ("SA.dat", "w");
    EBout = fopen ("EB.dat", "w");
    SBout = fopen ("SB.dat", "w");
    // do nothing for now
  }
};


#endif
