#ifndef BISECTION_SPHERE_BLOCK_CLASS_H
#define BISECTION_SPHERE_BLOCK_CLASS_H


#include "MoveBase.h"
#include "../PathDataClass.h"
#include "PermuteStage.h"
#include "CoupledPermuteStage.h"
#include "BisectionSphereStage.h"
#include "../Observables/ObservableVar.h"

/// This is the bisection move class inherited from ParticleMoveClass
/// Explanation of how bisection moves work is in  
/// Path Integrals in the theory of condensed helium
/// by D.M. Ceperley  (Review of Modern Physics 1995) section V.H
class BisectionSphereBlockClass : public MultiStageClass
{
private:
  int NumLevels;
  int StepsPerBlock;
  bool HaveRefslice;
  int SpeciesNum;
  void ChooseTimeSlices();
  //  void WriteRatio();
  //  ObservableDouble AcceptanceRatioVar;
public:


  /// Number of levels the bisection move works on 
  void Read(IOSectionClass &in);
  
  /// Override base class MakeMove to do a block of moves
  void MakeMove();

  BisectionSphereBlockClass(PathDataClass &pathData, IOSectionClass &out) : 
    MultiStageClass(pathData, out)

  { 
    // do nothing for now
  }
};


#endif
