#ifndef BISECTION_STAGE_CLASS_H
#define BISECTION_STAGE_CLASS_H

#include "MultiStage.h"
#include "../Observables/ObservableBase.h"

class BisectionStageClass : public LocalStageClass
{
private:
  int NumImage;

public:
  void WriteRatio();
  ObservableVarDouble AcceptRatioVar;
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  void Accept();
  void Reject();
  BisectionStageClass(PathDataClass &pathData, int level,
		      IOSectionClass outSection) : 
    LocalStageClass(pathData,outSection),
    NumImage(1),
    AcceptRatioVar("AcceptanceRatio",OutSection,myPathData.Path.Communicator) 
  { 
    //do nothing for now
    BisectionLevel = level;
  }
};

#endif
