#ifndef BISECTION_JOSEPHSON_STAGE_CLASS_H
#define BISECTION_JOSEPHSON_STAGE_CLASS_H

#include "MultiStage.h"
#include "../Observables/ObservableVar.h"

class BisectionJosephsonStageClass : public LocalStageClass

{
public:
  void WriteRatio();
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  void Accept();
  void Reject();
  int TotalLevels;
  BisectionJosephsonStageClass(PathDataClass &pathData, int level,
		      IOSectionClass outSection) : 
    LocalStageClass(pathData,outSection)
  { 
    //do nothing for now
    BisectionLevel = level;

  }
};

#endif
