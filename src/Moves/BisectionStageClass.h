#ifndef BISECTION_STAGE_CLASS_H
#define BISECTION_STAGE_CLASS_H

#include "MultiStage.h"


class BisectionStageClass : public StageClass
{
private:
  int NumImages;
  int Level;
public:
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  BisectionStageClass(PathDataClass &pathData, int level) : 
    StageClass(pathData), Level(level)
  { 
    //do nothing for now
  }
};

#endif
