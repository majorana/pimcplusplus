#ifndef BISECTION_STAGE_CLASS_H
#define BISECTION_STAGE_CLASS_H

#include "MultiStage.h"


class BisectionStageClass : public LocalStageClass
{
private:
  int NumImage;
  //  double LogSampleProb (int slice1, int slice2, 
  //			Array<int,1> particles);
//  double SamplePaths (int slice1, int slice2, Array<int,1> particles);
public:
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  BisectionStageClass(PathDataClass &pathData, int level) : 
    LocalStageClass(pathData), NumImage(1)
  { 
    //do nothing for now
    BisectionLevel = level;
  }
};

#endif
