#ifndef BISECTION_STAGE_CLASS_H
#define BISECTION_STAGE_CLASS_H

#include "MultiStage.h"


class BisectionStageClass : public LocalStageClass
{
private:
  int NumImage;
  VarClass *IOVar;
  IOSectionClass OutSection;
  bool FirstTime;

  //  double LogSampleProb (int slice1, int slice2, 
  //			Array<int,1> particles);
//  double SamplePaths (int slice1, int slice2, Array<int,1> particles);
public:
  void WriteRatio();
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  void Accept();
  void Reject();
  BisectionStageClass(PathDataClass &pathData, int level,
		      IOSectionClass outSection) : 
    LocalStageClass(pathData), NumImage(1), OutSection(outSection)
  { 
    //do nothing for now
    BisectionLevel = level;
    FirstTime=false;
    
  }
};

#endif
