#ifndef END_STAGE_CLASS_H
#define END_STAGE_CLASS_H

#include "MultiStage.h"

typedef enum {HEAD,TAIL} EndType;

class EndStageClass : public LocalStageClass
{
private:
  EndType Open;
  int NumLevels;
  void ChooseTimeSlices(int &slice1,int &slice2);
public:
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  EndStageClass(PathDataClass &pathData, int numLevels,
		IOSectionClass &outSection) : 
    LocalStageClass(pathData,outSection),
    NumLevels(numLevels)
  { 
    //do nothing for now

  }
};

#endif

