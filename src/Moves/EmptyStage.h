#ifndef EMPTY_STAGE_CLASS_H
#define EMPTY_STAGE_CLASS_H

#include "MultiStage.h"
#include "../Observables/ObservableBase.h"


class EmptyStageClass : public LocalStageClass
{
private:
  IOSectionClass OutSection;
  Array<int,1> AcceptRatio;
  ObservableVecDouble1 AcceptRatioVar;
  int EndAttempts;
public:
  void WriteRatio();
  void Accept();
  void Reject();
  void Read(IOSectionClass  &in);
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  EmptyStageClass(PathDataClass &pathData, 
		IOSectionClass &outSection) : 
    LocalStageClass(pathData,outSection),
    OutSection(outSection),
    AcceptRatioVar("Acceptance Ratio",OutSection,PathData.Path.Communicator)
  { 
    EndAttempts=0;
    AcceptRatio.resize(2);
    AcceptRatio=0;
    //do nothing for now

  }
};

#endif

