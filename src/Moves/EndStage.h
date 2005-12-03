#ifndef END_STAGE_CLASS_H
#define END_STAGE_CLASS_H

#include "MultiStage.h"
#include "../Observables/ObservableVar.h"

typedef enum {HEAD,TAIL} EndType;

class EndStageClass : public LocalStageClass
{
private:
  EndType Open;
  int NumLevels;
  void ChooseTimeSlices(int &slice1,int &slice2);
  IOSectionClass OutSection;
  Array<int,1> AcceptRatio;
  int EndAttempts;
  ObservableVecDouble1 AcceptRatioVar;
public:
  bool OnlyX;
  void WriteRatio();
  void Accept();
  void Reject();
  void Read(IOSectionClass  &in);
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  EndStageClass(PathDataClass &pathData, int numLevels,
		IOSectionClass &outSection) : 
    LocalStageClass(pathData,outSection),
    OutSection(outSection),
    AcceptRatioVar("Acceptance Ratio",OutSection,PathData.Path.Communicator),
    NumLevels(numLevels)
  { 
    AcceptRatio.resize(2);
    AcceptRatio=0;
    EndAttempts=0;
    //do nothing for now

  }
};

#endif

