#ifndef OPEN_STAGE_CLASS_H
#define OPEN_STAGE_CLASS_H

#include "MultiStage.h"
#include "../Observables/ObservableBase.h"
#include "PermuteStage.h"

class OpenStageClass : public PermuteStageClass
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
  bool Attempt (int &slice1, int &slice2, 
		Array<int,1> &activeParticles, double &prevActionChange);
  dVec RandomBoxLocation();
  //  int NumLevels;
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  OpenStageClass(PathDataClass &pathData, int speciesNum,
		 int numLevels,
		 IOSectionClass &outSection) : 
    //    LocalStageClass(pathData,outSection),
    PermuteStageClass(pathData, speciesNum, numLevels,outSection),
    OutSection(outSection),
    AcceptRatioVar("Acceptance Ratio",OutSection,PathData.Path.Communicator)
  { 
    EndAttempts=0;
    AcceptRatio.resize(2);
    AcceptRatio=0;
    ///    BisectionLevel=level;
    //do nothing for now

  }
};

#endif

