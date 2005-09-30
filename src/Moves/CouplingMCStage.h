#ifndef COUPLING_STAGE_MC_CLASS_H
#define COUPLING_STAGE_MC_CLASS_H

#include "MultiStage.h"
#include "../Observables/ObservableBase.h"


class CouplingMCStageClass : public LocalStageClass
{
private:
  void ChooseTimeSlices(int &slice1,int &slice2);
  IOSectionClass OutSection;
  Array<int,1> AcceptRatio;
  ObservableVecDouble1 AcceptRatioVar;
public:
  //  void WriteRatio();
  void Accept();
  void Reject();
  void Read(IOSectionClass  &in);
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  CouplingMCStageClass(PathDataClass &pathData, int numLevels,
		IOSectionClass &outSection) : 
    LocalStageClass(pathData,outSection),
    OutSection(outSection),
    AcceptRatioVar("Acceptance Ratio",OutSection,PathData.Path.Communicator)

  { 
    AcceptRatio.resize(2);
    AcceptRatio=0;

    //do nothing for now

  }
};

#endif

