#ifndef STRUCTURE_REJECT_STAGE_CLASS_H
#define STRUCTURE_REJECT_STAGE_CLASS_H

#include "MultiStage.h"
#include "../Observables/StructureFactor.h"

class StructureRejectStageClass : public LocalStageClass
{
private:
  int Species1;
  int Species2;
  double MaxValue;
  StructureFactorClass StructureFactor;
public:
  void Read(IOSectionClass &in);
  double Sample (int &slice1,int &slice2,
		 Array<int,1> &changedParticles); 
  bool Attempt(int &slice1, int &slice2,
	       Array<int,1> &activeParticles,
	       double &prevActionChange);

  StructureRejectStageClass (PathDataClass &pathData, IOSectionClass &in,
			     IOSectionClass &out) 
    : LocalStageClass(pathData,out),
    StructureFactor(pathData,in)
  {
    // do nothing for now
  }

};

#endif
