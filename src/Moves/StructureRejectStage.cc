#include "StructureRejectStage.h"


double StructureRejectStageClass::Sample(int &slice1,int &slice2,
				      Array<int,1> &changedParticles)
{
  //do nothing for now
}




void StructureRejectStageClass::Read (IOSectionClass &in)
{
  assert(in.ReadVar("MaxValue",MaxValue));
  StructureFactor.Read(in);
}

bool StructureRejectStageClass::Attempt (int &slice1, int &slice2, 
			       Array<int,1> &activeParticles, 
			       double &prevActionChange)
{
  bool toReturn=false;
  //  double tcounts=(double)StructureFactor.TotalCounts;
  //  double nslices=(double)PathData.Path.NumTimeSlices();
  double norm=PathData.Path.NumParticles();
  StructureFactor.Calculate();
  for (int counter=0;counter<StructureFactor.Sk.size();counter++){
    cerr<<StructureFactor.Sk(counter)/norm<<endl;
    if (StructureFactor.Sk(counter)/norm>4){
      toReturn=true;
    }
  }
  StructureFactor.Clear();
  return toReturn;
}
