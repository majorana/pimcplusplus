#include "Phik.h"



void PhiKClass::Accumulate()
{
  //Move the join to the end so we don't have to worry about permutations
  PathData.MoveJoin(PathData.NumTimeSlices()-1);
  
  NumSamples++;

  int slice1 = 0;
  int slice2 = PathData.Path.NumTimeSlices();
  int numTimeSlices=PathData.Path.NumTimeSlices();
  PathClass &Path = PathData.Path;
  for (int tau=0;tau<numTimeSlices;tau++)
    for (int slice=slice1;slice<slice2;slice++)
      PhiK(tau) += Path(slice,0)[0]*Path((slice+tau) % numTimeSlices,0)[0];
}


void PhiKClass::WriteBlock()
{
  int nslices=PathData.Path.TotalNumSlices;
  double norm = 1.0/((double)NumSamples*(double)nslices);
  PhiK=PhiK*norm;
  PhiKVar.Write(PhiK);  
  PhiK=0;
  NumSamples = 0;
}

void PhiKClass::Read(IOSectionClass &in)
{  
  ObservableClass::Read(in);
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
}

