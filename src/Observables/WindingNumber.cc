#include "WindingNumber.h"
#include <Common/MPI/Communication.h>



///This currently tells you the winding number of all the species. It
///might make sense to fix this so it tells only the winding number of
///a specific species
void WindingNumberClass::Accumulate()
{
  TotalDisp=0.0;
  int numLinks=PathData.Path.NumTimeSlices()-1;
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    for (int slice=0;slice<numLinks;slice++) {
      NumSamples++;
      dVec disp;
      disp=PathData.Path.Velocity(slice,slice+1,ptcl);
      TotalDisp(ptcl) =TotalDisp(ptcl)+ disp;
    }
  }
  PathData.Path.Communicator.Sum(TotalDisp,TempDisp);
  if (PathData.Path.Communicator.MyProc()==0) {  
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
      for (int dim=0;dim<NDIM;dim++){
	TotalW2[dim] =TotalW2[dim]+TempDisp(ptcl)[dim]*TempDisp(ptcl)[dim];
      }
    }
  }
  //  for (dim=0;dim<NDIM;dim++){
  //    allPtclDisp[dim]=(allPtclDisp[dim]*allPtclDisp[dim]);
  //  }
  //  TotalW2 += allPtclDispl;
}

void WindingNumberClass::Read(IOSectionClass& in)
{
  ObservableClass::Read(in);
}
void WindingNumberClass::WriteBlock()
{
  // Only processor 0 writes.
  if (PathData.Path.Communicator.MyProc()==0) 
    if (FirstTime) {
      FirstTime = false;
      WriteInfo();
      IOSection.WriteVar("Type",string("Vector"));
    }
  Array<double,1> dummy(3);
  for (int dim=0;dim<NDIM;dim++)
    dummy(dim)=TotalW2[dim];
  WNVar.Write(dummy);
  NumSamples = 0;
}

