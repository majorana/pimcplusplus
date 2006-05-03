#include "SuperfluidFraction.h"
#include <Common/MPI/Communication.h>

void 
SuperfluidFractionClass::Read(IOSectionClass& IO)
{
  WindingNumberClass::Read(IO);
  assert(SpeciesList.size()==1);
}

void
SuperfluidFractionClass::WriteBlock()
{
  int species=SpeciesList(0);
  double beta=PathData.Path.tau*PathData.Path.TotalNumSlices;
  int numParticles=
    PathData.Path.Species(species).LastPtcl-PathData.Path.Species(species).FirstPtcl;
  double factor=(2*PathData.Path.Species(species).lambda*
		 beta*numParticles);
  
  CalcWN2();
  // Only processor 0 writes.
  if (PathData.Path.Communicator.MyProc()==0) {
    if (FirstTime) {
      FirstTime = false;
      WriteInfo();
      IOSection.WriteVar("Type",string("Vector"));
    }
    WN2Array *=factor;
    SFVar.Write(WN2Array);
  }
}

