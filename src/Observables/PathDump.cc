#include "PathDump.h"


void PathDumpClass::Accumulate()
{
  TimesCalled++;
  if (TimesCalled % 5000==0){
    WriteBlock();
  }
}

void PathDumpClass::Read(IOSectionClass &in)
{
  ObservableClass::Read(in);
}

void PathDumpClass::WriteBlock()
{
  //Move the join to the end so we don't have to worry about permutations
  PathData.MoveJoin(PathData.NumTimeSlices()-1);

  int numPtcls = PathData.NumParticles();
  int numTimeSlices = PathData.NumTimeSlices();
  if (FirstTime){
    FirstTime=false;
    WriteInfo();
    Array<string,1> speciesNames(numPtcls);
    for (int speciesIndex=0;speciesIndex<PathData.NumSpecies();speciesIndex++){
      SpeciesClass &species = PathData.Path.Species(speciesIndex);
      for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++)
      	speciesNames(ptcl)=species.Name;
    }
    IOSection.WriteVar("SpeciesNames", speciesNames);
    IOSection.WriteVar("Type","Path");
    Array<double,4> pathArray(1,numPtcls,numTimeSlices,NDIM);
    for (int ptcl=0;ptcl<numPtcls;ptcl++){
      for (int slice=0;slice<numTimeSlices;slice++){
	for (int dim=0;dim<NDIM;dim++){
	  pathArray(0,ptcl,slice,dim)=PathData(slice,ptcl)[dim];
	}
      }
    }

    // Write the first path here
    IOSection.WriteVar("Path",pathArray);
    // Now get the pointer to it
    IOVar = IOSection.GetVarPtr("Path");
  }
  else {
    // Append the new path here.
    Array<double,3> pathArray(numPtcls,numTimeSlices,NDIM);
    for (int ptcl=0;ptcl<numPtcls;ptcl++)
      for (int slice=0;slice<numTimeSlices;slice++)
	for (int dim=0;dim<NDIM;dim++)
	  pathArray(ptcl,slice,dim)=PathData(slice,ptcl)[dim];
    IOVar->Append(pathArray);
    IOSection.FlushFile();
  }
}

