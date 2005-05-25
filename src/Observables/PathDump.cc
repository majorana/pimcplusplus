#include "PathDump.h"


void PathDumpClass::Accumulate()
{
  TimesCalled++;
  if (TimesCalled%DumpFreq == (DumpFreq-1)){
    WriteBlock();
  }
}

void PathDumpClass::Read(IOSectionClass &in)
{
  ObservableClass::Read(in);
  assert (in.ReadVar ("dumpFreq", DumpFreq));
}

void PathDumpClass::WriteBlock()
{
//   PathClass &Path = PathData.Path;
//   int first0, last0;
//   Path.SliceRange (0, first0, last0);
//   // First, let's get the reference slice to processor 0's first slice
//   while (Path.GetRefSlice() >= (last0-1)) {
//     PathData.MoveJoin((last0-1));
//     PathData.ShiftData (-(last0-1));
//     PathData.Join = 0;
//   }
//   PathData.MoveJoin(Path.GetRefSlice());
//   PathData.ShiftData (-Path.GetRefSlice());
//   assert (Path.GetRefSlice() == 0);
//   PathData.Join = 0;
//   // Now the reference slice should be at proc0's slice 0

//   //Move the join to the end so we don't have to worry about permutations
//   PathData.MoveJoin(PathData.NumTimeSlices()-1);
  PathClass &Path = PathData.Path;
  int refSave = PathData.Path.GetRefSlice();
  PathData.MoveRefSlice(0);

  if (PathData.Path.OpenPaths){
    Array<double,1> tailLoc(NDIM);
    OpenLinkVar.Write((int)PathData.Path.OpenLink);
    for (int dim=0;dim<NDIM;dim++){
      tailLoc(dim)=PathData.Path(PathData.Path.OpenLink,
				 PathData.Path.NumParticles())[dim];
    }
    TailLocVar.Write(tailLoc);
    RefLinkVar.Write(PathData.Path.RefSlice);
    OpenLinkPtclVar.Write(PathData.Path.OpenPtcl);
  }

  int start, end, numProcs, myProc;
  numProcs = Path.Communicator.NumProcs();
  myProc   = Path.Communicator.MyProc();

  Array<double,3> pathArray;
  int totalSlices = Path.TotalNumSlices;

  Path.SliceRange(numProcs-1, start,end);
  int maxShift = end-start;
  int slicesLeft = totalSlices;
  int offset = 0;
  int numPtcls = PathData.NumParticles();

  if (myProc == 0)
    pathArray.resize(numPtcls,totalSlices,NDIM);
  while (slicesLeft > maxShift) {
    // First copy
    PathData.MoveJoin(maxShift);
    if (myProc == 0)
      for (int i=0; i<maxShift; i++)
	for (int ptcl=0; ptcl < numPtcls; ptcl++)
	  for (int dim=0; dim<NDIM; dim++)
	    pathArray(ptcl, i+offset, dim) = Path(i,ptcl)[dim];
    // Now shift
    PathData.ShiftData(-maxShift);
    PathData.Join = 0;
    slicesLeft -= maxShift;
    offset += maxShift;
  }
  // Move join out of the way
  PathData.MoveJoin (Path.NumTimeSlices()-1);
  // Copy last part
  if (myProc == 0)
    for (int i=0; i<slicesLeft; i++)
      for (int ptcl=0; ptcl < numPtcls; ptcl++)
	for (int dim=0; dim<NDIM; dim++)
	  pathArray(ptcl, i+offset, dim) = Path(i,ptcl)[dim];
  
  // Reset path to original position
  PathData.MoveRefSlice (refSave);

  Array<int,1> permVec(numPtcls);
  Path.TotalPermutation(permVec);

  // Now write
  PathVar.Write (pathArray);
  PermVar.Write (permVec);


  if (FirstTime && (myProc == 0)) {
    WriteInfo();
    Array<string,1> speciesNames(numPtcls);
    for (int spIndex=0;spIndex<PathData.NumSpecies();spIndex++){
      SpeciesClass &species = PathData.Path.Species(spIndex);
      for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++)
	speciesNames(ptcl)=species.Name;
    }
    IOSection.WriteVar("SpeciesNames", speciesNames);
    IOSection.WriteVar("Type","Path");
  }
  FirstTime = false;
  IOSection.FlushFile();
//   int numPtcls = PathData.NumParticles();
//   int numTimeSlices = PathData.NumTimeSlices();
//   if (FirstTime){
//     FirstTime=false;
//     WriteInfo();
//     Array<string,1> speciesNames(numPtcls);
//     for (int speciesIndex=0;speciesIndex<PathData.NumSpecies();speciesIndex++){
//       SpeciesClass &species = PathData.Path.Species(speciesIndex);
//       for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++)
//       	speciesNames(ptcl)=species.Name;
//     }
//     IOSection.WriteVar("SpeciesNames", speciesNames);
//     IOSection.WriteVar("Type","Path");
//     Array<double,4> pathArray(1,numPtcls,numTimeSlices,NDIM);
//     for (int ptcl=0;ptcl<numPtcls;ptcl++){
//       for (int slice=0;slice<numTimeSlices;slice++){
// 	for (int dim=0;dim<NDIM;dim++){
// 	  pathArray(0,ptcl,slice,dim)=PathData(slice,ptcl)[dim];
// 	}
//       }
//     }

//     // Write the first path here
//     IOSection.WriteVar("Path",pathArray);
//     // Now get the pointer to it
//     IOVar = IOSection.GetVarPtr("Path");
//   }
//   else {
//     // Append the new path here.
//     Array<double,3> pathArray(numPtcls,numTimeSlices,NDIM);
//     for (int ptcl=0;ptcl<numPtcls;ptcl++)
//       for (int slice=0;slice<numTimeSlices;slice++)
// 	for (int dim=0;dim<NDIM;dim++)
// 	  pathArray(ptcl,slice,dim)=PathData(slice,ptcl)[dim];
//     IOVar->Append(pathArray);
//     IOSection.FlushFile();
//   }
}

