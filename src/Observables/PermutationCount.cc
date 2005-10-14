#include "PermutationCount.h"

///////////////////////////////////////////////////////

///////////////////////////////////////////////////////




void PermutationCountClass::WriteBlock()
{
  //  Array<double,1> CorSum(Correlation.size());
  //  Path.Communicator.Sum(Correlation, CorSum);

  CycleCount=CycleCount/TotalCounts;
  CycleCountVar.Write(CycleCount);
  CycleCountVar.Flush();
  CycleCount=0;
  TotalCounts=0;

}


void PermutationCountClass::Read(IOSectionClass &in)
{  

  ObservableClass::Read(in);
  assert(in.ReadVar("freq",Freq));
  assert(in.ReadVar("dumpFreq",DumpFreq));
  string speciesName;
  Species=-1;
//   assert(in.ReadVar("Species1",speciesName));
//   for (int spec=0;spec<PathData.NumSpecies();spec++){ //???what is Species ?
//     if (PathData.Species(spec).Name==speciesName){
//       Species=spec;
//     }
//   }
  if (PathData.Path.Communicator.MyProc()==0){
    IOSection.WriteVar("Type","PermutationCount");
    IOSection.WriteVar("Cumulative", false);
  }
  CycleCount.resize(PathData.Path.NumParticles());
  CycleCount=0;
  TotalCounts=0;
  /// Now write the one-time output variables
//   if (PathData.Path.Communicator.MyProc()==0)
//     WriteInfo();

}

void PermutationCountClass::Accumulate()
{
  TimesCalled++;
  if (TimesCalled % DumpFreq==0){
    WriteBlock();
  }
  if ((TimesCalled % Freq)!=0){
    return;
  }
  TotalCounts++;
  PathClass &Path= PathData.Path;
  Array<bool,1> countedAlready(PathData.Path.NumParticles());
  countedAlready =false;
  int ptcl=0;
  while (ptcl<PathData.Path.NumParticles()){
    if (!countedAlready(ptcl)){
      int startPtcl=ptcl;
      int roamingPtcl=ptcl;
      int cycleLength=0;
      roamingPtcl=PathData.Path.Permutation(roamingPtcl);
      while (roamingPtcl!=startPtcl){
	countedAlready(roamingPtcl)=true;
	cycleLength++;
	roamingPtcl=PathData.Path.Permutation(roamingPtcl);
      }
      CycleCount(cycleLength)++;
    }
    ptcl++;
  }
}


