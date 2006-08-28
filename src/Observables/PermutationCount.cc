/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

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
  TotalCounts++;
  PathClass &Path= PathData.Path;
  int N = PathData.Path.NumParticles();
  if (CountedAlready.size() != N) {
    CountedAlready.resize(N);
    TotalPerm.resize(N);
  }
  PathData.Path.TotalPermutation (TotalPerm);
  CountedAlready =false;
  int ptcl=0;
  while (ptcl < N) {
    if (!CountedAlready(ptcl)) {
      int startPtcl=ptcl;
      int roamingPtcl=ptcl;
      int cycleLength=0;
      roamingPtcl = TotalPerm(roamingPtcl);
      while (roamingPtcl!=startPtcl){
	CountedAlready(roamingPtcl)=true;
	cycleLength++;
	roamingPtcl=PathData.Path.Permutation(roamingPtcl);
      }
      CycleCount(cycleLength)++;
    }
    ptcl++;
  }
}


