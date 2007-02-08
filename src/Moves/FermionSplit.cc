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

#include "EmptyStage.h"
#include "FermionSplit.h"
#include "StructureRejectStage.h"
#include "CouplingStage.h"
//#include "WormPermuteStage.h"
#include "OpenStage.h"
#include "NoPermuteStage.h"


void FermionSplitClass::Read(IOSectionClass &in)
{
  string speciesString;
  assert(in.ReadVar("Species",speciesString));
  SpeciesNum=PathData.Path.SpeciesNum(speciesString);
  assert(in.ReadVar("NumLevels",NumLevels));
  LeviFlightFSPClass *leviFlight;
  leviFlight = new LeviFlightFSPClass(PathData,level,
				     IOSection);
  leviFlight->Actions.push_back(&PathData.Path.Kinetic);
  leviFlight->Actions.push_back(&PathData.Path.ShortRange);
}


void FermionSplitClass::ChooseTimeSlices()
{
#ifdef BUILD_DEV
  PathClassDev &Path = PathData.Path;
#else
  PathClass &Path = PathData.Path;
#endif
  int myProc = PathData.Path.Communicator.MyProc();
  int sliceSep = 1<<NumLevels;
  assert (sliceSep < PathData.Path.NumTimeSlices());
  int numLeft = PathData.Path.NumTimeSlices()-sliceSep;
  if (numLeft < 0) {
    cerr << "Not enough slices to bisect with " << NumLevels 
	 << " levels in FermionSplit.\n";
    abort();
  }
  Slice1 = PathData.Path.Random.LocalInt (numLeft);
  Slice2 = Slice1+sliceSep;
}


void FermionSplitClass::ChoosePermutations(CycleRep &myCycle)
{
  int ptcl1=PathData.Path.Random.LocalInt(PathData.Path.NumParticles());
  int ptcl2=PathData.Path.Random.LocalInt(PathData.Path.NumParticles());
  while (ptcl1==ptcl2){
    ptcl2=PathData.Path.Random.LocalInt(PathData.Path.NumParticles());
  }
  
  myCycle.Length=2;
  myCycle.CycleRep[0]=min(ptcl1,ptcl2);
  myCycle.CycleRep[1]=max(ptcl1,ptcl2);
}

void FermionSplitClass::MakeMove()
{
  CycleRep myCycle;
  ChooseTimeSlices(myCycle);
  PathData.MoveJoin(Slice2);
  ((PermuteStageClass*)PermuteStage)->InitBlock(Slice1,Slice2);
  
  ChoosePermutations();
  ActiveParticles.resize(1);
  for (int step=0; step<StepsPerBlock; step++) {
    ActiveParticles(0)=-1;
    MultiStageClass::MakeMove();
  }
  if (LowestLevel != 0)
    MakeStraightPaths();
}


void
FermionSplitClass::MakeStraightPaths()
{
#ifdef BUILD_DEV
  PathClassDev &Path = PathData.Path;
#else
  PathClass &Path = PathData.Path;
#endif
  SetMode(NEWMODE);
  int skip = 1<<LowestLevel;
  int first = Path.Species(SpeciesNum).FirstPtcl;
  int last = Path.Species(SpeciesNum).LastPtcl;
  double inc = 1.0/(double)skip;
  for (int slice=Slice1; slice < Slice2; slice += skip) 
    for (int ptcl=first; ptcl<=last; ptcl++) {
      dVec delta = Path.Velocity(slice, slice+skip, ptcl);
      double frac = inc;
      for (int s=1; s<skip; s++) {
	Path(slice+s,ptcl) = Path(slice,ptcl) + frac*delta;
	frac += inc;
      }
    }
  Array<int,1> ptcls(last-first+1);
  for (int ptcl=first; ptcl<=last; ptcl++)
    ptcls(ptcl-first) = ptcl;
  PathData.AcceptMove(Slice1, Slice2, ptcls);
}
  
      
