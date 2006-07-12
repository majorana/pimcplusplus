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
#include "BisectionBlock.h"
#include "StructureRejectStage.h"
#include "CouplingStage.h"
//#include "WormPermuteStage.h"
#include "OpenStage.h"
#include "NoPermuteStage.h"


void BisectionBlockClass::Read(IOSectionClass &in)
{
  //  bool orderNBosons;
  string permuteType, speciesName;
  //  StageClass *permuteStage;
  assert (in.ReadVar ("NumLevels", NumLevels));
  assert (NumLevels <= PathData.Actions.GetMaxLevels());
  LowestLevel = 0;
  if (!in.ReadVar ("LowestLevel", LowestLevel))
    LowestLevel = 0;
  assert (LowestLevel < NumLevels);
  assert (in.ReadVar ("Species", speciesName));
  assert (in.ReadVar ("StepsPerBlock", StepsPerBlock));
  Josephson=false;
  if (in.ReadVar("Josephson",Josephson))
    cerr << "Read in Josephson Data" << endl;
  //  in.ReadVar("OrderNBosons",orderNBosons);
  bool addStructureRejectStage=false;
  in.ReadVar("StructureReject",addStructureRejectStage);
  SpeciesNum = PathData.Path.SpeciesNum (speciesName);
  if    (PathData.Path.Species(SpeciesNum).GetParticleType() == FERMION){
    HaveRefslice=true;
  }
  HaveRefslice = 
    ((PathData.Path.Species(SpeciesNum).GetParticleType() == FERMION) &&
     (PathData.Actions.NodalActions(SpeciesNum) != NULL) &&
     (!PathData.Actions.NodalActions(SpeciesNum)->IsGroundState()));
//   cerr<<"I have a ref slice? "<<HaveRefslice<<" ";
//   cerr<<(PathData.Path.Species(SpeciesNum).GetParticleType() == FERMION); 
//   cerr<<(PathData.Actions.NodalActions(SpeciesNum) != NULL);
//   cerr<<(!PathData.Actions.NodalActions(SpeciesNum)->IsGroundState());
//   cerr<<(!(PathData.Actions.NodalActions(SpeciesNum)->IsGroundState()));
//   cerr<<endl;
  /// Set up permutation
  assert (in.ReadVar ("PermuteType", permuteType));
  if (permuteType == "TABLE") 
    PermuteStage = new TablePermuteStageClass(PathData, SpeciesNum, NumLevels,
					      IOSection);
  else if (permuteType=="COUPLE")
    PermuteStage= new CoupledPermuteStageClass(PathData,SpeciesNum,NumLevels,
					       IOSection);
  else if (permuteType=="WORMMOVE")
    PermuteStage=new OpenStageClass(PathData,SpeciesNum,NumLevels,
				    IOSection);
  else if (permuteType == "NONE") 
    PermuteStage = new NoPermuteStageClass(PathData, SpeciesNum, NumLevels,
					   IOSection);
//   else if (permuteType=="OPEN")
//     PermuteStage = new WormPermuteStageClass(PathData,SpeciesNum,NumLevels,
// 					     IOSection);
  else {
    cerr << "Unrecognized PermuteType, """ << permuteType << """\n";
    exit(EXIT_FAILURE);
  }
  PermuteStage->Read (in);
  Stages.push_back (PermuteStage);
  //   if (permuteType=="OPEN"){
//     EmptyStageClass *newStage=new EmptyStageClass(PathData,NumLevels-1,OutSection);
//     newStage->Read(in);
//     newStage->Actions.push_back(&PathData.Actions.Kinetic);
//     newStage->Actions.push_back(&PathData.Actions.OpenLoopImportance);
//     Stages.push_back(newStage);
//   }
//   if (PathData.Path.Random.Local()>0.5)
//     PathData.Path.ExistsCoupling=1;
//   else
//     PathData.Path.ExistsCoupling=0;
  for (int level=NumLevels-1; level>=LowestLevel; level--) {
    if (Josephson){
      BisectionJosephsonStageClass *newStage;
      //      BisectionStageClass *newStage;


      newStage = new BisectionJosephsonStageClass(PathData, level,
						  IOSection);
      if (level==LowestLevel)
	newStage->Actions.push_back(&PathData.Actions.DualHermele);
      //newStage->Actions.push_back(&PathData.Actions.Josephson);
      //      newStage = new BisectionStageClass(PathData, level,
      //					  IOSection);
      newStage->TotalLevels=NumLevels;
      newStage->BisectionLevel = level;
      Stages.push_back (newStage);
    

    }
    else{
      BisectionStageClass *newStage;
      newStage = new BisectionStageClass (PathData, level,
					  IOSection);
      newStage->TotalLevels=NumLevels;

      newStage->Actions.push_back(&PathData.Actions.Kinetic);
      //      newStage->Actions.push_back(&PathData.Actions.PairFixedPhase);
      
      if (PathData.Path.OpenPaths && level==LowestLevel){
	newStage->Actions.push_back(&PathData.Actions.OpenLoopImportance);
      }
      if (PathData.Path.OpenPaths && level>LowestLevel) // && permuteType=="OPEN")
	cerr<<"Don't look at short range"<<endl;
      else if (PathData.Path.OpenPaths && level>LowestLevel){
	newStage->Actions.push_back(&PathData.Actions.ShortRangeApproximate);
      }
      else if (PathData.Path.OrderN){
	newStage->Actions.push_back(&PathData.Actions.ShortRangeOn);
      }
      else{ // if (level==LowestLevel) //HACK HERE CURRENTLY 
	cerr<<"Adding short range values"<<endl;
	newStage->Actions.push_back(&PathData.Actions.ShortRange);
      }
      //      else
      //	int dummy=5;
      if (level == LowestLevel) {
	if (addStructureRejectStage){
	  newStage->Actions.push_back(&PathData.Actions.StructureReject);
	}
	///If it's David's long range class then do this
	if (PathData.Path.DavidLongRange){
	  newStage->Actions.push_back(&PathData.Actions.DavidLongRange);
	}
	else if (PathData.Actions.HaveLongRange()){
	  if (PathData.Actions.UseRPA)
	    newStage->Actions.push_back(&PathData.Actions.LongRangeRPA);
	  else
	    newStage->Actions.push_back(&PathData.Actions.LongRange);
	}
	if ((PathData.Actions.NodalActions(SpeciesNum)!=NULL)) {
	  cerr << "Adding fermion node action for species " 
	       << speciesName << endl;
	  newStage->Actions.push_back(PathData.Actions.NodalActions(SpeciesNum));
	}
      }
      newStage->BisectionLevel = level;
      Stages.push_back (newStage);

    }
  }
  // Add the second stage of the permutation step
  /// EVIL BAD ERROR!!!  Pushing onto the stack twice causes the stage
  /// to be accepted twice, which causes swapping the forward and
  // reverse tables twice!
  // Stages.push_back (PermuteStage);

//   ///HACK! Addding a stage that will reject the move if the structure
//   //factor gets too large
  bool useStructureRejectStage=false;
  in.ReadVar("StructureReject",useStructureRejectStage);
  if (useStructureRejectStage){
    StructureRejectStageClass* structureReject =
      new StructureRejectStageClass(PathData,in,IOSection);
    structureReject->Read(in);
    Stages.push_back(structureReject);
  }
}


void BisectionBlockClass::ChooseTimeSlices()
{
  //  if (PathData.Path.Communicator.MyProc()==0)
    //    cerr<<"Choosing time slices"<<endl;
  PathClass &Path = PathData.Path;
  int myProc = PathData.Path.Communicator.MyProc();
  // do something special to avoid moving reference slice
  if (HaveRefslice &&
      Path.SliceOwner(Path.GetRefSlice()) == myProc) {
    int bSlices = 1<<NumLevels;
    int myStart, myEnd;
    Path.SliceRange (myProc, myStart, myEnd);
    // refSlice is relative to my first slice
    int refSlice = Path.GetRefSlice() - myStart;
    int numSlices = Path.NumTimeSlices();
    assert(bSlices*2<numSlices);
    if (refSlice < bSlices) {  
      int numStarts = numSlices - bSlices - refSlice;
      Slice1 = refSlice + Path.Random.LocalInt(numStarts);
      Slice2 = Slice1 + bSlices;
      if (Slice2>=numSlices){
	cerr << "(Slice 2, numSlices): (" << Slice2 << ", " 
	     << numSlices << ")" << " "<< bSlices << " "
	     << numStarts << " " << refSlice << endl;
      }
      assert (Slice2 < numSlices);
    }
    else if (refSlice > (numSlices-1-bSlices) && refSlice<numSlices) {
      int numStarts = refSlice - bSlices + 1;
      Slice1 = Path.Random.LocalInt(numStarts);
      Slice2 = Slice1+bSlices;
      assert (Slice2 <= refSlice);
    }
    else {
      int numBefore = refSlice - bSlices +1;
      int numAfter = numSlices -bSlices - refSlice;
      int numStarts = numBefore+numAfter;
      int start = Path.Random.LocalInt (numStarts);
      if (start < numBefore)
	Slice1 = start;
      else
	Slice1 = start-numBefore+refSlice;
      Slice2 = Slice1+bSlices;
      assert ((refSlice <= Slice1) || (refSlice >= Slice2));
      assert (Slice2 < numSlices);
    }
  }
  // Bryan should fill this part in.
  //  else if (PathData.Path.OpenPtcl) {
    

  //  }
  else {
    int sliceSep = 1<<NumLevels;
    assert (sliceSep < PathData.Path.NumTimeSlices());
    int numLeft = PathData.Path.NumTimeSlices()-sliceSep;
    if (numLeft < 0) {
      cerr << "Not enough slices to bisect with " << NumLevels 
	   << " levels in BisectionBlock.\n";
      abort();
    }
    Slice1 = PathData.Path.Random.LocalInt (numLeft);
    Slice2 = Slice1+sliceSep;
  }
  //  if (PathData.Path.Communicator.MyProc()==0)
  //    cerr<<"Ending Choosing time slices"<<endl;
}

void BisectionBlockClass::MakeMove()
{
  ChooseTimeSlices();
  PathData.MoveJoin(Slice2);

//   if (PathData.Path.OrderN){
//     for (int slice=Slice1;slice<=Slice2;slice++)
//       PathData.Path.Cell.BinParticles(slice);
//   }

  ((PermuteStageClass*)PermuteStage)->InitBlock(Slice1,Slice2);
  ActiveParticles.resize(1);
  for (int step=0; step<StepsPerBlock; step++) {
    ActiveParticles(0)=-1;
    MultiStageClass::MakeMove();
  }
  if (LowestLevel != 0)
    MakeStraightPaths();
}


void
BisectionBlockClass::MakeStraightPaths()
{
  PathClass &Path = PathData.Path;
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
  
      
