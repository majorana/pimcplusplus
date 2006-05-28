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

#include "PathDataClass.h"


void PathDataClass::Next(int &slice, int &ptcl)
{
  if (Join==slice)
    ptcl=Path.Permutation(ptcl);
  slice=(slice+1) % Path.NumTimeSlices();
}


void 
PathDataClass::WormInfo(int &headSlice, int &headPtcl,
		    int &tailSlice, int &tailPtcl,
		    int &numEmpty, int &wormSize)
{
#if 1==2
  tailSlice=-1;tailPtcl=-1;headSlice=-1;headPtcl=-1;
  numEmpty=0;
  wormSize=0;
  for (int slice=0;slice<Path.NumTimeSlices();slice++)
    for (int ptcl=0;ptcl<Path.NumParticles();ptcl++){
      if (SliceFullandPreviousSliceEmpty(slice,ptcl)){
	headSlice=slice;
	headPtcl=ptcl;
      }
      if (SliceFullandNextSliceEmpty(slice,ptcl)){
	tailSlice=slice;
	tailPtcl=ptcl;
      }
      if (Path.ParticleExist(slice,ptcl)==0 && slice!=Path.NumTimeSlices()-1)
	numEmpty++;
    }
  int currSlice=headSlice;
  int currPtcl=headPtcl;
  while (currSlice!=tailSlice || currPtcl!=tailPtcl){
    if (currSlice!=Path.NumTimeSlices()-1)
      wormSize++;
    Next(currSlice,currPtcl);
  } 
#endif
}
void
PathDataClass::MoveTailToSliceZero()
{
  int tailPtcl,tailSlice;
  int lastSlice=NumTimeSlices()-1;
  MoveJoin(lastSlice);
  FindTail(tailSlice,tailPtcl);
  int needToShift=tailSlice;
  ShiftData(-needToShift);
  Join=lastSlice-needToShift;


}


///Will only work in serial
void PathDataClass::MoveLinkToEnd(int linkToMove)
{
  //  cerr<<"Being called "<<linkToMove<<endl;
  
  int endTimeSlice=Path.NumTimeSlices()-1; //las slice on processor
  int needToShift=endTimeSlice-linkToMove;
  MoveJoin(linkToMove);
  if (needToShift!=0){
    ShiftData(needToShift);
    Join=linkToMove+needToShift;
    assert(Join==Path.NumTimeSlices()-1);
  }
  return;


}


///Will only work in serial
void PathDataClass::MoveOpenLinkToEnd()
{
  //  cerr<<"Calling moveopenlinktoend"<<endl;
  int openLink=Path.OpenLink;
  int endTimeSlice=Path.NumTimeSlices()-1; //las slice on processor
  int needToShift=endTimeSlice-openLink;
  //  cerr<<openLink<<" "<<endTimeSlice<<" "<<needToShift<<endl;
    MoveJoin(openLink);
  if (needToShift!=0){
    ShiftData(needToShift);
    Join=Path.OpenLink;
  }
  //  cerr<<"Join: "<<Join<<" "<<Path.OpenLink<<" "<<Path.NumTimeSlices()-1<<endl;
  assert(Join==Path.NumTimeSlices()-1);
  return;


}

///Worm Moves////////
bool PathDataClass::SliceFullandNextSliceEmpty(int slice,int ptcl)
{
#if 1==2
  int nextSlice=(slice+1) % NumTimeSlices();
  int nextPtcl;
  if (Join==slice)
    nextPtcl=Path.Permutation(ptcl);
  else
    nextPtcl=ptcl;
  return (Path.ParticleExist(slice,ptcl)==1.0 && Path.ParticleExist(nextSlice,nextPtcl)==0.0);
#endif
}

bool PathDataClass::SliceFullandPreviousSliceEmpty(int slice,int ptcl)
{
#if 1==2
  int prevSlice=((slice-1)+Path.NumTimeSlices() ) % Path.NumTimeSlices();
  int prevPtcl=0;
  if (Join==prevSlice){
    while (Path.Permutation(prevPtcl)!=ptcl)
      prevPtcl++;
  }
  else
    prevPtcl=ptcl;

  return (Path.ParticleExist(slice,ptcl)==1.0 && Path.ParticleExist(prevSlice,prevPtcl)==0.0);

#endif
}
  

void PathDataClass::FindHead(int &headSlice,int &headPtcl)
{
  for (int slice=0;slice<Path.NumTimeSlices();slice++)
    for (int ptcl=0;ptcl<Path.NumParticles();ptcl++)
      if (SliceFullandPreviousSliceEmpty(slice,ptcl)){
	headSlice=slice;
	headPtcl=ptcl;
	return;
      }
}

void PathDataClass::FindTail(int &tailSlice,int &tailPtcl)
{
  for (int slice=0;slice<Path.NumTimeSlices();slice++)
    for (int ptcl=0;ptcl<Path.NumParticles();ptcl++)
      if (SliceFullandNextSliceEmpty(slice,ptcl)){
	tailSlice=slice;
	tailPtcl=ptcl;
	return;
      }
}



/////////////



void PathDataClass::Read (IOSectionClass &in)
{
  ///  //MINOR HACK!
  ////  Join=NumTimeSlices()-1;
  ///  //END MINOR HACK!
  int N = WorldComm.NumProcs();
  int procsPerClone = 1;
  if (N > 1) {
    assert (in.OpenSection ("Parallel"));
    assert (in.ReadVar("ProcsPerClone", procsPerClone));
    in.CloseSection();
  }
  
  // Setup Inter- and IntraComms
  assert ((N % procsPerClone) == 0);
  NumClones = N / procsPerClone;
  MyCloneNum = WorldComm.MyProc()/procsPerClone;
  // Create IntraComm
  WorldComm.Split(MyCloneNum, IntraComm);
  
  Array<int,1> ranks (NumClones);
  for (int clone=0; clone<NumClones; clone++)
    ranks(clone) = clone*procsPerClone;
  WorldComm.Subset (ranks, InterComm);
  
  int seed;
  bool haveSeed = in.ReadVar ("Seed", seed);
  // Now, set up random number generator


  //  int seed;
  if (in.ReadVar("Seed",Seed)){
    Random.Init (Seed, NumClones);
  }
  else {
    Seed=Random.InitWithRandomSeed(NumClones);
  }
  //    Random.Init (314159, numClones);
  
  Path.MyClone=IntraComm.MyProc()/procsPerClone;
}


void PathDataClass::MoveRefSlice (int absSlice)
{
  int first, last;
  int myProc   = Path.Communicator.MyProc();
  int numProcs = Path.Communicator.NumProcs();
  Path.SliceRange(numProcs-1, first, last);
  /// Min slices is the minimum number of slices owned by a processor
  /// It is the maximum we can shift the path at a time.
  int minSlices = last - first-1;
  Path.SliceRange(myProc, first, last);
  if (absSlice < Path.GetRefSlice()) { // Shift to left
    while ((Path.GetRefSlice()-absSlice) > minSlices) {
      MoveJoin(minSlices);
      ShiftData (-minSlices);
      Join = 0;
    }
    int shift = absSlice - Path.GetRefSlice();
    MoveJoin(-shift);
    ShiftData (shift);
    Join = 0;
  }
  else if (absSlice > Path.GetRefSlice()) { // Shift to right
    while ((absSlice-Path.GetRefSlice()) > minSlices) {
      MoveJoin(0);
      Path.ShiftData(minSlices);
      Join = minSlices;
    }
    int shift = absSlice - Path.GetRefSlice();
    MoveJoin(0);
    ShiftData(shift);
    Join = shift;
  }
  MoveJoin (Path.NumTimeSlices()-1);
  assert (Path.GetRefSlice() == absSlice);
}


#include <sys/time.h>

int
PathDataClass::GetWallTime()
{
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv, &tz);
  return (int)tv.tv_sec;
}


bool
PathDataClass::ExceededWallTime()
{
  if (MaxWallTime == -1)
    return false;
  return ((GetWallTime()-StartWallTime) > MaxWallTime);
}

void
PathDataClass::SetMaxWallTime(int maxWallTime)
{
  MaxWallTime = maxWallTime;
}

