#include "PathDataClass.h"

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

void PathDataClass::Read (IOSectionClass &in)
{
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

