#include "PathDataClass.h"

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
  int numClones = N / procsPerClone;
  MyCloneNum = WorldComm.MyProc()/procsPerClone;
  // Create IntraComm
  WorldComm.Split(MyCloneNum, IntraComm);
  
  Array<int,1> ranks (numClones);
  for (int clone=0; clone<numClones; clone++)
    ranks(clone) = clone*procsPerClone;
  WorldComm.Subset (ranks, InterComm);
  
  // Now, set up random number generator
  Random.Init (314159, numClones);
}


void PathDataClass::MoveRefSlice (int absSlice)
{
  int first, last;
  int myProc   = Path.Communicator.MyProc();
  int numProcs = Path.Communicator.NumProcs();
  Path.SliceRange(numProcs-1, first, last);
  /// Min slices is the minimum number of slices owned by a processor
  /// It is the maximum we can shift the path at a time.
  int minSlices = last - first;
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
