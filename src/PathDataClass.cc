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
  
  if (IntraComm.MyProc()==0) {    
    Array<int,1> ranks (numClones);
    for (int clone=0; clone<numClones; clone++)
      ranks(clone) = clone*procsPerClone;
    WorldComm.Subset (ranks, InterComm);
  }
  
  // Now, set up random number generator
  Random.Init (314159, numClones);
}
