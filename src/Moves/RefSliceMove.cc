#include "RefSliceMove.h"

void RefSliceMoveClass::Read(IOSectionClass &in)
{
  string permuteType, speciesName;
  StageClass *permuteStage;
  assert (in.ReadVar ("NumLevels", NumLevels));
  assert (in.ReadVar ("Species", speciesName));
  assert (in.ReadVar ("name", Name));
  SpeciesNum = PathData.Path.SpeciesNum (speciesName);

  /// Set up permutation
  assert (in.ReadVar ("PermuteType", permuteType));
  if (permuteType == "TABLE") 
    permuteStage = new TablePermuteStageClass (PathData,SpeciesNum,NumLevels);
  else if (permuteType == "NONE") 
    permuteStage = new NoPermuteStageClass(PathData, SpeciesNum, NumLevels);
  else {
    cerr << "Unrecognized PermuteType, """ << permuteType << """\n";
    exit(EXIT_FAILURE);
  }
  permuteStage->Read (in);
  Stages.push_back (permuteStage);
  
  for (int level=NumLevels-1; level>=0; level--) {
    BisectionStageClass *newStage = new BisectionStageClass (PathData, level);
    newStage->Actions.push_back(&PathData.Actions.Kinetic);
    newStage->Actions.push_back(&PathData.Actions.ShortRange);
    if (level == 0) {
      if (PathData.Path.LongRange) 
      if (PathData.Actions.UseRPA)
      newStage->Actions.push_back(&PathData.Actions.LongRangeRPA);
      else
      newStage->Actions.push_back(&PathData.Actions.LongRange);
      
      if ((PathData.Actions.NodalActions(SpeciesNum)!=NULL)) {
	cerr << "Adding fermion node action for species " 
	     << speciesName << endl;
	newStage->Actions.push_back(PathData.Actions.NodalActions(SpeciesNum));
      }
    }
    newStage->BisectionLevel = level;
    Stages.push_back (newStage);
  }

  // Add the second stage of the permutation step
  Stages.push_back (permuteStage);
}


bool RefSliceMoveClass::NodeCheck()
{
  if (PathData.Actions.NodalActions(SpeciesNum) != NULL) {
    PathClass &Path = PathData.Path;
    // Broadcast the new reference path to all the other processors
    PathData.Path.BroadcastRefPath();
    
    // Calculate local nodal action
    SetMode(OLDMODE);
    double oldLocalNode = 
      PathData.Actions.NodalActions(SpeciesNum)->Action
      (0, Path.NumTimeSlices()-1, ActiveParticles, 0);
    SetMode(NEWMODE);
    double newLocalNode = 
      PathData.Actions.NodalActions(SpeciesNum)->Action
      (0, Path.NumTimeSlices()-1, ActiveParticles, 0);
    
//     if (newLocalNode > 1.0e20) 
//       cerr << "Local node crossing in NodeCheck().\n";

    // Do global sum over processors
    double localChange = newLocalNode - oldLocalNode;
    double globalChange = PathData.Path.Communicator.AllSum (localChange);
    bool toAccept = (-globalChange)>=log(PathData.Path.Random.Common()); 

//     fprintf (stderr, "old = %1.12e\n", oldLocalNode);
//     fprintf (stderr, "new = %1.12e\n", newLocalNode);
//     fprintf (stderr, "globalChange = %1.12e\n", globalChange);
//     fprintf (stderr, "localchange = %1.12e\n", localChange);
//     fprintf (stderr, "globalchange = %1.12e\n", globalChange);
//     cerr << "toAccept = " << (toAccept ? "true\n" : "false\n");

    return toAccept;
  }
  else
    return true;
}


/// This version is for the processor with the reference slice
void RefSliceMoveClass::MakeMoveMaster()
{
  PathClass &Path = PathData.Path;
  int myProc = PathData.Path.Communicator.MyProc();

  int firstSlice, lastSlice;
  Path.SliceRange (myProc, firstSlice, lastSlice);
  // Choose time slices
  int localRef = Path.GetRefSlice() - firstSlice;
  int bisectSlices = (1 << NumLevels);
  int minStart = max(0, localRef - bisectSlices);
  int maxStart = 
    min (Path.NumTimeSlices()-1-bisectSlices, localRef);
  maxStart = max (0, maxStart);
  Slice1 = Path.Random.LocalInt(maxStart-minStart) + minStart;
  Slice2 = Slice1 + bisectSlices;
  assert (Slice1 >= 0);
  assert ((localRef >= Slice1) && (localRef <= Slice2));

  ActiveParticles.resize(1);
  ActiveParticles(0) = -1;
  // Go through local stages
  bool toAccept=true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;
  while (stageIter!=Stages.end() && toAccept){
    toAccept = (*stageIter)->Attempt(Slice1,Slice2,
				     ActiveParticles,prevActionChange);
    stageIter++;
  }
  // Broadcast acceptance or rejection 
  int accept = toAccept ? 1 : 0;
  PathData.Path.Communicator.Broadcast (myProc, accept);

  // Now, if we accept local stages, move on to global nodal
  // decision. 
  if (toAccept) {
//     if ((NodeAccept+NodeReject) % 1000 == 999)
//       fprintf (stderr, "Node accept ratio = %5.3f\n",
// 	       (double)NodeAccept/(double)(NodeAccept+NodeReject));
    if (NodeCheck()) {
      NodeAccept++;
      Accept();
      Path.RefPath.AcceptCopy();
    }
    else {
      NodeReject++;
      Reject();
      Path.RefPath.RejectCopy();
    }
  }
  // Otherwise, reject the whole move
  else {
    Reject();
    Path.RefPath.RejectCopy();
  }

}


/// This version is for processors that do no own the reference slice 
void RefSliceMoveClass::MakeMoveSlave()
{
  PathClass &Path = PathData.Path;
  int myProc = PathData.Path.Communicator.MyProc();
  int master = Path.SliceOwner (Path.GetRefSlice());

//   /// Choose time slices for local bisections
//   /// We don't have the reference slice, so we can be cavalier
//   int firstSlice, lastSlice;
//   Path.SliceRange (myProc, firstSlice, lastSlice);
//   // Choose time slices
//   int bisectSlices = (1 << NumLevels);
//   int maxStart = Path.NumTimeSlices()-bisectSlices;
//   int slice1 = Path.Random.LocalInt(maxSlice);
//   int slice2 = slice1 + bisectSlices;

  int accept;
  /// Receive broadcast from Master.
  cerr << "Before SetMode(NEWMODE), mode = " << 
  SetMode (NEWMODE);
  PathData.Path.Communicator.Broadcast (master, accept);
  if (accept==1) {
    if (NodeCheck()) {
      Path.RefPath.AcceptCopy();
      NodeAccept++;
    }
    else {
      Path.RefPath.RejectCopy();
      NodeReject++;
    }
  }
}


void RefSliceMoveClass::MakeMove()
{
  PathClass &Path = PathData.Path;
  MasterProc = Path.SliceOwner (Path.GetRefSlice());
  if (PathData.Path.Communicator.MyProc() == MasterProc)
    MakeMoveMaster();
  else
    MakeMoveSlave();
  if ((NodeAccept+NodeReject) % 1000 == 999)
    fprintf (stderr, "Node accept ratio = %5.3f\n",
	     (double)NodeAccept/(double)(NodeAccept+NodeReject));
}
