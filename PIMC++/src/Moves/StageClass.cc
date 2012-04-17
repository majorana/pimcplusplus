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

#include "StageClass.h"
#include "sys/time.h"

void StageClass::Read(IOSectionClass &in)
{
  ///Do nothing for now
}
void StageClass::WriteRatio()
{
  AcceptRatioVar.Write((double)NumAccepted/(double)NumAttempted);
  AcceptRatioVar.Flush();
}

//BUG: DOES NOT HAVE CORRECT SLICES!
void StageClass::Accept()
{
  NumAccepted++;
  NumAttempted++;
  for (list<ActionBaseClass*>::iterator actionIter = Actions.begin(); actionIter != Actions.end(); actionIter++)
    (*actionIter) -> AcceptCopy(0,0);
}

//BUG: DOES NOT HAVE CORRECT SLICES
void StageClass::Reject()
{
  NumAttempted++;
  for (list<ActionBaseClass*>::iterator actionIter = Actions.begin(); actionIter != Actions.end(); actionIter++) {
    (*actionIter) -> RejectCopy(0,0);
  }
}


bool LocalStageClass::Attempt(int &slice1, int &slice2, Array<int,1> &activeParticles, double &prevActionChange)
{
  struct timeval start, end;
  struct timezone tz;
  //cerr << BisectionLevel << " ";

  SetMode (NEWMODE);
  double sampleRatio = Sample(slice1,slice2,activeParticles);
  //cerr << sampleRatio << " ";

  SetMode (OLDMODE);
  gettimeofday(&start, &tz);
  double oldAction = StageAction(slice1,slice2,activeParticles);
  //cerr << oldAction << " ";

  SetMode(NEWMODE);
  double newAction = StageAction(slice1,slice2,activeParticles);
  //cerr << newAction << endl;

  gettimeofday(&end, &tz);
  double currActionChange = newAction - oldAction;
  double logAcceptProb = log(sampleRatio) - currActionChange + prevActionChange;
  double ran_number = PathData.Path.Random.Local();
  bool toAccept = logAcceptProb >= log(ran_number); // Accept condition

  if (abs(newAction) > 1e50 || abs(oldAction) > 1e50) { // Nodal Rejection
    if (toAccept) {
      cerr << "Broken Staging: " << BisectionLevel << " " << toAccept << " " << slice1 << " " << slice2 << " " << PathData.Path.GetRefSlice() << " " << PathData.Path.SliceOwner(PathData.Path.GetRefSlice()) << " " << PathData.Path.Communicator.MyProc() << " " << oldAction << " " << newAction << endl;
      toAccept = 0;
      assert(1==2);
    }
  }

  if (toAccept){
    NumAccepted++;
  }
  NumAttempted++;
  prevActionChange = currActionChange;

  TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);

  return toAccept;
}


bool CommonStageClass::Attempt (int &slice1, int &slice2, Array<int,1> &activeParticles, double &prevActionChange)
{
  assert (slice1 == 0);
  assert (slice2 == PathData.NumTimeSlices()-1);

  SetMode (NEWMODE);
  double sampleRatio=Sample(slice1,slice2,activeParticles);
  double logSampleRatio = log(sampleRatio);
  logSampleRatio = PathData.Path.Communicator.AllSum (logSampleRatio);
  SetMode(OLDMODE);
  double oldAction= GlobalStageAction(activeParticles);
  SetMode(NEWMODE);
  double newAction = GlobalStageAction(activeParticles);
  double currActionChange=newAction-oldAction;
  double logAcceptProb=logSampleRatio-currActionChange+prevActionChange;
  bool toAccept = logAcceptProb>=log(PathData.Path.Random.Common()); /// Accept condition
  if (toAccept) {
    NumAccepted++;
  }
  NumAttempted++;
  prevActionChange=currActionChange;

  return toAccept;
}


double StageClass::StageAction (int startSlice, int endSlice, const Array<int,1> &changedParticles)
{
  double TotalAction = 0.0;
  list<ActionBaseClass*>::iterator actionIter = Actions.begin();
  while (actionIter != Actions.end()) {
    double TempAction = ((*actionIter) -> Action(startSlice, endSlice, changedParticles, BisectionLevel));
    TotalAction += TempAction;
    //cerr << (*actionIter) -> GetName() << "(" << TempAction << ") ";
    actionIter++;
  }

  return TotalAction;
}


inline double StageClass::GlobalStageAction (const Array<int,1> &changedParticles)
{
  int slice1 = 0;
  int slice2 = PathData.Path.NumTimeSlices()-1;
  double localAction = 0.0;
  list<ActionBaseClass*>::iterator actionIter=Actions.begin();
  while (actionIter!=Actions.end()) {
    double action = ((*actionIter) -> Action(slice1, slice2, changedParticles, BisectionLevel));
    localAction += action;
    actionIter++;
  }
  double globalAction = PathData.Path.Communicator.AllSum (localAction);

  return globalAction;
}

