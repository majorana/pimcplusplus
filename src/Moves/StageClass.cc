#include "StageClass.h"


void StageClass::Read(IOSectionClass &in)
{
  ///Do nothing for now
}
void StageClass::WriteRatio()
{
  AcceptRatioVar.Write((double)NumAccepted/(double)NumAttempted);
  AcceptRatioVar.Flush();
}

void StageClass::Accept()
{
  NumAccepted++;
  NumAttempted++;
}

void StageClass::Reject()
{
  NumAttempted++;
}


bool LocalStageClass::Attempt(int &slice1, int &slice2, 
			      Array<int,1> &activeParticles,
			      double &prevActionChange)
{
  SetMode (NEWMODE);
  double sampleRatio=Sample(slice1,slice2,activeParticles);
  SetMode(OLDMODE);
  double oldAction=StageAction(slice1,slice2,activeParticles);
  SetMode(NEWMODE);
  double newAction =StageAction(slice1,slice2,activeParticles);
  double currActionChange=newAction-oldAction;
  double logAcceptProb=log(sampleRatio)-currActionChange+prevActionChange;
  bool toAccept = logAcceptProb>=log(PathData.Path.Random.Local()); /// Accept condition
  if (toAccept)
    NumAccepted++;
  NumAttempted++;
  prevActionChange=currActionChange;
  return toAccept;
}

bool CommonStageClass::Attempt(int &slice1, int &slice2, 
			       Array<int,1> &activeParticles,
			       double &prevActionChange)
{
  assert (slice1 == 0);
  assert (slice2 == PathData.NumTimeSlices()-1);

  SetMode (NEWMODE);
  double sampleRatio=Sample(slice1,slice2,activeParticles);
  SetMode(OLDMODE);
  double oldAction= GlobalStageAction(activeParticles);
  SetMode(NEWMODE);
  double newAction = GlobalStageAction(activeParticles);
  //  perr << "oldAction = " << oldAction
  //       << "newAction = " << newAction << endl;
  double currActionChange=newAction-oldAction;
  double logAcceptProb=log(sampleRatio)-currActionChange+prevActionChange;
  bool toAccept = logAcceptProb>=log(PathData.Path.Random.Common()); /// Accept condition
  if (toAccept)
    NumAccepted++;
  NumAttempted++;
  prevActionChange=currActionChange;
  return toAccept;
}

double StageClass::StageAction(int startSlice,int endSlice,
				      const Array<int,1> &changedParticles)
{
  double TotalAction=0.0;
  list<ActionBaseClass*>::iterator actionIter=Actions.begin();
  while (actionIter!=Actions.end()){
    TotalAction += 
      ((*actionIter)->Action(startSlice, endSlice, changedParticles,
			     BisectionLevel));
    actionIter++;
  }
  return TotalAction;
}
