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

#include "PreSampling.h"


void PreSamplingClass::Read(IOSectionClass& in)
{
  cerr << "Pre-Sampling Move read..." << endl;
  int stages = 1;
  in.ReadVar("NumPreStages",stages);
	Array<string,1> methodList;
	assert(in.ReadVar("MoveMethod",methodList));
  assert(methodList.size() == stages);
  Array<int,1> numActions(stages);
  numActions = 1;
  assert(in.ReadVar("NumActions", numActions));
  assert((numActions.size()-1) == stages);
  int startIndex = 0;
  for(int s=0; s<stages; s++){
    StageClass* MoveStage;
    string method = methodList(s);
    int actionsToRead = numActions(s);
    cerr << "  Init " << s+1 << " of " << stages << " stages: " << method << endl;
	  if(method == "Translate"){
	  	cerr << "Creating new Translate move...";
    	MoveStage = new MoleculeTranslate(PathData, IOSection, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "Rotate"){
	  	cerr << "Creating new Rotate move...";
    	MoveStage = new MoleculeRotate(PathData, IOSection, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "Multi"){
	  	cerr << "Creating new Multiple move...";
    	MoveStage = new MoleculeMulti(PathData, IOSection, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "Stretch"){
	  	cerr << "Creating new bond-stretching move...";
    	MoveStage = new BondStretch(PathData, IOSection, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "ForceBias"){
	  	cerr << "Creating new Force Bias move...";
    	MoveStage = new MoleculeForceBiasMove(PathData, IOSection, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else {
	  	cerr << "ERROR: method " << method << " is not supported." << endl;
	  }
    MoveStage->Read(in);
    PreStages.push_back(MoveStage);
    startIndex += numActions(s);
  }

  string finalMethod;
  assert(in.ReadVar("FinalStageMethod",finalMethod));
  int actionsToRead = numActions(numActions.size());
	if (finalMethod == "Dummy"){
	  cerr << "Creating new dummy stage to evaluate the specified action (should be preceded by an actual move stage which is evaluate with a \"cheap\" action...";
    FinalStage = new DummyEvaluate(PathData, IOSection, actionsToRead, startIndex);
    cerr << " done." << endl;
  } else {
    cerr << "Actually, I'm not supporting anything other than ``Dummy'' for the final move stage for this algorithm" << endl;
    assert(0);
  }

  // initialize stuff
  NumPtcls = PathData.Path.NumParticles();
  NumSlices = PathData.Path.NumTimeSlices();
  InitialPath.resize(NumSlices, NumPtcls);
  cerr << "PreSamplingClass read I have " << NumPtcls << " particles and " << NumSlices << " time slices" << endl;

}

void PreSamplingClass::WriteRatio()
{
  list<StageClass*>::iterator stageIter=PreStages.begin();
  while (stageIter!=PreStages.end()){
    (*stageIter)->WriteRatio();
    stageIter++;
  }
  // handle it here; don't use the default
  //MoveClass::WriteRatio();
  RatioVar.Write(double(NumFinalAccept)/NumSteps);
}

void PreSamplingClass::StoreInitialPath()
{
  for(int s=0; s<NumSlices; s++){
    for(int n=0; n<NumPtcls; n++){
      InitialPath(s,n) = PathData.Path(s,n);
    }
  }
}

void PreSamplingClass::AssignInitialPath()
{
  for(int s=0; s<NumSlices; s++){
    for(int n=0; n<NumPtcls; n++){
      PathData.Path.SetPos(s,n,InitialPath(s,n));
    }
  }
}

void PreSamplingClass::MakeMove()
{
  StoreInitialPath();  
  NumSteps++;
  NumPreAccept = 0;
  PreDeltaAction = 0.0;
  bool finalAccept, toAccept;
  double prevActionChange = 0.0;

  for(int i=0; i<TotalNumPreSteps; i++){
    list<StageClass*>::iterator stageIter=PreStages.begin();
    while (stageIter!=PreStages.end())
    {
      prevActionChange = 0.0;
      toAccept = (*stageIter)->Attempt(Slice1,Slice2,
	  			     ActiveParticles,prevActionChange);

      if(toAccept){
        NumPreAccept++;
        PreDeltaAction += prevActionChange;
        Accept();
        (*stageIter)->Accept();
      } else {
        Reject();
        (*stageIter)->Reject();
      }

      stageIter++;
    }

  }

  finalAccept = FinalStage->Attempt(Slice1, Slice2, ActiveParticles, PreDeltaAction);

  if(finalAccept){
    NumFinalAccept++;
    Accept();
    FinalStage->Accept();
  } else {
    AssignInitialPath();
    Accept();
    FinalStage->Reject();
  }

}


void PreSamplingClass::Accept()
{
  PathData.AcceptMove(Slice1,Slice2,ActiveParticles);
  //NumAccepted++;
}

void PreSamplingClass::Reject()
{
  PathData.RejectMove(Slice1,Slice2,ActiveParticles);
}

