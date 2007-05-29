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
  in.ReadVar("NumPreSteps",TotalNumPreSteps);
  int stages = 1;
  in.ReadVar("NumPreStages",stages);
	Array<string,1> methodList;
	assert(in.ReadVar("MoveMethod",methodList));
  assert(methodList.size() == stages);
  Array<int,1> numActions(stages);
  numActions = 1;
  assert(in.ReadVar("NumActions", numActions));
  //assert((numActions.size()-1) == stages);
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
  int actionsToRead = numActions(numActions.size()-1);
	if (finalMethod == "Dummy"){
	  cerr << "Creating new dummy stage to evaluate the specified action (should be preceded by an actual move stage which is evaluate with a \"cheap\" action...";
    FinalStage = new PreSampleDummy(PathData, IOSection, actionsToRead, startIndex);
    cerr << " done." << endl;
  } else {
    cerr << "Actually, I'm not supporting anything other than ``Dummy'' for the final move stage for this algorithm" << endl;
    assert(0);
  }
  FinalStage->Read(in);

  // initialize stuff
  NumPtcls = PathData.Path.NumParticles();
  NumSlices = PathData.Path.NumTimeSlices();
  InitialPath.resize(NumSlices, NumPtcls);
  cerr << "PreSamplingClass read I have " << NumPtcls << " particles and " << NumSlices-1 << " time slices" << endl;

}

void PreSamplingClass::WriteRatio()
{
  list<StageClass*>::iterator stageIter=PreStages.begin();
  while (stageIter!=PreStages.end()){
    (*stageIter)->WriteRatio();
    stageIter++;
  }
  FinalStage->WriteRatio();
  // handle it here; don't use the default
  //MoveClass::WriteRatio();
  RatioVar.Write(double(NumFinalAccept)/NumSteps);
}

void PrintPaths(int S, int N, PathDataClass& PD)
{
  for(int s=0; s<S-1; s++){
    for(int n=0; n<N; n++){
      cerr << s << " " << n << " ";
      SetMode(OLDMODE);
      cerr << PD.Path(s,n) << " ";
      SetMode(NEWMODE);
      cerr << PD.Path(s,n) << endl;
    }
  }
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
  //cerr << "STORED PATH IS" << endl;
  //for(int n=0; n<NumPtcls; n++)
  //  cerr << 0 << " " << n << " " << InitialPath(0,n) << endl;

  //cout << "INITIAL PATHS OLD NEW" << endl;
  //PrintPaths(NumSlices, NumPtcls, PathData);
  NumSteps++;
  NumPreAccept = 0;
  PreDeltaAction = 0.0;
  bool finalAccept, toAccept;
  double prevActionChange = 0.0;
  Array<bool,1> alreadyActive(NumPtcls);
  alreadyActive = false;
  Array<int,1> myActiveP(0);

  for(int i=0; i<TotalNumPreSteps; i++){
    //cout << i << " ";
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
        //PathData.AcceptMove(Slice1,Slice2,ActiveParticles);
        (*stageIter)->Accept();
        //cout << setw(12) << prevActionChange;

        for(int p=0; p<ActiveParticles.size(); p++){
          int ptcl = ActiveParticles(p);
          if(!alreadyActive(ptcl)){
            alreadyActive(ptcl) = true;
            int size = myActiveP.size();
            myActiveP.resizeAndPreserve(size+1);
            myActiveP(size) = ptcl;
          }
        }

      } else {
        //PathData.RejectMove(Slice1,Slice2,ActiveParticles);
        Reject();
        (*stageIter)->Reject();
        //cout << setw(12) << 0.0 << " ";
      }

      //cout << "now activeP is " << myActiveP << endl;
      stageIter++;
    }
    //cout << PreDeltaAction << endl;

  }

  //cout << "############################################################################################################" << endl;
  //cerr << "Presampled Paths are OLD NEW" << endl;
  //PrintPaths(NumSlices, NumPtcls, PathData);
  //cerr << "Now resetting old path" << endl;
  SetMode(OLDMODE);
  AssignInitialPath();
  //cout << NumSteps << " FINAL ATTEMPT with prevActChg " << PreDeltaAction << endl;
  //cout << "FINAL PATHS OLD NEW" << endl;
  //PrintPaths(NumSlices, NumPtcls, PathData);
  //ActiveParticles.resize(NumPtcls);
  //for(int p=0; p<ActiveParticles.size(); p++)
  //  ActiveParticles(p) = p;
  ActiveParticles.resize(myActiveP.size());
  for(int p=0; p<ActiveParticles.size(); p++)
    ActiveParticles(p) = myActiveP(p);

  finalAccept = FinalStage->Attempt(Slice1, Slice2, ActiveParticles, PreDeltaAction);
  //finalAccept = FinalStage->Attempt(Slice1, Slice2, ActiveParticles, PreDeltaAction);

  if(finalAccept){
    NumFinalAccept++;
    FinalStage->Accept();
    Accept();
  } else {
    //SetMode(NEWMODE);
    //AssignInitialPath();
    //Accept();
    FinalStage->Reject();
    Reject();
  }
  
  //cout << NumSteps << "; presample accept " << NumPreAccept << "; FINAL STAGE ACCEPT ";
  cout << "FINAL ACCEPT ";
  cout << finalAccept;
  cout << " ratio " << setw(10) << double(NumFinalAccept)/NumSteps << endl;
}


void PreSamplingClass::Accept()
{
  //cout << "ACCEPT SLICES " << Slice1 << " " << Slice2 << " PTCLS " << ActiveParticles << endl;
  PathData.AcceptMove(Slice1,Slice2+1,ActiveParticles);
  //int N = ActiveParticles.size();
  //Array<dVec, 2> P(NumSlices, N);
  //SetMode(NEWMODE);
  //for(int s=0; s<NumSlices; s++){
  //  for(int p=0; p<N; p++){
  //    P(s,p) = PathData.Path(s,ActiveParticles(p));
  //  }
  //}

  //SetMode(OLDMODE);
  //for(int s=0; s<NumSlices; s++){
  //  for(int p=0; p<N; p++){
  //    PathData.Path.SetPos(s,ActiveParticles(p),P(s,p));
  //  }
  //}

  //NumAccepted++;
}

void PreSamplingClass::Reject()
{
  PathData.RejectMove(Slice1,Slice2+1,ActiveParticles);
  //int N = ActiveParticles.size();
  //Array<dVec, 2> P(NumSlices, N);
  //SetMode(OLDMODE);
  //for(int s=0; s<NumSlices; s++){
  //  for(int p=0; p<N; p++){
  //    P(s,p) = PathData.Path(s,ActiveParticles(p));
  //  }
  //}

  //SetMode(NEWMODE);
  //for(int s=0; s<NumSlices; s++){
  //  for(int p=0; p<N; p++){
  //    PathData.Path.SetPos(s,ActiveParticles(p),P(s,p));
  //  }
  //}
}

PreSampleDummy::PreSampleDummy(PathDataClass& PathData, IOSectionClass& IO, int actionsToRead, int startIndex):
  DummyEvaluate(PathData, IO, actionsToRead, startIndex)
{
  toRead = actionsToRead;
  startI = startIndex;
}

bool PreSampleDummy::Attempt(int &slice1, int &slice2, 
			      Array<int,1> &activeParticles,
			      double &prevActionChange)
{
  SetMode (NEWMODE);
  double sampleRatio=Sample(slice1,slice2,activeParticles);
  SetMode(OLDMODE);
  double oldAction=StageAction(slice1,slice2,activeParticles);
  double oldPreAction=PreSampleAction(slice1,slice2,activeParticles);
  SetMode(NEWMODE);
  double newAction =StageAction(slice1,slice2,activeParticles);
  double newPreAction=PreSampleAction(slice1,slice2,activeParticles);
  cout << "FINAL DELTA-ACTION " << setw(12) << newAction-oldAction;
  cout << "PRE DELTA-ACTION " << setw(12) << newPreAction-oldPreAction;
  double currActionChange=newAction-oldAction;
  prevActionChange = newPreAction - oldPreAction;

  //cout << "ATTEMPT log(sample) " << log(sampleRatio) << " -currActChg " << -currActionChange << " prevActChg " << prevActionChange << endl;
  double logAcceptProb=log(sampleRatio)-currActionChange+prevActionChange;
  //cout << "  so logAcceptProb is " << logAcceptProb << endl;
  //  AcceptProb=exp(logAcceptProb);
  //  OldAcceptProb=exp(log(1/sampleRatio)+currActionChange);
  //  if (AcceptProb>1.0)
  //    AcceptProb=1.0;
  //  if (OldAcceptProb>1.0)
  //    OldAcceptProb=1.0;
  double logRand = log(PathData.Path.Random.Local());
  bool toAccept = logAcceptProb>= logRand; /// Accept condition
  //cout << "  and toAccept is " << toAccept << " based on rn " << logRand << endl;
  if (toAccept){
    NumAccepted++;
  }
  NumAttempted++;
  //cout<<"Curr action change is "<<currActionChange<<endl;
  prevActionChange=currActionChange;

  return toAccept;
}

double PreSampleDummy::PreSampleAction(int startSlice,int endSlice,
				      const Array<int,1> &changedParticles)
{
  double TotalAction=0.0;
  list<ActionBaseClass*>::iterator actionIter=PreActions.begin();
  while (actionIter!=PreActions.end()){
    TotalAction += 
      ((*actionIter)->Action(startSlice, endSlice, changedParticles,
			     BisectionLevel));
    actionIter++;
  }
  return TotalAction;
}

void PreSampleDummy::Read (IOSectionClass &in){
	/// Read in the molecule to move by string or int
	cerr << "MolMoveClass::Read" << endl;
	string setMolecule;
	int indexMolecule;
	if(in.ReadVar("Molecule",setMolecule)){
		cerr << "got Molecule " << setMolecule << endl;
		molIndex = PathData.Mol.Index(setMolecule);
		numMol = PathData.Mol.NumMol(molIndex);
	}
	else if (in.ReadVar("Molecule",indexMolecule)){
    molIndex = indexMolecule;
		numMol = PathData.Mol.NumMol(molIndex);
	}
	else{
		molIndex = 0;
		numMol = PathData.Mol.NumMol(molIndex);
    cerr << "BY DEFAULT, ";
	}
	cerr << "Selected molecule index " << molIndex << "(" << PathData.Mol.NameOf(molIndex) << ") with " << numMol << endl;
	
  /// Read in update mode: GLOBAL, SEQUENTIAL, SINGLE	
	string setMode;
	assert(in.ReadVar("Mode", setMode));
	if(setMode == "GLOBAL"){
		mode = GLOBAL;
    MoveList = PathData.Mol.MolOfType(molIndex);
    cerr << "GLOBAL move init MoveList is " << MoveList << endl;
		//MoveList.resize(numMol);
		//for(int m=0; m<numMol; m++)
		//	MoveList(m) = m + PathData.Path.offset[molIndex];
		//cerr << "Set for GLOBAL particle updates." << endl;
	}
	else if(setMode == "SEQUENTIAL"){
		mode = SEQUENTIAL;
		MoveList.resize(1);
		MoveList(0) = 0;
		cerr << "Set for SEQUENTIAL particle updates." << endl;
	}
	else{
		cerr << "Set for SINGLE particle updates." << endl;
		MoveList.resize(1);
	}
	
  cerr << "Looking for NumPreActions " << endl;
  assert(in.ReadVar("NumPreActions", numPre));
  cerr << "Looking for NumFinalActions " << endl;
  assert(in.ReadVar("NumFinalActions", numFinal));

  Array<string,1> ActionList;
  assert (in.ReadVar ("Actions", ActionList));
  cerr << "  Looking for " << toRead << " actions starting at index " << startI << endl;
  assert(toRead == numFinal);
  //assert ((toRead + startI) <= ActionList.size());
  assert ((numPre + numFinal + startI) == ActionList.size());
  cerr << "PRE looking for " << numPre << endl;
	for(int a=0; a<numPre; a++){
    cerr << "Read in PreSample Actions " << ActionList << endl;
		string setAction = ActionList(a+startI);
	  if(setAction == "MoleculeInteractions"){
			// read should be done in Actions now
			//PathData.Actions.MoleculeInteractions.Read(in);
  		PreActions.push_back(&PathData.Actions.MoleculeInteractions);
			cerr << "  Added Molecule PreActions" << endl;
		}else if(setAction == "ST2Water"){
			//ActionBaseClass* newAction(PathData.Actions.CEIMCAction);
  		PreActions.push_back(&PathData.Actions.ST2Water);
			cerr << "Added ST2Water action" << endl;
		}else if(setAction == "Kinetic"){
  		PreActions.push_back(&PathData.Actions.Kinetic);
			cerr << "Added Kinetic action" << endl;
#ifdef USE_QMC
		}else if(setAction == "CEIMCAction"){
			//ActionBaseClass* newAction(PathData.Actions.CEIMCAction);
			PathData.Actions.CEIMCAction.Read(in);
  		PreActions.push_back(&PathData.Actions.CEIMCAction);
			cerr << "Added CEIMC calculation of BO energy" << endl;
#endif
		}else if(setAction == "LongRangeCoulomb"){
  		PreActions.push_back(&PathData.Actions.LongRangeCoulomb);
			cerr << "Added long-range coulomb interaction" << endl;
		}else if(setAction == "IonInteraction"){
  		PreActions.push_back(&PathData.Actions.IonInteraction);
			cerr << "Added intermolecular ion-ion interaction" << endl;
		} else if(setAction == "QBoxAction"){
  		PreActions.push_back(&PathData.Actions.QBoxAction);
			cerr << "Computing action with QBox DFT code" << endl;
		} else
    	cerr << "You specified " << setAction << ", which is not supported for this type of move" << endl;
	}

  cerr << "FINAL looking for " << numFinal << endl;
	for(int a=0; a<numFinal; a++){
    cerr << "Read in actions " << ActionList << endl;
		string setAction = ActionList(a+startI+numPre);
	  if(setAction == "MoleculeInteractions"){
			// read should be done in actions now
			//PathData.Actions.MoleculeInteractions.Read(in);
  		Actions.push_back(&PathData.Actions.MoleculeInteractions);
			cerr << "  Added Molecule Actions" << endl;
		}else if(setAction == "ST2Water"){
			//ActionBaseClass* newAction(PathData.Actions.CEIMCAction);
  		Actions.push_back(&PathData.Actions.ST2Water);
			cerr << "Added ST2Water action" << endl;
		}else if(setAction == "Kinetic"){
  		Actions.push_back(&PathData.Actions.Kinetic);
			cerr << "Added Kinetic action" << endl;
#ifdef USE_QMC
		}else if(setAction == "CEIMCAction"){
			//ActionBaseClass* newAction(PathData.Actions.CEIMCAction);
			PathData.Actions.CEIMCAction.Read(in);
  		Actions.push_back(&PathData.Actions.CEIMCAction);
			cerr << "Added CEIMC calculation of BO energy" << endl;
#endif
		}else if(setAction == "LongRangeCoulomb"){
  		Actions.push_back(&PathData.Actions.LongRangeCoulomb);
			cerr << "Added long-range coulomb interaction" << endl;
		}else if(setAction == "IonInteraction"){
  		Actions.push_back(&PathData.Actions.IonInteraction);
			cerr << "Added intermolecular ion-ion interaction" << endl;
		} else if(setAction == "QBoxAction"){
  		Actions.push_back(&PathData.Actions.QBoxAction);
			cerr << "Computing action with QBox DFT code" << endl;
		} else
    	cerr << "You specified " << setAction << ", which is not supported for this type of move" << endl;
	}
}
