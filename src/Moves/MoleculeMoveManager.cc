#include "MoleculeMoveManager.h"

void MoleculeManagerClass::Read(IOSectionClass &in)
{
  cerr << "MolMoveMgr read..." << endl;
  int stages = 1;
  in.ReadVar("Stages",stages);
	Array<string,1> methodList;
	assert(in.ReadVar("MoveMethod",methodList));
  assert(methodList.size() == stages);
  Array<int,1> numActions(stages);
  numActions = 1;
  in.ReadVar("NumActions", numActions);
  assert(numActions.size() == stages);
  int startIndex = 0;
  for(int s=0; s<stages; s++){
    string method = methodList(s);
    int actionsToRead = numActions(s);
    cerr << "  Init " << s+1 << " of " << stages << " stages: " << method << endl;
	  if(method == "Translate"){
	  	cerr << "Creating new Translate move...";
    	MoveStage = new MoleculeTranslate(PathData, IOSection, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "ParticleTranslate"){
	  	cerr << "Creating new INTRAmolecular displacement move...";
    	MoveStage = new ParticleTranslate(PathData, IOSection, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "DimerMove"){
	  	cerr << "Creating new DimerMove to change separation of TWO molecules...";
    	MoveStage = new DimerMove(PathData, IOSection, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "Rotate"){
	  	cerr << "Creating new Rotate move...";
    	MoveStage = new MoleculeRotate(PathData, IOSection, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else if (method == "AVB"){
	  	cerr << "Creating new Aggregation Volume Bias move...";
    	MoveStage = new AVBMove(PathData, IOSection, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  	cerr << "See B. Chen and J.I. Siepmann, J. Phys. Chem. B 104, 8725 (2000)" << endl;
	  } else if (method == "Dummy"){
	  	cerr << "Creating new dummy stage to evaluate the specified action (should be preceded by an actual move stage which is evaluate with a \"cheap\" action...";
      assert(s>0); // meant to be the second stage after pre-rejection by a cheap action
    	MoveStage = new DummyEvaluate(PathData, IOSection, actionsToRead, startIndex);
	  	cerr << " done." << endl;
	  } else {
	  	cerr << "ERROR: method " << method << " is not supported." << endl;
	  }
    MoveStage->Read(in);
    Stages.push_back(MoveStage);
    startIndex += numActions(s);
  }
}
