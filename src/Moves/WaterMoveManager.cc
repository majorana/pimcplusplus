#include "WaterMoveManager.h"

void WaterMoveClass::Read(IOSectionClass &in)
{
	string method;
	assert(in.ReadVar("MoveMethod",method));
	if(method == "Translate"){
		cerr << "Creating new Translate move...";
  	MoveStage = new WaterTranslate(PathData, IOSection);
		cerr << " done." << endl;
	} else if (method == "Rotate"){
		cerr << "Creating new Rotate move...";
  	MoveStage = new WaterRotate(PathData, IOSection);
		cerr << " done." << endl;
	} else if (method == "AVB"){
		cerr << "Creating new Aggregation Volume Bias move...";
  	MoveStage = new AVBMove(PathData, IOSection);
		cerr << " done." << endl;
		cerr << "See B. Chen and J.I. Siepmann, J. Phys. Chem. B 104, 8725 (2000)" << endl;
	} else {
		cerr << "ERROR: method " << method << " is not supported." << endl;
	}
  MoveStage->Read(in);
  Stages.push_back(MoveStage);
}
