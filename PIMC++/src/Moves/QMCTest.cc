#include "QMCTest.h"

void QMCTestMove::AssignPtcl(int mol,Array<int,1>& activeParticles){
  for (int i = 0;i<5;i++)
    activeParticles(i) = mol + PathData.Mol.NumMol()*i;
}

dVec QMCTestMove::Translate(double epsilon)
{
	double x = (PathData.Path.Random.Local()-0.5)*2*epsilon;	
	double y = (PathData.Path.Random.Local()-0.5)*2*epsilon;	
	double z = (PathData.Path.Random.Local()-0.5)*2*epsilon;	
	dVec translate;
	translate[0] = x;
	translate[1] = y;
	translate[2] = z;
	return translate;
}

void QMCTestMove::MakeMove()
{
  double step = 0.3;
  //double step = Step; // Using Step from input file
  int speciesO = PathData.Path.SpeciesNum("O");
  int speciesp = PathData.Path.SpeciesNum("p");
  int speciese = PathData.Path.SpeciesNum("e");
  int  numWater=PathData.Path.Species(speciesO).LastPtcl-
    PathData.Path.Species(speciesO).FirstPtcl+1;
  int choosemol = (int)floor(PathData.Path.Random.Local()*numWater);
  
  // choose a time slice to move
  int numSlices = PathData.Path.TotalNumSlices;
  int slice=0;
  int endSlice = 0;
  int startSlice = 0;
  if(numSlices>1){
    int P_max = numSlices - 1;
    slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
    startSlice = slice-1;
    endSlice = slice+1;
  }

  AssignPtcl(choosemol, ActiveParticles);

  double oldAction = 0.0;
 oldAction += PathData.Actions.MoleculeInteractions.SingleAction(slice,slice,ActiveParticles,0);

  dVec move = Translate(step); 
  double move_mag = sqrt(move(0)*move(0) + move(1)*move(1) + move(2)*move(2));
  for(int i=0;i<ActiveParticles.size();i++){
    dVec old_coord = PathData.Path(slice,ActiveParticles(i));
    dVec new_coord = old_coord + move; 
    PathData.Path.SetPos(slice,ActiveParticles(i),new_coord);
  }
 
  double newAction = 0.0;
 newAction += PathData.Actions.MoleculeInteractions.SingleAction(slice,slice,ActiveParticles,0);
 //if(PathData.UsingQMC){
    //cerr << "calling QMCSampling from " << PathData.Actions.QMCSampling.myProc << endl;
    newAction += PathData.Actions.QMCSampling.SingleAction(slice,slice,ActiveParticles,0);
  //}
  if (-(newAction-oldAction)>=log(PathData.Path.Random.Local())){
    PathData.AcceptMove(startSlice,endSlice,ActiveParticles);
    NumAccepted++;
    total_r_squared += move_mag*move_mag;
  }
  else {
    PathData.RejectMove(startSlice,endSlice,ActiveParticles);
  }
  TimesCalled++;
  if (TimesCalled%10000 == 0){
    cerr << TimesCalled << " moves; current QMCTest ratio is " << double(NumAccepted)/TimesCalled << " with step size " << step << endl;
    double avg = sqrt(total_r_squared)/NumAccepted;
    double diff = total_r_squared/TimesCalled;
    cerr << "QMCTest diffusion value is " << diff << " avg step is " << avg << endl;
  }
}
