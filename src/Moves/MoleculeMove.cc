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

#include "MoleculeMove.h"

void MoleculeRotate::Read(IOSectionClass &moveInput) {
  cerr << "  Rotate read in" << endl;
  //assert(moveInput.ReadVar("SetTheta",Theta));
  assert(moveInput.ReadVar("SetAngle",Theta));
	MolMoveClass::Read(moveInput);
}

void MoleculeTranslate::Read(IOSectionClass &moveInput) {
  cerr << "  Molecule Translate read in" << endl;
  assert(moveInput.ReadVar("SetStep",Step));
	MolMoveClass::Read(moveInput);
}

void DimerMove::Read(IOSectionClass &moveInput) {
  cerr << "  Dimer Move read in" << endl;
  assert(moveInput.ReadVar("SetStep",Step));
	MolMoveClass::Read(moveInput);
  assert(mode == GLOBAL);
}

void ParticleTranslate::Read(IOSectionClass &moveInput) {
  cerr << "  Ptcl Trans read in" << endl;
  assert(moveInput.ReadVar("SetSigma",Sigma));
	MolMoveClass::Read(moveInput);
}

void DummyEvaluate::Read(IOSectionClass &moveInput) {
	MolMoveClass::Read(moveInput);
}

//void MoleculeTranslate::MakeMove() {
double MoleculeTranslate::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles) {
  //double step = 0.3;
  double step = Step; // Using Step from input file

	if(mode == SINGLE){
  	int choosemol = (int)floor(PathData.Path.Random.Local()*PathData.Path.numMol);
		MoveList(0) = choosemol;
		activeParticles.resize(PathData.Path.MolMembers(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = PathData.Path.MolMembers(MoveList(0))(i);
	}
	else if(mode == SEQUENTIAL){
		activeParticles.resize(PathData.Path.MolMembers(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = PathData.Path.MolMembers(MoveList(0))(i);
	}
	else if(mode == GLOBAL){
		activeParticles.resize(PathData.Path.NumParticles());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = i;
	}
//cerr << counter << ", ";
//  int choosemol = counter;
//	counter = (counter+1)%PathData.Path.numMol;
//	activeParticles.resize(MolMembers(choosemol).size());
//	for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = MolMembers(choosemol)(i);
//cerr << "  chose molecule " << choosemol << " of " << PathData.Path.numMol << endl;
//cerr << "  Getting info from " << MolMembers(choosemol) << endl;
//cerr << "  activeParticles is " << activeParticles << endl;

//	activeParticles.resize(PathData.Path.NumParticles());
//	for(int x=0; x<activeParticles.size(); x++) activeParticles(x) = x;

  // choose a time slice to move
  int numSlices = PathData.Path.TotalNumSlices;
	int slice =0;
  slice1 = slice;
  slice2 = slice;
  if(numSlices>1){
    int P_max = numSlices - 1;
    slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
    //slice1 = slice;
    //slice2 = slice;
    slice1 = slice-1;
    slice2 = slice+1;
  }

	double move_mag_sq = 0.0;
	for(int activeMol=0; activeMol<MoveList.size(); activeMol++){
		//cerr << "  Before move: chose slice " << slice << endl;
		//for(int i=0; i<activeParticles.size(); i++) cerr << "  " << activeParticles(i) << ": " << PathData.Path(slice,activeParticles(i)) << endl;
  	dVec move = TranslateMol(slice,PathData.Path.MolMembers(MoveList(activeMol)),step); 
		//cerr << "  After move: of " << move << endl;
		//for(int i=0; i<activeParticles.size(); i++) cerr << "  " << activeParticles(i) << ": " << PathData.Path(slice,activeParticles(i)) << endl;
  	move_mag_sq += move(0)*move(0) + move(1)*move(1) + move(2)*move(2);
	}

  // this block is a hack to output configurations
  // in a format for QBox
	int numMol = PathData.Path.numMol;
	int output1 = 0;
	int output2 = output1 + 1;
	int output3 = output2 + 1;
	if(numMoves == output1){
		ofstream before("H2O.64.rigid.global.1.dat");
		for(int i=0; i<numMol; i++){
			dVec RO = PathData.Path(0,i);
			dVec RH1 = PathData.Path(0,i+numMol);
			dVec RH2 = PathData.Path(0,i+2*numMol);
			before << "O" << " " << RO(0) << " " << RO(1) << " " << RO(2) << endl;
			before << "H" << " " << RH1(0) << " " << RH1(1) << " " << RH1(2) << endl;
			before << "H" << " " << RH2(0) << " " << RH2(1) << " " << RH2(2) << endl;
		}
		before.close();
	}
	else if(numMoves == output2){
		ofstream after("H2O.64.rigid.global.2.dat");
		for(int i=0; i<numMol; i++){
			dVec RO = PathData.Path(0,i);
			dVec RH1 = PathData.Path(0,i+numMol);
			dVec RH2 = PathData.Path(0,i+2*numMol);
			after << "O" << " " << RO(0) << " " << RO(1) << " " << RO(2) << endl;
			after << "H" << " " << RH1(0) << " " << RH1(1) << " " << RH1(2) << endl;
			after << "H" << " " << RH2(0) << " " << RH2(1) << " " << RH2(2) << endl;
		}
		after.close();
	}
	else if(numMoves == output3){
		ofstream third("H2O.64.rigid.global.3.dat");
		for(int i=0; i<numMol; i++){
			dVec RO = PathData.Path(0,i);
			dVec RH1 = PathData.Path(0,i+numMol);
			dVec RH2 = PathData.Path(0,i+2*numMol);
			third << "O" << " " << RO(0) << " " << RO(1) << " " << RO(2) << endl;
			third << "H" << " " << RH1(0) << " " << RH1(1) << " " << RH1(2) << endl;
			third << "H" << " " << RH2(0) << " " << RH2(1) << " " << RH2(2) << endl;
		}
		third.close();
	}

  if (numMoves%10000 == 0 && numMoves>0){
    cerr << numMoves << " moves; current translate ratio is " << double(numAccepted)/numMoves << " with step size " << step << endl;
    double avg = sqrt(move_mag_sq)/numAccepted;
    double diff = move_mag_sq/numMoves;
    cerr << "TRANSLATE diffusion value is " << diff << " avg step is " << avg << endl;
  }
  numMoves++;

//	if(counter == 0)
//		cerr << endl << "I have average action " << UAction << "/" << numMoves << " = " << UAction/numMoves << endl;
	//cerr << "+";

	if(mode == SEQUENTIAL) Advance();

	return 1;
}

double DimerMove::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles) {
  //cerr << " DimerMove::Sample ";
  // load activeParticles; this move is hard-wired to be a GLOBAL update
	activeParticles.resize(PathData.Path.NumParticles());
	for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = i;

  double step = Step; // Using Step from input file
  // this is only written to work for a dimer!
  assert(MoveList.size() == 2);

  int numSlices = PathData.Path.TotalNumSlices;
	int slice =0;
  slice1 = 0;
  slice2 = 0;
  if(numSlices>1){
    int P_max = numSlices - 1;
    slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
    //slice1 = slice;
    //slice2 = slice;
    slice1 = slice-1;
    slice2 = slice+1;
  }

  MoveDimerSeparation(slice,PathData.Path.MolMembers(0), PathData.Path.MolMembers(1),step); 

  if (numMoves%10000 == 0 && numMoves>0){
    cerr << numMoves << " moves; current Dimer Move ratio is " << double(numAccepted)/numMoves << " with step size " << step << " really numAcc is " << numAccepted << endl;
  }
  numMoves++;

	return 1;
}

//void MoleculeRotate::MakeMove()
double MoleculeRotate::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles) {
  double dtheta = M_PI*Theta; // Using Theta from input file
  //cerr << " MoleculeRotate::Sample ";
  //double dtheta = 2*M_PI*0.3;

	if(mode == SINGLE){
  	int choosemol = (int)floor(PathData.Path.Random.Local()*PathData.Path.numMol);
		MoveList(0) = choosemol;
		activeParticles.resize(PathData.Path.MolMembers(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = PathData.Path.MolMembers(MoveList(0))(i);
	}
	else if(mode == SEQUENTIAL){
		activeParticles.resize(PathData.Path.MolMembers(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = PathData.Path.MolMembers(MoveList(0))(i);
	}
	else if(mode == GLOBAL){
		activeParticles.resize(PathData.Path.NumParticles());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = i;
	}

	// choose a time slice to move
	int numSlices = PathData.Path.TotalNumSlices;
	int slice=0;
	slice1 = 0;
	slice2 = 0;
	if(numSlices>1){
	  int P_max = numSlices - 1;
	  slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
    //slice1 = slice;
    //slice2 = slice;
	  slice1 = slice-1;
	  slice2 = slice+1;
	}

	//cerr << "Rotate choosing slice " << slice << " and molecule " << MoveList << " wit members " << PathData.Path.MolMembers(MoveList(0)) << endl;
	for(int activeMol=0; activeMol<MoveList.size(); activeMol++){
		//activeParticles.resize(MolMembers(MoveList(activeMol)).size());
		//for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = MolMembers(MoveList(activeMol))(i);
		
		double theta = 2*(PathData.Path.Random.Local()-0.5)*dtheta;
		RotateMol(slice,PathData.Path.MolMembers(MoveList(activeMol)), theta);

	}

  if (numMoves%10000 == 0 && numMoves>0){
    cerr << numMoves << " moves; current rotate ratio is " << double(numAccepted)/numMoves << " with angle size " << dtheta << endl;
  }
  numMoves++;

	if(mode == SEQUENTIAL) Advance();

	return 1;
}

double ParticleTranslate::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles) {
  numMoves++;
	//cerr << " ParticleTranslate::Sample" << endl;

	if(mode == SINGLE){
  	int choosemol = (int)floor(PathData.Path.Random.Local()*PathData.Path.numMol);
		MoveList(0) = choosemol;
		activeParticles.resize(PathData.Path.MolMembers(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = PathData.Path.MolMembers(MoveList(0))(i);
	}
	else if(mode == SEQUENTIAL){
		activeParticles.resize(PathData.Path.MolMembers(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = PathData.Path.MolMembers(MoveList(0))(i);
	}
	else if(mode == GLOBAL){
		activeParticles.resize(PathData.Path.NumParticles());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = i;
	}

  // choose a time slice to move
  int numSlices = PathData.Path.TotalNumSlices;
	int slice =0;
  slice1 = 0;
  slice2 = 0;
  if(numSlices>1){
    int P_max = numSlices - 1;
    slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
    //slice1 = slice;
    //slice2 = slice;
    slice1 = slice-1;
    slice2 = slice+1;
  }
	//cerr << numMoves << ": ParticleTranslate moving " << activeParticles << endl;
	for(int activeMol=0; activeMol<MoveList.size(); activeMol++){
		for(int activeP = 0; activeP<PathData.Path.MolMembers(MoveList(activeMol)).size(); activeP++){
			int movePtcl = PathData.Path.MolMembers(MoveList(activeMol))(activeP);
      //if(Path.ParticleSpeciesNum(movePtcl) != 0){
			  //cerr << "Moving ptcl " <<  movePtcl << " at slice " << slice << " from " << PathData.Path(slice,movePtcl);
  		TranslatePtcl(slice, movePtcl, Sigma);
			  //cerr << " to " << PathData.Path(slice,movePtcl) << endl;
      //}
		}
	}

  if (numMoves%10000 == 0){
    cerr << numMoves << " moves; current PARTICLE translate ratio is " << double(numAccepted)/numMoves << " with step size " << Sigma << endl;
  }
//	if(counter == 0)
//		cerr << endl << "I have average action " << UAction << "/" << numMoves << " = " << UAction/numMoves << endl;
	//cerr << "+";

	if(mode == SEQUENTIAL) Advance();

	//cerr << "Goodbye Ptcl Translate Sample " << endl;
	return 1;
}

double DummyEvaluate::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles) {
  numMoves++;

  if (numMoves%10000 == 0){
    cerr << numMoves << " moves; current DUMMY accept ratio is " << double(numAccepted)/numMoves << endl;
  }

	return 1;
}
