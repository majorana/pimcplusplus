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
  //assert(moveInput.ReadVar("SetTheta",Theta));
  assert(moveInput.ReadVar("SetAngle",Theta));
	MolMoveClass::Read(moveInput);
}

void MoleculeTranslate::Read(IOSectionClass &moveInput) {
  assert(moveInput.ReadVar("SetStep",Step));
	MolMoveClass::Read(moveInput);
}

//void MoleculeTranslate::MakeMove() {
double MoleculeTranslate::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles) {
	//cerr << " MoleculeTranslate::Sample ";
  //double step = 0.3;
  double step = Step; // Using Step from input file

	if(mode == SINGLE){
  	int choosemol = (int)floor(PathData.Path.Random.Local()*PathData.Path.numMol);
		MoveList(0) = choosemol;
		activeParticles.resize(MolMembers(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = MolMembers(MoveList(0))(i);
	}
	else if(mode == SEQUENTIAL){
		activeParticles.resize(MolMembers(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = MolMembers(MoveList(0))(i);
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
  slice1 = 0;
  slice2 = 0;
  if(numSlices>1){
    int P_max = numSlices - 1;
    slice = (int)floor(P_max*PathData.Path.Random.Local()) + 1;
    slice1 = slice-1;
    slice2 = slice+1;
  }

	double move_mag_sq = 0.0;
	for(int activeMol=0; activeMol<MoveList.size(); activeMol++){
		//cerr << "  Before move: chose slice " << slice << endl;
		//for(int i=0; i<activeParticles.size(); i++) cerr << "  " << activeParticles(i) << ": " << PathData.Path(slice,activeParticles(i)) << endl;
  	dVec move = TranslateMol(slice,MolMembers(MoveList(activeMol)),step); 
		//cerr << "  After move: of " << move << endl;
		//for(int i=0; i<activeParticles.size(); i++) cerr << "  " << activeParticles(i) << ": " << PathData.Path(slice,activeParticles(i)) << endl;
  	move_mag_sq += move(0)*move(0) + move(1)*move(1) + move(2)*move(2);
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

//void MoleculeRotate::MakeMove()
double MoleculeRotate::Sample(int &slice1,int &slice2, Array<int,1> &activeParticles) {
  double dtheta = 2*M_PI*Theta; // Using Theta from input file
//cerr << " MoleculeRotate::Sample ";
  //double dtheta = 2*M_PI*0.3;

	if(mode == SINGLE){
  	int choosemol = (int)floor(PathData.Path.Random.Local()*PathData.Path.numMol);
		MoveList(0) = choosemol;
		activeParticles.resize(MolMembers(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = MolMembers(MoveList(0))(i);
	}
	else if(mode == SEQUENTIAL){
		activeParticles.resize(MolMembers(MoveList(0)).size());
		for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = MolMembers(MoveList(0))(i);
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
	  slice1 = slice-1;
	  slice2 = slice+1;
	}

	//cerr << "Rotate choosing slice " << slice << " and molecule " << MoveList << " wit members " << MolMembers(MoveList(0)) << endl;
	for(int activeMol=0; activeMol<MoveList.size(); activeMol++){
		//activeParticles.resize(MolMembers(MoveList(activeMol)).size());
		//for(int i=0; i<activeParticles.size(); i++) activeParticles(i) = MolMembers(MoveList(activeMol))(i);
		
		double theta = 2*(PathData.Path.Random.Local()-0.5)*dtheta;
		RotateMol(slice,MolMembers(MoveList(activeMol)),theta);

	}

  if (numMoves%10000 == 0 && numMoves>0){
    cerr << numMoves << " moves; current rotate ratio is " << double(numAccepted)/numMoves << " with angle size " << dtheta << endl;
  }
  numMoves++;

	if(mode == SEQUENTIAL) Advance();

	return 1;
}
