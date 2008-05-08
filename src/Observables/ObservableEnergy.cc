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

#include "ObservableEnergy.h"

// These are included for the new runtime 
// specification of energy observables 
// to compute; see below -John
#include "../Actions/MoleculeInteractionsClass.h"
#include "../Actions/EAMClass.h"
#include "../Actions/QBoxAction.h"
#include "../Actions/ST2WaterClass.h"
#include "../Actions/QMCSamplingClass.h"

// Fix to include final link between link M and 0
void EnergyClass::Accumulate()
{

  //Move the join to the end so we don't have to worry about
  //permutations
  PathData.MoveJoin(PathData.NumTimeSlices()-1);

  NumSamples++;

  	//map<double> Energies;
  double kinetic, dUShort, dULong, node, vShort, vLong, tip5p, dUNonlocal,
    residual;
  PathData.Actions.Energy (kinetic, dUShort, dULong, node, vShort, vLong,
			   dUNonlocal,residual);
  //PathData.Actions.Energy(Energies);
	//kinetic = Energies["kinetic"];
	//dUShort = Energies["dUShort"];
	//dULong = Energies["dULong"];
	//node = Energies["node"];
	//vShort = Energies["vShort"];
	//vLong = Energies["vLong"];
  int myGetPermNumber=GetPermNumber();


  TotalSum   += kinetic + dUShort + dULong + node + dUNonlocal;// + tip5p;
  KineticSum += kinetic;
  dUShortSum += dUShort;
  dULongSum  += dULong;
  NodeSum    += node;
  VShortSum  += vShort;
  VLongSum   += vLong;
  dUNonlocalSum += dUNonlocal;
  Residual += residual;
  //  cerr<<"My get perm number is "<<myGetPermNumber<<endl;
  EnergyVals(myGetPermNumber)+=kinetic;
  //  cerr<<"ENERGY VVLAS"<<EnergyVals<<"ASDF "<<kinetic<<endl;

  int slice1 = 0;
  int slice2 = PathData.Path.NumTimeSlices() - 1;
	for(int n=0; n<numEnergies; n++){
		//cerr << "ObsEnergy computing dBeta";
		double otherE = OtherActions[n]->d_dBeta(slice1,slice2,0);
		OtherSums[n] += otherE;
		TotalSum += otherE;
		//cerr << "  finished." << endl;
	}

//   double kAction, uShortAction, uLongAction, nodeAction;
//   PathData.Actions.GetActions(kAction, uShortAction, uLongAction, nodeAction);
//   double totalAction = (kAction + uShortAction + uLongAction + nodeAction);
//   TotalActionSum += totalAction;
//   PathData.Path.Communicator.Sum (totalAction);
//   ExpTotalActionSum += exp(-totalAction+4.8879e4);
/// Removing this total action stuff for now.
//   double kAction, uShortAction, uLongAction, nodeAction;
//   PathData.Actions.GetActions(kAction, uShortAction, uLongAction, nodeAction);
//   TotalActionSum += (kAction + uShortAction + uLongAction + nodeAction);

  //  TIP5PSum   += tip5p;
	//cerr << "ObservableEnergy leaving Accumulate" << endl;
}


// // Fix to include final link between link M and 0
// void EnergyClass::Accumulate()
// {
//   TimesCalled++;
//   if (TimesCalled % DumpFreq==0)
//     WriteBlock();

//   if ((TimesCalled % Freq)!=0){
//     return;
//   }
//   //Move the join to the end so we don't have to worry about permutations
//   PathData.MoveJoin(PathData.NumTimeSlices()-1);
//   // Loop over all links
//   int numLinks = PathData.NumTimeSlices()-1; 
  
//   NumSamples++;

//   /// CHECK code
//   double Echeck = 0.0;
//   double spring, dU, V = 0.0;
//   spring = dU = 0.0;
//   for (int slice=0; slice<numLinks; slice++) {
//     double sp, du,v;
//     PathData.Action.Energy(slice, 0, sp, du);
//     v=PathData.Action.PotentialEnergy(slice);
//     dU += du; //*PathData.Path.Weight;
//     spring += sp; //*PathData.Path.Weight;
//     V += v; //*PathData.Path.Weight;
//   }

//   double node = 0.0;
//   for (int species=0; species<PathData.Path.NumSpecies(); species++)
//     if (PathData.Actions.NodalActions(species) != NULL)
//       node += PathData.Actions.NodalActions(species)->d_dBeta(0,numLinks,0);
//   Echeck = spring + dU + node;
//   ESum += Echeck;
//   VSum += V;
//   SSum += spring;
//   FSum += dU;
//   NodeSum += node;
// }


void EnergyClass::ShiftData (int NumTimeSlices)
{
  // Do nothing
}

int EnergyClass::GetPermNumber()
{
  int totalPerms=0;
  PathClass &Path= PathData.Path;
  int N = PathData.Path.NumParticles();
  if (CountedAlready.size() != N) {
    CountedAlready.resize(N);
    TotalPerm.resize(N);
  }
  PathData.Path.TotalPermutation (TotalPerm);
  CountedAlready =false;
  int ptcl=0;
  /// Only proc 0 gets TotalPerm
  if (Path.Communicator.MyProc() == 0) 
    while (ptcl < N) {
      if (!CountedAlready(ptcl)) {
	int startPtcl=ptcl;
	int roamingPtcl=ptcl;
	int cycleLength=0;
	roamingPtcl = TotalPerm(roamingPtcl);
	while (roamingPtcl!=startPtcl){
	  CountedAlready(roamingPtcl)=true;
	  cycleLength++;
	  roamingPtcl=TotalPerm(roamingPtcl);
	}
	//	CycleCount(cycleLength)++;
	totalPerms+=cycleLength;
      }
      ptcl++;
    }
  return totalPerms;



}
void EnergyClass::WriteBlock()
{
  int nslices=PathData.Path.TotalNumSlices;
  double norm = 1.0/((double)NumSamples*(double)nslices);
  
  TotalVar.Write   (Prefactor*PathData.Path.Communicator.Sum(TotalSum)*norm);
  KineticVar.Write (Prefactor*PathData.Path.Communicator.Sum(KineticSum)*norm);
  dUShortVar.Write (Prefactor*PathData.Path.Communicator.Sum(dUShortSum)*norm);
  dULongVar.Write  (Prefactor*PathData.Path.Communicator.Sum(dULongSum)*norm);
  NodeVar.Write    (Prefactor*PathData.Path.Communicator.Sum(NodeSum)*norm);
  VShortVar.Write  (Prefactor*PathData.Path.Communicator.Sum(VShortSum)*norm);
  VLongVar.Write   (Prefactor*PathData.Path.Communicator.Sum(VLongSum)*norm);
  dUNonlocalVar.Write   (Prefactor*PathData.Path.Communicator.Sum(dUNonlocalSum)*norm);
  ResidualVar.Write(Prefactor*PathData.Path.Communicator.Sum(Residual)*norm);
  //  cerr<<"norm is "<<norm<<endl;
  //  cerr<<"norm"<<EnergyVals<<endl;
  EnergyVals=EnergyVals*norm;
  //  cerr<<"norm"<<EnergyVals<<endl;
  EnergyValsVar.Write(EnergyVals);
  EnergyVals=0.0;
	for(int n=0; n<numEnergies; n++){
		OtherVars[n]->Write(Prefactor*PathData.Path.Communicator.Sum(OtherSums[n])*norm);
		OtherSums[n] = 0.0;
	}
  //  TotalActionVar.Write(PathData.Path.Communicator.Sum(TotalActionSum)/(double)(NumSamples));
  //  ExpTotalActionVar.Write(ExpTotalActionSum/(double)NumSamples);
  //  TIP5PVar.Write(PathData.Path.Communicator.Sum(TIP5PSum)*norm);

  //  cerr << "Total " << TotalSum*norm << " and per molecule: " << TotalSum*norm/PathData.Path.numMol << endl;
  //  cerr << "Kinetic " << KineticSum*norm << " and per molecule: " << KineticSum*norm/PathData.Path.numMol << endl;
  //  cerr << "TIP5P " << TIP5PSum*norm << " and per molecule: " << TIP5PSum*norm/PathData.Path.numMol << endl;
  
  TotalSum       = 0.0;
  KineticSum     = 0.0;
  dUShortSum     = 0.0;
  dULongSum      = 0.0;
  NodeSum        = 0.0;
  VShortSum      = 0.0;
  VLongSum       = 0.0;
  dUNonlocalSum = 0.0;
  Residual=0.0;
  EnergyVals=0.0;
//   TotalActionSum = 0.0;
//   ExpTotalActionSum = 0.0;
  //  TIP5PSum   = 0.0;
  NumSamples = 0;
}

//   double totSum;
//   double totNumSamples;
  
//   double myAvg = ESum/(double)NumSamples; //everybody should have the same number of samples for this to be happy
//   double myVAvg= VSum/(double)NumSamples;
//   double mySAvg= SSum/(double)NumSamples;
//   double myFAvg= FSum/(double)NumSamples;
//   double myNodeAvg = NodeSum/(double)NumSamples;
//   double avg = PathData.Path.Communicator.Sum(myAvg);
//   double vavg =PathData.Path.Communicator.Sum(myVAvg);
//   double savg =PathData.Path.Communicator.Sum(mySAvg);
//   double favg =PathData.Path.Communicator.Sum(myFAvg);
//   double NodeAvg =PathData.Path.Communicator.Sum(myNodeAvg);
//   avg  = avg/(double)PathData.Path.TotalNumSlices;
//   vavg =vavg/(double)PathData.Path.TotalNumSlices;
//   savg =savg/(double)PathData.Path.TotalNumSlices;
//   favg =favg/(double)PathData.Path.TotalNumSlices;
//   NodeAvg = NodeAvg/(double)(PathData.Path.TotalNumSlices);
//   // Only processor 0 writes.
//   if (PathData.Path.Communicator.MyProc()==0) {
//     cerr << "myAvg = " << myAvg << endl;
//     cerr << "avg = " << avg << endl;
//     cerr << "Pot avg = " << vavg << endl;
//     cerr << "S avg = " << savg << endl;
//     cerr << "U avg = " <<favg <<endl;
//     cerr << "NodeAvg = " <<NodeAvg <<endl;
//     if (FirstTime) {
//       FirstTime = false;
//       WriteInfo();
//       IOSection.WriteVar("Type","Scalar");
//       Array<double,1> dummy(1);
//       dummy(0)=avg;
//       IOSection.WriteVar ("TotalEnergy", dummy);
//       dummy(0)=vavg;
//       IOSection.WriteVar ("PotentialEnergy",dummy);
//       dummy(0)=savg;
//       IOSection.WriteVar ("SpringEnergy",dummy);
//       dummy(0)=favg;
//       IOSection.WriteVar ("DBetaEnergy",dummy);
//       dummy(0)=NodeAvg;
//       IOSection.WriteVar ("NodeEnergy",dummy);
//       IOVar = IOSection.GetVarPtr("TotalEnergy");
//       IOVVar= IOSection.GetVarPtr("PotentialEnergy");
//       IOSVar= IOSection.GetVarPtr("SpringEnergy");
//       IOUVar= IOSection.GetVarPtr("DBetaEnergy");
//       IONodeVar= IOSection.GetVarPtr("NodeEnergy");
//     }
//     else {
//       IOVar->Append(avg);
//       IOVVar->Append(vavg);
//       IOSVar->Append(savg);
//       IOUVar->Append(favg);
//       IONodeVar->Append(NodeAvg);
//       IOSection.FlushFile();
//     }
//   }
//   ESum       = 0.0;
//   VSum       = 0.0;
//   SSum       = 0.0;
//   FSum       = 0.0;
//   NodeSum    = 0.0;
//  NumSamples = 0;
//}

void EnergyClass::Read(IOSectionClass &in)
{  
  ObservableClass::Read(in);
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
  EnergyVals.resize(PathData.Path.NumParticles()*2);
  EnergyVals=0.0;
	// New code to read in and accumulate
	// energies from specified action objects
	// Added by John on June 16 2006
  Array<string,1> EnergyStrings(0);
  in.ReadVar("ComputeEnergies", EnergyStrings);
  numEnergies = EnergyStrings.size();
	OtherActions.resize(numEnergies);
	OtherVars.resize(numEnergies);
	OtherSums.resize(numEnergies);	
	for(int n=0; n<numEnergies; n++){
    OtherActions[n] = PathData.Actions.GetAction(EnergyStrings(n));
    cerr << "Energy observable added action with label " << EnergyStrings(n) << endl;

    // old action read in is now deprecated
		//if(EnergyStrings(n) == "MoleculeInteractions"){
		//	OtherActions[n] = &PathData.Actions.MoleculeInteractions;
		//	// read should be done in action
		//	//PathData.Actions.MoleculeInteractions.Read(in);
		//} else if(EnergyStrings(n) == "ST2WaterClass"){
		//	OtherActions[n] = &PathData.Actions.ST2Water;
		//} else if(EnergyStrings(n) == "EAMClass"){
		//	OtherActions[n] = &PathData.Actions.EAM;
		//} else if(EnergyStrings(n) == "KineticRotorActionClass"){
		//	OtherActions[n] = &PathData.Actions.KineticRotor;
		//} else if(EnergyStrings(n) == "FixedAxisRotorActionClass"){
		//	OtherActions[n] = &PathData.Actions.FixedAxisRotor;
		//} else if(EnergyStrings(n) == "QBoxActionClass"){
		//	OtherActions[n] = &PathData.Actions.QBoxAction;
		//} else if(EnergyStrings(n) == "QMCSamplingClass"){
		//	OtherActions[n] = &PathData.Actions.QMCSampling;
		//} else if(EnergyStrings(n) == "IonIonActionClass"){
		//	OtherActions[n] = &PathData.Actions.IonInteraction;
		//} else if(EnergyStrings(n) == "LongRangeCoulomb"){
		//	OtherActions[n] = &PathData.Actions.LongRangeCoulomb;
		//}
		//// Other action objects can be specified here of course
//#ifdef USE_QMC
		//else if(EnergyStrings(n) == "CEIMCActionClass"){
  	//	OtherActions[n] = &PathData.Actions.CEIMCAction;
    //  PathData.Actions.CEIMCAction.Read(IOSection);
		//}
//#endif
		//else {
		//	cerr << "You specified " << EnergyStrings(n) << ", which is not supported for runtime inclusion as a computed energy observable." << endl;
		//}
		OtherVars[n] = new ObservableDouble(EnergyStrings(n), IOSection, PathData.Path.Communicator);
		OtherSums[n] = 0.0;
	}
	// End John's block of code

}


////Code for energy sign class


// Fix to include final link between link M and 0
void EnergySignClass::Accumulate()
{
  TimesCalled++;
  if (TimesCalled % DumpFreq==0)
    WriteBlock();

  if ((TimesCalled % Freq)!=0){
    return;
  }
  //Move the join to the end so we don't have to worry about permutations
  PathData.MoveJoin(PathData.NumTimeSlices()-1);
  
  NumSamples++;

  double kinetic, dUShort, dULong, node, vShort, vLong, dUNonlocal;
  PathData.Actions.Energy (kinetic, dUShort, dULong, node, vShort, vLong,
			   dUNonlocal);
  
  TotalSum   += (kinetic + dUShort + dULong + node)/* *PathData.Path.Weight*/;
  KineticSum += kinetic/* * PathData.Path.Weight*/;
  dUShortSum += dUShort/* * PathData.Path.Weight*/;
  dULongSum  += dULong/* * PathData.Path.Weight*/;
  NodeSum    += node/* * PathData.Path.Weight*/;
  VShortSum  += vShort/* * PathData.Path.Weight*/;
  VLongSum   += vLong/* * PathData.Path.Weight*/;
  dUNonlocalSum += dUNonlocal;
}



void EnergySignClass::ShiftData (int NumTimeSlices)
{
  // Do nothing
}

void EnergySignClass::WriteBlock()
{
  int nslices=PathData.Path.TotalNumSlices;
  double norm = 1.0/((double)NumSamples*(double)nslices);

  TotalVar.Write(PathData.Path.Communicator.Sum(TotalSum)*norm);
  KineticVar.Write(PathData.Path.Communicator.Sum(KineticSum)*norm);
  dUShortVar.Write(PathData.Path.Communicator.Sum(dUShortSum)*norm);
  dULongVar.Write(PathData.Path.Communicator.Sum(dULongSum)*norm);
  NodeVar.Write(PathData.Path.Communicator.Sum(NodeSum)*norm);
  VShortVar.Write(PathData.Path.Communicator.Sum(VShortSum)*norm);
  VLongVar.Write(PathData.Path.Communicator.Sum(VLongSum)*norm);
  dUNonlocalVar.Write(PathData.Path.Communicator.Sum(dUNonlocalSum)*norm);

  if (PathData.Path.Communicator.MyProc()==0)
    IOSection.FlushFile();
  
  TotalSum   = 0.0;
  KineticSum = 0.0;
  dUShortSum = 0.0;
  dULongSum  = 0.0;
  NodeSum    = 0.0;
  VShortSum  = 0.0;
  VLongSum   = 0.0;
  dUNonlocalSum = 0.0;
  NumSamples = 0;
}


void EnergySignClass::Read(IOSectionClass &in)
{  
  ObservableClass::Read(in);
  assert(in.ReadVar("freq",Freq));
  assert(in.ReadVar("dumpFreq",DumpFreq));
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
}
