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
#include "../Actions/DavidLongRangeClassYk.h"
#include "../Actions/ShortRangeOn_diagonal_displace_Class.h"

// Fix to include final link between link M and 0
void EnergyClass::Accumulate()
{

  // Move the join to the end so we don't have to worry about permutations
  PathData.MoveJoin(PathData.NumTimeSlices()-1);

  NumSamples++;

  double kinetic, dUShort, dULong, node, vShort, vLong, tip5p, dUNonlocal, residual;
  PathData.Actions.Energy (kinetic, dUShort, dULong, node, vShort, vLong, dUNonlocal,residual);

  if (abs(node) > 1e50) { // Check for odd node energies
    cerr<<"NODE BLOW UP, SETTING NODE ENERGY TO 0.0 :: ENERGIES: "<<kinetic<<" "<<dUShort<<" "<<dULong<<" "<<node<<" "<<vShort<<" "<<vLong<<" "<<dUNonlocal<<endl;
    //abort();
    node = 0.0;
  }

  Array<int,1> changedParticles(PathData.Path.NumParticles());
  for (int i=0; i<changedParticles.size(); i++)
    changedParticles(i) = i;
  int M=PathData.Path.NumTimeSlices();
  int myGetPermNumber = GetPermNumber();
  double localSum = 0.0;
  localSum += kinetic + dUShort + dULong + node + dUNonlocal;
  TotalSum   += localSum;
  KineticSum += kinetic;
  dUShortSum += dUShort;
  dULongSum  += dULong;
  NodeSum    += node;
  VShortSum  += vShort;
  VLongSum   += vLong;
  dUNonlocalSum += dUNonlocal;
  Residual += residual;
  EnergyVals(myGetPermNumber) += localSum;

  int slice1 = 0;
  int slice2 = PathData.Path.NumTimeSlices() - 1;
  for(int n=0; n<numEnergies; n++) {
    double otherE = OtherActions[n]->d_dBeta(slice1,slice2,0);
    OtherSums[n] += otherE;
    localSum += otherE;
    TotalSum += otherE;
  }


  double completeSum=PathData.Path.Communicator.Sum(localSum)/(double)PathData.Path.TotalNumSlices;
  if (DoHist)
    EnergyHistogram.add(PathData.Path.Communicator.Sum(localSum)/(double)PathData.Path.TotalNumSlices,1.0);
}


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
  if (Path.Communicator.MyProc() == 0) {
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
        //  CycleCount(cycleLength)++;
        totalPerms+=cycleLength;
      }
      ptcl++;
    }
  }

  return totalPerms;
}


void EnergyClass::WriteBlock()
{
  int nslices=PathData.Path.TotalNumSlices;
  double norm = 1.0/((double)NumSamples*(double)nslices);
  if (FirstTime){
    FirstTime=false;
    Array<double,1> vtail;
    vtail.resize(PathData.Actions.PairArray.size());
    double longrange_vtail;
    for (int i=0;i<PathData.Actions.PairArray.size();i++)
      vtail(i)=((DavidPAClass*)(PathData.Actions.PairArray(i)))->Vimage;
    if (PathData.Path.DavidLongRange){
      DavidLongRangeClassYk *lr = (DavidLongRangeClassYk*)(&(PathData.Actions.DavidLongRange));
      longrange_vtail=0.5*lr->yk_zero(0)*PathData.Path.NumParticles()/Path.GetVol();
      // cerr << "Writing David Long Range energy" << endl;
    }
    VTailSRVar.Write(vtail);
    VTailLRVar.Write(longrange_vtail);
    if (DoHist) {
      HistStart.Write(EnergyHistogram.startVal);
      HistEnd.Write(EnergyHistogram.endVal);
      NumPoints.Write(EnergyHistogram.NumPoints);
    }

  }
  if (DoHist) {
    for (int i=0;i<EnergyHistogram.histogram.size();i++)
      EnergyHistogram.histogram[i]=EnergyHistogram.histogram[i]*norm*nslices*Prefactor;
    Array<double,1> EnergyHistogramTemp(&(EnergyHistogram.histogram[0]),shape(EnergyHistogram.histogram.size()),neverDeleteData);
    EnergyHistogramVar.Write(EnergyHistogramTemp);
  }
  TotalVar.Write   (Prefactor*PathData.Path.Communicator.Sum(TotalSum)*norm);
  KineticVar.Write (Prefactor*PathData.Path.Communicator.Sum(KineticSum)*norm);
  dUShortVar.Write (Prefactor*PathData.Path.Communicator.Sum(dUShortSum)*norm);
  dULongVar.Write  (Prefactor*PathData.Path.Communicator.Sum(dULongSum)*norm);
  NodeVar.Write    (Prefactor*PathData.Path.Communicator.Sum(NodeSum)*norm);
  VShortVar.Write  (Prefactor*PathData.Path.Communicator.Sum(VShortSum)*norm);
  VLongVar.Write   (Prefactor*PathData.Path.Communicator.Sum(VLongSum)*norm);
  dUNonlocalVar.Write   (Prefactor*PathData.Path.Communicator.Sum(dUNonlocalSum)*norm);
  ResidualVar.Write(Prefactor*PathData.Path.Communicator.Sum(Residual)*norm);
  EnergyVals=EnergyVals*norm;
  EnergyValsVar.Write(EnergyVals);
  EnergyVals=0.0;
  for(int n=0; n<numEnergies; n++){
    OtherVars[n]->Write(Prefactor*PathData.Path.Communicator.Sum(OtherSums[n])*norm);
    OtherSums[n] = 0.0;
  }

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
  NumSamples = 0;
  if (DoHist)
    EnergyHistogram.Clear();
}

void EnergyClass::Read(IOSectionClass &in)
{
  ObservableClass::Read(in);
  if (PathData.Path.Communicator.MyProc()==0) {
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
  EnergyVals.resize(PathData.Path.NumParticles()*2);
  EnergyVals=0.0;
  Array<string,1> EnergyStrings(0);
  if (in.ReadVar("ComputeEnergies", EnergyStrings))
    numEnergies = EnergyStrings.size();
  else
    numEnergies = 0;
  OtherActions.resize(numEnergies);
  OtherVars.resize(numEnergies);
  OtherSums.resize(numEnergies);
  for(int n=0; n<numEnergies; n++) {
    OtherActions[n] = PathData.Actions.GetAction(EnergyStrings(n));
    cerr << "Energy observable added action with label " << EnergyStrings(n) << endl;
    OtherVars[n] = new ObservableDouble(EnergyStrings(n), IOSection, PathData.Path.Communicator);
    OtherSums[n] = 0.0;
  }
  double histStart=0.0;
  double histEnd=1.0;
  int histPoints=0;
  in.ReadVar("HistStart",histStart);
  in.ReadVar("HistEnd",histEnd);
  in.ReadVar("HistPoints",histPoints);
  if(!histPoints) {
    cerr << PathData.Path.Communicator.MyProc() << " Not Doing Energy Histogram" << endl;
    DoHist = false;
  }
  if(DoHist) {
    EnergyHistogram.Init(histPoints,histStart,histEnd);
    EnergyHistogramSum.resize(EnergyHistogram.histogram.size());
  }
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
  double FullWeight;
  double currWeight=PathData.Path.Weight;
  //PathData.Path.Communicator.GatherProd(currWeight,FullWeight,0);
  FullWeight = 1;
  NumSamples++;

  double kinetic, dUShort, dULong, node, vShort, vLong, dUNonlocal;
  PathData.Actions.Energy (kinetic, dUShort, dULong, node, vShort, vLong,
			   dUNonlocal);
  //  cerr<<"ENERGIES: "<<kinetic<<" "<<duShort<<" "<<dULong<<" "<<node<<" "<<vShort<<" "<<vLong<<" "<<dUNonlocal<<endl;
  TotalSum   += (kinetic + dUShort + dULong + node)*FullWeight;
  KineticSum += kinetic*FullWeight;/* * PathData.Path.Weight*/;
  dUShortSum += dUShort*FullWeight;/* * PathData.Path.Weight*/;
  dULongSum  += dULong*FullWeight;/* * PathData.Path.Weight*/;
  NodeSum    += node*FullWeight;/* * PathData.Path.Weight*/;
  VShortSum  += vShort*FullWeight;/* * PathData.Path.Weight*/;
  VLongSum   += vLong*FullWeight;/* * PathData.Path.Weight*/;
  dUNonlocalSum += dUNonlocal*FullWeight;
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
  assert(in.ReadVar("Frequency",Freq));
  //  assert(in.ReadVar("dumpFreq",DumpFreq));
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
}
