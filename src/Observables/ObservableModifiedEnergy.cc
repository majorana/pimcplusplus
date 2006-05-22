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

#include "ObservableModifiedEnergy.h"


// Fix to include final link between link M and 0
void ModifiedEnergyClass::Accumulate()
{
  //Move the join to the end so we don't have to worry about permutations
  PathData.MoveJoin(PathData.NumTimeSlices()-1);
  
  NumSamples++;

  double kinetic, dUShort, dULong, node, vShort, vLong, tip5p, st2, rotkin, p2rotkin;
  PathData.Actions.Energy (kinetic, dUShort, dULong, node, vShort, vLong);

  int slice1 = 0;
  int slice2 = PathData.Path.NumTimeSlices()-1;
  tip5p = PathData.Actions.TIP5PWater.d_dBeta(slice1,slice2,0);
  st2 = PathData.Actions.ST2Water.d_dBeta(slice1,slice2,0);
  p2rotkin = PathData.Actions.TIP5PWater.NewRotKinEnergy(slice1,slice2,0); // Calculates angles; should return 3/2 kT
  //rotkin = PathData.Actions.TIP5PWater.RotationalEnergy(slice1,slice2,0);  // Calculates phi and theta; should return kT
  rotkin = PathData.Actions.TIP5PWater.FixedAxisEnergy(slice1,slice2,0); //Fixed Axis
//  p2rotkin = PathData.Actions.TIP5PWater.SecondProtonKineticEnergy(slice1,slice2,0);

  TotalSum   += kinetic + dUShort + dULong + node + st2;// +tip5p + rotkin;
  KineticSum += kinetic;
  dUShortSum += dUShort;
  dULongSum  += dULong;
  NodeSum    += node;
  VShortSum  += vShort;
  VLongSum   += vLong;
  TIP5PSum   += tip5p;
  ST2Sum   += st2;
  RotKinSum  += rotkin;
  P2RotKinSum += p2rotkin;
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


void ModifiedEnergyClass::ShiftData (int NumTimeSlices)
{
  // Do nothing
}

void ModifiedEnergyClass::WriteBlock()
{
  int nslices=PathData.Path.TotalNumSlices;
cerr << "okay nslices is " << nslices << endl;
  double norm = 1.0/((double)NumSamples*(double)nslices);

// NORMALIZE BY NUMBER OF MOLECULES *******
  KineticSum = KineticSum/PathData.Path.numMol;
  TIP5PSum = TIP5PSum/PathData.Path.numMol;
  ST2Sum = ST2Sum/PathData.Path.numMol;
//  RotKinSum = RotKinSum/PathData.Path.numMol;
//  P2RotKinSum = P2RotKinSum/PathData.Path.numMol;
  TotalSum = TotalSum/PathData.Path.numMol;// + RotKinSum;
  
  TotalVar.Write(PathData.Path.Communicator.Sum(TotalSum)*norm);
  KineticVar.Write(PathData.Path.Communicator.Sum(KineticSum)*norm);
//  dUShortVar.Write(PathData.Path.Communicator.Sum(dUShortSum)*norm);
//  dULongVar.Write(PathData.Path.Communicator.Sum(dULongSum)*norm);
//  NodeVar.Write(PathData.Path.Communicator.Sum(NodeSum)*norm);
//  VShortVar.Write(PathData.Path.Communicator.Sum(VShortSum)*norm);
//  VLongVar.Write(PathData.Path.Communicator.Sum(VLongSum)*norm);
  TIP5PVar.Write(PathData.Path.Communicator.Sum(TIP5PSum)*norm);
  ST2Var.Write(PathData.Path.Communicator.Sum(ST2Sum)*norm);
  RotKinVar.Write(PathData.Path.Communicator.Sum(RotKinSum)*norm);
  P2RotKinVar.Write(PathData.Path.Communicator.Sum(P2RotKinSum)*norm);

  cerr << "Total " << TotalSum*norm << " per molecule." << endl;
  cerr << "Kinetic " << KineticSum*norm << " per molecule." << endl;
  cerr << "TIP5P " << TIP5PSum*norm << " per molecule." << endl;
  cerr << "ST2 " << ST2Sum*norm << " per molecule." << endl;
  cerr << "Rotational Kinetic " << RotKinSum*norm << " per molecule." << endl;
  cerr << "Second Proton Rotational Kinetic " << P2RotKinSum*norm << " per molecule." << endl;
  
  TotalSum   = 0.0;
  KineticSum = 0.0;
  dUShortSum = 0.0;
  dULongSum  = 0.0;
  NodeSum    = 0.0;
  VShortSum  = 0.0;
  VLongSum   = 0.0;
  TIP5PSum   = 0.0;
  ST2Sum   = 0.0;
  RotKinSum   = 0.0;
  P2RotKinSum   = 0.0;
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

void ModifiedEnergyClass::Read(IOSectionClass &in)
{  
  ObservableClass::Read(in);
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
}


////Code for energy sign class


// Fix to include final link between link M and 0
/*
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

  double kinetic, dUShort, dULong, node, vShort, vLong;
  PathData.Actions.Energy (kinetic, dUShort, dULong, node, vShort, vLong);
  
  TotalSum   += (kinetic + dUShort + dULong + node)*PathData.Path.Weight;
  KineticSum += kinetic * PathData.Path.Weight;
  dUShortSum += dUShort * PathData.Path.Weight;
  dULongSum  += dULong * PathData.Path.Weight;
  NodeSum    += node * PathData.Path.Weight;
  VShortSum  += vShort * PathData.Path.Weight;
  VLongSum   += vLong * PathData.Path.Weight;
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

  if (PathData.Path.Communicator.MyProc()==0)
    IOSection.FlushFile();
  
  TotalSum   = 0.0;
  KineticSum = 0.0;
  dUShortSum = 0.0;
  dULongSum  = 0.0;
  NodeSum    = 0.0;
  VShortSum  = 0.0;
  VLongSum   = 0.0;
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
*/
