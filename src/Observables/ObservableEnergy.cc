#include "ObservableEnergy.h"


// Fix to include final link between link M and 0
void EnergyClass::Accumulate()
{
  //Move the join to the end so we don't have to worry about
  //permutations
  PathData.MoveJoin(PathData.NumTimeSlices()-1);
  
  NumSamples++;

  double kinetic, dUShort, dULong, node, vShort, vLong, tip5p;
  PathData.Actions.Energy (kinetic, dUShort, dULong, node, vShort, vLong);

  int slice1 = 0;
  int slice2 = PathData.Path.TotalNumSlices;
  // tip5p = PathData.Actions.TIP5PWater.d_dBeta(slice1,slice2,0);
  
  TotalSum   += kinetic + dUShort + dULong + node;// + tip5p;
  KineticSum += kinetic;
  dUShortSum += dUShort;
  dULongSum  += dULong;
  NodeSum    += node;
  VShortSum  += vShort;
  VLongSum   += vLong;

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

void EnergyClass::WriteBlock()
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

  double kinetic, dUShort, dULong, node, vShort, vLong;
  PathData.Actions.Energy (kinetic, dUShort, dULong, node, vShort, vLong);
  
  TotalSum   += (kinetic + dUShort + dULong + node)/* *PathData.Path.Weight*/;
  KineticSum += kinetic/* * PathData.Path.Weight*/;
  dUShortSum += dUShort/* * PathData.Path.Weight*/;
  dULongSum  += dULong/* * PathData.Path.Weight*/;
  NodeSum    += node/* * PathData.Path.Weight*/;
  VShortSum  += vShort/* * PathData.Path.Weight*/;
  VLongSum   += vLong/* * PathData.Path.Weight*/;
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
