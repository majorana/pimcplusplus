#include "WindingNumber.h"
#include "../Common/MPI/Communication.h"



///This currently tells you the winding number of all the species. It
///might make sense to fix this so it tells only the winding number of
///a specific species
void WindingNumberClass::Accumulate()
{
  TimesCalled++;

  if (TimesCalled % DumpFreq==0){
    WriteBlock();
  }
  if ((TimesCalled % Freq)!=0){
    return;
  }
  TotalDisp=0.0;
  int numLinks=PathData.Path.NumTimeSlices()-1;
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    for (int slice=0;slice<numLinks;slice++) {
      NumSamples++;
      dVec disp;
      disp=PathData.Path.Velocity(slice,slice+1,ptcl);
      TotalDisp(ptcl) =TotalDisp(ptcl)+ disp;
    }
  }
  PathData.Path.Communicator.Sum(TotalDisp,TempDisp);
  if (PathData.Path.Communicator.MyProc()==0) {  
    for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
      for (int dim=0;dim<NDIM;dim++){
	TotalW2[dim] =TotalW2[dim]+TempDisp(ptcl)[dim]*TempDisp(ptcl)[dim];
      }
    }
  }
  //  for (dim=0;dim<NDIM;dim++){
  //    allPtclDisp[dim]=(allPtclDisp[dim]*allPtclDisp[dim]);
  //  }
  //  TotalW2 += allPtclDispl;
}

void WindingNumberClass::Read(IOSectionClass& in)
{
  ObservableClass::Read(in);
  assert(in.ReadVar("freq",Freq));
  assert(in.ReadVar("dumpFreq",DumpFreq));
}
void WindingNumberClass::WriteBlock()
{
//   double totSum;
//   double totNumSamples;
  
//   double myAvg = ESum/(double)NumSamples; //everybody should have the same number of samples for this to be happy
//   double avg = PathData.Path.Communicator.Sum(myAvg);
//   double vavg =PathData.Path.Communicator.Sum(myVAvg);
//   double savg =PathData.Path.Communicator.Sum(mySAvg);
//   double favg =PathData.Path.Communicator.Sum(myFAvg);
//   avg  = avg/(double)PathData.Path.TotalNumSlices;
//   vavg =vavg/(double)PathData.Path.TotalNumSlices;
//   savg =savg/(double)PathData.Path.TotalNumSlices;
//   favg =favg/(double)PathData.Path.TotalNumSlices;
//   // Only processor 0 writes.
  if (PathData.Path.Communicator.MyProc()==0) {
//     cerr << "myAvg = " << myAvg << endl;
//     cerr << "avg = " << avg << endl;
//     cerr << "Pot avg = " << vavg << endl;
//     cerr << "S avg = " << savg << endl;
//     cerr << "U avg = " <<favg <<endl;
    if (FirstTime) {
      FirstTime = false;
      WriteInfo();
      IOSection.WriteVar("Type","Vector");
       Array<double,2> dummy(1,3);
       for (int dim=0;dim<NDIM;dim++)
	 dummy(0,dim)=TotalW2[dim];
       IOSection.WriteVar ("WindingNumber", dummy);
//       dummy(0)=vavg;
//       IOSection.WriteVar ("PotentialEnergy",dummy);
//       dummy(0)=savg;
//       IOSection.WriteVar ("SpringEnergy",dummy);
//       dummy(0)=favg;
//       IOSection.WriteVar ("DBetaEnergy",dummy);
       IOVar = IOSection.GetVarPtr("WindingNumber");
//       IOVVar= IOSection.GetVarPtr("PotentialEnergy");
//       IOSVar= IOSection.GetVarPtr("SpringEnergy");
//       IOUVar= IOSection.GetVarPtr("DBetaEnergy");
    }
    else {
      Array<double,1> dummy(3);
      for (int dim=0;dim<NDIM;dim++)
	dummy(dim)=TotalW2[dim];
      IOVar->Append(dummy);
//       IOVVar->Append(vavg);
//       IOSVar->Append(savg);
//       IOUVar->Append(favg);
      IOSection.FlushFile();
    }
  }
//   ESum = 0.0;
//   VSum = 0.0;
//   SSum = 0.0;
//   FSum=0.0;
  NumSamples = 0;
}

