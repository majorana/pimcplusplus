#include "Weight.h"


// Fix to include final link between link M and 0
void WeightClass::Accumulate()
{
  //  cerr<<PathData.InterComm.MyProc()<<": I have been called "<<TimesCalled<<endl;
  TimesCalled++;
  if (TimesCalled % DumpFreq==0)
    WriteBlock();

  if ((TimesCalled % Freq)!=0){
    return;
  }
  NumSamples++;
  //  cerr<<"Accumulating the obseravble"<<endl;
  //  cerr<<"Attempting the coupled permute stage "<<PathData.InterComm.MyProc()<<endl;
  Array<int,1> coupledWeight(1);
  Array<int,1> coupledWeightSend(1);
  coupledWeightSend(0)=PathData.Path.Weight;
  int sendProc;
  int recvProc;
  int numProcs=PathData.InterComm.NumProcs();
  int myProc=PathData.InterComm.MyProc();
  sendProc=(myProc+1) % numProcs;
  recvProc=((myProc-1) + numProcs) % numProcs;
  PathData.InterComm.SendReceive (sendProc, coupledWeightSend,
				  recvProc, coupledWeight);
  //  cerr<<"I am "<<recvProc<<" and my friend is "<<sendProc<<endl;
  //  cerr<<myProc<<": My weight is "<<coupledWeightSend(0)<<" and my friends weight is "<<coupledWeight(0)<<endl;

//   if (abs(coupledWeight(0)*coupledWeightSend(0)-1)<=1e-10){
//     Weight(0)=Weight(0)+1.0;
//   }
//   else {
//     Weight(1)=Weight(1)+1.0;
//   }
  if (coupledWeight(0)==1 && coupledWeightSend(0)==1)
    Weight(0)=Weight(0)+1.0;
  else if (coupledWeight(0)==-1 && coupledWeightSend(0)==-1)
    Weight(1)=Weight(1)+1.0;
  else if (coupledWeight(0)==1 && coupledWeightSend(0)==-1)
    Weight(2)=Weight(2)+1.0;
  else if (coupledWeight(0)==-1 && coupledWeightSend(0)==1)
    Weight(3)=Weight(3)+1.0;

  //  cerr<<myProc<<" "<<Weight(0)/NumSamples<<" "<<Weight(1)/NumSamples<<" "<<NumSamples<<endl;
}

void WeightClass::ShiftData (int NumTimeSlices)
{
  // Do nothing
}

void WeightClass::WriteBlock()
{

  Array<double,1> weight(4);
  for (int counter=0;counter<4;counter++){
    weight(counter)=Weight(counter)/(double)NumSamples;
  }
    //  weight=Weight/(double)NumSamples;
  cerr<<PathData.InterComm.MyProc()<<" "<<weight<<endl;
  Tot.Write(weight);
  if (PathData.Path.Communicator.MyProc()==0)
    IOSection.FlushFile();
  Weight=0.0;
  NumSamples = 0; 


}

void WeightClass::Read(IOSectionClass &in)
{  
  ObservableClass::Read(in);
  assert(in.ReadVar("freq",Freq));
  assert(in.ReadVar("dumpFreq",DumpFreq));
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","CorrelationFunction");
  }
}
