#include "Time.h"


// Fix to include final link between link M and 0
void MCTimeClass::Accumulate()
{

  TimesCalled++;
  if (TimesCalled % DumpFreq==0)
    WriteBlock();

  if ((TimesCalled % Freq)!=0){
    return;
  }
  
}


void MCTimeClass::WriteBlock()
{
  if (FirstTime){
    FirstTime=false;
    MoveTime.resize(Moves.size());
    ObservableTime.resize(Observables.size());
    MoveTime=0;
    ObservableTime=0;
    TotalTime=0;
    cerr<<"My moves size is "<<MoveTime.size()<<endl;
    cerr<<"My observable size is "<<ObservableTime.size()<<endl;
  }
  //  cerr<<"This is being called"<<endl;
  TotalTime+=(double)(clock()-StartTime)/(double)CLOCKS_PER_SEC;
  StartTime=clock();
  TotalTimeVar.Write(TotalTime);
//   double moveNorm=0;
//   double observeNorm=0;
//   for (int counter=0;counter<Moves.size();counter++){
//     moveNorm+=Moves(counter)->SecondsInMove;
//   }
//   for (int counter=0;counter<Observables.size();counter++){
//     observeNorm+=Observables(counter)->SecondsInObservable;
//  }
  for (int counter=0;counter<Moves.size();counter++){
    MoveTime(counter)=(Moves(counter)->TimeSpent)/TotalTime;
  }
  for (int counter=0;counter<Observables.size();counter++){
    ObservableTime(counter)=(Observables(counter)->TimeSpent)/TotalTime;
  }

  MoveTimeVar.Write(MoveTime);
  ObservableTimeVar.Write(ObservableTime);


  if (PathData.Path.Communicator.MyProc()==0)
    IOSection.FlushFile();
  
}

void MCTimeClass::Read(IOSectionClass &in)
{  

  ObservableClass::Read(in);
  cerr<<"Hello I'm reading!"<<endl;
  assert(in.ReadVar("freq",Freq));
  assert(in.ReadVar("dumpFreq",DumpFreq));
  cerr<<"Dump freq is "<<DumpFreq<<endl;
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
  StartTime=0;
}



void MCTimeClass::WriteInfo()
{
//   Array<string,1> movesString;
//   Array<string,1> observeString;
//   movesString.resize(Moves.size());
//   observeString.resize(Observables.size());
//   for (int counter=0;counter<Moves.size();counter++){
//     movesString(counter)=Moves(counter)->Name;
//   }
//   for (int counter=0;counter<Observables.size();counter++){
//     observeString(counter)=Observables(counter)->Name;
//   }
  //  IOSection.WriteVar("Move Names",movesString);
  //  IOSection.WriteVar("Observable Names",observeString);
  
}
