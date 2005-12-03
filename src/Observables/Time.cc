#include "Time.h"



void MCTimeClass::Accumulate()
{

  
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
  TotalTime+=(double)(clock()-StartTime)/(double)CLOCKS_PER_SEC;
  StartTime=clock();
  TotalTimeVar.Write(TotalTime);
  list<MoveClass*>::iterator moveIter;
  int counter=0;
  for (moveIter=Moves.begin();moveIter!=Moves.end();moveIter++){
    MoveTime(counter)=((*moveIter)->TimeSpent)/TotalTime;
    counter++;
  }
  counter = 0;
  list<ObservableClass*>::iterator observeIter; 
  for (observeIter=Observables.begin();observeIter!=Observables.end();observeIter++){
    ObservableTime(counter)=((*observeIter)->TimeSpent)/TotalTime;
    counter++;
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
