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
    Array<string,1> moveNames       (Moves.size());
    Array<string,1> observableNames (Observables.size());
    list<MoveClass*>::iterator moveIter;
    int i=0;
    for (moveIter=Moves.begin();moveIter!=Moves.end();moveIter++) {
      moveNames(i) = ((*moveIter)->Name);
      i++;
    }
    i = 0;
    list<ObservableClass*>::iterator observableIter;
    for (observableIter=Observables.begin();
	 observableIter!=Observables.end();observableIter++) {
      observableNames(i) = ((*observableIter)->Name);
      i++;
    }
    if (PathData.Path.Communicator.MyProc()==0) {
      IOSection.WriteVar("MoveNames", moveNames);
      IOSection.WriteVar("ObserveNames", observableNames);
    }
  }
  TotalTime+=(double)(clock()-StartTime)/(double)CLOCKS_PER_SEC;
  StartTime=clock();
  TotalTimeVar.Write(TotalTime);
  list<MoveClass*>::iterator moveIter;
  int i=0;
  for (moveIter=Moves.begin();moveIter!=Moves.end();moveIter++){
    MoveTime(i)=((*moveIter)->TimeSpent)/TotalTime;
    i++;
  }
  i = 0;
  list<ObservableClass*>::iterator observeIter; 
  for (observeIter=Observables.begin();
       observeIter!=Observables.end();observeIter++){
    ObservableTime(i)=((*observeIter)->TimeSpent)/TotalTime;
    i++;
  }
  
  MoveTimeVar.Write(MoveTime);
  ObservableTimeVar.Write(ObservableTime);


  if (PathData.Path.Communicator.MyProc()==0)
    IOSection.FlushFile();
  
}

void MCTimeClass::Read(IOSectionClass &in)
{  

  ObservableClass::Read(in);
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
//   for (int i=0;i<Moves.size();i++){
//     movesString(i)=Moves(i)->Name;
//   }
//   for (int i=0;i<Observables.size();i++){
//     observeString(i)=Observables(i)->Name;
//   }
  //  IOSection.WriteVar("Move Names",movesString);
  //  IOSection.WriteVar("Observable Names",observeString);
  
}
