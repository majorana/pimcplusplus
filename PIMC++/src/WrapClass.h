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

#ifndef WRAP_CLASS_H
#define WRAP_CLASS_H

#include "time.h"
#include "EventClass.h"
#include "Moves/MoveClass.h"
#include "Observables/ObservableClass.h"

class MoveWrap : public EventClass
{
 public:
  MoveClass* Move;
  void DoEvent() 
  {
    if (Move->PathData.ExceededWallTime())
      return;
    int start=clock();
    Move->MakeMove();
    int end=clock();
    Move->SecondsInMove += (double)(end-start)/(double)CLOCKS_PER_SEC;
  }
  void Read(IOSectionClass &IO) {Move->Read(IO);}
};


class ObserveWrap : public EventClass
{
 public:
  ObservableClass* Observe;
  ///Need to implement the do event here
  void DoEvent() 
  {
    if (Observe->PathData.ExceededWallTime())
      return;
    int start=clock();
    Observe->Accumulate();
    int end=clock();
    Observe->SecondsInObservable+=
      (double)(end-start)/(double)CLOCKS_PER_SEC;
  }
  void Read(IOSectionClass &IO) {Observe->Read(IO);}
};

class LoopClass : public EventClass 
{
  int Steps;
  ///This is a pointer to the Move List. Will eventually need to have
  ///a pointer to the Observable List as well
  Array<MoveClass*,1> *MoveList;
  Array<ObservableClass* ,1> *ObserveList;
  Array<EventClass*,1> Events;
 public:
  void DoEvent() {
    if ((MoveList->size()>0))
      if ((*MoveList)(0)->PathData.ExceededWallTime())
	return;
    //    cerr<<"The number of events are "<<Events.size()<<" "<<Steps<<endl;
    //    cerr<<"Steps are "<<Steps<<endl;
    //    cerr<<"The line after"<<endl;
    for (int countSteps=0;countSteps<Steps;countSteps++){
      for (int counter=0;counter<Events.size();counter++){
	Events(counter)->DoEvent();
      }
    }

  };
  void Read(IOSectionClass &IO);
  void Read(IOSectionClass &IO, int steps);
  LoopClass(Array<MoveClass*,1> *tempMoveList, 
	    Array<ObservableClass*,1> *tempObserveList){
    MoveList=tempMoveList;
    ObserveList=tempObserveList;
  }
   
};

#endif
