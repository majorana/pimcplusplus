#ifndef WRAP_CLASS_H
#define WRAP_CLASS_H


#include "EventClass.h"
#include "Moves/MoveClass.h"
#include "Observables/ObservableClass.h"

class MoveWrap : public EventClass
{
 public:
  MoveClass* Move;
  void DoEvent() {Move->MakeMove();}
  void Read(IOSectionClass &IO) {Move->Read(IO);}
};


class ObserveWrap : public EventClass
{
 public:
  ObservableClass* Observe;
  ///Need to implement the do event here
  void DoEvent() {Observe->Accumulate();}
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
