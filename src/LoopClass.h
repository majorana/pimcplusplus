#ifndef LOOP_CLASS_H
#define LOOP_CLASS_H

#include "EventClass.h"


/// Note:  LoopClass is not allow to write any output to its out
/// IOsection. 
class LoopClass : public EventClass 
{
protected:
  int NumSteps;
  std::list<EventClass*> &Moves, &Observables;
  std::list<EventClass*> Events;
  EventClass *FindMove(string name);
  EventClass *FindObservable(string name);
 public:
  void DoEvent();
  void Read(IOSectionClass &IO);
  void Read(IOSectionClass &IO, int steps);
  LoopClass(PathDataClass &pathData, IOSectionClass &out,
	    list<EventClass*> &moves,
	    list<EventClass*> &observables) : 
    EventClass(pathData, out), NumSteps(1), 
    Moves(moves), Observables(observables)
  {
    
  }
};



#endif
