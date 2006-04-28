#ifndef LOOP_CLASS_H
#define LOOP_CLASS_H

#include "EventClass.h"

class ObservableClass;
class MoveClass;

/// Note:  LoopClass is not allow to write any output to its out
/// IOsection. 
class LoopClass : public EventClass 
{
protected:
  int NumSteps;
  std::list<MoveClass*> &Moves;
  std::list<ObservableClass*> &Observables;
  std::list<EventClass*> Events;
  MoveClass *FindMove(string name);
  ObservableClass *FindObservable(string name);
 public:
  void DoEvent();
  void Read(IOSectionClass &IO);
  void Read(IOSectionClass &IO, int steps);
  void SetOutfile (IOSectionClass &out)
  { IOSection = out; }
  LoopClass(PathDataClass &pathData, IOSectionClass &out,
	    list<MoveClass*> &moves,
	    list<ObservableClass*> &observables) : 
    EventClass(pathData, out), NumSteps(1), 
    Moves(moves), Observables(observables)
  {
    
  }
};



#endif
