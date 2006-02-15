#ifndef WRITE_DATA_CLASS_H
#define WRITE_DATA_CLASS_H

#include "EventClass.h"
class ObservableClass;
class MoveClass;


/// Note:  LoopClass is not allow to write any output to its out
/// IOsection. 
class WriteDataClass : public EventClass 
{
protected:
  list<MoveClass*> &Moves;
  list<ObservableClass*> &Observables;
 public:
  void DoEvent();
  void Read(IOSectionClass &IO);
  WriteDataClass(PathDataClass &pathData, IOSectionClass &out,
	    list<MoveClass*> &moves,
	    list<ObservableClass*> &observables) : 
    EventClass(pathData, out),
    Moves(moves), Observables(observables)
  {
    Name = "WriteData";
  }
};



#endif
