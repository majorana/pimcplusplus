#include "LoopClass.h"
#include "PathDataClass.h"

void
LoopClass::DoEvent()
{
  std::list<EventClass*>::iterator iter;
  for (int step=0; step<NumSteps; step++)
    for (iter=Events.begin(); iter!=Events.end(); iter++) {
      if (PathData.ExceededWallTime())
	return;
      else
	(*iter)->DoEvent();
    }
}




EventClass*
LoopClass::FindMove(string name)
{
  list<EventClass*>::iterator iter;
  for (iter=Moves.begin(); iter != Moves.end(); iter++)
    if ((*iter)->Name == name)
      return (*iter);
}

EventClass*
LoopClass::FindObservable(string name) 
{
  list<EventClass*>::iterator iter;
  for (iter=Observables.begin(); iter != Observables.end(); iter++)
    if ((*iter)->Name == name)
      return (*iter);
  return NULL;
}
  

void
LoopClass::Read(IOSectionClass &in)
{
  assert (in.ReadVar("Steps", NumSteps));
  Read(in, NumSteps);
}

void
LoopClass::Read(IOSectionClass &in, int steps)
{
  NumSteps = steps;
  int numSections = in.CountSections();
  for (int secNum=0; secNum<numSections; secNum++) {
    in.OpenSection(secNum);
    if ((in.GetName()=="Move")){
      string name;
      assert (in.ReadVar("Name", name));
      EventClass *event = FindMove(name);
      if (event == NULL) {
	cerr << "Unknown move """ << name << """.\n";
	abort();
      }
      Events.push_back(event);
    }
    else if(in.GetName() == "Observe") {
      string name;
      assert (in.ReadVar("Name", name));
      EventClass *event = FindMove(name);
      if (event == NULL) {
	cerr << "Unknown observable """ << name << """.\n";
	abort();
      }
      Events.push_back(event);
    }
    else if (in.GetName() == "Loop") {
      LoopClass *newLoop = 
	new LoopClass (PathData, IOSection, Moves, Observables);
      newLoop->Read(in);
      Events.push_back(newLoop);
    }
  }
}


    
