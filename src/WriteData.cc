#include "WriteData.h"
#include "PathDataClass.h"
#include "Moves/MoveBase.h"
#include "Observables/ObservableBase.h"

void
WriteDataClass::DoEvent()
{
  list<MoveClass*>::iterator moveIter;
  for (moveIter=Moves.begin(); moveIter!=Moves.end(); moveIter++) 
    (*moveIter)->WriteRatio();

  list<ObservableClass*>::iterator observeIter;
  for (observeIter=Observables.begin(); observeIter!=Observables.end(); 
       observeIter++) 
    (*observeIter)->WriteBlock();
  if (PathData.Path.Communicator.MyProc()==0)
    IOSection.FlushFile();
}


void
WriteDataClass::Read(IOSectionClass &in)
{
  //do nothing for now
}



    
