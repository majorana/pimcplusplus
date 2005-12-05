#include "WriteData.h"
#include "PathDataClass.h"
#include "Moves/MoveBase.h"
#include "Observables/ObservableBase.h"

void
WriteDataClass::DoEvent()
{
  cerr<<"I am now calling write data"<<endl;
  list<MoveClass*>::iterator moveIter;
  for (moveIter=Moves.begin(); moveIter!=Moves.end(); moveIter++) 
    (*moveIter)->WriteRatio();

  list<ObservableClass*>::iterator observeIter;
  for (observeIter=Observables.begin(); observeIter!=Observables.end(); 
       observeIter++) 
    (*observeIter)->WriteBlock();
    


}


void
WriteDataClass::Read(IOSectionClass &in)
{
  //do nothing for now
}



    
