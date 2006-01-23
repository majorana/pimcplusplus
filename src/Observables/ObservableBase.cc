#include "ObservableBase.h"
#include "time.h"

void 
ObservableClass::WriteInfo()
{
  if (PathData.Path.Communicator.MyProc()==0)
    IOSection.WriteVar("Description",Description);
}

void 
ObservableClass::Read(IOSectionClass &in)
{
  assert(in.ReadVar("Frequency",Frequency));
  assert(in.ReadVar("Name",Name));
  if(!(in.ReadVar("Description",Description))){
    Description="No description available";
  }
}



void 
ObservableClass::DoEvent()
{
  int start=clock();
  if ((TimesCalled % Frequency) == 0)
    Accumulate();
  TimesCalled++;
  int end=clock();
  SecondsInObservable+=
    (double)(end-start)/(double)CLOCKS_PER_SEC;
  
}
