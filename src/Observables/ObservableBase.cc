#include "ObservableBase.h"
#include "time.h"


void
ObservableVar::Flush()
{
  if (Comm.MyProc() == 0)
    Out.FlushFile();
}


void 
ObservableClass::WriteInfo()
{
  if (PathData.Path.Communicator.MyProc()==0)
    IOSection.WriteVar("Description",Description);
}

void 
ObservableClass::Read(IOSectionClass &in)
{
  in.ReadVar("Prefactor", Prefactor);
  assert(in.ReadVar("Frequency",Frequency));
  assert(in.ReadVar("Name",Name));
  if(!(in.ReadVar("Description",Description))){
    Description="No description available";
  }
}


void 
ObservableClass::DoEvent()
{
  if ((TimesCalled % Frequency) == 0)
    Accumulate();
  TimesCalled++;
}
