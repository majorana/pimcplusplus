#include "ObservableBase.h"


void ObservableClass::WriteInfo()
{
  if (PathData.Path.Communicator.MyProc()==0)
    IOSection.WriteVar("Description",Description);
}

void ObservableClass::Read(IOSectionClass &in)
{
  assert(in.ReadVar("Name",Name));
  if(!(in.ReadVar("Description",Description))){
    Description="No description available";
  }
  
}

void ObservableVar::Flush()
{
  if (Comm.MyProc()==0)
    Out.FlushFile();

}
