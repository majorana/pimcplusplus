#include "ObservableBase.h"


void ObservableClass::WriteInfo()
{
  IOSection.WriteVar("Description",Description);
}

void ObservableClass::Read(IOSectionClass &in)
{
  assert(in.ReadVar("Name",Name));
  if(!(in.ReadVar("Description",Description))){
    Description="No description available";
  }
  
}
