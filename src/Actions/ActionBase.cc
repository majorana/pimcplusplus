#include "ActionBase.h"
#include "../PathDataClass.h"

void ActionBaseClass::Read (IOSectionClass &in)
{
  // Do nothing for now
}

ActionBaseClass::ActionBaseClass(PathDataClass &pathData) : 
  PathData(pathData), Path(pathData.Path)
{
  /* Do nothing */
}

PotentialBaseClass::PotentialBaseClass(PathDataClass &pathData) : 
  PathData(pathData), Path(pathData.Path)
{
  /* Do nothing */
}
