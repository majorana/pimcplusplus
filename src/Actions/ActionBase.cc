#include "ActionBase.h"
#include "../PathDataClass.h"

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
