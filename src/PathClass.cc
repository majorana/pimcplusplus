#include "PathClass.h"


void PathClass::shiftData(int slicesToShift,CommClass &Communicator)
{
  Positions.shiftData(slicesToShift,Communicator);
  TimeStamp.shiftData(slicesToShift,Communicator);


}
