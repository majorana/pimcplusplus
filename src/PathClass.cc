#include "PathClass.h"


void PathClass::shiftData(int slicesToShift,CommunicatorClass &Communicator)
{
  Positions.shiftData(slicesToShift,Communicator);
  TimeStamp.shiftData(slicesToShift,Communicator);


}
