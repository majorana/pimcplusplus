#include "PathClass.h"


void PathClass::ShiftData(int slicesToShift,CommunicatorClass &Communicator)
{
  Positions.ShiftData(slicesToShift,Communicator);
  TimeStamp.ShiftData(slicesToShift,Communicator);


}
