#include "MetaMoves.h"



void PrintMoveClass::Read(IOSectionClass &IO)
{
  string typeCheck;
  assert(IO.ReadVar("type",typeCheck));
  assert(typeCheck=="PrintMove");
  assert(IO.ReadVar("name",Name));
  assert(IO.ReadVar("toprint",MyString));
}

//////////////////////////////////////////////////////////////

void ShiftMoveClass::Read(IOSectionClass &theInput)
{
  string typeCheck;
  assert(theInput.ReadVar("type",typeCheck));
  assert(typeCheck=="ShiftMove");
  assert(theInput.ReadVar("name",Name));

}


void ShiftMoveClass::MakeMove()
{
  int numTimeSlicesToShift =(int)floor((PathData.NumTimeSlices()-1)*PathData.Path.Random.Common()); 
  PathData.MoveJoin(0);
  PathData.ShiftData(numTimeSlicesToShift);
  PathData.Join=numTimeSlicesToShift;
}
