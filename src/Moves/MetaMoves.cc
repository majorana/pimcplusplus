#include "MetaMoves.h"



void PrintMoveClass::Read(IOSectionClass &IO)
{
  string typeCheck;
  assert(IO.ReadVar("Type",typeCheck));
  assert(typeCheck=="PrintMove");
  assert(IO.ReadVar("toprint",MyString));
}

//////////////////////////////////////////////////////////////

void ShiftMoveClass::Read(IOSectionClass &theInput)
{
  string typeCheck;
  assert(theInput.ReadVar("Type",typeCheck));
  assert(typeCheck=="ShiftMove");
}


void ShiftMoveClass::MakeMove()
{
  int slice1, slice2;
  // The last processor will have the least number of slices
  // possible.  Use that for the maximum shift.
  PathData.Path.SliceRange(PathData.Path.Communicator.NumProcs()-1,
			   slice1, slice2);
  int maxSlices = slice2-slice1;

  int numTimeSlicesToShift = PathData.Path.Random.CommonInt(maxSlices);

  //  HACK! HACK! HACK!!NEED THESE LINES

  PathData.MoveJoin(0);
  PathData.ShiftData(numTimeSlicesToShift);
  PathData.Join=numTimeSlicesToShift;

 


  ///END HACK!!!
//   Array<int,1> changedParticles(PathData.Path.NumParticles());
//   for (int counter=0;counter<changedParticles.size();counter++){
//     changedParticles(counter)=counter;
//   }
   
//   cerr.precision(10);
//   cerr<<"My current kinetic action is "
//       <<  PathData.Actions.Kinetic.Action(0,PathData.Path.TotalNumSlices,
//  					  changedParticles,0)<<endl;
//   cerr<<"My current shortrange action is "
//       <<  PathData.Actions.ShortRange.Action(0,PathData.Path.TotalNumSlices,
// 					     changedParticles,0)<<endl;
//   cerr<<"My current approximate shortrange action is "
//       <<  PathData.Actions.ShortRangeApproximate.Action(0,PathData.Path.TotalNumSlices,changedParticles,0)<<endl;

}
