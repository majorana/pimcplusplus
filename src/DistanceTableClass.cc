#include "DistanceTableClass.h"

template class MirroredAntiSymmetricMatrixClass<ImageNumClass>;


void DistanceTableClass::ShiftData(int slicesToShift,
				   CommunicatorClass &Communicator)
{

  DistTable.ShiftData(slicesToShift,Communicator);
  DispTable.ShiftData(slicesToShift,Communicator);
  ImageNumTable.ShiftData(slicesToShift,Communicator);

  /*
  if (slicesToShift>0){
    for (int slice=0;slice<slicesToShift;slice++){
      UpdateAll(slice);
    }
  }
  else if (slicesToShift<0){
    for (int slice=Path.NumTimeSlices()-1;
	 slice>Path.NumTimeSlices()+slicesToShift;slice--){
      UpdateAll(slice);
    }
    }*/

}

