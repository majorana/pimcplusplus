#include "MirroredArrayClass.h"

void MirroredArrayClass::shiftData(int slicesToShift, CommunicatorClass &Communicator)
{
  int numProcs=Communicator.NumProcs();
  int myProc=Communicator.MyProc();
  int recvProc;
  int sendProc;
  int NumPtcls=AB.extent(1);
  int NumSlices=AB.extent(2);

  sendProc=(myProc+1) % numProcs;
  recvProc=((myProc-1) + numProcs) % numProcs;
  if (slicesToShift<0){
    int tempProc=sendProc;
    sendProc=recvProc;
  }

  ///First shifts the data in the A copy left or right by the appropriate amount
  if (slicesToShift>0){
    for (int ptclCounter=0;ptclCounter<NumPtcls;ptclCounter++){
      for (int sliceCounter=NumSlices-1;sliceCounter>=slicesToShift;sliceCounter--){
	AB(0,ptclCounter,sliceCounter)=AB(0,ptclCounter,sliceCounter-slicesToShift);	
      }
    }
  }
  else{
    for (int ptclCounter=0;ptclCounter<NumPtcls;ptclCounter++){
      for (int sliceCounter=0;sliceCounter<NumSlices+slicesToShift;sliceCounter++){
	AB(0,ptclCounter,sliceCounter)=AB(0,ptclCounter,sliceCounter-slicesToShift);
      }
    }
  }
  
  bufferSize=slicesToShift*NumPtcls;
  Array<T,1> sendBuffer(bufferSize), receiveBuffer(bufferSize);
  int startTimeSlice;
  if (slicesToShift>0)
    startTimeSlice=NumSlices-slicesToShift;
  else 
    startTimeSlice=0;

  int BufferCounter=0;
  for (int ptclCounter=0;ptclCounter<NumPtcls;ptclCounter++){
    for (int sliceCounter=startTimeSlice;sliceCounter<startTimeSlice+abs(slicesToShift);sliceCounter++){
     sendBuffer(BufferCounter)=AB(1,ptclCounter,sliceCounter);
     BufferCounter++;
    }
  }
  
  Communicator.SendReceive(sendProc, sendBuffer,receiveProc, receiveBuffer);
  
  if (slicesToShift>0)
    startTimeSlice=0;
  else 
    startTimeSlice=NumSlices+slicesToShift;

  BufferCounter=0;
  for (int ptclCounter=0;ptclCounter<NumPtcls;ptclCounter++){
    for (int sliceCounter=startTimeSlice;sliceCounter<startTimeSlice+abs(slicesToShift);sliceCounter++){
      AB(0,ptclCounter,sliceCounter)=receiveBuffer(BufferCounter);
      BufferCounter++;
    }
  }
  
  // Now copy A into B, since A has all the good, shifted data now.
  for (int ptclCounter=0;ptclCounter<NumPtcls;ptclCounter++)
    for (int sliceCounter=0; sliceCounter<NumSlices; sliceCounter++)
      AB(1,ptclCounter,sliceCounter) = AB(0,ptclCounter,sliceCounter);

  // And we're done!

}
