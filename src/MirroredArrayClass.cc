
#include "MirroredArrayClass.h"


template class MirroredArrayClass<int>;
template class MirroredArrayClass<dVec>;
template class MirroredArrayClass1D<int>;


int Write1=0;
int Write2=1;



template<class T>
void MirroredArrayClass1D<T>::Resize(int numParticles)
{
  AB.resize(2,numParticles);
}


template <class T>
void MirroredArrayClass1D<T>::Print()
{
  for (int counter=0;counter<AB.extent(0);counter++){
      cout<<AB(0,counter)<<" ";
    }
    cout<<endl;
}

//template <class T>
//void MirroredArrayClass1D<T>::ShiftData(int slicesToShift, CommunicatorClass &Communicator)
//{
//  int numProcs=Communicator.NumProcs();
//  int myProc=Communicator.MyProc();
//  int recvProc;
//  int sendProc;
//  int NumPtcls=1;


//  sendProc=(myProc+1) % numProcs;
//  recvProc=((myProc-1) + numProcs) % numProcs;
//  if (slicesToShift<0){
//    int tempProc=sendProc;
//    sendProc=recvProc;
//    recvProc=tempProc;
//  }
//  ///First shifts the data in the A copy left or right by the appropriate amount
//  if (slicesToShift>0){
//      for (int sliceCounter=NumSlices-1;sliceCounter>=slicesToShift;sliceCounter--){
//	AB(0,sliceCounter)=AB(0,sliceCounter-slicesToShift);	
//     }
//    }
//  else{

//      for (int sliceCounter=0;sliceCounter<NumSlices+slicesToShift;sliceCounter++){
//	AB(0,sliceCounter)=AB(0,sliceCounter-slicesToShift);
//      }
//    }

  
//  int bufferSize=abs(slicesToShift)*NumPtcls;
//  Array<T,1> sendBuffer(bufferSize), receiveBuffer(bufferSize);
//  int startTimeSlice;
//  if (slicesToShift>0)
//    startTimeSlice=NumSlices-slicesToShift;
//  else 
//    startTimeSlice=0;

//  int BufferCounter=0;
//  for (int ptclCounter=0;ptclCounter<NumPtcls;ptclCounter++){
//    for (int sliceCounter=startTimeSlice;sliceCounter<startTimeSlice+abs(slicesToShift);sliceCounter++){
//     sendBuffer(BufferCounter)=AB(1,ptclCounter,sliceCounter);
//     BufferCounter++;
//    }
//  }
  
//  Communicator.SendReceive(sendProc, sendBuffer,recvProc, receiveBuffer);
  
//  if (slicesToShift>0)
//    startTimeSlice=0;
//  else 
//    startTimeSlice=NumSlices+slicesToShift;

//  BufferCounter=0;
//  for (int ptclCounter=0;ptclCounter<NumPtcls;ptclCounter++){
//    for (int sliceCounter=startTimeSlice;sliceCounter<startTimeSlice+abs(slicesToShift);sliceCounter++){
//      AB(0,ptclCounter,sliceCounter)=receiveBuffer(BufferCounter);
//      BufferCounter++;
//    }
//  }
  
//  // Now copy A into B, since A has all the good, shifted data now.
//  for (int ptclCounter=0;ptclCounter<NumPtcls;ptclCounter++)
//    for (int sliceCounter=0; sliceCounter<NumSlices; sliceCounter++)
//      AB(1,ptclCounter,sliceCounter) = AB(0,ptclCounter,sliceCounter);

//  // And we're done!

//}



template<class T>
void MirroredArrayClass<T>::Resize(int numParticles,int numTimeSlices)
{
  AB.resize(2,numParticles,numTimeSlices);
}


template <class T>
void MirroredArrayClass<T>::MoveJoin(MirroredArrayClass1D<int> &PermMatrix,int oldJoin, int newJoin)
{
  if (newJoin>oldJoin){
    for (int counterP=0;counterP<NumParticles();counterP++){
      for (int counter=oldJoin;counter<newJoin;counter++){
	AB(0,counterP,counter)=AB(1,PermMatrix(counterP),counter);
      }
    }
  }
  else if (oldJoin>=newJoin){
    for (int counterP=0;counterP<NumParticles();counterP++){
      for (int counter=newJoin;counter<oldJoin;counter++){
	AB(0,PermMatrix(counterP),counter)=AB(1,counterP,counter);
      }
    }
  }

  //Now that we've copied the data from B into A, we need to copy the 
  //information into a

  if (newJoin>oldJoin){
    for (int counterP=0;counterP<NumParticles();counterP++){
      for (int counter=oldJoin;counter<newJoin;counter++){
	AB(0,counterP,counter)=AB(1,counterP,counter);
      }
    }
  }
  else if (oldJoin>=newJoin){
    for (int counterP=0;counterP<NumParticles();counterP++){
      for (int counter=newJoin;counter<oldJoin;counter++){
	AB(0,PermMatrix(counterP),counter)=AB(1,PermMatrix(counterP),counter);
      }
    }
  }



  return;

}


template <class T>
void MirroredArrayClass<T>::Print()
{
  for (int counter=0;counter<AB.extent(1);counter++){
    for (int counter2=0;counter2<AB.extent(2);counter2++){
      cout<<AB(0,counter,counter2)<<" ";
    }
    cout<<endl;
  }
}

template <class T>
void MirroredArrayClass<T>::ShiftData(int slicesToShift, CommunicatorClass &Communicator)
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
    recvProc=tempProc;
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
  
  int bufferSize=abs(slicesToShift)*NumPtcls;
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
  
  Communicator.SendReceive(sendProc, sendBuffer,recvProc, receiveBuffer);
  
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


