// #pragma implementation "MirroredArrayClass.h"
#include "MirroredArrayClass.h"


template class MirroredArrayClass<int>;
template class MirroredArrayClass<dVec>;


// #pragma implementation "MirroredArrayClass.h"

//template <class T>
//MirroredArrayClass<T>::MirroredArrayClass(int dim1, int dim2)
//{


//}

//template <class T>
//MirroredArrayClass<T>::MirroredArrayClass()
//{
//  cerr<<"You've called a not very useful constructor!!! WARNING!! WARNING!!!";
//  AB.resize(2,0,0); 

//}

int Write1=0;
int Write2=1;

template<class T>
void MirroredArrayClass<T>::resize(int numParticles,int numTimeSlices)
{
  AB.resize(2,numParticles,numTimeSlices);

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
void MirroredArrayClass<T>::shiftData(int slicesToShift, CommunicatorClass &Communicator)
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

//void foo()
//{
//  MirroredArrayClass<int> myArray;
//  CommunicatorClass myCommunicator;
//  myArray.Print();
//  myArray.shiftData(3,myCommunicator);
//  MirroredArrayClass<dVec> myArray2;
//  myArray2.resize(1,1);
//  myArray.resize(1,1);
//  myArray2.shiftData(3,myCommunicator);
//}
