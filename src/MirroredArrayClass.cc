
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
void MirroredArrayClass<T>::Resize(int numTimeSlices,int numParticles)
{
  AB.resize(2,numTimeSlices,numParticles);
}


///Moves the join. The join si defined so that the permutation happens
///between the time slice the join is on and the time slice after
///the join (i.e. a(join+1)=a(p(join)), but a(join) = a(join)
template <class T>
void MirroredArrayClass<T>::MoveJoin(MirroredArrayClass1D<int> &PermMatrix,
				     int oldJoin, int newJoin)
{
  if (newJoin>oldJoin){
    for (int timeSlice=oldJoin+1;timeSlice<=newJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
	AB(0,timeSlice,ptcl)=AB(1,timeSlice,PermMatrix(ptcl));
      }
    }
    //Now that we've copied the data from B into A, we need to copy the 
    //information into a
    for (int timeSlice=oldJoin+1;timeSlice<=newJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
	AB(0,timeSlice,ptcl)=AB(1,timeSlice,ptcl);
      }
    }
  }
  else if (oldJoin>=newJoin){
    for (int timeSlice=newJoin+1;timeSlice<=oldJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
	AB(0,timeSlice,PermMatrix(ptcl))=AB(1,timeSlice,ptcl);
      }
    }
    //Now that we've copied the data from B into A, we need to copy the 
    //information into a
    for (int timeSlice=newJoin+1;timeSlice<=oldJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
	AB(0,timeSlice,ptcl)=AB(1,timeSlice,ptcl);
      }
    }
  }
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
  int NumPtcls=AB.extent(2);
  int NumSlices=AB.extent(1);

  sendProc=(myProc+1) % numProcs;
  recvProc=((myProc-1) + numProcs) % numProcs;
  if (slicesToShift<0){
    int tempProc=sendProc;
    sendProc=recvProc;
    recvProc=tempProc;
  }

  ///First shifts the data in the A copy left 
  ///or right by the appropriate amount   
  if (slicesToShift>0){
    for (int ptclCounter=0;ptclCounter<NumPtcls;ptclCounter++){
      for (int sliceCounter=NumSlices-1;
	   sliceCounter>=slicesToShift;sliceCounter--){
	AB(0,sliceCounter,ptclCounter)=
	  AB(0,sliceCounter-slicesToShift,ptclCounter);
      }
    }
  }
  else{
    for (int ptclCounter=0;ptclCounter<NumPtcls;ptclCounter++){
      for (int sliceCounter=0;
	   sliceCounter<NumSlices+slicesToShift;sliceCounter++){
	AB(0,sliceCounter,ptclCounter)=
	  AB(0,sliceCounter-slicesToShift,ptclCounter);
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
    for (int sliceCounter=startTimeSlice;
	 sliceCounter<startTimeSlice+abs(slicesToShift);sliceCounter++){
      sendBuffer(BufferCounter)=AB(1,sliceCounter,ptclCounter);
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
    for (int sliceCounter=startTimeSlice;
	 sliceCounter<startTimeSlice+abs(slicesToShift);sliceCounter++){
      AB(0,sliceCounter,ptclCounter)=receiveBuffer(BufferCounter);
      BufferCounter++;
    }
  }
  
  // Now copy A into B, since A has all the good, shifted data now.
  for (int ptclCounter=0;ptclCounter<NumPtcls;ptclCounter++)
    for (int sliceCounter=0; sliceCounter<NumSlices; sliceCounter++)
      AB(1,sliceCounter,ptclCounter) = AB(0,sliceCounter,ptclCounter);

  // And we're done!

}


