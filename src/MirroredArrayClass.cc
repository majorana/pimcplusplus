#include "Common.h"
#include "MirroredArrayClass.h"


template class MirroredArrayClass<int>;
template class MirroredArrayClass<double>;
template class MirroredArrayClass<dVec>;
template class MirroredArrayClass1D<int>;
template class MirroredSymmetricMatrixClass<double>;
template class MirroredAntiSymmetricMatrixClass<dVec>;
template class MirroredAntiSymmetricMatrixClass<ImageNumClass>; 
template class MirroredAntiSymmetricMatrixClass<int>; 
template class MirroredSymmetricMatrixClass<int>; 

int Write1=0;
int Write2=1;



/****************************************************************/
/*                     MirroredArrayClass1D                     */            
/****************************************************************/
template <class T>
void MirroredArrayClass1D<T>::Print()
{
  for (int counter=0;counter<AB.extent(0);counter++){
      cout<<AB(0,counter)<<" ";
    }
    cout<<endl;
}


/****************************************************************/
/*                      MirroredArrayClass                      */            
/****************************************************************/
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
    //information into B
    for (int timeSlice=oldJoin+1;timeSlice<=newJoin;timeSlice++){ 
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
	AB(1,timeSlice,ptcl)=AB(0,timeSlice,ptcl);
      }
    }
  }
  else if (oldJoin>newJoin){
    //  else if (oldJoin>=newJoin){//CHANGED!
    for (int timeSlice=newJoin+1;timeSlice<=oldJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
	AB(0,timeSlice,PermMatrix(ptcl))=AB(1,timeSlice,ptcl);
      }
    }
    //Now that we've copied the data from B into A, we need to copy the 
    //information into B
    for (int timeSlice=newJoin+1;timeSlice<=oldJoin;timeSlice++){
      for (int ptcl=0;ptcl<NumParticles();ptcl++){
	AB(1,timeSlice,ptcl)=AB(0,timeSlice,ptcl);
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
void MirroredArrayClass<T>::ShiftData(int slicesToShift, PIMCCommunicatorClass &Communicator)
{
  int numProcs=Communicator.NumProcs();
  int myProc=Communicator.MyProc();
  int recvProc;
  int sendProc;
  int NumPtcls=AB.extent(2);
  int NumSlices=AB.extent(1);
  assert(abs(slicesToShift)<NumSlices);
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
  int BufferCounter=0;
  if (slicesToShift>0){
    startTimeSlice=NumSlices-slicesToShift;

    for (int ptclCounter=0;ptclCounter<NumPtcls;ptclCounter++){
      for (int sliceCounter=startTimeSlice;
	   sliceCounter<startTimeSlice+abs(slicesToShift);sliceCounter++){
	///If shifting forward, don't send the last time slice (so always)
	///send sliceCounter-1
	sendBuffer(BufferCounter)=AB(1,sliceCounter-1,ptclCounter);
	BufferCounter++;
      }
    }
  }
  else {
    startTimeSlice=0;
    for (int ptclCounter=0;ptclCounter<NumPtcls;ptclCounter++){
      for (int sliceCounter=startTimeSlice;
	   sliceCounter<startTimeSlice+abs(slicesToShift);sliceCounter++){
	///If shifting backward, don't send the first time slice (so always)
	///send sliceCounter+1
	sendBuffer(BufferCounter)=AB(1,sliceCounter+1,ptclCounter);
	BufferCounter++;
      }
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


/****************************************************************/
/*                 MirroredSymmetricMatrixClass                 */            
/****************************************************************/

template <class T>
void MirroredSymmetricMatrixClass<T>::Print()
{
  for (int counter=0;counter<AB.extent(1);counter++){
    for (int counter2=0;counter2<AB.extent(2);counter2++){
      cout<<AB(0,counter,counter2)<<" ";
    }
    cout<<endl;
  }
}


template <class T>
void MirroredSymmetricMatrixClass<T>::ShiftData(int slicesToShift, PIMCCommunicatorClass &Communicator)
{
  int numProcs=Communicator.NumProcs();
  int myProc=Communicator.MyProc();
  int recvProc;
  int sendProc;
  int NumIndices=AB.extent(2);
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
    for (int pairIndex=0;pairIndex<NumIndices;pairIndex++){
      for (int sliceCounter=NumSlices-1;
	   sliceCounter>=slicesToShift;sliceCounter--){
	AB(0,sliceCounter,pairIndex)=
	  AB(0,sliceCounter-slicesToShift,pairIndex);
      }
    }
  }
  else{
    for (int pairIndex=0;pairIndex<NumIndices;pairIndex++){
      for (int sliceCounter=0;
	   sliceCounter<NumSlices+slicesToShift;sliceCounter++){
	AB(0,sliceCounter,pairIndex)=
	  AB(0,sliceCounter-slicesToShift,pairIndex);
      }
    }
  }
  
  int bufferSize=abs(slicesToShift)*NumIndices;
  Array<T,1> sendBuffer(bufferSize), receiveBuffer(bufferSize);
  int BufferCounter=0;
  int startTimeSlice;
  if (slicesToShift>0){
    startTimeSlice=NumSlices-slicesToShift;
    for (int pairIndex=0;pairIndex<NumIndices;pairIndex++){
      for (int sliceCounter=startTimeSlice;
	   sliceCounter<startTimeSlice+abs(slicesToShift);sliceCounter++){
	sendBuffer(BufferCounter)=AB(1,sliceCounter-1,pairIndex);
	BufferCounter++;
      }
    }

  }
  else {
    startTimeSlice=0;
    for (int pairIndex=0;pairIndex<NumIndices;pairIndex++){
      for (int sliceCounter=startTimeSlice;
	   sliceCounter<startTimeSlice+abs(slicesToShift);sliceCounter++){
	sendBuffer(BufferCounter)=AB(1,sliceCounter+1,pairIndex);
	BufferCounter++;
      }
    }
    
  }
  
  Communicator.SendReceive(sendProc, sendBuffer,recvProc, receiveBuffer);
  
  if (slicesToShift>0)
    startTimeSlice=0;
  else 
    startTimeSlice=NumSlices+slicesToShift;

  BufferCounter=0;
  for (int pairIndex=0;pairIndex<NumIndices;pairIndex++){
    for (int sliceCounter=startTimeSlice;
	 sliceCounter<startTimeSlice+abs(slicesToShift);sliceCounter++){
      AB(0,sliceCounter,pairIndex)=receiveBuffer(BufferCounter);
      BufferCounter++;
    }
  }
  
  // Now copy A into B, since A has all the good, shifted data now.
  for (int pairIndex=0;pairIndex<NumIndices;pairIndex++)
    for (int sliceCounter=0; sliceCounter<NumSlices; sliceCounter++)
      AB(1,sliceCounter,pairIndex) = AB(0,sliceCounter,pairIndex);

  // And we're done!

}

///Moves the join. The join si defined so that the permutation happens
///between the time slice the join is on and the time slice after
///the join (i.e. a(join+1)=a(p(join)), but a(join) = a(join)
template <class T>
void MirroredSymmetricMatrixClass<T>::MoveJoin(MirroredArrayClass1D<int> &PermMatrix,
				     int oldJoin, int newJoin)
{
  
  Array <int,1> PermPairIndex(AB.extent(2));
  ///Here we construct the permutations on the pairs
  int pairIndex=0;
  for (int ptcl1=0;ptcl1<NumPtcls;ptcl1++){
    int permPtcl1=PermMatrix(ptcl1);
    for (int ptcl2=0;ptcl2<=ptcl1;ptcl2++){
      int permPtcl2=PermMatrix(ptcl2);
      PermPairIndex(pairIndex)=PairIndex(permPtcl1,permPtcl2);
      pairIndex++;
    }
  }
  if (newJoin>oldJoin){
    for (int timeSlice=oldJoin+1;timeSlice<=newJoin;timeSlice++){
      for (int index=0;index<AB.extent(2);index++){
	AB(0,timeSlice,index)=AB(1,timeSlice,PermPairIndex(index));
      }
    }
    //Now that we've copied the data from B into A, we need to copy the 
    //information into B
    for (int timeSlice=oldJoin+1;timeSlice<=newJoin;timeSlice++){
      for (int index=0;index<AB.extent(2);index++){
	AB(1,timeSlice,index)=AB(0,timeSlice,index);
      }
    }
  }
  else if (oldJoin>newJoin){ //CHANGED!!
    for (int timeSlice=newJoin+1;timeSlice<=oldJoin;timeSlice++){
      for (int index=0;index<AB.extent(2);index++){
	AB(0,timeSlice,PermPairIndex(index))=AB(1,timeSlice,index);
      }
    }
    //Now that we've copied the data from B into A, we need to copy the 
    //information into B
    for (int timeSlice=newJoin+1;timeSlice<=oldJoin;timeSlice++){
      for (int index=0;index<NumParticles();index++){
	AB(1,timeSlice,index)=AB(0,timeSlice,index);
      }
    }
  }
}





/****************************************************************/
/*               MirroredAntiSymmetricMatrixClass               */            
/****************************************************************/


template <class T>
void MirroredAntiSymmetricMatrixClass<T>::Print()
{
  for (int counter=0;counter<AB.extent(1);counter++){
    for (int counter2=0;counter2<AB.extent(2);counter2++){
      cout<<AB(0,counter,counter2)<<" ";
    }
    cout<<endl;
  }
}

template <class T>
void MirroredAntiSymmetricMatrixClass<T>::ShiftData(int slicesToShift, 
						    PIMCCommunicatorClass 
						    &Communicator)
{
  int numProcs=Communicator.NumProcs();
  int myProc=Communicator.MyProc();
  int recvProc;
  int sendProc;
  int NumIndices=AB.extent(2);
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
    for (int pairIndex=0;pairIndex<NumIndices;pairIndex++){
      for (int sliceCounter=NumSlices-1;
	   sliceCounter>=slicesToShift;sliceCounter--){
	AB(0,sliceCounter,pairIndex)=
	  AB(0,sliceCounter-slicesToShift,pairIndex);
      }
    }
  }
  else{
    for (int pairIndex=0;pairIndex<NumIndices;pairIndex++){
      for (int sliceCounter=0;
	   sliceCounter<NumSlices+slicesToShift;sliceCounter++){
	AB(0,sliceCounter,pairIndex)=
	  AB(0,sliceCounter-slicesToShift,pairIndex);
      }
    }
  }
  
  int bufferSize=abs(slicesToShift)*NumIndices;
  Array<T,1> sendBuffer(bufferSize), receiveBuffer(bufferSize);
  int startTimeSlice;
  int BufferCounter=0;
  if (slicesToShift>0){
    startTimeSlice=NumSlices-slicesToShift;
    for (int pairIndex=0;pairIndex<NumIndices;pairIndex++){
      for (int sliceCounter=startTimeSlice;
	   sliceCounter<startTimeSlice+abs(slicesToShift);sliceCounter++){
	sendBuffer(BufferCounter)=AB(1,sliceCounter-1,pairIndex);
	BufferCounter++;
      }
    }
    
  }
  else {
    startTimeSlice=0;
    for (int pairIndex=0;pairIndex<NumIndices;pairIndex++){
      for (int sliceCounter=startTimeSlice;
	   sliceCounter<startTimeSlice+abs(slicesToShift);sliceCounter++){
	sendBuffer(BufferCounter)=AB(1,sliceCounter+1,pairIndex);
	BufferCounter++;
      }
    }
    
  }


  
  Communicator.SendReceive(sendProc, sendBuffer,recvProc, receiveBuffer);
  
  if (slicesToShift>0)
    startTimeSlice=0;
  else 
    startTimeSlice=NumSlices+slicesToShift;

  BufferCounter=0;
  for (int pairIndex=0;pairIndex<NumIndices;pairIndex++){
    for (int sliceCounter=startTimeSlice;
	 sliceCounter<startTimeSlice+abs(slicesToShift);sliceCounter++){
      AB(0,sliceCounter,pairIndex)=receiveBuffer(BufferCounter);
      BufferCounter++;
    }
  }
  
  // Now copy A into B, since A has all the good, shifted data now.
  for (int pairIndex=0;pairIndex<NumIndices;pairIndex++)
    for (int sliceCounter=0; sliceCounter<NumSlices; sliceCounter++)
      AB(1,sliceCounter,pairIndex) = AB(0,sliceCounter,pairIndex);

  // And we're done!

}

///Moves the join. The join si defined so that the permutation happens
///between the time slice the join is on and the time slice after
///the join (i.e. a(join+1)=a(p(join)), but a(join) = a(join)
template <class T>
void MirroredAntiSymmetricMatrixClass<T>::MoveJoin
(MirroredArrayClass1D<int> &PermMatrix, int oldJoin, int newJoin)
{
  
  Array <int,1> PermPairIndex(AB.extent(2));
  // Must flip sign if true
  Array <bool,1> FlipSign(AB.extent(2));
  ///Here we construct the permutations on the pairs
  int pairIndex=0;
  for (int ptcl1=0;ptcl1<NumPtcls;ptcl1++){
    int permPtcl1=PermMatrix(ptcl1);
    for (int ptcl2=0;ptcl2<=ptcl1;ptcl2++){
      int permPtcl2=PermMatrix(ptcl2);
      PermPairIndex(pairIndex)=PairIndex(permPtcl1,permPtcl2);
      FlipSign(pairIndex) = permPtcl2 > permPtcl1;
      pairIndex++;
    }
  }
  if (newJoin>oldJoin){
    for (int timeSlice=oldJoin+1;timeSlice<=newJoin;timeSlice++){
      for (int index=0;index<AB.extent(2);index++){
	if (FlipSign(index))
	  AB(0,timeSlice,index) = -AB(1,timeSlice,PermPairIndex(index));
	else
	  AB(0,timeSlice,index) =  AB(1,timeSlice,PermPairIndex(index));
      }
    }
    //Now that we've copied the data from B into A, we need to copy the 
    //information into B
    for (int timeSlice=oldJoin+1;timeSlice<=newJoin;timeSlice++){
      for (int index=0;index<AB.extent(2);index++){
	AB(1,timeSlice,index)=AB(0,timeSlice,index);
      }
    }
  }
  else if (oldJoin>newJoin){ //CHANGED!
    for (int timeSlice=newJoin+1;timeSlice<=oldJoin;timeSlice++){
      for (int index=0;index<AB.extent(2);index++){
	if (FlipSign(index))
	  AB(0,timeSlice,PermPairIndex(index))=-AB(1,timeSlice,index);
	else
	  AB(0,timeSlice,PermPairIndex(index))=AB(1,timeSlice,index);
      }
    }
    //Now that we've copied the data from B into A, we need to copy the 
    //information into B
    for (int timeSlice=newJoin+1;timeSlice<=oldJoin;timeSlice++){
      for (int index=0;index<NumParticles();index++){
	AB(1,timeSlice,index)=AB(0,timeSlice,index);
      }
    }
  }
}
