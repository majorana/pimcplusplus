#include "CommunicatorClass.h"


// Parallel MPI Communicator functions definitions
#ifdef USE_MPI


void PIMCCommunicatorClass::SendReceive(int SendProc, 
					const Array<ImageNumClass,1> &SendBuff,
					int RecvProc, 
					Array<ImageNumClass,1> &RecvBuff)
{
  int *SendPtr = (int *) SendBuff.data();
  int *RecvPtr = (int*) RecvBuff.data();
  int NumSend = SendBuff.size();
  int NumRecv = RecvBuff.size();
   MPI_Status status;
   MPI_Sendrecv (SendPtr, NumSend, MPI_INT, SendProc, 3, 
		 RecvPtr, NumRecv, MPI_INT, RecvProc, 3,
		 MPIComm, &status);
}



// Serial definition
#else

//// GRRRRRRRRRRRRRRRRRRRR!!!!!!!!!!!!!!111
void PIMCCommunicatorClass::SendReceive(int SendProc, 
					const Array<ImageNumClass,1> &SendBuff,
					int RecvProc,       
					Array<ImageNumClass,1> &RecvBuff)
{
  RecvBuff = SendBuff;
}



#endif
