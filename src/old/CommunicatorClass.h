#ifndef COMM_CLASS_H
#define COMM_CLASS_H

#include "Common.h"
#include "Common/MPI/Communication.h"




///Here we include mpi if we are doing parallel runs

///Communicator class. Communicates with other processors.
class PIMCCommunicatorClass : public CommunicatorClass
{
 public:

  ///Sends and receives an array of ImageNumClass
  using CommunicatorClass::SendReceive;
  void SendReceive (int sendProc, const Array<ImageNumClass,1> &sendBuff,
  		    int recvProc,       Array<ImageNumClass,1> &recvBuff);

  using CommunicatorClass::Sum;
  void Sum (Array<dVec,1> &sendBuff, Array<dVec,1> &recvBuff);
};


#endif
