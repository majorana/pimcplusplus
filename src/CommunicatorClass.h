#ifndef COMM_CLASS_H
#define COMM_CLASS_H

#include "Common.h"


///Here we include mpi if we are doing parallel runs
#ifdef PARALLEL
#include <mpi.h>
#endif

///Communicator class. Communicates with other processors.
class CommunicatorClass
{
 public:
#ifdef PARALLEL
  ///If we are in parallel mode, we need an MPI communicator
  MPI_Comm my_mpi_comm;  
#endif

  ///How many processors are there.  In serial mode there is 1 
  int NumProcs() const;

  // Number of my processor. In serial mode, youa re number 0.
  int MyProc() const;

  
  void ReceiveFromAll(const Array<dVec,1> &recvBuff);
  void Send(int sendProc,const Array<dVec,1> &sendBuff);

  ///Sends and receives an array of dVec
  void SendReceive (int sendProc, const Array<dVec,1> &sendBuff,
		    int recvProc,       Array<dVec,1> &recvBuff);

  ///Sends and receives an array of double
  void SendReceive (int sendProc, const Array<double,1> &sendBuff,
		    int recvProc,       Array<double,1> &recvBuff);

  ///Sends and receives an array of int
  void SendReceive (int sendProc, const Array<int,1> &sendBuff,
		    int recvProc,       Array<int,1> &recvBuff);

  ///Sends and receives an array of ImageNumClass
  void SendReceive (int sendProc, const Array<ImageNumClass,1> &sendBuff,
		    int recvProc,       Array<ImageNumClass,1> &recvBuff);

  ///Sums up the vectors in sendBuff.  Processor 0 only gets the
  ///resulting sum.
  void Sum (Array<double,1> &sendBuff, Array<double,1> &recvBuff);
  ///Sums up all values a.  Only processor 0 gets the result.  All
  ///other processors return 0;
  double Sum (double a);
};


#endif
