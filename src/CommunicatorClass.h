#ifndef COMM_CLASS_H
#define COMM_CLASS_H

#include "Common.h"

#ifdef PARALLEL
#include <mpi.h>
#endif


class CommunicatorClass
{
 public:
#ifdef PARALLEL
  MPI_Comm my_mpi_comm;  
#endif

  int NumProcs() const;
  int MyProc() const;

  void SendReceive (int SendProc, const Array<dVec,1> &SendBuff,
		    int RecvProc,       Array<dVec,1> &RecvBuff);

  void SendReceive (int SendProc, const Array<double,1> &SendBuff,
		    int RecvProc,       Array<double,1> &RecvBuff);

  void SendReceive (int SendProc, const Array<int,1> &SendBuff,
		    int RecvProc,       Array<int,1> &RecvBuff);
};


#endif
