#include "CommunicatorClass.h"

// Parallel MPI Communicator functions definitions
#ifdef PARALLEL




// Serial definition
#else
int Communicator::NumProcs() const
{
  return(1);
}

int Communicator::MyProc() const
{
  return (0);
}

void Communicator::SendReceive(int SendProc, const Array<dVec,1> &SendBuff,
			       int RecvProc,       Array<dVec,1> &RecvBuff)
{
  RecvBuff = SendBuff;
}

void Communicator::SendReceive(int SendProc, const Array<double,1> &SendBuff,
			       int RecvProc,       Array<double,1> &RecvBuff)
{
  RecvBuff = SendBuff;
}

void Communicator::SendReceive(int SendProc, const Array<int,1> &SendBuff,
			       int RecvProc,       Array<int,1> &RecvBuff)
{
  RecvBuff = SendBuff;
}

#endif
