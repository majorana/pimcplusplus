#include "CommunicatorClass.h"

// Parallel MPI Communicator functions definitions
#ifdef PARALLEL

int CommClass::NumProcs() const
{
  int size;
  MPI_Comm_size(my_mpi_comm, &size);
  return (size);
}

int CommClass::MyProc() const
{
  int proc;
  MPI_Comm_rank(my_mpi_comm, &proc);
  return (proc);
}

void CommClass::SendReceive(int SendProc, const Array<dVec,1> &SendBuff,
			    int RecvProc,  Array<dVec,1> &RecvBuff)
{
  double *SendPtr = (double *)SendBuff.data();
  double *RecvPtr = (double *)RecvBuff.data();
  int NumSend = SendBuff.size()*NDIM;
  int NumRecv = RecvBuff.size()*NDIM;
  MPI_Status status;
  MPI_Sendrecv (SendPtr, NumSend, MPI_DOUBLE, SendProc, 1, 
		RecvPtr, NumRecv, MPI_DOUBLE, RecvProc, 1,
		my_mpi_comm, &status);
}



void CommClass::SendReceive(int SendProc, const Array<double,1> &SendBuff,
			    int RecvProc, Array<double,1> &RecvBuff)
{
  double *SendPtr = (double *)SendBuff.data();
  double *RecvPtr = RecvBuff.data();
  int NumSend = SendBuff.size();
  int NumRecv = RecvBuff.size();
  MPI_Status status;
  MPI_Sendrecv (SendPtr, NumSend, MPI_DOUBLE, SendProc, 2, 
		RecvPtr, NumRecv, MPI_DOUBLE, RecvProc, 2,
		my_mpi_comm, &status);
}



void CommClass::SendReceive(int SendProc, const Array<int,1> &SendBuff,
			    int RecvProc, Array<int,1> &RecvBuff)
{
  int *SendPtr = (int *) SendBuff.data();
  int *RecvPtr = RecvBuff.data();
  int NumSend = SendBuff.size();
  int NumRecv = RecvBuff.size();
  MPI_Status status;
  MPI_Sendrecv (SendPtr, NumSend, MPI_INT, SendProc, 3, 
		RecvPtr, NumRecv, MPI_INT, RecvProc, 3,
		my_mpi_comm, &status);
}





// Serial definition
#else
int CommClass::NumProcs() const
{
  return(1);
}

int CommClass::MyProc() const
{
  return (0);
}

void CommClass::SendReceive(int SendProc, const Array<dVec,1> &SendBuff,
			       int RecvProc,       Array<dVec,1> &RecvBuff)
{
  RecvBuff = SendBuff;
}

void CommClass::SendReceive(int SendProc, const Array<double,1> &SendBuff,
			       int RecvProc,       Array<double,1> &RecvBuff)
{
  RecvBuff = SendBuff;
}

void CommClass::SendReceive(int SendProc, const Array<int,1> &SendBuff,
			       int RecvProc,       Array<int,1> &RecvBuff)
{
  RecvBuff = SendBuff;
}

#endif
