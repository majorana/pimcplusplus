#ifndef COMMUNICATION_H
#define COMMUNICATION_H
#include "../../Common/Blitz.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace COMM
{
#ifdef USE_MPI
  inline void Init (int argc, char **argv)
    {
      MPI_Init(&argc, &argv);
    }
  inline void Finalize ()
    {
      MPI_Finalize();
    }
  inline int WorldProc()
  {
    int proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    return (proc);
  }
#else // Serial version
  inline void Init (int argc, char **argv)
    {
    }
  inline void Finalize ()
    {
    }
  inline int WorldProc()
  {
    return (0);
  }

#endif
}  


class CommStatusClass
{
public:
  int Source;
  int Tag;
  int Error;
  int Length;
};


class CommunicatorClass
{
public:
#ifdef USE_MPI
  ///If we are in parallel mode, we need an MPI communicator
  MPI_Comm MPIComm;  
  
  /// Sets this communicator to be that of all the processes
  /// (i.e. MPI_WORLD)
  void SetWorld();
  int MyProc();
  int NumProcs();
  void Send (void *sendBuf, int count, MPI_Datatype datatype,
	     int dest, int tag);
  void Send (int toProc, Array<double,1> &buff);
  void BroadCast (int root, Array<double,1> &buff);
  void Receive (void *recvBuf, int count, MPI_Datatype datatype,
		int source, int tag);
  void Receive (int toProc, Array<double,1> &buff);
  bool Probe(int source, int tag, CommStatusClass &status);
  void AllGather(void *sendbuf, int sendcount, 
		 MPI_Datatype sendtype, 
		       void* recvbuf, int recvcount,
		 MPI_Datatype recvtype);
  void AllGather (Array<double,1> &SendVec, 
		  Array<double,1> &RecvVec);
  void Split (int color, CommunicatorClass &newComm);
  void Subset (Array<int,1> ranks, CommunicatorClass &newComm);

  ///Sends and receives an array of dVec
  void SendReceive (int sendProc, const Array<Vec3,1> &sendBuff,
		    int recvProc,       Array<Vec3,1> &recvBuff);

  ///Sends and receives an array of dVec
  void SendReceive (int sendProc, const Array<Vec2,1> &sendBuff,
		    int recvProc,       Array<Vec2,1> &recvBuff);


  ///Sends and receives an array of double
  void SendReceive (int sendProc, const Array<double,1> &sendBuff,
		    int recvProc,       Array<double,1> &recvBuff);

  ///Sends and receives an array of complex
  void SendReceive (int sendProc, const Array<complex<double>,1> &sendBuff,
		    int recvProc,       Array<complex<double>,1> &recvBuff);
  
  ///Sends and receives an array of int
  void SendReceive (int sendProc, const Array<int,1> &sendBuff,
		    int recvProc,       Array<int,1> &recvBuff);



  ///Sums up the vectors in sendBuff.  Processor 0 only gets the
  ///resulting sum.
  void Sum (Array<int,1> &sendBuff, Array<int,1> &recvBuff);

  ///Sums up the vectors in sendBuff.  Processor 0 only gets the
  ///resulting sum.
  void Sum (Array<double,1> &sendBuff, Array<double,1> &recvBuff);
  ///Sums up all values a.  Only processor 0 gets the result.  All
  ///other processors return 0;
  double Sum (double a);

  CommunicatorClass()
  {
    SetWorld();
  }
  

#else   // Serial version
  inline void SetWorld()
  {
    // Do nothing
  }
  inline int MyProc()
  {
    return 0;
  }
  inline int NumProcs()
  {
    return 1;
  }
  inline void AllGather (Array<double,1> &SendVec, 
			 Array<double,1> &RecvVec)
  {
    RecvVec = SendVec;
  }
  inline void Split (int color, CommunicatorClass &newComm)
  {
    // Do nothing
  }
  inline void Subset (Array<int,1> ranks, CommunicatorClass &newComm)
  {
    if (ranks.size() !=1) {
      cerr << "Serial verion of code does not support nontrivial "
	   << "subsets.  Exiting.\n";
      exit(1);
    }
  }
  inline void Send (int toProc, Array<double,1> &buff)
  {
    cerr << "Sends not supported in serial mode.\n";
    exit(1);
  }
  inline void BroadCast(int root, Array<double,1> &buff)
  { /* Do nothing in serial mode */ }
  inline void Receive (int toProc, Array<double,1> &buff)
  {
    cerr << "Receives not supported in serial mode.\n";
    exit(1);
  }

  ///Sends and receives an array of dVec
  void SendReceive (int sendProc, const Array<Vec3,1> &sendBuff,
		    int recvProc,       Array<Vec3,1> &recvBuff)
  {
    recvBuff=sendBuff;
  }
  
  ///Sends and receives an array of dVec
  void SendReceive (int sendProc, const Array<Vec2,1> &sendBuff,
		    int recvProc,       Array<Vec2,1> &recvBuff)
  {
    recvBuff=sendBuff;
  }


  ///Sends and receives an array of double
  void SendReceive (int sendProc, const Array<double,1> &sendBuff,
		    int recvProc,       Array<double,1> &recvBuff)
  {
    recvBuff=sendBuff;
  }
  
  ///Sends and receives an array of complex
  void SendReceive (int sendProc, const Array<complex<double>,1> &sendBuff,
		    int recvProc,       Array<complex<double>,1> &recvBuff)
  {
    recvBuff=sendBuff;
  }
  ///Sends and receives an array of int
  void SendReceive (int sendProc, const Array<int,1> &sendBuff,
		    int recvProc,       Array<int,1> &recvBuff)
  {
    recvBuff=sendBuff;
  }



  ///Sums up all values a.  Only processor 0 gets the result.  All
  ///other processors return 0;
  double Sum (double a)
  {
    return a;
  }

  ///Sums up the vectors in sendBuff.  Processor 0 only gets the
  ///resulting sum.
  void Sum (Array<int,1> &sendBuff, Array<int,1> &recvBuff)
  {
    recvBuff=sendBuff;
  }
  void Sum (Array<double,1> &sendBuff, Array<double,1> &recvBuff)
  {
    recvBuff=sendBuff;
  }
  


#endif
};  // End class CommunicatorClass
#endif // End ifdef COMMUNICATION_H

