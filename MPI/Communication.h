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


class CommunicatorClass
{
public:
#ifdef USE_MPI
  ///If we are in parallel mode, we need an MPI communicator
  MPI_Comm MPIComm;  
  
  /// Sets this communicator to be that of all the processes
  /// (i.e. MPI_WORLD)
  inline void SetWorld()
  {
    MPIComm = MPI_COMM_WORLD;
  }

  inline int MyProc()
    {
      int MyRank;
      MPI_Comm_rank(MPIComm, &MyRank);
      return MyRank;
    }
  inline int NumProcs()
    {
      int NumProcs;
      MPI_Comm_size(MPIComm, &NumProcs);
      return NumProcs;
    }
  inline void Send (void *sendBuf, int count, MPI_Datatype datatype,
		    int dest, int tag)
  {
    MPI_Send (sendBuf, count, datatype, dest, tag, MPIComm);
  }

  inline void Send (int toProc, Array<double,1> buff)
  {
    Send(buff.data(), buff.size(), MPI_DOUBLE, toProc, 1);
  }

  inline void Receive (void *recvBuf, int count, MPI_Datatype datatype,
		       int source, int tag)
  {
    MPI_Status status;
    MPI_Recv (recvBuf, count, datatype, source, tag, MPIComm, &status);
  }

  inline void Receive (int toProc, Array<double,1> buff)
  {
    Receive(buff.data(), buff.size(), MPI_DOUBLE, toProc, 1);
  }




  inline void AllGather(void *sendbuf, int sendcount, 
		       MPI_Datatype sendtype, 
		       void* recvbuf, int recvcount,
		       MPI_Datatype recvtype)
    {
      MPI_Allgather (sendbuf, sendcount, sendtype, 
		     recvbuf, recvcount, recvtype, MPIComm);
    }

  inline void AllGather (Array<double,1> &SendVec, 
			 Array<double,1> &RecvVec)
    {
      int numProcs = NumProcs();
      assert ((RecvVec.size()/SendVec.size())==numProcs);	
      int count = SendVec.size();
      double *sendbuf = SendVec.data();
      double *recvbuf = RecvVec.data();

      AllGather(sendbuf, count, MPI_DOUBLE,
		recvbuf, count, MPI_DOUBLE);
    }

  inline void Split (int color, CommunicatorClass &newComm)
  {
    MPI_Comm_split(MPIComm, color, 0, &(newComm.MPIComm));
  }
  inline void Subset (Array<int,1> ranks, CommunicatorClass &newComm)
  {
    MPI_Group myGroup, newGroup;
    MPI_Comm_group (MPIComm, &myGroup);
    MPI_Group_incl(myGroup, ranks.size(), ranks.data(), &newGroup);
    MPI_Comm_create(MPIComm, newGroup, &(newComm.MPIComm));
  }

#else   // Serial version
  inline int MyNode()
    {
      return 0;
    }
  inline int NumNodes()
    {
      return 1;
    }
  inline void AllGatherVector (Array<double,1> &SendVec, 
			       Array<double,1> &RecvVec)
    {
      RecvVec = SendVec;
    }

#endif
};  // End class CommunicatorClass
#endif // End ifdef COMMUNICATION_H

