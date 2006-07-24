/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef COMMUNICATION_H
#define COMMUNICATION_H
#include "../Blitz.h"
#include <fstream>

#ifdef USE_MPI
#include <mpi.h>
#endif

extern ostream perr;

namespace COMM
{
#ifdef USE_MPI
  inline void Init (int argc, char **argv)
  {
    MPI_Init(&argc, &argv);
    int proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    if (proc == 0)
      perr.rdbuf(cerr.rdbuf());
//      else 
//        perr.open("/dev/null", ios::out); 
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
    perr.rdbuf(cerr.rdbuf());
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
  void Send (int toProc, Array<int,1> &buff);
  void Broadcast (int root, int &val);
  void Broadcast (int root, bool &val);
  void Broadcast (int root, double &val);
  void Broadcast (int root, Array<int,1> &buff)
  void Broadcast (int root, Array<double,1> &buff);
  void Broadcast (int root, Array<double,2> &buff);
  void Broadcast (int root, Array<double,3> &buff);
  void Broadcast (int root, Array<complex<double>,2> &buff);
  void Broadcast (int root, Array<Vec2,1> &buff);
  void Broadcast (int root, Array<Vec3,1> &buff);
  void Broadcast (int root, Array<double,4> &buff);
  void Broadcast (int root, Array<complex<double>,4> &buff);
  void Receive (void *recvBuf, int count, MPI_Datatype datatype,
		int source, int tag);
  void Receive (int toProc, Array<double,1> &buff);
  void Receive (int toProc, Array<int,1> &buff);
  bool Probe(int source, int tag, CommStatusClass &status);
  void Gather (Array<complex<double>,1> &sendVec, 
	       Array<complex<double>,1> &recvVec, 
	       Array<int,1> &recvCounts, int root=0);
  void Gather (Array<TinyVector<double,3>,1> &sendVec, 
	       Array<TinyVector<double,3>,1> &recvVec, 
	       Array<int,1> &recvCounts, int root=0);
  void Gather (Array<double,1> &sendVec, Array<double,2> &recvMat,
	       int root=0);
  void Gather (Array<int,1> &sendVec, Array<int,2> &recvMat,
	       int root=0);
  void AllGather(void *sendbuf, int sendcount, 
		 MPI_Datatype sendtype, 
		       void* recvbuf, int recvcount,
		 MPI_Datatype recvtype);
  void AllGather (Array<double,1> &SendVec, 
		  Array<double,1> &RecvVec);
  void AllGather (Array<int,1> &SendVec, 
		  Array<int,1> &RecvVec);
  /// This function gathers all rows of a matrix to all processors.
  /// It assumes that each processor receives an equal number of rows,
  /// with any left over being distributed to the low number
  /// processors. E.g. If there are 8 rows and 3 processors, procs 0
  /// and 1 would get 3 rows and proc 2 would get 2 rows.
  void AllGatherRows (Array<complex<double>,2> &mat);
  void AllGatherRows (Array<double,2> &mat);

  /// This function uses the same division stragegy as above, but
  /// gathers single elements, instead.
  void AllGatherVec (Array<double,1> &vec);

  void Split (int color, CommunicatorClass &newComm);
  void Subset (Array<int,1> &ranks, CommunicatorClass &newComm);

//   ///Sends and receives an array of dVec
//   void SendReceive (int sendProc, const Array<Vec3,1> &sendBuff,
// 		    int recvProc,       Array<Vec3,1> &recvBuff);

//   ///Sends and receives an array of dVec
//   void SendReceive (int sendProc, const Array<Vec2,1> &sendBuff,
//		    int recvProc,       Array<Vec2,1> &recvBuff);

  template<int N>
  void SendReceive (int sendProc, const Array<double,N> &sendBuff,
		    int recvProc,       Array<double,N> &recvBuff)
  {
    MPI_Status status;
    MPI_Sendrecv((void*)sendBuff.data(), sendBuff.size(), MPI_DOUBLE, sendProc, 2,
		 (void*)recvBuff.data(), recvBuff.size(), MPI_DOUBLE, recvProc, 2,
		 MPIComm, &status);
  }

  template<int N, int M>
  void SendReceive (int sendProc, const Array<TinyVector<double,M>,N> &sendBuff,
		    int recvProc,       Array<TinyVector<double,M>,N> &recvBuff)
  {
    MPI_Status status;
    MPI_Sendrecv((void*)sendBuff.data(), M*sendBuff.size(), MPI_DOUBLE, sendProc, 2,
		 (void*)recvBuff.data(), M*recvBuff.size(), MPI_DOUBLE, recvProc, 2,
		 MPIComm, &status);
  }
  
  template<int N>
  void SendReceive (int sendProc, const Array<complex<double>,N> &sendBuff,
		    int recvProc,       Array<complex<double>,N> &recvBuff)
  {
    MPI_Status status;
    MPI_Sendrecv((void*)sendBuff.data(), 2*sendBuff.size(), MPI_DOUBLE, sendProc, 2,
		 (void*)recvBuff.data(), 2*recvBuff.size(), MPI_DOUBLE, recvProc, 2,
		 MPIComm, &status);
  }

  template<int N, int M>
  void SendReceive (int sendProc, const Array<TinyVector<complex<double>,M>,N> &sendBuff,
		    int recvProc,       Array<TinyVector<complex<double>,M>,N> &recvBuff)
  {
    MPI_Status status;
    MPI_Sendrecv((void*)sendBuff.data(), 2*M*sendBuff.size(), MPI_DOUBLE, sendProc, 2,
		 (void*)recvBuff.data(), 2*M*recvBuff.size(), MPI_DOUBLE, recvProc, 2,
		 MPIComm, &status);
  }

  template<int N>
  void SendReceive (int sendProc, const Array<int,N> &sendBuff,
		    int recvProc,       Array<int,N> &recvBuff)
  {
    MPI_Status status;
    MPI_Sendrecv((void*)sendBuff.data(), sendBuff.size(), MPI_INT, sendProc, 3,
		 (void*)recvBuff.data(), recvBuff.size(), MPI_INT, recvProc, 3,
		 MPIComm, &status);
  }

  template<int N, int M>
  void SendReceive (int sendProc, const Array<TinyVector<int,M>,N> &sendBuff,
		    int recvProc,       Array<TinyVector<int,M>,N> &recvBuff)
  {
    MPI_Status status;
    MPI_Sendrecv((void*)sendBuff.data(), M*sendBuff.size(), MPI_INT, sendProc, 3,
		 (void*)recvBuff.data(), M*recvBuff.size(), MPI_INT, recvProc, 3,
		 MPIComm, &status);
  }

//   ///Sends and receives an array of double
//   void SendReceive (int sendProc, const Array<double,1> &sendBuff,
// 		    int recvProc,       Array<double,1> &recvBuff);

//   ///Sends and receives an array of double
//   void SendReceive (int sendProc, const Array<double,2> &sendBuff,
// 		    int recvProc,       Array<double,2> &recvBuff);

//   ///Sends and receives an array of double
//   void SendReceive (int sendProc, const Array<double,3> &sendBuff,
// 		    int recvProc,       Array<double,3> &recvBuff);


//   ///Sends and receives an array of complex
//   void SendReceive (int sendProc, const Array<complex<double>,1> &sendBuff,
// 		    int recvProc,       Array<complex<double>,1> &recvBuff);
  
//   ///Sends and receives an array of int
//   void SendReceive (int sendProc, const Array<int,1> &sendBuff,
// 		    int recvProc,       Array<int,1> &recvBuff);



  ///Sums up the vectors in sendBuff.  Processor 0 only gets the
  ///resulting sum.
  void Sum (Array<int,1> &sendBuff, Array<int,1> &recvBuff);

  ///Sums up the vectors in sendBuff.  Processor 0 only gets the
  ///resulting sum.
  void Sum (Array<double,1> &sendBuff, Array<double,1> &recvBuff);


  void Sum (Array<Vec2,1> &sendBuff, Array<Vec2,1> &recvBuff);
  void Sum (Array<Vec3,1> &sendBuff, Array<Vec3,1> &recvBuff);
  ///Sums up all values a.  Only processor 0 gets the result.  All
  ///other processors return 0;
  double Sum (double a);

  /// Sums up all values of a on all processors.  All processors
  ///  get result.
  double AllSum (double a);
  void AllSum (Array<double,1> &in, Array<double,1> &out);
  void AllSum (Array<double,2> &in, Array<double,2> &out);
  void AllSum (Array<double,3> &in, Array<double,3> &out);
  void AllSum (Array<TinyVector<double,2>,1> &in, 
	       Array<TinyVector<double,2>,1> &out);
  void AllSum (Array<TinyVector<double,3>,1> &in, 
	       Array<TinyVector<double,3>,1> &out);
  void AllAnd (bool &TorF);
  template<int N> 
  inline void AllMax (TinyVector<int,N> &vec) {
    TinyVector<int, N> outVec;
    MPI_Allreduce(&(vec[0]), &(outVec[0]), N, MPI_INT, MPI_MAX, MPIComm);
    vec = outVec;
  }
  void BarrierSync();
  void PrintSync();

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

  inline void Gather (Array<complex<double>,1> &sendVec, 
		      Array<complex<double>,1> &recvVec, 
		      Array<int,1>& recvCounts, int root=0)
  {
    recvVec = sendVec;
  }

  inline void Gather (Array<TinyVector<double,3>,1> &sendVec, 
		      Array<TinyVector<double,3>,1> &recvVec, 
		      Array<int,1>& recvCounts, int root=0)
  {
    recvVec = sendVec;
  }
  inline void Gather (Array<double,1> &sendVec, 
		      Array<double,2> &recvMat,
		      int root = 0)
  {
    assert (recvMat.rows() == NumProcs());
    assert (recvMat.cols() == sendVec.size());
    recvMat(0,Range::all()) = sendVec;
  }
  inline void Gather (Array<int,1> &sendVec, 
		      Array<int,2> &recvMat,
		      int root = 0)
  {
    assert (recvMat.rows() == NumProcs());
    assert (recvMat.cols() == sendVec.size());
    recvMat(0,Range::all()) = sendVec;
  }

  inline void AllGather (Array<double,1> &SendVec, 
			 Array<double,1> &RecvVec)
  {
    RecvVec = SendVec;
  }

  inline void AllGather (Array<int,1> &SendVec, 
			 Array<int,1> &RecvVec)
  {
    RecvVec = SendVec;
  }

  /// This function gathers all rows of a matrix to all processors.
  /// It assumes that each processor receives an equal number of rows,
  /// with any left over being distributed to the low number
  /// processors. E.g. If there are 8 rows and 3 processors, procs 0
  /// and 1 would get 3 rows and proc 2 would get 2 rows.
  inline void AllGatherRows (Array<double,2> &mat) 
  {
    // Do nothing
  }
  inline void AllGatherRows (Array<complex<double>,2> &mat) 
  {
    // Do nothing
  }

  /// This function uses the same division stragegy as above, but
  /// gathers single elements, instead.
  void AllGatherVec (Array<double,1> &vec) 
  {
    // do nothing
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
      abort();
    }
  }
  inline void Send (int toProc, Array<double,1> &buff)
  {
    cerr << "Sends not supported in serial mode.\n";
    abort();
  }

  inline void Send (int toProc, Array<int,1> &buff)
  {
    cerr << "Sends not supported in serial mode.\n";
    abort();
  }

  template<typename T>
  void Broadcast(int root, T &val) { }

  inline void Receive (int toProc, Array<double,1> &buff)
  { cerr << "Receives not supported in serial mode.\n";  abort(); }
  inline void Receive (int toProc, Array<int,1> &buff)
  { cerr << "Receives not supported in serial mode.\n";    abort(); }

  template<typename T>
  void SendReceive(int sendProc, const T &sendBuff,
		   int recvProc,       T &recvBuff)
  { recvBuff = sendBuff; }

  ///Sums up all values a.  Only processor 0 gets the result.  All
  ///other processors return 0;
  double Sum (double a)
  { return a; }

  template<class T, int N>
  inline void Sum(Array<T,N> &sendBuff, Array<T,N> &recvBuff)
  { recvBuff = sendBuff; }
  
  /// Sums up all values of a on all processors.  All processors get result.
  inline double AllSum (double a)
  {  return a; } 

  template<class T, int N>
  inline void AllSum (Array<T, N> &in, Array<T, N> &out)
  { out = in; }

  inline void AllAnd (bool &TorF)
  { /* do nothing */ }

  template<int N> 
  inline void
  AllMax (TinyVector<int,N> &vec) {
    // do nothing
  }


  inline void BarrierSync() {}
  inline void PrintSync() {}

#endif
};  // End class CommunicatorClass
#endif // End ifdef COMMUNICATION_H

