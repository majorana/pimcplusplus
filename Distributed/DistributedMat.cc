#include "DistributedMat.h"

void DistributedMat::AllGather()
{
  int maxElements = NumElements(0);
  int numProcs = MyComm.NumProcs();
  int row, col;

  // Allocate send an receive buffers
  Array<double,1> SendBuffer(maxElements),
    RecvBuffer(maxElements*numProcs);
  
  // Now fill send buffer with my elements
  for (int i=0; i<MyNumElements(); i++)
    {
      MyElement(i, row, col);
      SendBuffer(i) = Mat(row,col);
    }

  // Do collective communication
  MyComm.AllGather(SendBuffer, RecvBuffer);

  // Now, put elements where they should be
  for (int proc=0; proc<numProcs; proc++)
    {
      for (int i=0; i<NumElements(proc); i++)
	{
	  Element(proc, i, row, col);
	  Mat(row,col) = RecvBuffer(i+maxElements*proc);
	}
    }
}



void DistributedSymmMat::AllGather()
{
  int maxElements = NumElements(0);
  int numProcs = MyComm.NumProcs();
  int row, col;

  // Allocate send an receive buffers
  Array<double,1> SendBuffer(maxElements),
    RecvBuffer(maxElements*numProcs);
  
  // Now fill send buffer with my elements
  for (int i=0; i<MyNumElements(); i++)
    {
      MyElement(i, row, col);
      SendBuffer(i) = (*this)(row,col);
    }

  // Do collective communication
  MyComm.AllGather(SendBuffer, RecvBuffer);

  // Now, put elements where they should be
  for (int proc=0; proc<numProcs; proc++)
    {
      for (int i=0; i<NumElements(proc); i++)
	{
	  Element(proc, i, row, col);
	  (*this)(row,col) = RecvBuffer(i+maxElements*proc);
	}
    }
}



void DistributedArray3::AllGather()
{
  int maxElements = NumElements(0);
  int numProcs = MyComm.NumProcs();
  int row, col;

  int depth = Mat.extent(2);
  // Allocate send an receive buffers
  Array<double,2> SendBuffer(maxElements*depth),
    RecvBuffer(maxElements*numProcs*depth);
  
  // Now fill send buffer with my elements
  for (int i=0; i<MyNumElements(); i++) {
    MyElement(i, row, col);
    for (int k=0; k<depth; k++)
      SendBuffer(i*depth+k) = Mat(row,col,k);
  }

  // Do collective communication
  MyComm.AllGather(SendBuffer, RecvBuffer);

  // Now, put elements where they should be
  for (int proc=0; proc<numProcs; proc++)
    {
      for (int i=0; i<NumElements(proc); i++)
	{
	  Element(proc, i, row, col);
	  for (int k=0; k<depth; k++)
	    Mat(row,col,k) = RecvBuffer((i+maxElements*proc)*depth+k);
	}
    }
}



void TestDistributedMat()
{
  CommunicatorClass comm;
  comm.SetWorld();
  
  DistributedSymmMat  dmat(5,comm);
  int row, col;

  for (int i=0; i<dmat.MyNumElements(); i++)
    {
      dmat.MyElement(i,row,col);
      dmat(row,col) = 10.0*row+col;
    }
  dmat.Print();
  dmat.AllGather();
  dmat.Print();
}


/* main(int argc, char **argv)
{
  COMM::Init(argc, argv);
  TestDistributedMat();
  COMM::Finalize();
}
*/
