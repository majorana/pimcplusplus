#include "DistributedMat.h"

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

void TestDistributedArray3()
{
  const int L=5;
  const int M=7;
  const int N=3;
  CommunicatorClass comm;
  comm.SetWorld();
  DistributedArray3 d3(L,M,N,comm);
  int row,col;
  for (int i=0; i<d3.MyNumElements(); i++) {
    d3.MyElement(i,row,col);
    for (int j=0; j<N; j++){
      d3(row,col,j) = 100.0*row + 10.0*col + (double)(j+1);
    }
  }
  d3.AllGather();
  bool passed = true;
  for (int i=0; i<L; i++)
    for (int j=0; j<M; j++)
      for (int k=0; k<N; k++) 
	if (d3(i,j,k) != (100.0*i + 10.0*j + (double)(k+1))) {
	  cerr << "Error in TestDistributedArray3 at i="
	       << i << " j=" << j << " k=" << k << endl;
	  cerr << d3(i,j,k) << endl;
	  passed = false;
	}
	
}


main(int argc, char **argv)
{
  COMM::Init(argc, argv);
  TestDistributedMat();
  TestDistributedArray3();
  COMM::Finalize();
}
