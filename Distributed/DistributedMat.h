#ifndef DISTRIBUTED_MAT_H
#define DISTRIBUTED_MAT_H

#include "Common/MPI/Communication.h"

class DistributedMat
{
private:
  CommunicatorClass MyComm;
  Array<double,2> Mat;
public:
  inline double operator()(int row, int col) const
  {
    return Mat(row,col);
  }
  inline double& operator()(int row, int col)
  {
    return Mat(row,col);
  }
  inline void Resize(int rows, int cols)
  { 
    Mat.resize(rows,cols); 
  }
  /// Returns the number of elements this processor is responsible for
  inline int NumElements(int proc)
  {
    int totalElements = Mat.size();
    int numProcs = MyComm.NumProcs();
    int Elements  = totalElements/numProcs;
    if ((totalElements%numProcs)>proc)
      Elements++;
    return (Elements);
  }
  inline int MyNumElements()
  {
    return (NumElements(MyComm.MyProc()));
  }
  
  /// Sets the row and column of the i'th element process proc
  /// is responsible for.
  inline void Element(int proc, int index, int &row, int &col)
  {
    int NumProcs = MyComm.NumProcs();
    int elem = index*NumProcs + proc;
    row = elem/Mat.cols();
    col = elem%Mat.cols();
  }    
      
  /// Returns the row and column of the ith element this processor is
  /// responsible for. 
  inline void MyElement(int i, int &row, int &col)
  {
    int MyProc = MyComm.MyProc();
    Element(MyProc, i, row, col);
  }
  
  inline void Print()
  {
    cerr << "MyProc = " << MyComm.MyProc() << endl;
    cerr << Mat << endl;
  }

  /// Gathers the elements from all the processors to all the
  /// processors 
  void AllGather();

  DistributedMat(CommunicatorClass comm)
  {
    MyComm = comm;

  }
  DistributedMat(int rows, int cols, CommunicatorClass comm)
  {
    MyComm = comm;
    Mat.resize(rows,cols);
  }
};





/// Symmetric version of previous class
class DistributedSymmMat
{
private:
  CommunicatorClass MyComm;
  /// Only store lower triangular part
  Array<double,1> LowerTri;
  inline int index(int row, int col) const
  {
    if (col > row)
      return (((col*(col+1))>>1) + row);
    else
      return (((row*(row+1))>>1) + col);
  }
  int NumRows;
public:
  inline double operator()(int row, int col) const
  {
    return LowerTri(index(row,col));
  }
  inline double& operator()(int row, int col)
  {
    return LowerTri(index(row,col));
  }
  inline void Resize(int rows)
  { 
    NumRows = rows;
    LowerTri.resize((rows*(rows+1))/2);
  }
  /// Returns the number of elements this processor is responsible for
  inline int NumElements(int proc)
  {
    int totalElements = LowerTri.size();
    int numProcs = MyComm.NumProcs();
    int Elements  = totalElements/numProcs;
    if ((totalElements%numProcs)>proc)
      Elements++;
    return (Elements);
  }
  inline int MyNumElements()
  {
    return (NumElements(MyComm.MyProc()));
  }
  
  /// Sets the row and column of the i'th element process proc
  /// is responsible for.
  inline void Element(int proc, int i, int &row, int &col)
  {
    int NumProcs = MyComm.NumProcs();
    int elem = i*NumProcs + proc;
    // Add the 0.5 to elem to avoid tiny errors in sqrt
    row = (int)floor(0.5*(-1.0+sqrt(1.0+8*(elem+0.5))));
    col = elem - (row*(row+1)>>1);
  }    
      
  /// Returns the row and column of the ith element this processor is
  /// responsible for. 
  inline void MyElement(int i, int &row, int &col)
  {
    int MyProc = MyComm.MyProc();
    Element(MyProc, i, row, col);
  }
  
  inline void Print()
  {
    // Quick hack -- copy into a full matrix
    Array<double,2> Mat(NumRows, NumRows);
    for (int row=0; row<NumRows; row++)
      for (int col=0; col<NumRows; col++)
	Mat(row,col) = (*this)(row,col);
    cerr << "MyProc = " << MyComm.MyProc() << endl;
    cerr << Mat << endl;
  }

  /// Gathers the elements from all the processors to all the
  /// processors 
  void AllGather();

  DistributedSymmMat(CommunicatorClass comm)
  {
    MyComm = comm;
  }
  DistributedSymmMat(int rows, CommunicatorClass comm)
  {
    MyComm = comm;
    Resize(rows);
  }
};


class DistributedArray3
{
private:
  CommunicatorClass MyComm;
  Array<double,3> Mat;
public:
  inline double operator()(int i, int j, int k) const
  {
    return Mat(i,j,k);
  }
  inline double& operator()(int i, int j, int k)
  {
    return Mat(i,j,k);
  }
  inline void Resize(int rows, int cols, int depth)
  { 
    Mat.resize(rows,cols,depth); 
  }
  /// Returns the number of elements this processor is responsible for
  inline int NumElements(int proc)
  {
    int totalElements = Mat.extent(0)*Mat.extent(1);
    int numProcs = MyComm.NumProcs();
    int Elements  = totalElements/numProcs;
    if ((totalElements%numProcs)>proc)
      Elements++;
    return (Elements);
  }
  inline int MyNumElements()
  {
    return (NumElements(MyComm.MyProc()));
  }
  
  /// Sets the row and column of the i'th element process proc
  /// is responsible for.
  inline void Element(int proc, int index, int &row, int &col)
  {
    int NumProcs = MyComm.NumProcs();
    int elem = index*NumProcs + proc;
    row = elem/Mat.cols();
    col = elem%Mat.cols();
  }    
      
  /// Returns the row and column of the ith element this processor is
  /// responsible for. 
  inline void MyElement(int i, int &row, int &col)
  {
    int MyProc = MyComm.MyProc();
    Element(MyProc, i, row, col);
  }
  
  inline void Print()
  {
    cerr << "MyProc = " << MyComm.MyProc() << endl;
    cerr << Mat << endl;
  }

  /// Gathers the elements from all the processors to all the
  /// processors 
  void AllGather();

  DistributedMat(CommunicatorClass comm)
  {
    MyComm = comm;

  }
  DistributedArray3(int rows, int cols, int depth, CommunicatorClass comm)
  {
    MyComm = comm;
    Mat.resize(rows,cols,depth);
  }
};


#endif
