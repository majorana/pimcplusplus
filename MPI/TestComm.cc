#include "Communication.h"

bool TestAllGatherRows()
{
  CommunicatorClass Comm;
  Comm.SetWorld();

  int cols = 3;
  int rows = 2*Comm.NumProcs()+1;
  Array<complex<double>,2> mat(rows,cols);
  mat = 0.0;
  int currRow = 0;
  for (int proc=0; proc<Comm.NumProcs(); proc++) {
    int procRows = rows/Comm.NumProcs() + ((rows%Comm.NumProcs())>proc);
    if (proc == Comm.MyProc()) {
      for (int row=currRow; row<(currRow+procRows); row++)
	for (int col=0; col<cols; col++)
	  mat(row,col) = (double)row + (double)col/10.0;
    }
    currRow += procRows;
  }
//   if (Comm.MyProc() == 0) 
//     cerr << "mat = " << mat << endl;

  Comm.AllGatherRows(mat);
  bool passed = true;
  for (int row=0; row<rows; row++)
    for (int col=0; col<cols; col++)
      passed = passed && 
	(fabs(mat(row,col).real()-((double)row + (double)col/10.0)) < 1.0e-12);
//   if (Comm.MyProc() == 0) 
//     cerr << "mat = " << mat << endl;
  return passed;
}


main(int argc, char *argv[])
{
  COMM::Init(argc, argv);
  CommunicatorClass comm;
  comm.SetWorld();

  bool passed = TestAllGatherRows();  
  if (comm.MyProc()==0)
    cerr << "AllGatherRows() check:  " 
	 << (passed ? "passed." : "failed.") << endl;
}
