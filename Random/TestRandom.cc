#include "Random.h"

int main(int argc, char** argv)
{
  MPI_Init (&argc, &argv);
  CommunicatorClass comm;

  comm.SetWorld();
  RandomClass rand(comm);
  rand.Init();
  cerr << "Proc " << comm.MyProc() << "  local random num = " << rand.Local() << endl;
  cerr << "Proc " << comm.MyProc() << " common random num = " << rand.Common() << endl;
  MPI_Finalize();
}

  

  
