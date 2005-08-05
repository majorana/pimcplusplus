#include "Random.h"

int main(int argc, char** argv)
{
#ifdef USE_MPI
  MPI_Init (&argc, &argv);
#endif
  CommunicatorClass comm;

  comm.SetWorld();
  RandomClass rand(comm);
  rand.Init(314159);
  int numLocal;
  if (comm.MyProc()==0) 
    numLocal = 4;
  else
    numLocal = 6;
  for (int i=0; i<numLocal; i++)
      cerr << "Proc " << comm.MyProc() << "  local random num = " 
	   << rand.Local() << endl;

  cerr << "Proc " << comm.MyProc() << " common random num = " 
       << rand.Common() << endl;

#ifdef USE_MPI
  MPI_Finalize();
#endif
}

  

  
