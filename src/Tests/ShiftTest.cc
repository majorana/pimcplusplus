#include "../MirroredArrayClass.h"
#include <unistd.h>







void ShiftTest()
{
  const int numSlices=10;
  const int numParticles=1;
  CommunicatorClass myCommunicator;
#ifdef PARALLEL
  myCommunicator.my_mpi_comm = MPI_COMM_WORLD;
#endif
  cerr<<"The number of processors is "<<myCommunicator.NumProcs()<<endl;
  int myProc = myCommunicator.MyProc();
  int numProcs=myCommunicator.NumProcs();

  MirroredArrayClass<int> myArray(numSlices,numParticles);
  for (int slice=0;slice<numSlices;slice++){
    for (int ptcl=0;ptcl<numParticles;ptcl++){
      myArray.Set(slice,ptcl,myProc*(numSlices-1)+slice);
    }
  }
  if (myProc==myCommunicator.NumProcs()-1){
    for (int ptcl=0;ptcl<numParticles;ptcl++){
      myArray.Set(numSlices-1,ptcl,0);
    }
  }

  for (int proc=0;proc<numProcs;proc++){
    if (myProc == proc){
      cout<<"Proc number "<<myProc<<" printing A:"<<endl; 
      myArray.Print();
    }
    cout<<endl<<endl;
    sleep(2);
  }

  myArray.ShiftData(3,myCommunicator);
  cout<<endl<<endl;

  for (int proc=0;proc<numProcs;proc++){
    if (myProc == proc){
      cout<<"Proc number "<<myProc<<" printing B:"<<endl; 
      myArray.Print();
    }
    cout<<endl<<endl;
    sleep(2);
  }
}



void ShiftTestSymmetric()
{
  const int numSlices=10;
  const int numParticles=1;
  CommunicatorClass myCommunicator;
#ifdef PARALLEL
  myCommunicator.my_mpi_comm = MPI_COMM_WORLD;
#endif
  cerr<<"The number of processors is "<<myCommunicator.NumProcs()<<endl;
  int myProc = myCommunicator.MyProc();
  int numProcs=myCommunicator.NumProcs();

  MirroredSymmetricMatrixClass<int> myArray(numSlices,numParticles);
  for (int slice=0;slice<numSlices;slice++){
    for (int ptcl=0;ptcl<numParticles;ptcl++){
      myArray.Set(slice,ptcl,myProc*(numSlices-1)+slice);
    }
  }
  if (myProc==myCommunicator.NumProcs()-1){
    for (int ptcl=0;ptcl<numParticles;ptcl++){
      myArray.Set(numSlices-1,ptcl,0);
    }
  }

  for (int proc=0;proc<numProcs;proc++){
    if (myProc == proc){
      cout<<"Proc number "<<myProc<<" printing A:"<<endl; 
      myArray.Print();
    }
    cout<<endl<<endl;
    sleep(2);
  }

  myArray.ShiftData(3,myCommunicator);
  cout<<endl<<endl;

  for (int proc=0;proc<numProcs;proc++){
    if (myProc == proc){
      cout<<"Proc number "<<myProc<<" printing B:"<<endl; 
      myArray.Print();
    }
    cout<<endl<<endl;
    sleep(2);
  }
}



void ShiftTestAntiSymmetric()
{
  const int numSlices=10;
  const int numParticles=1;
  CommunicatorClass myCommunicator;
#ifdef PARALLEL
  myCommunicator.my_mpi_comm = MPI_COMM_WORLD;
#endif
  cerr<<"The number of processors is "<<myCommunicator.NumProcs()<<endl;
  int myProc = myCommunicator.MyProc();
  int numProcs=myCommunicator.NumProcs();

  MirroredAntiSymmetricMatrixClass<int> myArray(numSlices,numParticles);
  for (int slice=0;slice<numSlices;slice++){
    for (int ptcl=0;ptcl<numParticles;ptcl++){
      myArray.Set(slice,ptcl,myProc*(numSlices-1)+slice);
    }
  }
  if (myProc==myCommunicator.NumProcs()-1){
    for (int ptcl=0;ptcl<numParticles;ptcl++){
      myArray.Set(numSlices-1,ptcl,0);
    }
  }

  for (int proc=0;proc<numProcs;proc++){
    if (myProc == proc){
      cout<<"Proc number "<<myProc<<" printing A:"<<endl; 
      myArray.Print();
    }
    cout<<endl<<endl;
    sleep(2);
  }

  myArray.ShiftData(3,myCommunicator);
  cout<<endl<<endl;

  for (int proc=0;proc<numProcs;proc++){
    if (myProc == proc){
      cout<<"Proc number "<<myProc<<" printing B:"<<endl; 
      myArray.Print();
    }
    cout<<endl<<endl;
    sleep(2);
  }
}





int main(int argc, char** argv )
{

  MPI_Init(&argc, &argv);
  cerr<<"Normal"<<endl;
  ShiftTest();
  cerr<<"Symmetric"<<endl;
  ShiftTestSymmetric();
  cerr<<"Anti-Symmetric"<<endl;
  ShiftTestAntiSymmetric();

  MPI_Finalize();






}
