#include "Common.h"
#include "PathClass.h"
#include "IdenticalParticlesClass.h"
#include "ActionClass.h"
#include "PathDataClass.h"
#include "BisectionMoveClass.h"
#include "MirroredArrayClass.h"
#include "ObservableClass.h"


void setupAction(ActionClass &myActionClass,PathDataClass &myPathData)


{
  myActionClass.PairActionVector.resize(1);
  myActionClass.PairActionVector(0).ReadDavidSquarerFile("../inputs/ep_beta1.0.dm");
  myActionClass.PairMatrix.resize(2,2);
  myActionClass.PairMatrix=0;
  myActionClass.tau=myActionClass.PairActionVector(0).tau;
  cerr << "Tau = " << myActionClass.tau << endl;
  myActionClass.mySpeciesArray=&(myPathData.SpeciesArray);

     
}


void setupIDParticleArray(PathDataClass &myPathData)
{
  SpeciesArrayClass &myIDParticles=myPathData.SpeciesArray;
  double tau=1.0; //This better be the same as in the squarer file! UGly!
  int NumTimeSlices=50;
  ElectronsClass *myElectronptr = new ElectronsClass;
  ProtonsClass *myProtonptr = new ProtonsClass;
  ElectronsClass &myElectrons = *myElectronptr;
  ProtonsClass &myProtons = *myProtonptr;
  myElectrons.NumParticles=1;
  myProtons.NumParticles=1;
  myElectrons.lambda=0.5;
  myProtons.lambda=0;
   
  myElectrons.Path.resize(1,NumTimeSlices);
  myProtons.Path.resize(1,NumTimeSlices);
  setMode(BOTHMODE);
  dVec zeroVector=0;
  for  (int counter=0;counter<myProtons.NumParticles;counter++){
    for (int counter2=0;counter2<NumTimeSlices;counter2++){
      myProtons.Path.SetPos(counter,counter2,zeroVector);
    }
  }
  double sigma=sqrt(2*0.5*tau);
  dVec electronVector;
  for (int counter=0;counter<myElectrons.NumParticles;counter++){
    for (int counter2=0;counter2<NumTimeSlices;counter2++){
      electronVector=GaussianRandomVec(sigma);
      myElectrons.Path.SetPos(counter,counter2,electronVector);
    }
  }

  myPathData.SpeciesArray.resize(2);
  myPathData.SpeciesArray.Set(0,myElectrons);
  myPathData.SpeciesArray.Set(1,myProtons);
  myPathData.NumTimeSlices=NumTimeSlices;
}
  

void setupMove(BisectionMoveClass &myBisectionMove,ShiftMove &myShiftMove, PathDataClass &thePathData)
{

  Array<int,1> ActiveSpecies(1);
  ActiveSpecies(0) = 0;
  myBisectionMove.PathData=&thePathData;
  cerr << "ActiveSpecies = " << ActiveSpecies << endl;
  myBisectionMove.SetActiveSpecies(ActiveSpecies);
  myBisectionMove.SetNumParticlesToMove(1);
  myBisectionMove.StartTimeSlice=0;
  myBisectionMove.NumLevels=4;
  myShiftMove.PathData=&thePathData;

}






  
   


#include <mpi.h>

#include <unistd.h>


void TestShift()
{

 CommunicatorClass myCommunicator;
  myCommunicator.my_mpi_comm = MPI_COMM_WORLD;
  int MyProc = myCommunicator.MyProc();

  MirroredArrayClass<int> myArray(5,5);
  for (int counter=0;counter<5;counter++){
    for (int counter2=0;counter2<5;counter2++){
      myArray.Set(counter,counter2,counter*10+counter2+MyProc*5);
    }
  }
  if (myCommunicator.MyProc() == 0)
    myArray.Print();
  cout<<endl<<endl;
  sleep(2);
  if (myCommunicator.MyProc() == 1)
    myArray.Print();

  myArray.shiftData(-3,myCommunicator);
  cout<<endl<<endl;

  sleep(2);
  if (myCommunicator.MyProc() == 0)
    myArray.Print();
  
  sleep(2);
  
  if (myCommunicator.MyProc() == 1)
    myArray.Print();


}


int main(int argc, char **argv)

{
  MPI_Init(&argc, &argv);

 


  
  PathDataClass myPathData;
  PairCorrelation PC;
  PC.PathData = &myPathData;
  PC.Species1 = 0;
  PC.Species2 = 1;
  PC.Initialize();
  //  ActionClass myActionClass;
  setupIDParticleArray(myPathData);
  myPathData.Communicator.my_mpi_comm = MPI_COMM_WORLD;
  setupAction(myPathData.Action,myPathData);
  BisectionMoveClass myBisectionMove;
  ShiftMove myShiftMove;
  setupMove(myBisectionMove,myShiftMove,myPathData);
  //  cerr<<"The size of the SpeciesArray is ";
  //  cerr << (myBisectionMove.PathData)->SpeciesArray.size()<<endl;
  //  cerr<<"What the action class thinks the size is: ";
  //  cerr<<  myActionClass.mySpeciesArray->size()<<endl;
  for (int counter=0;counter<1000000;counter++){
    if ((counter % 1000) == 0)
      cerr << "Step #" << counter << ":\n";
    for (int counter2=0;counter2<2;counter2++){
      //cerr << "Doing step " << counter << endl;
      
      myBisectionMove.makeMove();
      if (counter > 100000)
      PC.Accumulate();
    }
    myShiftMove.makeMove();
  }
  PC.Print();
  //cerr<<"done! done!"<<endl;
  MPI_Finalize();
}
  
