#include <fstream>
#include "Common.h"
#include "PathClass.h"
#include "SpeciesClass.h"
#include "ActionClass.h"
#include "PathDataClass.h"
#include "BisectionMoveClass.h"
#include "MirroredArrayClass.h"
#include "ObservableClass.h"
#include "DistanceTablePBCClass.h"
#include "DistanceTableFreeClass.h"
#include "Common/IO/InputOutput.h"
// #include "InputOutputASCII.h"

void setupAction(ActionClass &myActionClass)


{
  myActionClass.PairActionVector.resize(1);
  //  myActionClass.PairActionVector(0).ReadDavidSquarerFile("../inputs/ep_beta1.0.dm");
  myActionClass.PairMatrix.resize(2,2);
  myActionClass.PairMatrix=0;
  //myActionClass.tau=myActionClass.PairActionVector(0).tau;
  cerr << "Tau = " << myActionClass.tau << endl;
  //  myActionClass.mySpeciesArray =&(myPathData.SpeciesArray);

     
}






void setupIDParticleArray(PathDataClass &myPathData)
{

  //  SpeciesArrayClass &myIDParticles=myPathData.SpeciesArray;
  double tau=1.0; //This better be the same as in the squarer file! UGly!
  int NumTimeSlices=50;

  FermionClass *myElectronptr = new FermionClass;
  FermionClass *myProtonptr = new FermionClass;
  myPathData.Path.SetTimeSlices(NumTimeSlices);
  //  (*myElectronptr).NumParticles=1;
  //  (*myProtonptr).NumParticles=1;
  (*myElectronptr).lambda=0.5;
  (*myProtonptr).lambda=0;
  (*myElectronptr).NumParticles=1;
  (*myProtonptr).NumParticles=1;
  myPathData.Path.AddSpecies(myElectronptr);
  myPathData.Path.AddSpecies(myProtonptr);
  myPathData.Path.Allocate();
    
  //  (*myElectronptr).Path.Resize(1,NumTimeSlices);
  //  (*myProtonptr).Path.Resize(1,NumTimeSlices);
  SetMode(BOTHMODE);
  dVec zeroVector=0;
  for (int counter=(*myProtonptr).FirstPtcl;counter<=(*myProtonptr).LastPtcl;counter++){
    for (int counter2=0;counter2<NumTimeSlices;counter2++){
      myPathData.Path.SetPos(counter2,counter,zeroVector);
    }
  }
  //  for  (int counter=0;counter<(*myProtonptr).NumParticles();counter++){
  //    for (int counter2=0;counter2<NumTimeSlices;counter2++){
  //      (*myProtonptr).Path.SetPos(counter,counter2,zeroVector);
  //    }
  //  }
  double sigma=sqrt(2*0.5*tau);
  dVec electronVector;
  for (int counter=(*myElectronptr).FirstPtcl;counter<=(*myElectronptr).LastPtcl;counter++){
    for (int counter2=0;counter2<NumTimeSlices;counter2++){
      electronVector=GaussianRandomVec(sigma);
      myPathData.Path.SetPos(counter2,counter,electronVector);
    }
  }
  //  for (int counter=0;counter<(*myElectronptr).NumParticles();counter++){
  //    for (int counter2=0;counter2<NumTimeSlices;counter2++){
  //a      electronVector=GaussianRandomVec(sigma);
  //      (*myElectronptr).Path.SetPos(counter,counter2,electronVector);
  //    }
  //  }
  //  (*myElectronptr).InitPermMatrix();
  //  (*myProtonptr).InitPermMatrix();


  //  ProtonsClass &myProtons = *myProtonptr;
  //  ElectronsClass &myElectrons= *myElectronptr;
  //  myPathData.Path.SpeciesArray.Resize(2);
  //  myPathData.Path.SpeciesArray.Set(0,myElectrons);
  //  myPathData.Path.SpeciesArray.Set(1,myProtons);
  //  myPathData.NumTimeSlices=NumTimeSlices;
}
  

void setupMove(BisectionMoveClass &myBisectionMove,ShiftMoveClass &myShiftMove, PathDataClass &thePathData)
{

  Array<int,1> ActiveSpecies(1);
  ActiveSpecies(0) = 0;
  //  myBisectionMove.PathData=&thePathData;
  cerr << "ActiveSpecies = " << ActiveSpecies << endl;
  myBisectionMove.SetActiveSpecies(ActiveSpecies);
  myBisectionMove.SetNumParticlesToMove(1);
  myBisectionMove.StartTimeSlice=0;
  myBisectionMove.NumLevels=3;
  //  myShiftMove.PathData=&thePathData;

}






  
   


#include <mpi.h>

#include <unistd.h>


void TestShift()
{

 CommunicatorClass myCommunicator;
#ifdef PARALLEL
  myCommunicator.my_mpi_comm = MPI_COMM_WORLD;
#endif
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

  myArray.ShiftData(-29,myCommunicator);
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
  IOSectionClass inSection;
  assert(inSection.OpenFile("hydrogen.in"));  
  //  inSection.PrintTree();
  cerr<<"I've opened the file\n";
  assert(inSection.OpenSection("System"));
  myPathData.Path.Read(inSection);
  inSection.CloseSection(); // "System"
  cerr<<"I'm right before the action\n";
  inSection.OpenSection("Action");
  myPathData.Action.Read(inSection);
  inSection.CloseSection(); //"Action"

  // Let's test the action
  LinearGrid qgrid(0.0, 3.0, 31);
  for (int i=0; i<qgrid.NumPoints; i++) {
    double q = qgrid(i);
    double U = myPathData.Action.PairActionVector(0)->U(q,0.0, 4.0*q*q, 3);
    cerr << q << " " << U << endl;
  }

  IOSectionClass out;
  out.NewFile ("Observables2.h5");
  out.NewSection("Energies");
  TotalEnergyClass TotE(myPathData, out);
  out.CloseSection();
  out.NewSection("gofr");
  PairCorrelationClass gofr(myPathData,out);
  out.CloseSection();
  //Observable setup Hack!
  //PairCorrelation PC(myPathData);
  //  PC.PathData = &myPathData;
  //PC.Species1 = 0;
  //PC.Species2 = 1;
  //PC.Initialize();
  //Observable Setup Done

  //Move Setup Hack
  BisectionMoveClass myBisectionMove(myPathData);
  ShiftMoveClass myShiftMove(myPathData);
  setupMove(myBisectionMove,myShiftMove,myPathData);
  //Move Setup Done
  
  ///Here we are setting up distance table
  DistanceTablePBCClass *myDistTable=
    new DistanceTablePBCClass(myPathData.Path);
  myPathData.DistanceTable=myDistTable;
  myPathData.Action.DistanceTable=myDistTable;
  myPathData.DistanceTable->UpdateAll();
  ///Done setting up distance table

#ifdef PARALLEL
  myPathData.Communicator.my_mpi_comm = MPI_COMM_WORLD;
#endif
  

//   //  ActionClass myActionClass;
//   //  setupIDParticleArray(myPathData);
//   IOSectionClass *theInput=new InputSectionHDF5Class();
//   theInput->OpenFile("inputFile","nonameyet",NULL);
//   IOSectionClass *pathInput;
//   theInput->FindSection("PathInfo",pathInput,true);
//   pathInput->Rewind();
//   myPathData.Path.Read(pathInput);
//   //#ifdef PARALLEL
//   //  myPathData.Communicator.my_mpi_comm = MPI_COMM_WORLD;
//   //#endif

//   setupAction(myPathData.Action);
//   BisectionMoveClass myBisectionMove(myPathData);
//   ShiftMoveClass myShiftMove(myPathData);
//   setupMove(myBisectionMove,myShiftMove,myPathData);
//   myPathData.Path.Box[0]=0.25*3.14159265358979323846;
//   myPathData.Path.Box[1]=0.25*3.14159265358979323846;
//   myPathData.Path.Box[2]=0.25*3.14159265358979323846;
//   DistanceTablePBCClass *myDistTable=
//     new DistanceTablePBCClass(myPathData.Path);
//   myPathData.DistanceTable=myDistTable;

//   myPathData.Action.DistanceTable=myDistTable;
//   myPathData.DistanceTable->UpdateAll();
//   //  cerr<<"The size of the SpeciesArray is ";
//   //  cerr << (myBisectionMove.PathData)->SpeciesArray.size()<<endl;
//   //  cerr<<"What the action class thinks the size is: ";
//   //  cerr<<  myActionClass.mySpeciesArray->size()<<endl;
//   //  PrintConfigClass myPrintConfig(myPathData);
  ofstream outfile;
  outfile.open("ourPath.dat");
  int steps=2000000;
  for (int counter=0;counter<steps;counter++){
    if (counter>steps/8 && (counter % 100)==0){
      TotE.Accumulate();
      gofr.Accumulate();
    }
    if (counter>steps/8 && (counter % 1000) == 0){
      TotE.WriteBlock();
      cerr << "Step #" << counter << ":\n";
//       for (int slice=0;slice<myPathData.Path.NumTimeSlices();slice++){
// 	outfile<<myPathData.Path(slice,0)[0]<<" ";
// 	outfile<<myPathData.Path(slice,0)[1]<<" ";
// 	outfile<<myPathData.Path(slice,0)[2]<<" ";
// 	outfile<<endl;
//       }
//       outfile<<myPathData.Path(0,0)[0]<<" ";
//       outfile<<myPathData.Path(0,0)[1]<<" ";
//       outfile<<myPathData.Path(0,0)[2]<<" ";
//       outfile<<endl;
    }
      
    for (int counter2=0;counter2<2;counter2++){
      //cerr << "Doing step " << counter << endl;
      
      myBisectionMove.MakeMove();
      if (counter >= steps/3)
	;//PC.Accumulate();
      //      myPrintConfig.Print();
    }
    myShiftMove.MakeMove();
  }
  gofr.WriteBlock();
  out.CloseFile();
  //PC.Print();
  cout<<"My acceptance ratio is "<<myBisectionMove.AcceptanceRatio()<<endl;
  //cerr<<"done! done!"<<endl;
  MPI_Finalize();
  outfile.close();
  
}
  
