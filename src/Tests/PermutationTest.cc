#include <fstream>
#include "../Common.h"
#include "../PathClass.h"
#include "../SpeciesClass.h"
#include "../ActionClass.h"
#include "../PathDataClass.h"
#include "../BisectionMoveClass.h"
//#include "../MirroredArrayClass.h"
#include "../Observables/ObservableClass.h"
#include "../Common/IO/InputOutput.h"




int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);



  PathDataClass myPathData;
#ifdef PARALLEL
  myPathData.Communicator.my_mpi_comm = MPI_COMM_WORLD;
#endif

  IOSectionClass inSection;
  assert(inSection.OpenFile("Permutation.in"));  
  //  inSection.PrintTree();
  cerr<<"I've opened the file\n";
  assert(inSection.OpenSection("System"));
  myPathData.Path.Read(inSection);
  inSection.CloseSection(); // "System"
  
  



  //Move Setup Hack
  BisectionMoveClass myBisectionMove(myPathData);
  ShiftMoveClass myShiftMove(myPathData);

  //Move Setup Done
  
  ///Here we are setting up distance table
  DistanceTablePBCClass *myDistTable=
    new DistanceTablePBCClass(myPathData.Path);
  myPathData.DistanceTable=myDistTable;
  myPathData.Action.DistanceTable=myDistTable;
  myPathData.DistanceTable->UpdateAll();
  ///Done setting up distance table
  myPathData.Path.Print();
  cerr<<"My join is at "<<myPathData.Join<<endl;

  SetMode(BOTHMODE);
  myPathData.Path.Permutation.Set(0,1);
  myPathData.Path.Permutation.Set(1,0);
  myPathData.Path.Permutation.Print();
  myPathData.MoveJoin(2);
  myPathData.Path.Print();
  cerr<<"My join is at "<<myPathData.Join<<endl;

}

