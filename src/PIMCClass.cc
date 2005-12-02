#include "PIMCClass.h"
#include "Moves/MoveClass.h"
#include "Observables/ObservableClass.h"
#include <sstream>
#include <fstream>
#include <Common/Blitz.h>
#include <Common/IO/FileExpand.h>

void PIMCClass::Read(IOSectionClass &in)
{
  // Read the parallelization strategy
  PathData.Read (in);
  perr << "Finished PathData Read.\n";

  // Read in the system information and allocate the path
  assert(in.OpenSection("System"));
  PathData.Path.Read(in);
  in.CloseSection();
  perr << "Finished Path read.\n";

  // Read in the action information
  assert(in.OpenSection("Action"));
  PathData.Actions.Read(in);
  in.CloseSection();
  perr << "Finished Actions read.\n";

  // Now actually initialize the paths
  perr << "Before InitPaths.\n";
  assert(in.OpenSection("System"));
  PathData.Path.InitPaths(in);
  in.CloseSection();
  perr << "Done InitPaths.\n";
  if (PathData.Path.UseCorrelatedSampling())
    PathData.Path.SetIonConfig(0);
  
  perr << "Initializing Actions caches.\n";
  PathData.Actions.Init();
  perr << "done.\n";

  // Read in the Observables
  assert(in.OpenSection("Observables"));
  ReadObservables(in);
  perr << "Finished Observables Read.\n";
  in.CloseSection();

  bool iAmRootProc = (PathData.Path.Communicator.MyProc()==0);

  if (iAmRootProc)
    OutFile.NewSection("Actions");
  if (PathData.Actions.HaveLongRange()) {
    assert (in.OpenSection ("Action"));
    PathData.Actions.LongRange.Init (in, OutFile);
    if (PathData.Actions.UseRPA)
      PathData.Actions.LongRangeRPA.Init(in);
    in.CloseSection();
  }

  if (iAmRootProc) {
    PathData.Actions.WriteInfo(OutFile);
    OutFile.CloseSection(); // "Actions"
  }

  // Read in the Moves
  assert(in.OpenSection("Moves"));
  ReadMoves(in);
  perr << "Finished Moves Read.\n";
  in.CloseSection();

  // Read in the Algorithm
  assert(in.OpenSection("Algorithm"));
  ReadAlgorithm(in);
  in.CloseSection();
}




void PIMCClass::ReadObservables(IOSectionClass &in)
{
  int myProc=PathData.Path.Communicator.MyProc();
  bool iAmRoot= myProc==0;
  if (iAmRoot) {
    string outFileBase;
    assert(in.ReadVar("OutFileBase",outFileBase));
    int fileStart;
    if (!in.ReadVar("FileStart", fileStart))
      fileStart = 0;
    // Allow for tilde-expansion in these files
    outFileBase = ExpandFileName (outFileBase);
    ostringstream cloneNum;
    cloneNum << (PathData.GetCloneNum() + fileStart);
    OutFileName = 
      outFileBase+ "." + cloneNum.str() + ".h5";
    OutFile.NewFile(OutFileName);

    /////////////////////////////////////
    // Write input file to output file //
    /////////////////////////////////////
    ifstream infile;
    infile.open(in.GetFileName().c_str());
    infile.seekg(0,ios::end);
    int length=infile.tellg();
    infile.seekg(0,ios::beg);
    char *buffer=new char[length+1];
    infile.read(buffer,length);
    buffer[length] = '\0';
    infile.close(); 
    string fileCopy(buffer);   
    delete buffer;
    OutFile.WriteVar("InputFile",fileCopy);


    OutFile.NewSection("RunInfo");
    RunInfo.Write(OutFile);
    OutFile.CloseSection();
    OutFile.NewSection("System");
    WriteSystemInfo();
    OutFile.CloseSection(); // "System" 
    OutFile.NewSection ("Observables");
    Array<double,1> weights;
    if (in.ReadVar("Weights", weights)) {
      double myWeight = weights(PathData.GetCloneNum());
      OutFile.WriteVar("Weight", myWeight);
    }
  }
  int numOfObservables=in.CountSections("Observable");
  Observables.resize(numOfObservables);
  
  for (int counter=0;counter<numOfObservables;counter++){
    in.OpenSection("Observable",counter);
    string observeType;
    assert(in.ReadVar("type",observeType));
    if (iAmRoot)
      OutFile.NewSection(observeType);
    ObservableClass* tempObs;
    switch (observeType) { 
      case "PairCorrelation": 
	tempObs = new PairCorrelationClass(PathData,OutFile);
      case "nofr":
	tempObs = new nofrClass(PathData,OutFile);
      case "AngularMomentumCorrelation":
	tempObs = new AngularMomCor(PathData,OutFile);
      case "DropletSuperfluidity":
	tempObs = new SuperfluiDrop(PathData,OutFile);
      case "Vacancy":
	tempObs = new VacancyLocClass(PathData,OutFile);
      case "Coupling":
	tempObs = new CouplingClass(PathData,OutFile);
      case "Energy":
	tempObs = new EnergyClass(PathData,OutFile);
      case "EnergySign":
	tempObs = new EnergySignClass(PathData,OutFile);
      case "ModifiedEnergy":
	tempObs = new ModifiedEnergyClass(PathData,OutFile);
      case "AutoCorr":
	tempObs = new AutoCorrClass(PathData,OutFile);
      case "DistanceToOpen":
	tempObs = new HeadLocClass(PathData,OutFile);
      case "VacancyLocation":
	tempObs = new VacancyLocClass(PathData,OutFile);
      case "TimeAnalysis":
	tempObs = new MCTimeClass(PathData,OutFile,Moves,Observables);
      case "Angular":
	tempObs =  new AngularClass(PathData,OutFile);
      case "PathDump":
	tempObs = new PathDumpClass(PathData,OutFile);
      case "WindingNumber":
	tempObs = new WindingNumberClass(PathData,OutFile);
      case "Vacancy":
	tempObs = new VacancyLocClass(PathData,OutFile);
      case "CycleCount":
	tempObs = new PermutationCountClass(PathData,OutFile);
      case "StructureFactor":
	tempObs = new StructureFactorClass(PathData,OutFile);
      case "Sign":
	tempObs = new WeightClass(PathData,OutFile);
      case "Forces":
	tempObs = new ForcesClass(PathData,OutFile);
	//  case "OpenOrientation":
	//   tempObs = new OpenOrientationClass(PathData,OutFile);
      default:
	perr << "We do not recognize the observable " << observeType << endl;
	abort();
    }
    tempObs->Read(in);
    Observables(counter)=tempObs;
    if (iAmRoot)
      OutFile.CloseSection();
    in.CloseSection();//Observable
  }
  if (iAmRoot)
    OutFile.CloseSection(); // "Observables"
}




void PIMCClass::ReadMoves(IOSectionClass &in)
{

  int numOfMoves=in.CountSections("Move");
  Moves.resize(numOfMoves);
  int steps;
  int myProc=PathData.Path.Communicator.MyProc();
  bool iAmRoot = (myProc == 0);
  
  if (iAmRoot)
    OutFile.NewSection("Moves");
  for (int counter=0;counter<numOfMoves;counter++){
    string moveName;
    in.OpenSection("Move",counter);
    string MoveType;
    assert(in.ReadVar("type",MoveType));
    if (iAmRoot)
      OutFile.NewSection(MoveType);
//     if (MoveType=="Bisection") {
//       moveName = "BisectionMove";
//       Moves(counter)=new BisectionMoveClass(PathData, OutFile);
//       Moves(counter)->Read(in);
//     }
//     else 
/*    if (MoveType=="OpenBisection") {
      moveName = "OpenBisectionMove";
      Moves(counter)=new OpenBisectionMoveClass(PathData, OutFile);
      Moves(counter)->Read(in);
    }
    else*/
    if (MoveType=="ShiftMove") {
      Moves(counter)=new ShiftMoveClass(PathData, OutFile);
      Moves(counter)->Read(in);
    }
    //    else if (MoveType=="EndSTage"){
    //      Moves(counter)=new OpenBisectionMoveClass(PathData,OutFile);
    //      Moves(counter)->Read(in);
    //    }
    else if (MoveType=="PrintMove") {
      Moves(counter)=new PrintMoveClass(PathData, OutFile);
      Moves(counter)->Read(in);
    }
//     else if (MoveType=="CycleBlock") {
//       moveName = "CycleBlockMove";
//       Moves(counter)=new CycleBlockMoveClass(PathData, OutFile);
//       Moves(counter)->Read(in);
//     }
//     else if (MoveType=="PermMove") {
//       moveName= "PermMove";
//       Moves(counter)=new PermMove(PathData, OutFile);
//       Moves(counter)->Read(in);
//     }
    else if (MoveType=="BisectionBlock"){
      moveName="BisectionBlock";
      Moves(counter)=new BisectionBlockClass(PathData,OutFile);
      Moves(counter)->Read(in); 
    }
    else if (MoveType=="CorrelatedBisectionBlock"){
      moveName="CorrelatedBisectionBlock";
      Moves(counter)=new CorrelatedBisectionBlockClass(PathData,OutFile);
      Moves(counter)->Read(in); 
    }
    else if (MoveType=="CouplingMove"){
      moveName="CouplingMove";
      Moves(counter)=new CouplingMoveClass(PathData,OutFile);
      Moves(counter)->Read(in);
    }
    else if (MoveType=="BisectionSphereBlock"){
      moveName="BisectionSphereBlock";
      Moves(counter)=new BisectionSphereBlockClass(PathData,OutFile);
      Moves(counter)->Read(in);
    }
    else if (MoveType=="CenterDroplet"){
      moveName="CenterDroplet";
      Moves(counter)=new CenterDropletClass(PathData,OutFile);
      Moves(counter)->Read(in);
    }
    //    else if (MoveType=="StructureReject"){
    //      moveName="StructureReject";
    //      Moves(counter)=new StructureRejectClass(PathData,OutFile);
    //      Moves(counter)->Read(in);
    //    }
    else if (MoveType=="OpenEnd"){
      moveName="OpenEnd";
      Moves(counter)=new OpenEndMoveClass(PathData,OutFile);
      Moves(counter)->Read(in);
    }
    else if (MoveType=="RefSlice"){
      moveName="RefSlice";
      Moves(counter)=new RefSliceMoveClass(PathData,OutFile);
      Moves(counter)->Read(in); 
    }
    else if (MoveType=="Displace"){
      moveName="Displace";
      Moves(counter)=new DisplaceMoveClass(PathData,OutFile);
      Moves(counter)->Read(in); 
    }
    else if (MoveType=="WaterRotate") {
      Moves(counter)=new WaterRotate(PathData, OutFile);
      Moves(counter)->Read(in);
    }
    else if (MoveType=="WaterTranslate") {
      Moves(counter)=new WaterTranslate(PathData, OutFile);
      Moves(counter)->Read(in);
    }
    else if (MoveType=="WaterTranslateRing") {
      Moves(counter)=new WaterTranslateRing(PathData, OutFile);
      Moves(counter)->Read(in);
    }
    else if (MoveType=="LocalFlip") {
      Moves(counter)=new LocalFlip(PathData, OutFile);
      Moves(counter)->Read(in);
    }
    else if (MoveType=="GlobalFlip") {
      Moves(counter)=new GlobalFlip(PathData, OutFile);
      Moves(counter)->Read(in);
    }
    else if (MoveType=="Langevin") {
      Moves(counter) = new LangevinMoveClass(PathData, OutFile);
      Moves(counter)->Read(in);
    }
    else {
      perr<<"This type of move is not recognized: "<< MoveType <<endl;
      abort();
    }
    if (iAmRoot)
      OutFile.CloseSection();
    in.CloseSection();
  }
  if (iAmRoot)
    OutFile.CloseSection (); // "Moves"
  
}



void PIMCClass::ReadAlgorithm(IOSectionClass &in)
{
  int maxWallTime;
  if (in.ReadVar("MaxWallTime", maxWallTime)) {
    PathData.SetMaxWallTime(maxWallTime);
    int hours = maxWallTime/3600;
    int minutes = (maxWallTime-3600*hours)/60;
    int seconds = maxWallTime%60;
    perr << "Maximum wall time is " << hours 
	 << ((hours != 1) ? " hours, " : " hour, ") << minutes
	 << ((minutes != 1) ? " minutes, and " : " minute, and ") << seconds 
	 << ((seconds != 1) ? " seconds.\n" : " second.\n");
  }
  Algorithm.Read(in,1);
}


void PIMCClass::Run()
{
  Algorithm.DoEvent();
  //  Array<MoveClass*,1> Moves;
  for (int counter=0;counter<Moves.size();counter++){
    cout<<"My name is "<<((MoveClass*)Moves(counter))->Name<<endl;
    cout<<"My acceptance ratio is "<<((MoveClass*)Moves(counter))->AcceptanceRatio()<<endl;
  }
  
}

void PIMCClass::WriteSystemInfo()
{
  dVec box = PathData.Path.GetBox();
  Array<double,1> boxArray(3);
  boxArray(0) = box[0];   boxArray(1) = box[1];   boxArray(2) = box[2];
  OutFile.WriteVar ("Box", boxArray);
  OutFile.WriteVar("tau",PathData.Path.tau);
  OutFile.WriteVar("NumTimeSlices",PathData.Path.TotalNumSlices);
  OutFile.WriteVar("seed",PathData.Seed);
  for (int speciesIndex=0; speciesIndex < PathData.Path.NumSpecies(); 
       speciesIndex++) {
    SpeciesClass &species = PathData.Path.Species(speciesIndex);
    OutFile.NewSection("Species");
    OutFile.WriteVar ("Name", species.Name);
    OutFile.WriteVar ("NumParticles", species.NumParticles);
    OutFile.WriteVar ("lambda", species.lambda);
    ParticleType type = species.GetParticleType();
    if (type == FERMION)
      OutFile.WriteVar ("ParticleType", "Fermion");
    if (type == BOSON)
      OutFile.WriteVar ("ParticleType", "Boson");
    if (type == BOLTZMANNON)
      OutFile.WriteVar ("ParticleType", "Boltzmannon");
    OutFile.CloseSection(); //"Species"
  }
}
