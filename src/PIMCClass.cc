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
  

  if (PathData.Path.ExistsCoupling){
    int myProc=PathData.InterComm.MyProc();
    PathData.Path.ExistsCoupling=(double)(myProc)/100;
  }
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
  
  for (int counter=0;counter<numOfObservables;counter++){
    in.OpenSection("Observable",counter);
    string observeType, observeName;
    assert(in.ReadVar("Type",observeType));
    assert(in.ReadVar("Name",observeName));
    if (iAmRoot)
      OutFile.NewSection(observeType);
    ObservableClass* tempObs;
    if (observeType=="PairCorrelation") 
	tempObs = new PairCorrelationClass(PathData,OutFile);
    else if (observeType=="nofr")
      tempObs=new nofrClass(PathData,OutFile);
    else if (observeType=="AngularMomentumCorrelation")
      tempObs= new AngularMomCor(PathData,OutFile);
    else if (observeType=="DropletSuperfluidity")
      tempObs = new SuperfluiDrop(PathData,OutFile);
    else if (observeType=="Vacancy")
      tempObs = new VacancyLocClass(PathData,OutFile);
    else if (observeType=="Coupling")
      tempObs = new CouplingClass(PathData,OutFile);
    else if (observeType=="Energy")
      tempObs = new EnergyClass(PathData,OutFile);
    else if (observeType=="EnergySign")
      tempObs = new EnergySignClass(PathData,OutFile);
    else if (observeType=="ModifiedEnergy")
      tempObs = new ModifiedEnergyClass(PathData,OutFile);
    else if (observeType=="AutoCorr")
      tempObs = new AutoCorrClass(PathData,OutFile);
    else if (observeType=="DistanceToOpen")
      tempObs = new HeadLocClass(PathData,OutFile);
    else if (observeType=="PhiK")
      tempObs = new PhiKClass(PathData,OutFile);
    else if (observeType=="JosephsonPathDump")
      tempObs = new JosephsonPathDumpClass(PathData,OutFile);
    else if (observeType=="VacancyLocation")
      tempObs = new VacancyLocClass(PathData,OutFile);
    else if (observeType=="TimeAnalysis")
      tempObs = new MCTimeClass(PathData,OutFile,Moves,Observables);
    else if (observeType=="Angular")
      tempObs =  new AngularClass(PathData,OutFile);
    else if (observeType=="PathDump")
      tempObs = new PathDumpClass(PathData,OutFile);
    else if (observeType=="WindingNumber")
      tempObs = new WindingNumberClass(PathData,OutFile);
    else if (observeType=="Vacancy")
      tempObs = new VacancyLocClass(PathData,OutFile);
    else if (observeType=="CycleCount")
      tempObs = new PermutationCountClass(PathData,OutFile);
    else if (observeType=="StructureFactor")
      tempObs = new StructureFactorClass(PathData,OutFile);
    else if (observeType=="Sign")
      tempObs = new WeightClass(PathData,OutFile);
    else if (observeType=="Forces")
      tempObs = new ForcesClass(PathData,OutFile);
    //  else if ( "OpenOrientation")
    //   tempObs = new OpenOrientationClass(PathData,OutFile);
    else {
      perr << "We do not recognize the observable " << observeType << endl;
      abort();
    }
    tempObs->Name = observeName;
    tempObs->Read(in);
    Observables.push_back(tempObs);
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
  int steps;
  int myProc=PathData.Path.Communicator.MyProc();
  bool iAmRoot = (myProc == 0);
  MoveClass* move;  
  if (iAmRoot)
    OutFile.NewSection("Moves");
  for (int counter=0;counter<numOfMoves;counter++){
    in.OpenSection("Move",counter);
    string moveType, moveName;
    assert(in.ReadVar("Type",moveType));
    assert(in.ReadVar("Name",moveName));
    if (iAmRoot)
      OutFile.NewSection(moveType);
    if (moveType=="ShiftMove")
      move = new ShiftMoveClass(PathData, OutFile);
    else if (moveType=="PrintMove")
      move = new PrintMoveClass(PathData, OutFile);
    else if (moveType=="BisectionBlock")
      move = new BisectionBlockClass(PathData,OutFile);
    else if (moveType=="CorrelatedBisectionBlock")
      move = new CorrelatedBisectionBlockClass(PathData,OutFile);
    else if (moveType=="CouplingMove")
      move = new CouplingMoveClass(PathData,OutFile);
    else if (moveType=="CenterOfMass")
      move = new CenterOfMassMOveClass(PathData,OutFile);
    else if (moveType=="BisectionSphereBlock")
      move = new BisectionSphereBlockClass(PathData,OutFile);
    else if (moveType=="CenterDroplet")
      move = new CenterDropletClass(PathData,OutFile);
    else if (moveType=="GrowWorm")
      move = new WormGrowMoveClass(PathData,OutFile);
    else if (moveType=="OpenEnd")
      move = new OpenEndMoveClass(PathData,OutFile);
    else if (moveType=="RefSlice")
      move = new RefSliceMoveClass(PathData,OutFile);
    else if (moveType=="Displace")
      move = new DisplaceMoveClass(PathData,OutFile);
    else if (moveType=="WaterRotate")
      move = new WaterRotate(PathData, OutFile);
    else if (moveType=="WaterTranslate")
      move = new WaterTranslate(PathData, OutFile);
    else if (moveType=="WaterTranslateRing")
      move = new WaterTranslateRing(PathData, OutFile);
    else if (MoveType=="LocalFlip")
      move =new LocalFlip(PathData,OutFile);
    else if (moveType=="GlobalJosephson")
      move = new GlobalJosephsonMove(PathData, OutFile);
    else if (moveType=="LocalFlip")
      move = new LocalFlip(PathData, OutFile);
    else if (moveType=="GlobalFlip")
      move = new GlobalFlip(PathData, OutFile);
    else if (moveType=="Langevin")
      move = new LangevinMoveClass(PathData, OutFile);
    else{
      perr<<"This type of move is not recognized: "<< moveType <<endl;
      abort();
    }
    move->Name = moveName;
    move->Read(in);
    Moves.push_back(move);
    if (iAmRoot)
      OutFile.CloseSection();
    in.CloseSection();
  }
  if (iAmRoot) {
    OutFile.CloseSection (); // "Moves"
    OutFile.FlushFile();
  }
  
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
//   for (int counter=0;counter<Moves.size();counter++){
//     cout<<"My name is "<<((MoveClass*)Moves(counter))->Name<<endl;
//     cout<<"My acceptance ratio is "<<((MoveClass*)Moves(counter))->AcceptanceRatio()<<endl;
//   }
  
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
