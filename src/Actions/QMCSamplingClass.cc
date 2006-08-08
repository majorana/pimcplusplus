#include "QMCSamplingClass.h"
#include "../PathDataClass.h"


#ifdef USE_QMC
	#include "Message/Communicate.h"
	#include "Utilities/OhmmsInfo.h"
	#include <QMCDrivers/QMCDriver.h>
	//#include <QMCApp/QMCInterface.h>

	std::string
	CEIMCActionClass::GetName()
	{
	  return "CEIMC";
	}
	
	CEIMCActionClass::CEIMCActionClass(PathDataClass &pathData) : ActionBaseClass (pathData){
	  //cerr << "CEIMCActionClass Constructor: " << endl;
	}
	
	double CEIMCActionClass::SingleAction(int slice1,int slice2,const Array<int,1> &activeParticles,int level){
	  double Utotal = 0.0;
	  double U;
		int slice = 0;
		// not sure about looping over slices; not done right now
	  //for(int slice=slice1; slice<slice2; slice++) left bracket
		string setPtclSet = ptclSet0;
		if((GetMode() == NEWMODE) && correlated)
			setPtclSet = ptclSet1;
		// set ion positions for either mode
		for(int i=0; i<activeParticles.size(); i++){
		  int ptcl = activeParticles(i);
		  PathData.qmc->SetPtclPos(setPtclSet, ptcl, PathData.Path(slice,ptcl).data());
		}

		if((GetMode() == NEWMODE) && correlated){
				//cerr << "  QMCAction: Setting up VMCMultiple run...";
				PathData.qmc->SetVMCMultiple(50);
				//cerr << " done." << endl;
	
	  		PathData.qmc->process();
	
				EnergyDiffIndex = PathData.qmc->qmcDriver->addObservable("DiffS0S1");
				EnergyIndex0 = PathData.qmc->qmcDriver->addObservable("LE0");
				EnergyIndex1 = PathData.qmc->qmcDriver->addObservable("LE1");
		}
		else {
			PathData.qmc->SetVMC(50);
			PathData.qmc->process();
			EnergyIndex0 = PathData.qmc->qmcDriver->addObservable("LocalEnergy");
		}

	  PathData.qmc->execute();
	
		vector<double> Uvalues;
		if((GetMode() == NEWMODE) && correlated){
			PathData.qmc->qmcDriver->Estimators->getData(EnergyDiffIndex,Uvalues); // get energy difference
			double deltaE, variance;
			QuickAvg(&Uvalues,deltaE,variance);
			cerr << "	computed average " << deltaE << " with variance " << variance << " from " << Uvalues.size() << " entries." << endl;
  		Utotal += (deltaE + PathData.Path.tau*variance/2); // penalty included
  		//Utotal += deltaE; // no penalty; biased estimator
		}
		else{
			PathData.qmc->qmcDriver->Estimators->getData(EnergyIndex0,Uvalues); // get energy
			double E, variance;
			QuickAvg(&Uvalues,E,variance);
			cerr << "	computed average " << E << " with variance " << variance << " from " << Uvalues.size() << " entries." << endl;
	  	Utotal += E;
		}
	  // end slice for loop
		return Utotal*PathData.Path.tau;
	}
	
	void QuickAvg(std::vector<double>* values, double& avg, double& var){
	  double total = 0.0;
		double totalSq = 0.0;
	  int T = values->size();
	  for(int t=0; t<T; t++){
			double v = (*values)[t];
	    total += v;
	    totalSq += v*v;
	  }
		avg = total/T;
		var = totalSq/T - avg*avg;
		var *= 1.0/T;
	}
	
	double CEIMCActionClass::d_dBeta (int slice1, int slice2, int level){
	  double Utotal = 0.0;
	  for(int slice=slice1; slice<slice2; slice++){
	  	PathData.qmc->SetVMC(100);
	
	  	PathData.qmc->process();
	
			EnergyIndex0 = PathData.qmc->qmcDriver->addObservable("LocalEnergy");
	  	PathData.qmc->execute();
	
			vector<double> Evalues;
			PathData.qmc->qmcDriver->Estimators->getData(EnergyIndex0,Evalues); // get energy difference
			double Uavg, variance;
			QuickAvg(&Evalues,Uavg,variance);
			cerr << "	ENERGY computed average " << Uavg << " with variance " << variance << " from " << Evalues.size() << " entries." << endl;
	    Utotal += Uavg;
	  }
	  return Utotal;
	}
	
	void CEIMCActionClass::Read (IOSectionClass &in){
		//cerr << "In CEIMCActionClass::Read.  Got indices" << endl;
		//EnergyIndex0 = PathData.qmc->qmcDriver->addObservable("LE0");
		//assert(EnergyIndex0 != -1);
		//EnergyIndex1 = PathData.qmc->qmcDriver->addObservable("LE1");
		//assert(EnergyIndex1 != -1);
		//EnergyDiffIndex = PathData.qmc->qmcDriver->addObservable("DiffS0S1");
		//assert(EnergyDiffIndex != -1);
		//cerr << "EnergyIndex0 at " << EnergyIndex0 << endl;
		//cerr << "EnergyIndex1 at " << EnergyIndex1 << endl;
		//cerr << "EnergyDiffIndex at " << EnergyDiffIndex << endl;
		assert(in.ReadVar("Correlated",correlated));
		if(correlated){
			cerr << "Using correlated samping to compute energy differences." << endl;
		}
		else{
			cerr << "NOT using correlated sampling for energy differences." << endl;
		}
		if(in.ReadVar("ParticleSet0", ptclSet0) && in.ReadVar("ParticleSet1", ptclSet1))
			cerr << "  Got qmcPACK particlesets " << ptclSet0 << " and " << ptclSet1 << " for correlated sampling." << endl;
		else{
			cerr << "  Didn't read in particleset names; using default ion particleset 'i'" << endl;
			ptclSet0 = ptclSet1 = "i";
		}
	}
#endif

QMCSamplingClass::QMCSamplingClass(PathDataClass &pathData) : 
  ActionBaseClass (pathData)
{
  manager = PathData.InterComm.MyProc();
  myProc = PathData.InterComm.MyProc();
  workerSize = PathData.InterComm.NumProcs() - 1;
  int index = 0;
  //cerr << "In constructor.  I'm " << myProc << " of " << workerSize+1 << endl;
}

IonIonActionClass::IonIonActionClass(PathDataClass &pathData): 
  ActionBaseClass(pathData)
{
  prefactor = 1.0;
  cutoff = 1e10;
}

double IonIonActionClass::SingleAction(int slice1,int slice2,
				       const Array<int,1> &activeParticles,
				       int level)
{
  //cerr << "Calculating Ion-Ion interaction" << endl;
  double Utotal = 0.0;
  
  for (int c=0; c<Path.DoPtcl.size(); c++) Path.DoPtcl(c)=true;
  int speciese=Path.SpeciesNum("e");
  
  for(int i=0; i<activeParticles.size(); i++){
    int ptcl1 = activeParticles(i);
    Path.DoPtcl(ptcl1) = false;
    
    for(int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++){
      int species1 = Path.ParticleSpeciesNum(ptcl1);
      int species2 = Path.ParticleSpeciesNum(ptcl2);
      // make sure we haven't already done it AND that it's not an electron
      if(Path.DoPtcl(ptcl2) && (species2 != speciese)){
	//cerr << "  " << ptcl1 << ", " << ptcl2 << ": ";
				for (int slice=slice1;slice<slice2;slice++){
				  double r;
				  dVec Rvec;
				  PathData.Path.DistDisp(slice, ptcl1, ptcl2, r, Rvec);
				  Utotal += prefactor*
				    PathData.Species(species1).Charge*
				    PathData.Species(species2).Charge*(1.0/r - 1.0/cutoff);
				  //cerr << "r=" << r << "; Utotal = " << Utotal << endl;
				  //cerr << "  assembling prefactor " << prefactor << "* charge1 " << PathData.Species(species1).Charge << "* charge2 " << PathData.Species(species2).Charge << " *(1/r " << 1.0/r << " - 1.0/cutoff " << -1.0/cutoff << endl;
				}
      }
    }
  }
  return Utotal*PathData.Path.tau;
}

double IonIonActionClass::d_dBeta (int slice1, int slice2, int level){
  double Utotal = 0.0;
  
  for (int c=0; c<Path.DoPtcl.size(); c++) Path.DoPtcl(c)=true;
  int speciese=Path.SpeciesNum("e");
  
  for(int ptcl1=0; ptcl1<PathData.Path.NumParticles(); ptcl1++){
    Path.DoPtcl(ptcl1) = false;
    
    for(int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++){
      int species1 = Path.ParticleSpeciesNum(ptcl1);
      int species2 = Path.ParticleSpeciesNum(ptcl2);
      // make sure we haven't already done it AND that it's not an electron
      if(Path.DoPtcl(ptcl2) && (species2 != speciese)){
	for (int slice=slice1;slice<slice2;slice++){
	  double r;
	  dVec Rvec;
	  PathData.Path.DistDisp(slice, ptcl1, ptcl2, r, Rvec);
	  Utotal += prefactor*PathData.Species(species1).Charge*PathData.Species(species2).Charge*(1.0/r);
	}
      }
    }
  }
  return Utotal;
}

void IonIonActionClass::Read (IOSectionClass &in){
  in.ReadVar ("SetCutoff", cutoff);
  in.ReadVar ("SetPrefactor", prefactor);
}

double QMCSamplingClass::SingleAction(int slice1,int slice2,const Array<int,1> &activeParticles,int level){
  /*
    int num;
    if(myProc == 0) num = 1;
    cerr << "PIMC++ " << myProc << ": Broadcasting " << num << endl;
    PathData.InterComm.Broadcast(0,num);
    double total = 0.0;
    double result = 0.0;
    double random;
    if(myProc!=0){
    random = PathData.Random.Local()-0.5;
    PathData.InterComm.Send(0, &random);
    cerr << "  PIMC++ " << myProc << ": Sent " << random << "... " << endl;
    }
    else{
    for(int p = 1; p<=workerSize; p++){
    //  cerr << p << " ";
    PathData.InterComm.Receive(p, &result);
    cerr << "  PIMC++ " << myProc << ": Received " << result << "... " << endl;
    total += result;
    }
    }*/
  return 0.0;
}

double QMCSamplingClass::d_dBeta (int slice1, int slice2, int level){
  
}

void QMCSamplingClass::Read (IOSectionClass &in){
  
}

std::string 
QMCSamplingClass::GetName()
{
  return "QMCSampling";
}

std::string
IonIonActionClass::GetName()
{
  return "IonIon";
}

