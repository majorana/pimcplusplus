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
	  out << "CEIMCActionClass Constructor: " << endl;
		out.open("QMCManager.out");
	}
	
	double CEIMCActionClass::SingleAction(int slice1,int slice2,const Array<int,1> &activeParticles,int level){
	  double Utotal = 0.0;
	  double U;
		int slice = 0;
		bool newmode = false;
		if(GetMode() == NEWMODE)
			newmode = true;
		out << PathData.QMCComm.MyProc() << " broadcasting newmode " << newmode << "... ";
		PathData.QMCComm.Broadcast(0,newmode);
		out << "done" << endl;
		// not sure about looping over slices; not done right now
	  //for(int slice=slice1; slice<slice2; slice++) left bracket
		Array<string,1> setPtclSet(PathData.Path.NumSpecies());
		for(int s=0; s<setPtclSet.size(); s++) setPtclSet(s) = PathData.ptclSet0(s);
		if(newmode && correlated)
			for(int s=0; s<setPtclSet.size(); s++) setPtclSet(s) = PathData.ptclSet1(s);
		// set ion positions for either mode
		out << "SetPtclSet contains..." << endl;
		for(int s=0; s<setPtclSet.size(); s++){
			out << s << " " << setPtclSet(s) << endl;
		}
		out << "PtclSets contain..." << endl;
		for(int s=0; s<setPtclSet.size(); s++){
			out << s << " " << PathData.ptclSet0(s) << ", " << PathData.ptclSet1(s) << endl;
		}
		int ptclSize = activeParticles.size();
		out << "activeParticles.size is " << ptclSize << endl;
		out << PathData.QMCComm.MyProc() << " broadcasting ptclsize " << ptclSize << "... ";
		PathData.QMCComm.Broadcast(0,ptclSize);
		out << "done" << endl;
		Array<int,1> SpeciesList(ptclSize);
		Array<int,1> OffsetList(ptclSize);
		Array<Vec3,1> CoordList(ptclSize);
		for(int i=0; i<ptclSize; i++){
		  int ptcl = activeParticles(i);
			int mySpecies = PathData.Path.ParticleSpeciesNum(ptcl);
			int offset = PathData.Path.Species(mySpecies).FirstPtcl;
			SpeciesList(i) = mySpecies;
			OffsetList(i) = ptcl - offset;
			CoordList(i) = PathData.Path(slice,ptcl);
			out << i << "[" << ptcl << ", " << mySpecies << ", " << offset << "], ";
			out << endl << "going to update the position of ptcl " << ptcl << " of set " << setPtclSet(mySpecies) << endl;
		  PathData.qmc->SetPtclPos(setPtclSet(mySpecies), ptcl - offset, PathData.Path(slice,ptcl).data());
		}
		out << PathData.QMCComm.MyProc() << " broadcasting SpeciesList " << SpeciesList << "... ";
		PathData.QMCComm.Broadcast(0, SpeciesList);
		out << PathData.QMCComm.MyProc() << " broadcasting OffsetList " << OffsetList << "... ";
		PathData.QMCComm.Broadcast(0, OffsetList);
		out << PathData.QMCComm.MyProc() << " broadcasting CoordList " << CoordList << "... ";
		PathData.QMCComm.Broadcast(0, CoordList);

		if(correlated){
			if(newmode){
				out << "  QMCAction: Setting up VMCMultiple run...";
	  		PathData.qmc->SetVMCMultiple(dt, walkers, steps, blocks);
				out << " done." << endl;
	
	  		PathData.qmc->process();
	
				EnergyDiffIndex = PathData.qmc->qmcDriver->addObservable("DiffS0S1");
				EnergyIndex0 = PathData.qmc->qmcDriver->addObservable("LE0");
				EnergyIndex1 = PathData.qmc->qmcDriver->addObservable("LE1");
				out << "  QMCAction: calling execute...";
	  		PathData.qmc->execute();
				out << " done." << endl;
			}
		}
		else {
	  	PathData.qmc->SetVMC(dt, walkers, steps, blocks);
			PathData.qmc->process();
			EnergyIndex0 = PathData.qmc->qmcDriver->addObservable("LocalEnergy");
	  	PathData.qmc->execute();
		}

		vector<double> Uvalues(0);
	
		if(correlated){
			if(newmode){
				out << "going to collect energy differences" << endl;
				PathData.qmc->qmcDriver->Estimators->getData(EnergyDiffIndex,Uvalues); // get energy difference
				out << "got " << Uvalues.size() << " for myself" << endl;
				out << "Now collect from " << PathData.QMCComm.NumProcs() << " processes" << endl;
				for(int q=1; q<PathData.QMCComm.NumProcs(); q++){
					//vector<double>* DataBuff;
					out << "receiving from " << q << ": expecting " << blocks << " entries" << endl;
					//PathData.QMCComm.Receive (DataBuff, blocks, MPI_DOUBLE, q, 9);
					Array<double, 1> DataBuff(blocks);
					PathData.QMCComm.Receive (q,DataBuff);
					out << " complete.  Got " << DataBuff.size() << " entries: " << DataBuff << endl;
					//for(int f=0; f<DataBuff.size(); f++) out << DataBuff[f] << ", ";
					//out << "]]]" << endl;
					for(int e=0; e<DataBuff.size(); e++){
						Uvalues.push_back(DataBuff(e));
					}
					out << "now Uvalues has " << Uvalues.size() << " entries" << endl;
				}
				out << "completed loop over receives" << endl;
				double deltaE, variance;
				QuickAvg(&Uvalues,deltaE,variance);
				out << "	computed average " << deltaE << " with variance " << variance << " from " << Uvalues.size() << " entries." << endl;
  			Utotal += (deltaE + PathData.Path.tau*variance/2); // penalty included
  			//Utotal += deltaE; // no penalty; biased estimator
			}
		}
		else{
			PathData.qmc->qmcDriver->Estimators->getData(EnergyIndex0,Uvalues); // get energy
			for(int q=1; q<PathData.QMCComm.NumProcs(); q++){
				Array<double, 1> DataBuff;
				PathData.QMCComm.Receive (q, DataBuff);
				for(int e=0; e<DataBuff.size(); e++){
					Uvalues.push_back(DataBuff(e));
				}
			}
			double E, variance;
			QuickAvg(&Uvalues,E,variance);
			out << "	computed average " << E << " with variance " << variance << " from " << Uvalues.size() << " entries." << endl;
	  	Utotal += E; // no penalty
		}
	  // end slice for loop
		out << "QMCAction returning..." << endl;
		return Utotal*PathData.Path.tau;
	}

	/*
	// I think I can scrap this, but just in case...
	double CEIMCActionClass::SingleAction(int slice1,int slice2,const Array<int,1> &activeParticles,int level){
	  double Utotal = 0.0;
	  double U;
		int slice = 0;
		// not sure about looping over slices; not done right now
	  //for(int slice=slice1; slice<slice2; slice++) left bracket
		Array<string,1> setPtclSet(PathData.Path.NumSpecies());
		for(int s=0; s<setPtclSet.size(); s++) setPtclSet(s) = ptclSet0(s);
		if((GetMode() == NEWMODE) && correlated)
			for(int s=0; s<setPtclSet.size(); s++) setPtclSet(s) = ptclSet1(s);
		// set ion positions for either mode
		out << "SetPtclSet contains..." << endl;
		for(int s=0; s<setPtclSet.size(); s++){
			out << s << " " << setPtclSet(s) << endl;
		}
		out << "PtclSets contain..." << endl;
		for(int s=0; s<setPtclSet.size(); s++){
			out << s << " " << ptclSet0(s) << ", " << ptclSet1(s) << endl;
		}
		out << "activeParticles.size is " << activeParticles.size() << endl;
		for(int i=0; i<activeParticles.size(); i++){
		  int ptcl = activeParticles(i);
			int mySpecies = PathData.Path.ParticleSpeciesNum(ptcl);
			int offset = PathData.Path.Species(mySpecies).FirstPtcl;
			out << i << "[" << ptcl << ", " << mySpecies << ", " << offset << "], ";
			out << endl << "going to update the position of ptcl " << ptcl << " of set " << setPtclSet(mySpecies) << endl;
		  PathData.qmc->SetPtclPos(setPtclSet(mySpecies), ptcl - offset, PathData.Path(slice,ptcl).data());
		}

		if(correlated){
			if(GetMode() == NEWMODE){
				//out << "  QMCAction: Setting up VMCMultiple run...";
	  		PathData.qmc->SetVMCMultiple(dt, walkers, steps, blocks);
				//out << " done." << endl;
	
	  		PathData.qmc->process();
	
				EnergyDiffIndex = PathData.qmc->qmcDriver->addObservable("DiffS0S1");
				EnergyIndex0 = PathData.qmc->qmcDriver->addObservable("LE0");
				EnergyIndex1 = PathData.qmc->qmcDriver->addObservable("LE1");
	  		PathData.qmc->execute();
			}
		}
		else {
	  	PathData.qmc->SetVMC(dt, walkers, steps, blocks);
			PathData.qmc->process();
			EnergyIndex0 = PathData.qmc->qmcDriver->addObservable("LocalEnergy");
	  	PathData.qmc->execute();
		}

		vector<double> Uvalues;
	
		if(correlated){
			if(GetMode() == NEWMODE){
				PathData.qmc->qmcDriver->Estimators->getData(EnergyDiffIndex,Uvalues); // get energy difference
				double deltaE, variance;
				QuickAvg(&Uvalues,deltaE,variance);
				out << "	computed average " << deltaE << " with variance " << variance << " from " << Uvalues.size() << " entries." << endl;
  			Utotal += (deltaE + PathData.Path.tau*variance/2); // penalty included
  			//Utotal += deltaE; // no penalty; biased estimator
			}
		}
		else{
			PathData.qmc->qmcDriver->Estimators->getData(EnergyIndex0,Uvalues); // get energy
			double E, variance;
			QuickAvg(&Uvalues,E,variance);
			out << "	computed average " << E << " with variance " << variance << " from " << Uvalues.size() << " entries." << endl;
	  	Utotal += E;
		}
	  // end slice for loop
		return Utotal*PathData.Path.tau;
	}
	*/
	
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
	  	PathData.qmc->SetVMC(dt, walkers, steps, blocks);
	
	  	PathData.qmc->process();
	
			EnergyIndex0 = PathData.qmc->qmcDriver->addObservable("LocalEnergy");
	  	PathData.qmc->execute();
	
			vector<double> Evalues;
			PathData.qmc->qmcDriver->Estimators->getData(EnergyIndex0,Evalues); // get energy difference
			double Uavg, variance;
			QuickAvg(&Evalues,Uavg,variance);
			out << "	ENERGY computed average " << Uavg << " with variance " << variance << " from " << Evalues.size() << " entries." << endl;
	    Utotal += Uavg;
	  }
	  return Utotal;
	}
	
	void CEIMCActionClass::Read (IOSectionClass &in){
		//out << "In CEIMCActionClass::Read.  Got indices" << endl;
		//EnergyIndex0 = PathData.qmc->qmcDriver->addObservable("LE0");
		//assert(EnergyIndex0 != -1);
		//EnergyIndex1 = PathData.qmc->qmcDriver->addObservable("LE1");
		//assert(EnergyIndex1 != -1);
		//EnergyDiffIndex = PathData.qmc->qmcDriver->addObservable("DiffS0S1");
		//assert(EnergyDiffIndex != -1);
		//out << "EnergyIndex0 at " << EnergyIndex0 << endl;
		//out << "EnergyIndex1 at " << EnergyIndex1 << endl;
		//out << "EnergyDiffIndex at " << EnergyDiffIndex << endl;

		dt = PathData.dt;
		walkers = PathData.walkers;
		steps = PathData.steps;
		blocks = PathData.blocks;
		correlated = PathData.correlated;
	}
#endif

QMCSamplingClass::QMCSamplingClass(PathDataClass &pathData) : 
  ActionBaseClass (pathData)
{
  manager = PathData.InterComm.MyProc();
  myProc = PathData.InterComm.MyProc();
  workerSize = PathData.InterComm.NumProcs() - 1;
  int index = 0;
  //out << "In constructor.  I'm " << myProc << " of " << workerSize+1 << endl;
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
  //out << "Calculating Ion-Ion interaction" << endl;
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
      //if(Path.DoPtcl(ptcl2) && (species2 != speciese) && Path.MolRef(ptcl1)!=Path.MolRef(ptcl2)){
      if(Path.DoPtcl(ptcl2)){
	//out << "  " << ptcl1 << ", " << ptcl2 << ": ";
				for (int slice=slice1;slice<slice2;slice++){
				  double r;
				  dVec Rvec;
				  PathData.Path.DistDisp(slice, ptcl1, ptcl2, r, Rvec);
				  Utotal += prefactor*
				    PathData.Species(species1).Charge*
				    PathData.Species(species2).Charge*(1.0/r - 1.0/cutoff);
				  //out << "r=" << r << "; Utotal = " << Utotal << endl;
				  //out << "  assembling prefactor " << prefactor << "* charge1 " << PathData.Species(species1).Charge << "* charge2 " << PathData.Species(species2).Charge << " *(1/r " << 1.0/r << " - 1.0/cutoff " << -1.0/cutoff << endl;
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
      //if(Path.DoPtcl(ptcl2) && (species2 != speciese) && Path.MolRef(ptcl1)!=Path.MolRef(ptcl2)){
      if(Path.DoPtcl(ptcl2)){
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
    out << "PIMC++ " << myProc << ": Broadcasting " << num << endl;
    PathData.InterComm.Broadcast(0,num);
    double total = 0.0;
    double result = 0.0;
    double random;
    if(myProc!=0){
    random = PathData.Random.Local()-0.5;
    PathData.InterComm.Send(0, &random);
    out << "  PIMC++ " << myProc << ": Sent " << random << "... " << endl;
    }
    else{
    for(int p = 1; p<=workerSize; p++){
    //  out << p << " ";
    PathData.InterComm.Receive(p, &result);
    out << "  PIMC++ " << myProc << ": Received " << result << "... " << endl;
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

