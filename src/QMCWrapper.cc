#include "QMCWrapper.h"
#include <sstream>

void QMCWrapperClass::DummyInit(PathDataClass& PathData){
	ostringstream filename;
	filename << "QMCDummy." << PathData.MetaWorldComm.MyProc() << "." << PathData.QMCComm.MyProc() << ".out";
	out.open(filename.str().c_str());
	initialized = true;
}

void QMCWrapperClass::QMCDummy(PathDataClass& PathData){

#if USE_QMC
	if(!initialized)
		DummyInit(PathData);

	out << "proc" << PathData.QMCComm.MyProc() << "): In QMCDummy correlated is " << correlated << " but PathData.corr is " << PathData.correlated << endl;
	// SLICE SHOULD BE SENT, RIGHT??
	int slice = 0;

	bool newmode;
	out << "(proc" << PathData.QMCComm.MyProc() << "): Broadcast newmode received ";
	PathData.QMCComm.Broadcast(0, newmode);
	out << "(proc" << PathData.QMCComm.MyProc() << ") " << newmode << endl;
	// not sure about looping over slices; not done right now
  //for(int slice=slice1; slice<slice2; slice++) left bracket
	Array<string,1> setPtclSet(PathData.ptclSet0.size());
	out << "size of setptclset is " << setPtclSet.size() << " which should match " << PathData.ptclSet0.size() << endl;
	for(int s=0; s<setPtclSet.size(); s++) setPtclSet(s) = PathData.ptclSet0(s);
	if(newmode && correlated)
		for(int s=0; s<setPtclSet.size(); s++) setPtclSet(s) = PathData.ptclSet1(s);
	// set ion positions for either mode
	out << "proc" << PathData.QMCComm.MyProc() << "SetPtclSet contains..." << endl;
	for(int s=0; s<setPtclSet.size(); s++){
		out << "(proc" << PathData.QMCComm.MyProc() << ") " << s << " " << setPtclSet(s) << endl;
	}
	//out << "proc" << PathData.QMCComm.MyProc() << "PtclSets contain..." << endl;
	for(int s=0; s<setPtclSet.size(); s++){
		out << "(proc" << PathData.QMCComm.MyProc() << ") " << s << " " << PathData.ptclSet0(s) << ", " << PathData.ptclSet1(s) << endl;
	}
	int ptclSize;
	out << "(proc" << PathData.QMCComm.MyProc() << "): Broadcast ptclSize received ";
	PathData.QMCComm.Broadcast(0,ptclSize);
	out << "(proc" << PathData.QMCComm.MyProc() << ") " << ptclSize << endl;
	Array<int,1> SpeciesIndex(ptclSize);
	Array<int,1> OffsetList(ptclSize);
	Array<Vec3,1> CoordList(ptclSize);
	out << "(proc" << PathData.QMCComm.MyProc() << "): Broadcast SpeciesIndex received ";
	PathData.QMCComm.Broadcast(0, SpeciesIndex);
	out << "(proc" << PathData.QMCComm.MyProc() << ") " << SpeciesIndex << endl;
	out << "(proc" << PathData.QMCComm.MyProc() << "): Broadcast OffsetList received ";
	PathData.QMCComm.Broadcast(0, OffsetList);
	out << "(proc" << PathData.QMCComm.MyProc() << ") " << OffsetList << endl;
	out << "(proc" << PathData.QMCComm.MyProc() << "): Broadcast CooordList received ";
	PathData.QMCComm.Broadcast(0, CoordList);
	out << "(proc" << PathData.QMCComm.MyProc() << ") " << CoordList << endl;
	out << "ok now at for loop for " << ptclSize << " iterations..." << endl;
	//out << "I have setPtclSet " << setPtclSet << endl;
	//out << "I have OffsetList " << OffsetList << endl;
	//out << "I have CoordList " << CoordList << endl;
	//PathData.qmc->SetPtclPos(setPtclSet(SpeciesIndex(0)), OffsetList(0), CoordList(0).data());
	//out << "just did a single call to SetPtclPos." << endl;
	for(int i=0; i<ptclSize; i++){
		out << i << " of " << ptclSize << ": ";
		out << "going to update the position of ptcl " << OffsetList(i) << " of set " << setPtclSet(SpeciesIndex(i)) << endl;
	  PathData.qmc->SetPtclPos(setPtclSet(SpeciesIndex(i)), OffsetList(i), CoordList(i).data());
	}
	out << "on other side of for loop" << endl;

	if(correlated){
		if(newmode){
			out << "  QMCWrapper: Setting up VMCMultiple run...";
  		PathData.qmc->SetVMCMultiple(dt, walkers, steps, blocks);
			out << " done." << endl;

  		PathData.qmc->process();

			EnergyDiffIndex = PathData.qmc->qmcDriver->addObservable("DiffS0S1");
			EnergyIndex0 = PathData.qmc->qmcDriver->addObservable("LE0");
			EnergyIndex1 = PathData.qmc->qmcDriver->addObservable("LE1");
			out << "  QMCWrapper: calling execute...";
  		PathData.qmc->execute();
			out << "done" << endl;
		}
	}
	else {
		out << "NOT CORRELATED!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  	PathData.qmc->SetVMC(dt, walkers, steps, blocks);
		PathData.qmc->process();
		EnergyIndex0 = PathData.qmc->qmcDriver->addObservable("LocalEnergy");
  	PathData.qmc->execute();
	}

	vector<double> Uvalues(0);
	//vector<double> DataBuff(0);

	if(correlated){
		if(newmode){
			PathData.qmc->qmcDriver->Estimators->getData(EnergyDiffIndex,Uvalues); // get energy difference
			Array<double, 1> DataBuff(Uvalues.size());
			for(int u=0; u<DataBuff.size(); u++)
				DataBuff(u) = Uvalues[u];
			out << " Sending " << DataBuff.size() << " entries:" << endl;
			for(int f=0; f<DataBuff.size(); f++) out << DataBuff(f) << ", ";
			out << endl;
			assert(blocks == DataBuff.size());
			//PathData.QMCComm.Send (&DataBuff, DataBuff.size(), MPI_DOUBLE, 0, 1);
			PathData.QMCComm.Send (0, DataBuff);
			out << "Send complete" << endl;
		}
	}
	else{
		out << "NOT CORRELATED GETTING ENERGY!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		PathData.qmc->qmcDriver->Estimators->getData(EnergyIndex0,Uvalues); // get energy
		Array<double, 1> DataBuff(Uvalues.size());
		for(int u=0; u<DataBuff.size(); u++)
			DataBuff(u) = Uvalues[u];
		PathData.QMCComm.Send (0, DataBuff);
	}

	out << "(proc" << PathData.QMCComm.MyProc() << ") Leaving QMCWrapper" << endl;

#endif
}
