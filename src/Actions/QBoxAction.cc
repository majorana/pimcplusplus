#include "../PathDataClass.h"
#include "QBoxAction.h"
#include <sstream>

int qcount = 0;

std::string
QBoxActionClass::GetName()
{
  return "QBox_DFT_Action";
}
	
QBoxActionClass::QBoxActionClass(PathDataClass &pathData) : ActionBaseClass (pathData)
{
  qb = new qbox::qbLib(PathData.WorldComm.MPIComm);
  isAction = false;
}
	
double QBoxActionClass::SingleAction(int slice1,int slice2,const Array<int,1> &activeParticles,int level){
  if(slice1 == 0 && slice2 == 0){
    slice1 -= 1;
    slice2 += 1;
  }
  else if (slice1 == 0 && slice2 == PathData.Path.TotalNumSlices-1) {
    slice1 -= 1;
    slice2 += 1;
  }
  slice1 += 1;
  slice2 -= 1;
	int myAge = PathData.moveClock;
	double Utotal = 0.0;
  isAction = true;
	Utotal = ComputeEnergy(slice1, slice2, activeParticles, level);
	newEnergy = Utotal;
	return Utotal*PathData.Path.tau;
}

double QBoxActionClass::ComputeEnergy(int slice1,int slice2,const Array<int,1> &activeParticles,int level){
	//int slice = 0;
  double Utotal = 0.0;
    for(int slice=slice1; slice<=slice2; slice++) {
		  SetPtclPos(slice, activeParticles);
      ostringstream runCmd;
      runCmd << "run 0 " << steps;
      qcount ++;
      qb->issueCommand(runCmd.str());
      double getEnergy = qb->getEnergy();
      cerr << qcount <<" " << PathData.WorldComm.MyProc() << " " << getEnergy << endl;
		  Utotal += getEnergy;
  }
  //cerr << qcount <<" " << PathData.WorldComm.MyProc() << " " << Utotal << endl;
  return Utotal;
}


double QBoxActionClass::d_dBeta (int slice1, int slice2, int level){
	//int slice = 0;
  double Utotal = 0.0;
    for(int slice=slice1; slice<=slice2; slice++) {
      isAction = false;
	    Array<int,1> activeParticles;
	    activeParticles.resize(PathData.Path.NumParticles());
	    for(int a=0; a<activeParticles.size(); a++){
	    	activeParticles(a) = a;
	    }
	    SetPtclPos(slice, activeParticles);
	    //toqbox << "run 1 " << steps << endl;
      ostringstream runCmd;
      runCmd << "run 0 " << steps;
      qcount ++;
      qb->issueCommand(runCmd.str());
      double getEnergy = qb->getEnergy();
	    Utotal += getEnergy;
    }
	return Utotal;
}

void QBoxActionClass::SetPtclPos(int slice, Array<int,1> activeParticles){
	for(int a=0; a<activeParticles.size(); a++){
		int ptcl = activeParticles(a);
		dVec R = PathData.Path(slice,ptcl);
    PathData.Path.PutInBox(R);
    ostringstream runCmd;
		//toqbox << "move " << ptclID(ptcl) << " to " << R(0) << " " << R(1) << " " << R(2) << endl;
		runCmd << "move " << ptclID(ptcl) << " to " << R(0) << " " << R(1) << " " << R(2);
		cerr << PathData.WorldComm.MyProc() << " moving ptcl " << ptcl << "(" << ptclID(ptcl) << ") to " << R << endl;
	  //myWrite(towrite, runCmd.str());
    qb->issueCommand(runCmd.str());
	}
}
		
void QBoxActionClass::Read (IOSectionClass &in){
	  cerr.precision(16);
	  cerr << "QBoxAction Read" << endl;
	  int numPtcls;
	  assert(in.ReadVar("NumParticles",numPtcls));
	  assert(in.ReadVar("ParticleIDs",ptclID));
	  assert(numPtcls == ptclID.size());
	  assert(in.ReadVar("Steps", steps));
	  feedback = false;
	  //if(in.ReadVar("AutoConverge",feedback)){
	  //	if(feedback){
	  //		assert(in.ReadVar("Tolerance",tol));
	  //		assert(in.ReadVar("Strength",w));
	  //	}
	  //}
    
	  //string filename, fromqboxString, toqboxString;
	  //assert(in.ReadVar("InitFile", filename));
	  //assert(in.ReadVar("OutputStream", toqboxString));
	  //assert(in.ReadVar("InputStream", fromqboxString));
	  //cerr << "Opening pipes " << fromqboxString << " and " << toqboxString << "...";
	  //toqbox.open(toqboxString.c_str());
	  //fromqbox.open(fromqboxString.c_str());

    //assert(in.ReadVar("Command", command));
    string filename;
	  assert(in.ReadVar("InitFile", filename));
    cerr << "  Launching Qbox...";
    qb->launch(filename.c_str());
    //popen2(command.c_str(), &towrite, &toread);
    cerr << " SUCCESSFUL, ready to run" << endl;
    //cerr << " trying init run" << endl;
    //ostringstream runCmd;
    //runCmd << "run 0 10";
    //qb->issueCommand(runCmd.str());
    cerr << "QBoxAction Read finished." << endl;
}
  
void QBoxActionClass::AcceptCopy (int slice1, int slice2){
	age = PathData.moveClock;
	prevEnergy = newEnergy;
	//cerr << "QBOX accept: caching " << prevEnergy << " = " << newEnergy << endl;
}

void QBoxActionClass::RejectCopy (int slice1, int slice2){
	age = PathData.moveClock;
	//prevEnergy = oldEnergy;
	//cerr << "QBOX reject: caching " << prevEnergy << endl;
}
