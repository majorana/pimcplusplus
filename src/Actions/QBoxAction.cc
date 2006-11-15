#include "QBoxAction.h"
#include "../PathDataClass.h"
#include <sstream>

std::string
QBoxActionClass::GetName()
{
  return "QBox_DFT_Action";
}
	
QBoxActionClass::QBoxActionClass(PathDataClass &pathData) : ActionBaseClass (pathData){
  cerr << "QBoxActionClass Constructor: " << endl;
}
	
double QBoxActionClass::SingleAction(int slice1,int slice2,const Array<int,1> &activeParticles,int level){
	int myAge = PathData.moveClock;
	double Utotal = 0.0;
	if(GetMode()==OLDMODE){
		if((myAge-1)==age){
			//Utotal = prevEnergy;
			Utotal = ComputeEnergy(slice1, slice2, activeParticles, level);
			cerr << "Cached energy " << prevEnergy << " from move " << age << " at move " << myAge << " compares with computed energy " << Utotal << endl;
		}
		else{
			Utotal = ComputeEnergy(slice1, slice2, activeParticles, level);
			//oldEnergy = Utotal;
		}
	}
	else {
		Utotal = ComputeEnergy(slice1, slice2, activeParticles, level);
		newEnergy = Utotal;
	}
	return Utotal*PathData.Path.tau;
}

double QBoxActionClass::ComputeEnergy(int slice1,int slice2,const Array<int,1> &activeParticles,int level){
		int slice = 0;
		double Utotal = 0.0;
		SetPtclPos(slice, activeParticles);
		toqbox << "run 1 " << steps << endl;
		Utotal = Collect();
		cerr << "Qbox action received energy " << Utotal << endl;
	return Utotal;
}


double QBoxActionClass::d_dBeta (int slice1, int slice2, int level){
	int slice = 0;
	Array<int,1> activeParticles;
	activeParticles.resize(PathData.Path.NumParticles());
	for(int a=0; a<activeParticles.size(); a++){
		activeParticles(a) = a;
	}
	SetPtclPos(slice, activeParticles);
	toqbox << "run 1 " << steps << endl;
	double Utotal = Collect();
	return Utotal;
}

void QBoxActionClass::SetPtclPos(int slice, Array<int,1> activeParticles){
	for(int a=0; a<activeParticles.size(); a++){
		int ptcl = activeParticles(a);
		dVec R = PathData.Path(slice,ptcl);
		//cerr << "moving ptcl " << ptcl << "(" << ptclID(ptcl) << ") to " << R << endl;
		toqbox << "move " << ptclID(ptcl) << " to " << R(0) << " " << R(1) << " " << R(2) << endl;
	}
}

double QBoxActionClass::Collect(){
	double e0, e1, value;
	e0 = 0.0;
	string tag;
	while(fromqbox){
		fromqbox >> tag;
		if(tag == "run"){
			cerr << "QBox running for " << steps << " iterations";
		}
		else if(tag == "<etotal_int>"){
			e1 = e0;
			fromqbox >> tag;
			e0 = atof(tag.c_str());
			cerr << ".";
		}
		else if(tag == "<etotal>"){
			fromqbox >> tag;
			value = atof(tag.c_str());
			if(feedback){
				double diff = abs(e0 - e1);
				cerr << " converged to " << diff << endl;
				if(diff > tol){
					steps = int(steps*(1+w))+1;
				}else{
					steps = int(steps/(1+w));
					if(steps<4)
						steps = 4;
				}
			}
			cerr << "QBOX DFT RETURNING ETOTAL " << value << endl;
			return value;
		}
	}
	return 0.0;
}
			
		
void QBoxActionClass::Read (IOSectionClass &in){
	cerr.precision(8);
	cerr << "QBoxAction Read" << endl;
	int numPtcls;
	assert(in.ReadVar("NumParticles",numPtcls));
	assert(in.ReadVar("ParticleIDs",ptclID));
	assert(numPtcls == ptclID.size());
	assert(in.ReadVar("Steps", steps));
	feedback = false;
	if(in.ReadVar("AutoConverge",feedback)){
		if(feedback){
			assert(in.ReadVar("Tolerance",tol));
			assert(in.ReadVar("Strength",w));
		}
	}
	string filename, fromqboxString, toqboxString;
	assert(in.ReadVar("InitFile", filename));
	assert(in.ReadVar("OutputStream", fromqboxString));
	assert(in.ReadVar("InputStream", toqboxString));
	cerr << "Opening pipes " << fromqboxString << " and " << toqboxString << endl;
	toqbox.open(toqboxString.c_str());
	fromqbox.open(fromqboxString.c_str());
	toqbox << filename << endl;
	prevEnergy = Collect();
	age = PathData.moveClock;
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
