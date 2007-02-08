#include "QBoxAction.h"
#include "../PathDataClass.h"
#include <sstream>

bool Extract(string s, string find, string& data);

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
		//if((myAge-1)==age){
		//	Utotal = prevEnergy;
		//	cerr << "Cached energy " << prevEnergy << " from move " << age << " at move " << myAge << " compares with computed energy " << Utotal << endl;
		//}
		//else{
		Utotal = ComputeEnergy(slice1, slice2, activeParticles, level);
    //oldEnergy = Utotal;
		//}
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
		//toqbox << "run 1 " << steps << endl;
    ostringstream runCmd;
    runCmd << "run 1 " << steps;
		myWrite(towrite, runCmd.str());
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
	//toqbox << "run 1 " << steps << endl;
  ostringstream runCmd;
  runCmd << "run 1 " << steps;
	myWrite(towrite, runCmd.str());
	double Utotal = Collect();
	return Utotal;
}

void QBoxActionClass::SetPtclPos(int slice, Array<int,1> activeParticles){
	for(int a=0; a<activeParticles.size(); a++){
		int ptcl = activeParticles(a);
		dVec R = PathData.Path(slice,ptcl);
		//cerr << "moving ptcl " << ptcl << "(" << ptclID(ptcl) << ") to " << R << endl;
    ostringstream runCmd;
		//toqbox << "move " << ptclID(ptcl) << " to " << R(0) << " " << R(1) << " " << R(2) << endl;
		runCmd << "move " << ptclID(ptcl) << " to " << R(0) << " " << R(1) << " " << R(2);
	  myWrite(towrite, runCmd.str());
	}
}

double QBoxActionClass::Collect(){
  cerr << "  In QBoxAction::Collect for " << steps << " iterations: ";
	double e0, e1, value;
	e0 = 0.0;
	string tag, data;
  bool found = false;
	while(!found){
    tag = myRead(toread);
    //cerr << tag << endl;
		//if(tag == "run"){
		//	cerr << "QBox running for " << steps << " iterations";
		//}
    if(Extract(tag, "<etotal_int>", data)){
			e1 = e0;
			//tag = myRead(toread);
			e0 = atof(data.c_str());
			cerr << ".";
		}
    else if(Extract(tag, "<etotal>", data)){
      found = true;
      cerr << "GOT TOTAL ENERGY " << data << endl;
    }
  }
	value = atof(data.c_str());
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

/*
double QBoxActionClass::Collect(){
  cerr << "  In QBoxAction::Collect... ";
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
  cerr << "DEAD" << endl;
	return 0.0;
}

*/
			
		
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
	//string filename, fromqboxString, toqboxString;
	//assert(in.ReadVar("InitFile", filename));
	//assert(in.ReadVar("OutputStream", toqboxString));
	//assert(in.ReadVar("InputStream", fromqboxString));
	//cerr << "Opening pipes " << fromqboxString << " and " << toqboxString << "...";
	//toqbox.open(toqboxString.c_str());
	//fromqbox.open(fromqboxString.c_str());

  assert(in.ReadVar("Command", command));
  string filename;
	assert(in.ReadVar("InitFile", filename));
  cerr << "  Launching Qbox...";
  popen2(command.c_str(), &towrite, &toread);
  cerr << " SUCCESSFUL" << endl;

  cerr << "  Sending input file " << filename << "... ";
	//toqbox << filename << endl;
  int byteSent = myWrite(towrite, filename);
  cerr << "done: sent " << byteSent << " bytes" << endl;
	prevEnergy = Collect();
	age = PathData.moveClock;
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

bool Extract(string s, string find, string& data){
  //cerr << "XT looking for " << find << " in " << s << endl;
  int size = s.size();
  int l = find.size();
  int i = 0;
  string full = "";
  string tag = "";
  data = "";

  char c;
  while(i<size){
    c = s[i];
    while(tag != find && i<size){
      if(c != ' '){
        tag += c;
      }
      full += c;
      i++;
      c = s[i];
    }
    //cerr << "  XT at " << i << " of " << size << " assembled tag " << tag << endl;

    bool done = false;
    if(tag == find){
      //cerr << "  XT matched " << tag << " and " << find << endl;
      int a = 0;
      while(i<size && !done){
        if(c != ' '){
          data += c;
          a++;
        }
        else{
          if(a > 0)
            done = true;
        }
        full += c;
        i++;
        c = s[i];
      }
      //cerr << "  XT index at " << i << " of " << size << " and data is " << data << endl;
      if(i == size)
        done = true;
      //cerr << "  XT assembled data " << data << " and done is " << done << endl;
      if(done)
        return done;
    }

    //cerr << i << " " << c << endl;
    full += c;
    i++;
  }
  //cerr << "  XT leaving after searching full string " << full << endl;
  return false;
}
