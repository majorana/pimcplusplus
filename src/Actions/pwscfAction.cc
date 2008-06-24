#include "../PathDataClass.h"
#include "pwscfAction.h"
#include <sstream>

std::string
pwscfActionClass::GetName()
{
  return "pwscf_DFT_Action";
}
	
pwscfActionClass::pwscfActionClass(PathDataClass &pathData) : ActionBaseClass (pathData){
  ////  cerr << "pwscfActionClass Constructor: " << endl;
  //out.open("convergence.dat");
  //out.precision(16);
  //out << "## Rx Ry Rz oldS(10) oldS(20) oldS(25) oldS(30) oldS(full) oldConverged Rx Ry Rz newS(10) newS(20) newS(25) newS(30) newS(full) newConverged" << endl;
  //out << "## Rxnew Rynew Rznew newS(5) newS(10) newS(15) newS(20) newS(full) newConverged" << endl;
  qCount = 0;
}
	
double pwscfActionClass::SingleAction(int slice1,int slice2,const Array<int,1> &activeParticles,int level){
  //cerr << "SingleAction " << slice1 << " " << slice2 << endl;
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
  qCount++;
	double Utotal = 0.0;
  Utotal = ComputeEnergy(slice1, slice2, activeParticles, level);
	return Utotal*PathData.Path.tau;
}

double pwscfActionClass::ComputeEnergy(int slice1,int slice2,const Array<int,1> &activeParticles,int level){
  //int slice = 0;
  double Utotal = 0.0;
  //cerr << "pwscf Compute Energy " << slice1 << " " << slice2 << endl;
  for(int slice=slice1; slice<=slice2; slice++) {
    ostringstream coordname;
    //coordname << "coords." << qCount << "." << slice;
    coordname << "coords.tmp.dat";
    //if(GetMode() == NEWMODE)
    //  coordname << ".new.dat";
    //else
    //  coordname << ".old.dat";
    SetPtclPos(slice, coordname.str().c_str());
    ostringstream command;
    command << "./exec_pwscf.py " << coordname.str() << " energy.out";
    cerr << "Issuing command " << command.str() << endl;
    int exitcode = system(command.str().c_str());
    if(exitcode != 0) {
      cerr << "Execution error " << exitcode << endl;
      exit(1);
    }
    std::ifstream in("energy.out");
    std::string line;
    in >> line;
    double e;
    e = atof(line.c_str());
    cerr << std::setprecision(12) << "Got energy " << e << std::endl;
    in.close();
    Utotal += e;
  }
  //cerr << "pwscf total energy " << Utotal << endl;
  return Utotal;
}


double pwscfActionClass::d_dBeta (int slice1, int slice2, int level){
	//int slice = 0;
  cerr << "pwscf dBeta " << slice1 << " " << slice2 << endl;
  double Utotal = 0.0;
  for(int slice=slice1; slice<=slice2; slice++) {
	  SetPtclPos(slice, "dBeta_coords.dat");
    ostringstream command;
    command << "./exec_pwscf.py dBeta_coords.dat energy.out";
    cerr << "Issuing command " << command.str() << endl;
    int exitcode = system(command.str().c_str());
    if(exitcode != 0) {
      cerr << "Execution error " << exitcode << endl;
      exit(1);
    }
    std::ifstream in("energy.out");
    std::string line;
    in >> line;
    double e;
    e = atof(line.c_str());
    cerr << std::setprecision(12) << "Got energy " << e << std::endl;
    in.close();
    Utotal += e;
  }
	return Utotal;
}

void pwscfActionClass::SetPtclPos(int slice, const char* filename){
  dVec box = PathData.Path.GetBox();
  ofstream out(filename);
	for(int a=0; a<PathData.Path.NumParticles(); a++){
		int ptcl = a;
    int speciesNum = PathData.Path.ParticleSpeciesNum(ptcl);
		string spec = PathData.Path.Species(speciesNum).Name;
    out << spec;
		dVec R = PathData.Path(slice,ptcl);
    PathData.Path.PutInBox(R);
    for(int x=0; x<NDIM; x++) {
      R(x) = (R(x) + box(x)/2)/box(x);
      out << " " << R(x);
    }
    out << endl;
	}
  out.close();
}

void pwscfActionClass::Read (IOSectionClass &in){
	cerr.precision(16);
	cerr << "pwscfAction Read" << endl;
  cerr << "QBoxAction Read finished." << endl;
}
