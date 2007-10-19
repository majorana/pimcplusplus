////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "../PathDataClass.h"
#include "KineticRotorClass.h"
#include "../Moves/MoveUtils.h"

KineticRotorClass::KineticRotorClass(PathDataClass &pathData ) : 
  ActionBaseClass (pathData)
{
}

///This has to be called after pathdata knows how many
///particles it has
void KineticRotorClass::Read(IOSectionClass& in)
{
  // read in moments of inertia
  // these are the defaults for water in Hartree(?)
  Ia = 4381.1018999999997;
  Ib = 8739.8981000000003;
  Ic = 13121.0;
  in.ReadVar("Ia",Ia);
  in.ReadVar("Ib",Ib);
  in.ReadVar("Ic",Ic);
  A = 1.0/(2*Ia);//1.1411e-4;
  B = 1.0/(2*Ib);//5.7207e-5;
  C = 1.0/(2*Ic);//3.8105e-5;

  int setJmax;
  assert(in.ReadVar("SetJmax",setJmax));
  int setPrec = 100;
  assert(in.ReadVar("SetPrecision",setPrec));
  //rho = new RotorRhoClass(A, B, C, setJmax, setPrec); 
  //rho->SetTau(PathData.Path.tau);

  //cerr << "WARNING: FOR TESTING I AM NOT READING IN GRIDPOINTS: THIS WILL BE SLOW" << endl;
  ReadGridPoints(in);
  ReadEnergyGridPoints(in);
  cerr << "LEAVING KINETIC ROTOR READ" << endl;
}

void KineticRotorClass::ReadGridPoints(IOSectionClass& in)
{
  // read in num of grid points
  int Ntheta, Nphi, Nchi;
  assert(in.ReadVar("NumThetaPoints",Ntheta));
  assert(in.ReadVar("NumPhiPoints",Nphi));
  //assert(in.ReadVar("NumChiPoints",Nchi));
  //double last = M_PI * (1. - 1./(Ntheta+1));
  //ThetaGrid = new LinearGrid(0., last, Ntheta);
  //last = 2*M_PI * (1. - 1./(Nphi+1));
  //PhiGrid = new LinearGrid(0., last, Nphi);
  //last = 2*M_PI * (1. - 1./(Nchi+1));
  //ChiGrid = new LinearGrid(0., last, Nchi);

  // hack testing spline init
  ThetaGrid = new LinearGrid(0., M_PI, Ntheta);
  //PhiGrid = new LinearGrid(0., 2*M_PI, Nphi);
  PhiGrid = new LinearGrid(0., 4*M_PI, Nphi);
  //ChiGrid = new LinearGrid(0., 2*M_PI, Nchi);

  // read in spline gridpoints
  string filename;
  assert(in.ReadVar("File",filename));
  //Array<double,3> values;
  Array<double,2> values;
  values.resize(Ntheta, Nphi);
  double vol = Ntheta * Nphi;// * Nchi;
  int count = 0;
  ifstream infile(filename.c_str());
  cerr << "Reading in grid points from" << filename << endl;
  string line;
  cerr << "Getting Jmax: ";
  infile >> line;
  cerr << line << endl;
  cerr << "Getting tau: ";
  infile >> line;
  cerr << line << endl;
  //while(infile) {}
  int countT=0;
  for(int nt=0; nt<Ntheta; nt++) {
    for(int np=0; np<Nphi; np++) {
      //for(int nc=0; nc<Nchi; nc++) {
      double theta, phi, chi, logRho;
      infile >> line;
      theta = atof(line.c_str());
      infile >> line;
      phi = atof(line.c_str());
      //infile >> line;
      //chi = atof(line.c_str());
      infile >> line;
      double checkRho = atof(line.c_str());
      infile >> line;
      infile >> line;

      // checking readin of angles
      double checkT, checkP;//, checkC;
      checkT = (*ThetaGrid)(nt);
      checkP = (*PhiGrid)(np);
      //checkC = (*ChiGrid)(nc);
      //cerr << "Comparing angles at " << nt << " " << np << endl;//" " << nc << endl;
      //cerr << "theta " << theta << " " << checkT << endl;
      //cerr << "phi " << phi << " " << checkP << endl;
      //cerr << "chi " << chi << " " << checkC << endl;
      assert(abs(checkT - theta) < 1e-5);
      assert(abs(checkP - phi) < 1e-4);
      //assert(abs(checkC - chi) < 1e-5);

      logRho = atof(line.c_str());
      //values(nt, np, nc) = logRho;
      values(nt, np) = logRho;
      if(checkP == 0) {//) && checkC == 0) {
        countT++;
        cout << countT << " " << (*ThetaGrid)(nt) << " " << (*PhiGrid)(np) << " " << values(nt,np) << " " << exp(values(nt,np)) << endl;
      }
      count ++;
      //}
    }
  }

  cerr << "Read in " << count << " line; expected " << vol << endl;
  cerr << "Array extents " << values.extent(0) << " " << values.extent(1) << endl;
  cerr << "Grid sizes " << ThetaGrid->NumPoints << " " << PhiGrid->NumPoints << endl;
  spline.Init(ThetaGrid, PhiGrid, values);
  cerr << "Spline initialized" << endl;
}

void KineticRotorClass::ReadEnergyGridPoints(IOSectionClass& in)
{
  // read in num of grid points
  int Ntheta, Nphi, Nchi;
  assert(in.ReadVar("NumEnergyThetaPoints",Ntheta));
  assert(in.ReadVar("NumEnergyPhiPoints",Nphi));
  //assert(in.ReadVar("NumEnergyChiPoints",Nchi));
  ThetaEnergyGrid = new LinearGrid(0., M_PI, Ntheta);
  PhiEnergyGrid = new LinearGrid(0., 4*M_PI, Nphi);
  //ChiEnergyGrid = new LinearGrid(0., 2*M_PI, Nchi);

  // read in spline gridpoints
  string filename;
  assert(in.ReadVar("EnergyFile",filename));
  Array<double,2> values;
  values.resize(Ntheta, Nphi);
  double vol = Ntheta * Nphi;
  int count = 0;
  ifstream infile(filename.c_str());
  cerr << "Reading in Energy grid points from" << filename << endl;
  string line;
  cerr << "Getting Jmax: ";
  infile >> line;
  cerr << line << endl;
  cerr << "Getting tau: ";
  infile >> line;
  cerr << line << endl;
  for(int nt=0; nt<Ntheta; nt++) {
    for(int np=0; np<Nphi; np++) {
      //for(int nc=0; nc<Nchi; nc++) {
      double theta, phi, chi, E;
      infile >> line;
      theta = atof(line.c_str());
      infile >> line;
      phi = atof(line.c_str());
      infile >> line;
      chi = atof(line.c_str());
      infile >> line;

      // checking readin of angles
      double checkT, checkP, checkC;
      checkT = (*ThetaEnergyGrid)(nt);
      checkP = (*PhiEnergyGrid)(np);
      //checkC = (*ChiEnergyGrid)(nc);
      //cerr << "Comparing angles" << endl;
      //cerr << "theta " << theta << " " << checkT << endl;
      //cerr << "phi " << phi << " " << checkP << endl;
      //cerr << "chi " << chi << " " << checkC << endl;
      assert(abs(checkT - theta) < 1e-5);
      assert(abs(checkP - phi) < 1e-4);
      //assert(abs(checkC - chi) < 1e-5);

      E = atof(line.c_str());
      //values(nt, np, nc) = E;
      values(nt, np) = E;
      //cerr << (*ThetaEnergyGrid)(nt) << " " << (*PhiEnergyGrid)(np) << " " << (*ChiEnergyGrid)(nc) << values(nt,np,nc) << endl;
      count ++;
      //}
    }
  }

  cerr << "Read in " << count << " line; expected " << vol << endl;
  //EnergySpline.Init(ThetaEnergyGrid, PhiEnergyGrid, ChiEnergyGrid, values);
  EnergySpline.Init(ThetaEnergyGrid, PhiEnergyGrid, values);
  cerr << "Energy Spline initialized" << endl;
}

double 
KineticRotorClass::SingleAction (int slice1, int slice2,
			    const Array<int,1> &changedParticles, int level)
{
  //cerr << "KineticRotor testing..." << endl;
  //int numP = 1000;
  ////double dt = M_PI/numP;
  //double dt = 0.00314;
  ////double t0 = dt*numP;
  //double t0 = 0.;
  //for(int t=0; t<numP; t++) {
  //  double theta, phi, chi;
  //  theta = dt*t + t0;
  //  phi = 0.;
  //  chi = 0.;
  //  //for(int p=0; p<numP; p++) {
  //  //  phi = 2*(dt*p + t0);
  //  //  for(int c=0; c<numP; c++) {
  //  //    chi = 2*(dt*c + t0);
  //      double KE = -1 * EnergySpline(theta, phi, chi);
  //      double logR = spline(theta, phi, chi);
  //      cout << theta << " " << phi << " " << chi << " " << logR << " " << exp(logR) << " " << KE << endl;
  //  //  }
  //  //}
  //}
  //exit(1);
  //////////////////////////////////

  double TotalK = 0.0;
  // need to map activePtcls to active molecules - probably could be done more cleanly
  Array<bool,1> activeMol;
  activeMol.resize(PathData.Mol.NumMol());
  activeMol = false;
  for(int pindex=0; pindex<changedParticles.size(); pindex++) {
    int myMol = PathData.Mol(changedParticles(pindex));
    activeMol(myMol) = true;
  }
  vector<int> molRoster(0);
  for(int mindex=0; mindex<activeMol.size(); mindex++) {
    if(activeMol(mindex))
      molRoster.push_back(mindex);
  }

  int numChangedMol = molRoster.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int molIndex=0; molIndex<numChangedMol; molIndex++){
    int mol = molRoster[molIndex];
    //double FourLambdaTauInv=1.0/(4.0*Path.Species(species).lambda*levelTau);
    for (int slice=slice1; slice < slice2;slice+=skip) {
      double theta, phi, chi;
      getAngles(slice, mol, slice+skip, mol, theta, phi, chi);
      // explicit computation:
      //cout << "Calling CalcRho with angles " << phi << " " << theta << " " << chi << " between slices " << slice << " " << slice+skip << endl;
      //double r = rho->CalcRhoAllJExact(phi, theta, chi);
      //cerr << '\n';
      //// need to check for negative rho; use noisy correction
      //if(r<0) {
      //  cout << "Got NEGATIVE rho " << r;
      //  double a = -1 * PathData.Path.Random.Local() * r;
      //  cout << " generated noise " << a << endl;
      //  r = a;
      //}
      //double logRho = log(r);
      //cout << "Computed rho " << r << " and log " << logRho << endl;
      //cerr << "Accessed spline at indices " << theta << " " << phi << " " << chi << " value ";
      //double logRho = spline(theta, phi, chi);
      double logRho = spline(theta, phi + chi);
      //cerr << logRho << endl; 
      TotalK += logRho;
    }
  }
  //cout << "returning Action " << TotalK << endl << endl;
  return (TotalK);
}

double KineticRotorClass::d_dBeta (int slice1, int slice2,
			      int level)
{
  int M = PathData.Path.NumTimeSlices()-1;
  cerr << "Calculating energy over slices " << slice1 << " " << slice2 << endl;
  double TotalK = 0.0;
  double Z = 0.0;
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int molIndex=0; molIndex<PathData.Mol.NumMol(); molIndex++){
    int mol = molIndex;
    for (int slice=slice1; slice < slice2;slice+=skip) {
      double theta, phi, chi;
      getAngles(slice, mol, slice+skip, mol, theta, phi, chi);
      //double KE = rho->CalcEnergyExact(phi, theta, chi);
      // typo error: negatiave sign in tabulated energy, corrected here!!
      double KE = -1 * EnergySpline(theta, phi + chi);
      double z = exp(spline(theta, phi + chi));
      //cerr << slice << " " << mol << " " << phi << " " << theta << " " << chi << " " << KE << endl;
      TotalK += KE;
      Z += z;
    }
  }
  //cerr << TotalK << " " << Z << " " << TotalK/Z << endl;
  return (TotalK/(M*Z));
}


string
KineticRotorClass::GetName()
{
  return "KineticRotor";
}

const double TOL = 1e-6;

bool isEqual(dVec u, dVec v)
{
  if (abs(u(0) - v(0))< TOL) {
    if (abs(u(1) - v(1)) < TOL) {
      if (abs(u(2) - v(2)) < TOL) {
    return true;
      }
    }
  }
  else
    return false;
}

bool isEqualOpp(dVec u, dVec v)
{
  if (u(0) + v(0) < TOL) {
    if (u(1) + v(1) < TOL) {
      if (u(2) + v(2) < TOL) {
        return true;
      }
    }
  }
  else
    return false;
}

// map cartesian coordinates of two WATER molecules
// to an angular displacement expressed in terms of Euler angles
// Note this is specialized for a 3D rotor such as WATER!!
void KineticRotorClass::getAngles(int slice1, int mol1, int slice2, int mol2, double& theta, double& phi, double& chi)
{
  // get molecule coordinates WRT molecule center
  Array<int,1> mol1Members, mol2Members;
  PathData.Mol.MembersOf(mol1Members, mol1);
  PathData.Mol.MembersOf(mol2Members, mol2);
  dVec p1a = PathData.Path(slice1, mol1Members(1)) - PathData.Path(slice1, mol1Members(0));
  dVec p1b = PathData.Path(slice1, mol1Members(2)) - PathData.Path(slice1, mol1Members(0));
  dVec p2a = PathData.Path(slice2, mol2Members(1)) - PathData.Path(slice2, mol2Members(0));
  dVec p2b = PathData.Path(slice2, mol2Members(2)) - PathData.Path(slice2, mol2Members(0));

  // compute bisector, norm, and in-plane unit vectors
  dVec b1 = GetBisector(p1a, p1b);
  dVec n1 = Normalize(cross(p1a, p1b));
  dVec u1 = Normalize(cross(b1, n1));
  dVec b2 = GetBisector(p2a, p2b);
  dVec n2 = Normalize(cross(p2a, p2b));
  dVec u2 = Normalize(cross(b2, n2));

  // theta is angle between in-plane vectors u
  // (rotation about norm to molecule plane)
  theta = GetAngle(u2, u1);

  // axis of rotation for theta is cross of in-plane vectors u
  dVec A;
  if(isEqual(u1,u2))
    A = n1;
  else if(isEqualOpp(u1, u2))
    A = n1;
  else
    A = Normalize(cross(u1,u2));

  // phi is angle between norm n1 and A
  if(isEqual(n1,A))
    phi = 0.;
  else if(isEqualOpp(A, n1))
    phi = M_PI;
  else
    phi = GetAngle(n1, A);
  
  // chi is angle between norm A and n2
  if(isEqual(A, n2))
    chi = 0.;
  else if(isEqualOpp(A, n2))
    chi = M_PI;
  else
    chi = GetAngle(A, n2);

  //if(isnan(chi)) {
  //  cerr << "Crap chi=NAN A " << A << " " << n1 << " " << n2 << endl;
  //  cerr << "isEqual gives " << isEqual(A,n2) << " and negative " << isEqualOpp(A, n2) << " between " << A << " " << (-1*n2) << endl;
  //  cerr << "p1a " << p1a << " p1b " << p1b << " p2a " << p2a << " p2b " << p2b << endl;
  //}
  //if(isnan(phi)) {
  //  cerr << "Crap phi=NAN A " << A << " " << n1 << " " << n2 << endl;
  //  cerr << "isEqual gives " << isEqual(A,n1) << " and negative " << isEqualOpp(A, n1) << " between " << A << " " << (-1*n1) << endl;
  //  cerr << "p1a " << p1a << " p1b " << p1b << " p2a " << p2a << " p2b " << p2b << endl;
  //}
  if(phi<0){
    cerr << "NEGATIVE phi " << phi << endl;
    exit(1);
  }
  if(chi<0) {
    cerr << "NEGATIVE chi " << chi << endl;
    exit(1);
  }
}

FixedAxisRotorClass::FixedAxisRotorClass(PathDataClass &pathData ) : 
  ActionBaseClass (pathData), rotor(pathData)
{
  // hard-wired values
  Ia = 4381.1018999999997;
  Ib = 8739.8981000000003;
  Ic = 13121.0;
  C = 0.5*2.65146*0.0493089;
}

double FixedAxisRotorClass::SingleAction (int slice1, int slice2, 
		       const Array<int,1> &changedParticles, int level)
{
  double TotalK = 0.0;
  // need to map activePtcls to active molecules - probably could be done more cleanly
  Array<bool,1> activeMol;
  activeMol.resize(PathData.Mol.NumMol());
  activeMol = false;
  for(int pindex=0; pindex<changedParticles.size(); pindex++) {
    int myMol = PathData.Mol(changedParticles(pindex));
    activeMol(myMol) = true;
  }
  vector<int> molRoster(0);
  for(int mindex=0; mindex<activeMol.size(); mindex++) {
    if(activeMol(mindex))
      molRoster.push_back(mindex);
  }

  int numChangedMol = molRoster.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int molIndex=0; molIndex<numChangedMol; molIndex++){
    int mol = molRoster[molIndex];
    for (int slice=slice1; slice < slice2;slice+=skip) {
      double theta, phi, chi;
      rotor.getAngles(slice, mol, slice+skip, mol, theta, phi, chi);
      double r = CalcRho(phi, theta, chi);
      // need to check for negative rho; use noisy correction
      //if(r<0) {
      //  cout << "Got NEGATIVE rho " << r;
      //  double a = -1 * PathData.Path.Random.Local() * r;
      //  cout << " generated noise " << a << endl;
      //  r = a;
      //}

      double logRho = log(r);
      //cout << "Computed FIXED AXIS rho " << r << " and log " << logRho << endl;
      //double logRho = spline(theta, phi, chi);
      TotalK += logRho;
    }
  }
  //cout << "returning FIXED AXIS Action " << TotalK << endl << endl;
  return (TotalK);
}

double FixedAxisRotorClass::d_dBeta (int slice1, int slice2, int level)
{
  int M = PathData.Path.NumTimeSlices()-1;
  //cerr << "Calculating FIXED AXIS energy over slices " << slice1 << " " << slice2 << " with " << M << " total slices" << endl;
  double TotalK = 0.0;
  double Z = 0.0;
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int molIndex=0; molIndex<PathData.Mol.NumMol(); molIndex++){
    int mol = molIndex;
    for (int slice=slice1; slice < slice2;slice+=skip) {
      double theta, phi, chi;
      rotor.getAngles(slice, mol, slice+skip, mol, theta, phi, chi);
      double r = CalcRho(phi, theta, chi);
      TotalK += r * CalcFAEnergy(phi, theta, chi);
      Z += r;
    }
  }
  return (TotalK/(M*Z));

}

void FixedAxisRotorClass::Read (IOSectionClass &in)
{
  // read in moments of inertia
  // these are the defaults for water
  Ia = 4381.1018999999997;
  Ib = 8739.8981000000003;
  Ic = 13121.0;
  in.ReadVar("Ia",Ia);
  in.ReadVar("Ib",Ib);
  in.ReadVar("Ic",Ic);

  // arbitrary constant (?)
  //C = 0.5*2.65146*0.0493089;
  C = 1.0;
}

string FixedAxisRotorClass::GetName()
{
  return "FixedAxisRotor";
}

double FixedAxisRotorClass::CalcRho(double phi, double theta, double psi)
{
  double tau = PathData.Path.tau;
  double chi = cos(theta/2) * cos((psi + phi)/2);
  double eta = sin(theta/2) * cos((psi - phi)/2);
  double xi = sin(theta/2) * sin((psi - phi)/2);
  double zeta = cos(theta/2) * sin((psi + phi)/2);

  double gamma = 2 * acos(chi);
  double sinG2 = sin(gamma/2);
  double gammaOverSinG2 = gamma/sinG2;
  //if(phi==0. && psi==0. && theta==0.)
  if((phi==0. && psi==0.) || theta==0.)
    gammaOverSinG2 = 2.;

  double nDotI = eta*eta*Ib + xi*xi*Ic + zeta*zeta*Ia;

  double SFA;
  SFA = 0.5 * gammaOverSinG2 * gammaOverSinG2 * nDotI/tau;

  double D = sqrt((Ia * Ib * Ic)/(tau*tau*tau))*gammaOverSinG2;

  double Q = C * D * exp(-SFA);
  //cerr << "C " << C << " D " << D << " SFA " << SFA << endl;
  return Q;
}

double FixedAxisRotorClass::CalcFAEnergy(double phi, double theta, double psi)
{
  double tau = PathData.Path.tau;
  double chi = cos(theta/2) * cos((psi + phi)/2);
  double eta = sin(theta/2) * cos((psi - phi)/2);
  double xi = sin(theta/2) * sin((psi - phi)/2);
  double zeta = cos(theta/2) * sin((psi + phi)/2);

  double gamma = 2 * acos(chi);
  double sinG2 = sin(gamma/2);
  double gammaOverSinG2 = gamma/sinG2;
  //if(gammaOverSinG2 > 1e2) {
  //  cerr << phi << " " << theta << " " << psi << " " << gamma << " " << sinG2 << " " << gammaOverSinG2 << endl;
  //}
  if((phi==0. && psi==0.) || theta==0.)
    gammaOverSinG2 = 2.;

  //double nDotI = eta*eta*Ib + xi*xi*Ic + zeta*zeta*Ia;

  //double SFA;
  //SFA = 0.5 * gammaOverSinG2 * gammaOverSinG2 * nDotI/tau;

  //double D = sqrt((Ia * Ib * Ic)/(tau*tau*tau))*gammaOverSinG2;

  double Q = 1.0/(2 * tau * tau) * gammaOverSinG2 * gammaOverSinG2
    + 1.5/tau;
  return Q;
}
