#include "RadialWF.h"
#include "../Integration/RungeKutta.h"

int RadialWF::TurningIndex()
{
  // Start in classically forbidden region and search inward toward
  // origin looking for classically allowed region
  Grid &grid = *u.grid;
  int index = grid.NumPoints-1;
  bool done = false;

  while (!done && (index >=0))
    {
      double r = grid(index);
      double A = pot->A(r);
      double B = pot->V(r);
      double V = pot->V(r);
      double dAdr = pot->dVdr(r);
      double E = Energy;
      double sl = (double) l;
      
      double Deriv2 = dAdr/r + sl*(sl+1.0)*B/(r*r) + 2.0*(V-E);
      
      if (Deriv2 < 0.0)
	done = true;
      else
	index--;
    }
  // If we didn't fine a turning point, select the middle grid point;
  if (index == -1)
    index = grid.ReverseMap(0.5*grid(grid.NumPoints-1));

  return (index);
}
void 
RadialWF::OriginBC(double r0, double &u0, double &du0)
{
  if (!pot->IsPH())
    {  // Use relativistic equations
      const double alpha = 1.0/137.036;
      double Z = -pot->V(r0)*r0;
      cerr << "Z = " << Z << endl;
      double a2Z2 = alpha*alpha*(double)Z*(double)Z;

      double sl = (double) l;
      double gamma = 
	(sl+1.0)*sqrt((sl+1.0)*(sl+1.0)-a2Z2) /
	(2.0*sl + 1.0); 
      if (l!=0)
	gamma += sl*sqrt(sl*sl-a2Z2)/(2.0*sl+1.0);
      
      u0 = pow(r0, gamma);
      du0 = gamma * pow(r0,gamma-1.0);
    }
  else
    {  // Use pseudoHamiltonian equations
      if (l == 0)
	{
	  double E = Energy;
	  double A = pot->A(r0);
	  double B = pot->B(r0);
	  double V = pot->V(r0);
	  if (V > Energy)
	    {
	      double kappa = sqrt(2.0 * (V-E)/A);
	      u0 = sinh(kappa * r0);
	      du0 = kappa * cosh (kappa * r0);
	    }
	  else
	    {
	      double k = sqrt(2.0 * (E-V)/A);
	      u0 = sin(k * r0);
	      du0 = k * cos(k * r0);
	    }
	}
      else {  // l != 0
	u0 = 1.0;
	du0 = (double)(l+1)/r0;
      }
    }  
  // Flip sign for odd number of nodes
  if (((n-l-1) % 2) == 1)
    {
      u0 *= -1.0;
      du0 *= -1.0;
    }

}


void
RadialWF::IntegrateOut()
{
  Grid &grid = *u.grid;
  Array<Vec2,1> Temp(grid.NumPoints);
  // Set up initial conditions:
  double r0 = grid(0); 
  OriginBC(r0, Temp(0)[0], Temp(0)[1]);

  // Now do integration
  if (pot->IsPH()) {
    PHDerivs derivs(*this);
    RungeKutta<PHDerivs,Vec2> integrator(derivs);
    integrator.Integrate(grid, 0, grid.NumPoints-1, Temp);
  }
  else {
    NonPHDerivs derivs(*this);
    RungeKutta<NonPHDerivs,Vec2> integrator(derivs);
    integrator.Integrate(grid, 0, grid.NumPoints-1, Temp);    
  }

  // Copy results of integration into RadialWF
  for (int i=0; i < grid.NumPoints; i++) {
    u(i) = Temp(i)[0];
    dudr(i) = Temp(i)[1];
  }
  Normalize();
}


double RadialWF::IntegrateInOut (int &tindex)
{
  Grid &grid = *u.grid;
  // Find classical turning point
  tindex = TurningIndex();
  // Compute starting value and derivative at origin
  double r0 = grid(0);
  OriginBC (r0, uduVec(0)[0], uduVec(0)[1]);

  // Do integration from the origin to the turning point
  if (pot->IsPH()) {
    PHDerivs derivs(*this);
    RungeKutta<PHDerivs,Vec2> integrator(derivs);
    integrator.Integrate(grid, 0, tindex, uduVec);
  }
  else {
    NonPHDerivs derivs(*this);
    RungeKutta<NonPHDerivs,Vec2> integrator(derivs);
    integrator.Integrate(grid, 0, tindex, uduVec);    
  } 

  for (int j=0; j<=tindex; j++) {
    u(j) = uduVec(j)[0];
    dudr(j) = uduVec(j)[1];
  }

  // Now initilize value and derivative at rmax
  int endPoint = grid.NumPoints-1;
  double rend = grid(endPoint);
  InfinityBC(rend, uduVec(endPoint)[0], uduVec(endPoint)[1]);
  
  // Do integration from the right to the turning point
  if (pot->IsPH()) {
    PHDerivs derivs(*this);
    RungeKutta<PHDerivs,Vec2> integrator(derivs);
    integrator.Integrate(grid, endPoint, tindex, uduVec);
  }
  else {
    NonPHDerivs derivs(*this);
    RungeKutta<NonPHDerivs,Vec2> integrator(derivs);
    integrator.Integrate(grid, endPoint, tindex, uduVec);    
  } 
  if ((u(tindex)*uduVec(tindex)[0]) < 0.0)
    for (int i=tindex; i<=endPoint; i++)
      uduVec(i) *= -1.0;

  double cuspValue = uduVec(tindex)[1]/uduVec(tindex)[0] -
    dudr(tindex)/u(tindex);

  // Copy values in, normalizing to make R continuous
  double factor = u(tindex) / uduVec(tindex)[0];
  for (int i=tindex+1; i<grid.NumPoints; i++) {
    u(i)    = factor*uduVec(i)[0];
    dudr(i) = factor*uduVec(i)[1];
  }
  return (cuspValue);
}



void 
RadialWF::InfinityBC(double rend, double &uend, double &duend)
{
  // Infinity should be outside core, so we don't need to worry about
  // A and B.
  double Vend = pot->V(rend);
  double dVdrend = pot->dVdr(rend);
  double E = Energy;
  double k = sqrt(l*(l+1.0)/(rend*rend) + (Vend-E));
  uend = 1.0;
  duend = -uend * (k + 0.5*rend/k*(dVdrend - 2.0*l*(l+1.0)/(rend*rend*rend)));
}



void RadialWF::Solve(double tolerance)
{
  int tindex;
  IntegrateInOut(tindex);
  int NumNodes = CountNodes();
  double Ehigh, Elow, Etrial;
  Ehigh = 0.0;
  Grid &grid = *u.grid;
  int TotalNodes = n-l-1;

  char fname[100];
  snprintf ("WFn%dl%d.h5", 100, n, l);
  IOSectionClass out;
  out.NewFile(fname);    
  out.WriteVar ("x", u.grid->Points());
  Array<double,2> tmp(u.grid->NumPoints,1);
  tmp(Range::all,0) = u.Data();
  out.WriteVar ("u", tmp);
  VarClass *varPtr = out.GetVarPtr("u");

  if (!pot->IsPH()){
    double N = n;
    double Z = -pot->V(grid(0))*grid(0);
    Elow = -1.5*Z*Z/(N*N);
  } 
  else {
    // Eigenenergy can't be lower than lowest potential -- otherwise
    // we'd have no turning point.
    Elow = pot->V(grid(0));
    for (int i=1; i<grid.NumPoints; i++) {
      double r = grid(i);
      if (pot->V(r) < Elow)
	Elow = pot->V(r);
    }
  }

  double Eold;
  Eold = Etrial = Energy;
  bool done = false;
  while (!done) {
    cerr << "Etrial = " << Etrial << endl;
    //cerr << "Ehigh = " << Ehigh << " Elow = " << Elow << endl;
    double CuspValue = IntegrateInOut(tindex);
    //cerr << "Cusp value = " << CuspValue << "\n";
    
    NumNodes = CountNodes();
    //cerr << "NumNodes = " << NumNodes << "\n";
    
    if (NumNodes > TotalNodes)
      Ehigh = Etrial;
    else if (NumNodes < TotalNodes)
      Elow = Etrial;
    else if (CuspValue < 0.0)
      Elow = Etrial;
    else if (CuspValue > 0.0)
      Ehigh = Etrial;
    else if (isnan(CuspValue))
      Elow = Etrial;
    
    Normalize();
    varPtr->Append(u.Data());
    double A = pot->A(grid(tindex));
    double C = 0.5 * A * CuspValue;
    double u0 = u(tindex);
    double Etest = Eold - 4.0*M_PI*C*u0*u0;
    if ((Etest > Elow) && (Etest < Ehigh) && (NumNodes == TotalNodes))
      Etrial = Etest;
    else 
      Etrial = 0.5 * (Ehigh + Elow);
    
    Energy = Etrial;
    
    if ((NumNodes == TotalNodes) && (fabs(Etrial - Eold) < tolerance))
      done = true;
    Eold = Etrial;       
  }
  IntegrateInOut(tindex);
  Normalize();
  NumNodes = CountNodes();
  if (NumNodes != TotalNodes) {
    cerr << "Node number error!  We have " << NumNodes 
	 << " nodes and want " << TotalNodes << ".\n";
//     IOSectionClass out;
//     out.NewFile ("BadWF.h5");
//     out.WriteVar ("u", u.Data());
//     out.CloseFile();
    cerr << "Energy = " << Energy << endl;
  }
  out.CloseFile();
}


void RadialWF::Normalize()
{
  normVec(0) = 0.0;
  NormalizeDeriv normDeriv(*this);
  RungeKutta<NormalizeDeriv,double> integrator(normDeriv);
  integrator.Integrate(*u.grid, 0, u.grid->NumPoints-1, normVec);
  
  double norm = sqrt(1.0/(4.0*M_PI*normVec(u.grid->NumPoints-1)));
  for (int i=0; i<u.grid->NumPoints; i++) {
    u(i) *= norm;
    dudr(i) *= norm;
  }
}


int RadialWF::CountNodes()
{
  int nodes=0;
  for (int i=0; i<(u.grid->NumPoints-1); i++)
    if (u(i)*u(i+1) < 0.0)
      nodes++;
  return (nodes);
}

void RadialWF::SetGrid(Grid *grid)
{
  int N = grid->NumPoints;
  uduVec.resize(N);
  normVec.resize(N);
  uduVec = 0.0;
  normVec = 0.0;
  u.Init (grid, normVec);
  dudr.Init(grid, normVec);
}

void RadialWF::SetPotential(Potential *newPot)
{
  pot = newPot;
}

// This will write all the information except the potential and the
// grid.  We assume that those are store somewhere else
void RadialWF::Write(IOSectionClass &out)
{
  out.WriteVar ("u", u.Data());
  out.WriteVar ("dudr", dudr.Data());
  out.WriteVar ("n", n);
  out.WriteVar ("l", l);
  out.WriteVar ("CoreNodes", CoreNodes);
  out.WriteVar ("Occupancy", Occupancy);
  out.WriteVar ("Weight", Weight);
  out.WriteVar ("Energy", Weight);
  out.WriteVar ("PartialNorm", PartialNorm);
  out.WriteVar ("Label", Label);
}

// This function assumes the grid has already been set.
void RadialWF::Read (IOSectionClass &in)
{
  bool succ;
  if (u.grid == NULL) {
    cerr << "Grid not set prior to calling RadialWF::Read.\n";
    exit(1);
  }
  in.ReadVar ("u", u.Data());
  in.ReadVar ("dudr", dudr.Data());
  assert(in.ReadVar("n", n));
  assert(in.ReadVar("l", l));
  assert(in.ReadVar("Occupancy", Occupancy));
  in.ReadVar("Energy", Energy, -1.0);
  in.ReadVar("Weight", Weight, 1.0);
  in.ReadVar("PartialNorm",PartialNorm);
  assert(in.ReadVar("Label", Label));
}
