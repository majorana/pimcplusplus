#include "NewAtom.h"
#include "../Integration/RungeKutta.h"
#include "../DFT/Functionals.h"

void DFTAtom::UpdateVHXC()
{
  UpdateChargeDensity();
  UpdateHartree();
  UpdateExCorr();

  for (int i=0; i < grid->NumPoints; i++)
    V.HXC(i) = Hartree(i) + ExCorr(i);
}

/// Radial WFs must be normalized before calling this function
void DFTAtom::UpdateChargeDensity()
{
  for (int i=0; i<grid->NumPoints; i++) {
    ChargeDensity(i) = 0.0;
    double r = (*grid)(i);
    double rinv2 = 1.0/(r*r);
    for (int j=0; j<RadialWFs.size(); j++) {
      double u = RadialWFs(j).u(i);
      ChargeDensity(i) += RadialWFs(j).Occupancy * u*u * rinv2;
    }
  }
}


class HartreeDeriv1
{
  DFTAtom &atom;
public:
  inline double operator()(double r, double sum)
  { return atom.Hartree1(r,sum); }

  HartreeDeriv1(DFTAtom &newAtom) : atom(newAtom) 
  { /* Do nothing */ }
};

class HartreeDeriv2
{
  DFTAtom &atom;
public:
  inline double operator()(double r, double sum)
  { return atom.Hartree2(r,sum); }

  HartreeDeriv2(DFTAtom &newAtom) : atom(newAtom) 
  { /* Do nothing */ }
};


void DFTAtom::UpdateHartree()
{
  int N = grid->NumPoints;
  
  HartreeDeriv1 H1(*this);
  HartreeDeriv2 H2(*this);
  RungeKutta<HartreeDeriv1,double> integrator1(H1);
  RungeKutta<HartreeDeriv2,double> integrator2(H2);
  temp(0) = 0.0;  temp2(0) = 0.0;
  integrator1.Integrate(*grid, 0, N-1, temp);
  integrator2.Integrate(*grid, 0, N-1, temp2);

  for (int i=0; i<N; i++) {
    double r = (*grid)(i);
    Hartree(i) = 4.0*M_PI*(temp(i)/r + (temp2(N-1)-temp2(i)));
  }
}

void DFTAtom::UpdateExCorr()
{
  int N = grid->NumPoints;
  for (int i=0; i<N; i++) {
    double Vxcup, Vxcdown;
    double upCharge = 0.5*ChargeDensity(i);
    double downCharge = upCharge;
    FortranExCorr(upCharge, downCharge, Vxcup, Vxcdown);
    ExCorr(i) = Vxcup;
  }
}


void DFTAtom::SetGrid(Grid *newgrid)
{
  grid = newgrid;
  int N = grid->NumPoints;
  temp.resize(N);
  temp2.resize(N);
  temp = 0.0;
  V.HXC.Init(grid, temp);
  ChargeDensity.Init(grid,temp);
  Hartree.Init(grid, temp);
  ExCorr.Init(grid,temp);
  for (int i=0; i<RadialWFs.size(); i++)
    RadialWFs(i).SetGrid(grid);
}

void DFTAtom::SetBarePot(Potential *newPot)
{
  BarePot = newPot;
  V.BarePot = BarePot;
  for (int i=0; i<RadialWFs.size(); i++)
    RadialWFs(i).SetPotential(&V);
}


void DFTAtom::SolveSC()
{
  int N = grid->NumPoints; 
  // Just rename temp and temp2 for clarity
  Array<double,1> &oldCharge = temp;
  Array<double,1> &newCharge = temp2;

  // First, zero out screening
  for (int i=0; i<N; i++)
    V.HXC(i) = 0.0;
  
  // Now solve radial equations
  for (int i=0; i<RadialWFs.size(); i++)
    RadialWFs(i).Solve();

  UpdateChargeDensity();
  oldCharge = 0.0;
  newCharge = ChargeDensity.Data();

  bool done = false;
  while (!done) {
    for (int i=0; i<N; i++)
      ChargeDensity(i) = newMix*(newCharge(i)) + (1.0-newMix)*oldCharge(i);
    UpdateHartree();
    UpdateExCorr();
    for (int i=0; i < grid->NumPoints; i++)
      V.HXC(i) = Hartree(i) + ExCorr(i);

    for (int i=0; i<RadialWFs.size(); i++) {
      RadialWFs(i).Solve();
      fprintf (stderr, "Energy(%d) = %1.16f\n", i, RadialWFs(i).Energy);
    }
    oldCharge = ChargeDensity.Data();
    UpdateChargeDensity();
    newCharge = ChargeDensity.Data();
  }
  

}


void DFTAtom::Read(IOSectionClass &in) 
{

}

void DFTAtom::Write(IOSectionClass &out)
{


}

void DFTAtom::CalcEnergies(double &kinetic, double &potential,
			   double &hartree, double &XC)
{

}
