#ifndef NEW_ATOM_H
#define NEW_ATOM_H

#include "RadialWF.h"

class Atom
{
protected:
  Grid *grid;
public:
  Array<RadialWF,1> RadialWFs;
  virtual void CalcVHXC() = 0;
  virtual void CalcEnergies (double &kinetic, double &potential,
			     double &exCorr);
  virtual void Write(IOSectionClass &out) = 0;
  virtual void Read(IOSectionClass &in) = 0;
  virtual void SetGrid(Grid *newGrid) = 0;
  virtual void SetBarePot (Potential *pot);
};

class DFTAtom : public Atom
{
private:
  CubicSpline ChargeDensity;
  CubicSpline Hartree;
  CubicSpline ExCorr;
  void UpdateChargeDensity();
  void UpdateHartree();
  void UpdateExCorr();
  Array<double,1> temp;
  Potential *BarePot;
public:
  XCPot xcpot;
  
  /// This function calculates the charge density, hartree and exchange
  /// potentials and places them in pot.
  void CalcVHXC();
  void Write (IOSectionClass &out);
  void Read  (IOSectionClass &in);
  void SetGrid (Grid *newGrid);
  void SetBarePot (Potential *pot);
  
  inline double Hartree1 (double r, double sum);
  inline double Hartree2 (double r, double sum);
};

inline double DFTAtom::Hartree1 (double r, double sum)
{
  double rho = ChargeDensity(r);
  return (rho*r*r);
}

inline double DFTAtom::Hartree2 (double r, double sum)
{
  double rho = ChargeDensity(r);
  return (rho*r);
}



#endif
