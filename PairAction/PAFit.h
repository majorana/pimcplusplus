#ifndef PA_FIT_H
#define PA_FIT_H

#include "Common/IO/InputOutput.h"
#include "Particle.h"
#ifdef MAKE_FIT
#include "Rho.h"
#endif

class PairActionFitClass
{
public:
  ParticleClass Particle1, Particle2;
  bool UsePBC;
  PseudoHamiltonian *PH;
  
#ifdef MAKE_FIT
  virtual void AddFit(Rho &rho) = 0;
#endif
  virtual bool ReadParams (IOSectionClass &inSection);
  virtual void Write (IOSectionClass &outSection);
  virtual bool Read (IOSectionClass &inSection);
  virtual double U(double r, double rp, double costheta) = 0;
};


class PA1DFitClass : public PairActionFitClass
{
private:
  bool GridIsMine;
  Array<double, 1> Coefs;
public:
  int Order;
  double SmallestBeta;
  TinyVector<double,3> Box;
  Array<MultiCubicSpline,1> Ukj;
  Array<MultiCubicSpline,1> dUkj;

#ifdef MAKE_FIT
  void AddFit (Rho &rho);
#endif
  bool ReadParams (IOSectionClass &inSection);
  void Write (IOSectionClass &outSection);
  bool Read  (IOSectionClass &inSection);
  double U(double r, double rp, double costheta);
}


class PA2DFitClass : public PairActionFitClass
{
private:
  bool GridIsMine;
  Array<double, 1> Coefs;
public:
  bool UsePBC;
  int Order, NumBetas;
  double SmallestBeta;
  TinyVector<double,3> Box;

  // Real space parameters
  /// The array is over different values of beta
  Array<MultiBicubicSpline,1> Uk;
  Array<MultiBicubicSpline,1> dUk;

  /// k-space paramters

  // Member functions
#ifdef MAKE_FIT
  void AddFit (Rho &rho);
#endif
  bool ReadParams (IOSectionClass &inSection);
  void Write (IOSectionClass &outSection);
  bool Read  (IOSectionClass &inSection);
  double U(double r, double rp, double costheta);
}

#endif
