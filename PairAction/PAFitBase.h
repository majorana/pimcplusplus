#ifndef PA_FIT_BASE_H
#define PA_FIT_BASE_H

#include "../IO/InputOutput.h"
#include "../PH/Potential.h"
#include "../Splines/CubicSpline.h"
#include "Particle.h"
#ifdef MAKE_FIT
#include "../../Rho.h"
#endif

class PairActionFitClass
{
protected:
  void ReadHeader(IOSectionClass &inSection);
public:
  ParticleClass Particle1, Particle2;
  bool UsePBC;
  int NumBetas;
  Array<double,1> Box;
  Potential *Pot;
  double SmallestBeta;
  double lambda;



#ifdef MAKE_FIT
  virtual void ReadParams  (IOSectionClass &inSection)   = 0;
  virtual void WriteBetaIndependentInfo (IOSectionClass &outSection) = 0;
  virtual void AddFit (Rho &rho) = 0;
  virtual void Error (Rho &rho, double &Uerror, double &dUerror)=0;
  virtual void WriteFits(IOSectionClass &outSection) = 0;
#endif
  virtual bool Read (IOSectionClass &inSection,
		     double smallestBeta, int NumBetas) = 0;
  virtual double U(double q, double z, double s2, int level) = 0;
  virtual double dU(double q, double z, double s2, int level) = 0;
  virtual double V(double r) {return -3.14159;}
};




// class PA2DFitClass : public PairActionFitClass
// {
// private:
//   bool GridIsMine;
//   Array<double, 1> Coefs;
// public:
//   bool UsePBC;
//   int Order;
//   TinyVector<double,3> Box;

//   // Real space parameters
//   /// The array is over different values of beta
//   Array<MultiBicubicSpline,1> Uk;
//   Array<MultiBicubicSpline,1> dUk;

//   /// k-space paramters

//   // Member functions
// #ifdef MAKE_FIT
//   void AddFit (Rho &rho);
// #endif
//   bool ReadParams (IOSectionClass &inSection);
//   void Write (IOSectionClass &outSection);
//   bool Read  (IOSectionClass &inSection);
//   double U(double r, double rp, double costheta);
// }



// class UfitClass
// {
// public:
//   int Order;
//   MultiCubicSpline Ukj;
//   Array<double,1> Coefs;
//   inline double operator() (double q, double s, double z)
//   {
//     Ukj(q, Coefs);
//     double r = q + 0.5*z;
//     double rp = q - 0.5*z;
//     // Diagonal contribution
//     double Usum = 0.5 * (Ukj(0, r) + Ukj(0,rp));
//     // Off-diagonal contribution
//     double z2 = z*z;
//     double s2 = s*s;
//     double s2inv = (s2==0.0) ? 0.0 : 1.0/s2;
//     double sto2k = s2;
//     for (int k=1; k<=Order; k++) {
//       double zto2j = 1.0;
//       double currS = sto2k;
//       for (int j=0; j<=k; j++) {
// 	double Ucoef = Coefs(k*(k+1)/2+j);
// 	Usum+=Ucoef*zto2j*currS;
// 	zto2j*=z2;
// 	currS*=s2inv;
//       }
//       sto2k*=s2;
//     }
//     return (Usum);
//   }
//   void Init (Grid *qgrid, int order, Array<double,2> Ukjvals)
//   {
//     Order = order;
//     Ukj.Init(qgrid, Ukjvals);
//     Coefs.resize(Ukjvals.cols());
//   }
// };





// #include "Common/SpecialFunctions/HermitePoly.h"
// class Ufit2Class
// {
// public:
//   int Order;
//   MultiBicubicSpline Uk, dUk;
//   Array<double,1> Coefs;
//   Array<double,1> Hn;
//   double sigmainv;

//   inline double operator() (double q, double s, double z)
//   {
//     double r = q + 0.5*z;
//     double rp = q - 0.5*z;
//     Uk(r, rp, Coefs);
//     double Usum = 0.0;
// //     double smin = fabs (r-rp);
// //     double delta = s - smin;
// //     double stok = 1.0;
// //     for (int k=0; k<=Order; k++) {
// //       Usum += stok * Coefs(k);
// //       stok *= s*s;
// //     }
//     double t = s * sigmainv;
//     HermitePoly (t, Hn);
//     for (int n=0; n<Coefs.size(); n++)
//       Usum += Coefs(n) * Hn(2*n);
//     return (Usum);
//   }

//   inline double dU (double q, double s, double z)
//     {
//     double r = q + 0.5*z;
//     double rp = q - 0.5*z;
//     dUk(r, rp, Coefs);
//     double dUsum = 0.0;
//     double t = s * sigmainv;
//     HermitePoly (t, Hn);
//     for (int n=0; n<Coefs.size(); n++)
//       dUsum += Coefs(n) * Hn(2*n);
//     return (dUsum);
//   }


//   void Init (Grid *rgrid, Array<double,3> &Ukvals, Array<double,3> &dUkvals)
//   {
//     Order = Ukvals.extent(2)-1;
//     Uk.Init(rgrid, rgrid, Ukvals);
//     dUk.Init(rgrid, rgrid, dUkvals);
//     Coefs.resize(Ukvals.extent(2));
//     Hn.resize(2*Ukvals.extent(2));
//   }
// };



#endif
