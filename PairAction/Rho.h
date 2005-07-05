#ifndef RHO_H
#define RHO_H

#include "U_l.h"
#include "Particle.h"


class Rho
{
private:
  double beta;
  CommunicatorClass WorldComm, GroupComm;
  int NumGroups;
  int ProcsPerGroup;
  int NumChannels;
  int ChannelsPerGroup;
  int MyGroup;

  inline int Index(int row, int col) const
  {
    if (row > col)
      return ((row*(row+1))/2+col);
    else
      return ((col*(col+1))/2+row);
  }
  void SetSCscale();
  void UdU_PH (double r, double rp, double costheta, 
	       Array<double,1> &Ulvec, Array<double,1> &dUlvec,
	       double &U, double &dU);
  void UdU_local (double r, double rp, double costheta, 
		  Array<double,1> &Ulvec, Array<double,1> &dUlvec,
		  double &U, double &dU);


  FILE *DebugFile;
public:
  Array<U_l, 1> U_ls;
  Grid *grid;
  Potential *Pot;
  CoreTransform Transform;
  double lambda;
  bool Extrapolate_ls;
  bool WriteDebug;

  void Initialize(int lmax, double this_lambda, 
		  double FinalBeta, int NumSquares,
		  Grid *grid_, Potential *pot_, double TailPower,
		  double AbsTol, double RelTol);
  void Initialize(int lmax, double this_lambda, 
		  Grid *grid_, Potential *pot);
  /// Gathers all the U_ls from all the processor groups
  void GatherU_ls();
  /// Distributes all the U_ls from processor 0 to all other
  /// processors.
  void BroadcastU_ls();
  /// Do a squaring, doubling beta
  void Square();
  void U_lArray(double r, double rp, 
		Array<double,1> &Ulvec, Array<double,1> &dUlvec);
  double U (double r, double rp, double costheta, Array<double,1> &Uvec);
  double U (double r, double rp, double costheta);
  void UdU (double r, double rp, double costheta, 
	    double &U, double &dU);
  void UdU (double r, double rp, double costheta, 
	    Array<double,1> &Ulvec, Array<double,1> &dUlvec,
	    double &U, double &dU);
  void UdU_Coulomb (double r, double rp, double costheta, 
		    double &U, double &dU);
  void Write(IOSectionClass &outSection);
  void WriteU_ls (IOSectionClass &outSection);
  void ReadU_ls(IOSectionClass &inSection);
  void Read (IOSectionClass &inSection);
  inline double Beta() {return beta;}
  Rho()
  {
    lambda = 0.5;
    Extrapolate_ls=true;
  }
};




#endif
