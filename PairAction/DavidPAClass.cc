#include "DavidPAClass.h"

#ifdef MAKE_FIT

void DavidPAClass::ReadParams(IOSectionClass &in)
{

}

void DavidPAClass::WriteBetaIndependentInfo (IOSectionClass &out)
{

}

void DavidPAClass::Error (Rho &rho, double &Uerror, double &dUerror)
{

}

void DavidPAClass::AddFit (Rho &rho)
{

}

void DavidPAClass::WriteFits(IOSectionClass &out)
{


}
#endif

double DavidPAClass::U (double q, double z, double s2, int level)
{
  double s=sqrt(s2);
  double uTemp;
  double duTemp;
  double vTemp;
  calcUsqz(s,q,z,level,uTemp,duTemp,vTemp);
  return uTemp;

}

double DavidPAClass::VV(double q, double z, double s2, int level)
{
  double s=sqrt(s2);
  double uTemp;
  double duTemp;
  double vTemp;
  calcUsqz(s,q,z,level,uTemp,duTemp,vTemp);
  return vTemp;

}

double DavidPAClass::dU(double q, double z, double s2, int level)
{
  double s=sqrt(s2);
  double uTemp;
  double duTemp;
  double vTemp;
  calcUsqz(s,q,z,level,uTemp,duTemp,vTemp);
  //  cerr<<"my vtemp is "<<vTemp<<endl;
  return duTemp;


}
