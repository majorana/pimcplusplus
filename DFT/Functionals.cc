#include "Functionals.h"

/////////////////////////////////////////////////////////////////
//                    LSDA Exchange Functions                  //
/////////////////////////////////////////////////////////////////

inline scalar f(scalar zeta)
{
  const scalar FourThirds = 4.0/3.0;
  const scalar TwoToOneThird = 1.25992104989487;//pow(2.0,1.0/3.0);
  return ((pow(1.0+zeta, FourThirds) + pow(1.0-zeta, FourThirds)
	   - 2.0) / (2.0 * (TwoToOneThird-1.0)));
}

inline scalar df_dzeta(scalar zeta)
{
  scalar numer = (4.0/3.0)*(pow(1.0+zeta,1.0/3.0)-pow(1.0-zeta,1.0/3.0));
  scalar denom = 2.0*(pow(2.0,1.0/3.0)-1.0);
  return (numer/denom);
}


void ExchangePotential (scalar nup, scalar ndown,
			scalar &Vup, scalar &Vdown)
{
  scalar n = nup+ndown;
  scalar zeta = (nup-ndown)/n;
  
  const scalar Third = 1.0/3.0;
  const scalar TwoToOneThird = 1.25992104989487;//pow(2.0,Third);
  
    scalar ExP = -3.0*pow(3.0*n/(8.0*M_PI), Third);
  //scalar ExP = -3.0/2.0*pow(3.0*n/(8.0*M_PI), Third);
  scalar ExF = TwoToOneThird * ExP;

  const scalar FourThirds = 4.0/3.0;
  scalar f = (pow(1.0+zeta, FourThirds) + pow(1.0-zeta, FourThirds)
	      - 2.0) / (2.0 * (TwoToOneThird-1.0));

  scalar dExP_dn = -pow(3.0/(8.0*M_PI*n*n), Third);
  scalar dExF_dn = TwoToOneThird * dExP_dn;

  scalar dEx_dn = dExP_dn + (dExF_dn - dExP_dn) * f;
  scalar df_dzeta =
    4.0*Third*(pow(1.0+zeta,Third)-pow(1.0-zeta,Third))/
    (2.0*(TwoToOneThird-1.0));

  scalar Ex = ExP + (ExF - ExP)*f;
  
  //fprintf (stderr, "C++ f = %1.12f\n", f);
  //fprintf (stderr, "C++ zeta = %1.12f\n", zeta);
  //fprintf (stderr, "C++ Ex = %1.12f\n", Ex);

  
  Vup   = Ex + n*dEx_dn - (zeta-1.0)*(ExF - ExP)*df_dzeta;
  Vdown = Ex + n*dEx_dn - (zeta+1.0)*(ExF - ExP)*df_dzeta;
  // Now the original VWN papers were in Rydbergs.  Multiply
  // everything by 0.5 to get in Hartrees
  Vup *= 0.5;
  Vdown *= 0.5;
  Ex *= 0.5;
  if (isnan(Vup))
    Vup = 0.0;
  if (isnan(Vdown))
    Vdown = 0.0;
  if (isnan(Ex))
    Ex = 0.0;
}
  

inline scalar F(scalar n, scalar A, scalar x0, scalar b, scalar c)
{
  const scalar sixth = 1.0/6.0;
  scalar x = pow(3.0/(4.0*M_PI*n),sixth);
  scalar X = x*(x+b) + c;
  scalar X0 = x0*(x0+b) + c;
  scalar Q = sqrt(4.0*c-b*b);
  scalar atanQ = atan(Q/(2.0*x+b));

  scalar term1 = log(x*x/X);
  scalar term2 = (2.0*b/Q)*atanQ;
  scalar term3 = -(b*x0/X0)*log((x-x0)*(x-x0)/X);
  scalar term4 = -(b*x0/X0)*(2.0*(b+2.0*x0)/Q)*atanQ;

  return (A*(term1+term2+term3+term4));
}


scalar dFterms(scalar n, scalar A, scalar x0, scalar b, scalar c)
{
  scalar eps = 1.0e-6;
  scalar np = n+eps;
  scalar nm = n-eps;
  const scalar sixth = 1.0/6.0;
  scalar x = pow(3.0/(4.0*M_PI*np),sixth);
  scalar X = x*(x+b) + c;
  scalar X0 = x0*(x0+b) + c;
  scalar Q = sqrt(4.0*c-b*b);
  scalar atanQ = atan(Q/(2.0*x+b));

  scalar term1p = log(x*x/X);
  scalar term2p = (2.0*b/Q)*atanQ;
  scalar term3p = -(b*x0/X0)*log((x-x0)*(x-x0)/X);
  scalar term4p = -(b*x0/X0)*(2.0*(b+2.0*x0)/Q)*atanQ;

  x = pow(3.0/(4.0*M_PI*nm),sixth);
  X = x*(x+b) + c;
  X0 = x0*(x0+b) + c;
  Q = sqrt(4.0*c-b*b);
  atanQ = atan(Q/(2.0*x+b));

  scalar term1m = log(x*x/X);
  scalar term2m = (2.0*b/Q)*atanQ;
  scalar term3m = -(b*x0/X0)*log((x-x0)*(x-x0)/X);
  scalar term4m = -(b*x0/X0)*(2.0*(b+2.0*x0)/Q)*atanQ;

  //fprintf (stderr, "dterm1 = %1.12f\n", A*(term1p-term1m)/(2.0*eps));
  //fprintf (stderr, "dterm2 = %1.12f\n", A*(term2p-term2m)/(2.0*eps));
  //fprintf (stderr, "dterm3 = %1.12f\n", A*(term3p-term3m)/(2.0*eps));
  //fprintf (stderr, "dterm4 = %1.12f\n", A*(term4p-term4m)/(2.0*eps));

  return (A*(term1p+term2p+term3p+term4p));
}



scalar dF_dn(scalar n, scalar A, scalar x0, scalar b, scalar c)
{
  const scalar sixth = 1.0/6.0;
  scalar x = pow(3.0/(4.0*M_PI*n),sixth);
  scalar X = x*(x+b) + c;
  scalar X0 = x0*(x0+b) + c;
  scalar Q = sqrt(4.0*c-b*b);
 
  //fprintf (stderr, "C++ x = %1.12f\n", x);
  //fprintf (stderr, "Q = %1.12f\n", Q);
 
  scalar n3 = n*n*n;
  scalar n7 = n3*n3*n;
  scalar prefactor = -(A/6.0) * pow(3.0/(4.0*M_PI*n7),sixth);
  scalar bp2x = 2.0*x + b;
  scalar term1 = 2.0/x - (bp2x)/X;
  scalar Q2m_bp2x_2_inv = 1.0/(Q*Q + (bp2x*bp2x));
  scalar term2 = -4.0 * b * Q2m_bp2x_2_inv;
  scalar term34pre = -b*x0/X0;
  scalar term3 = 2.0/(x-x0) - bp2x/X;
  scalar term4 = -4.0*(b+2.0*x0)*Q2m_bp2x_2_inv;

  //dFterms (n, A, x0, b, c);
  //fprintf (stderr, "term1 = %1.12f\n", prefactor*term1);
  //fprintf (stderr, "term2 = %1.12f\n", prefactor*term2);
  //fprintf (stderr, "term3 = %1.12f\n", term34pre*prefactor*term3);
  //fprintf (stderr, "term4 = %1.12f\n", term34pre*prefactor*term4);

  return (prefactor*(term1+term2+term34pre*(term3+term4)));
}


scalar dF_dn_FD (scalar n, scalar A, scalar x0, scalar b, scalar c)
{
  const scalar eps = 1.0e-6;
  scalar Fp = F(n+eps, A, x0, b, c);
  scalar Fm = F(n-eps, A, x0, b, c);
  return ((Fp-Fm)/(2.0*eps));
}



scalar Ec(scalar nup, scalar ndown)
{
  scalar n = nup + ndown;
  scalar zeta = (nup-ndown)/n;
  
  scalar EcP =    F(n, 0.0310907, -0.10498, 3.72744, 12.9352);
  //fprintf (stderr, "EcP = %1.12f\n", EcP);
  scalar EcF =    F(n, 0.01554535, -0.325, 7.06042, 18.0578);
  //fprintf (stderr, "EcF = %1.12f\n", EcF);
  scalar alphac = F(n, -1.0/(6.0*M_PI*M_PI), -0.00475840,
		    1.13107, 13.0045);
  //fprintf (stderr, "alphac = %1.12f\n", alphac);

  scalar f_zeta = f(zeta);
  scalar f_doubleprime = 4.0/(9*(pow(2,1.0/3.0)-1.0));
  scalar beta = f_doubleprime*(EcF-EcP)/alphac -1.0;
  scalar zeta2 = zeta*zeta;
  scalar zeta4 = zeta2*zeta2;
  scalar deltaEc = (alphac*f_zeta/f_doubleprime) * (1.0 + beta * zeta4);
  
  return (EcP + deltaEc);
}


void CheckCorrelationPotential (scalar  nup, scalar ndown)
{
  scalar eps = 1.0e-4;
  
  scalar n = nup + ndown;
  scalar zeta = (nup - ndown) / n;

  scalar np = n+eps;
  scalar nm = n-eps;
  scalar zetap = zeta + eps;
  scalar zetam = zeta - eps;

  scalar nupp   = 0.5*(np + np * zeta);
  scalar ndownp = 0.5*(np - np*zeta);
  scalar nupm   = 0.5*(nm + nm*zeta);
  scalar ndownm = 0.5*(nm - nm*zeta); 

  scalar Ecplus =  Ec(nupp, ndownp);
  scalar Ecminus = Ec(nupm, ndownm);

  scalar dEc_dn = (Ecplus - Ecminus)/(2.0*eps);

  scalar dEcP_dn = dF_dn(n, 0.0310907, -0.10498, 3.72744, 12.9352); 
  scalar ddeltaEc_dn = dEc_dn - dEcP_dn;
  
  //fprintf (stderr, "FD: ddeltaEc_dn = %1.12f\n", ddeltaEc_dn);




  scalar EcPp =    F(np, 0.0310907, -0.10498, 3.72744, 12.9352);
  scalar EcFp =    F(np, 0.01554535, -0.325, 7.06042, 18.0578);
  scalar alphacp = F(np, -1.0/(6.0*M_PI*M_PI), -0.00475840,
		     1.13107, 13.0045);

  scalar EcPm =    F(nm, 0.0310907, -0.10498, 3.72744, 12.9352);
  scalar EcFm =    F(nm, 0.01554535, -0.325, 7.06042, 18.0578);
  scalar alphacm = F(nm, -1.0/(6.0*M_PI*M_PI), -0.00475840,
		     1.13107, 13.0045);

  scalar f_doubleprime = 4.0/(9*(pow(2.0,1.0/3.0)-1.0));

  scalar betap = f_doubleprime*(EcFp-EcPp)/alphacp -1.0;
  scalar betam = f_doubleprime*(EcFm-EcPm)/alphacm -1.0;

  scalar zeta2 = zeta*zeta;
  scalar zeta4 = zeta2*zeta2;
  scalar deltaEcp = alphacp*f(zeta)/f_doubleprime * (1.0+betap*zeta4);
  scalar deltaEcm = alphacm*f(zeta)/f_doubleprime * (1.0+betam*zeta4);

  fprintf (stderr, "FD2: ddeltaEc_dn = %1.12f\n",
	   (deltaEcp-deltaEcm)/(2.0*eps)); 

  fprintf (stderr, "FD: dbeta_dn = %1.12f\n", (betap-betam)/(2.0*eps));

}



void CorrelationPotential(scalar  nup, scalar ndown,
			  scalar &Vup, scalar &Vdown)
{
  scalar EC = Ec(nup, ndown);
  //fprintf (stderr, "C++ EC = %1.12f\n", EC);
  scalar n = nup + ndown;
  scalar zeta = (nup - ndown)/n;
  scalar zeta2 = zeta*zeta;
  scalar zeta3 = zeta * zeta2;
  scalar zeta4 = zeta2*zeta2;

  scalar EcP =    F(n, 0.0310907, -0.10498, 3.72744, 12.9352);
  scalar EcF =    F(n, 0.01554535, -0.325, 7.06042, 18.0578);
  scalar alphac = F(n, -1.0/(6.0*M_PI*M_PI), -0.00475840,
		    1.13107, 13.0045);
  scalar dEcP_dn = dF_dn(n, 0.0310907, -0.10498, 3.72744, 12.9352); 
  //scalar dEcP_dn = dF_dn_FD(n, 0.0310907, -0.10498, 3.72744, 12.9352); 
  scalar dEcF_dn = dF_dn(n, 0.01554535, -0.325, 7.06042, 18.0578);
  scalar dalphac_dn = dF_dn(n, -1.0/(6.0*M_PI*M_PI), -0.00475840,
			    1.13107, 13.0045);  

  scalar f_zeta = f(zeta);
  scalar f_prime = df_dzeta(zeta);
  //fprintf (stderr, "f_prime = %1.12f\n", f_prime);
  scalar f_doubleprime = 4.0/(9.0*(pow(2.0,1.0/3.0)-1.0));
  scalar beta = f_doubleprime*(EcF - EcP)/alphac -1.0;
  scalar dbeta_dn = f_doubleprime * ((dEcF_dn -dEcP_dn)/alphac -
				     (EcF-EcP)*dalphac_dn/(alphac*alphac));

  //fprintf (stderr, "dbeta_dn     = %1.12f\n", dbeta_dn);
  scalar ddeltaEc_dn = 
    dalphac_dn * f_zeta/f_doubleprime *  (1.0+beta*zeta4) +
    alphac * (f_zeta/f_doubleprime) * dbeta_dn * zeta4;

  //fprintf (stderr, "ddeltaEc_dn =     %1.12f\n", ddeltaEc_dn);
  //CheckCorrelationPotential(nup, ndown);

  scalar ddeltaEc_dzeta =
    alphac/f_doubleprime*((1.0+beta*zeta4)*f_prime + 4.0*beta*f_zeta*zeta3);
  
  //fprintf (stderr, "ddeltaEc_dzeta = %1.12f\n", ddeltaEc_dzeta);

  Vup   = EC + n * (dEcP_dn + ddeltaEc_dn) - (zeta-1.0) * ddeltaEc_dzeta;
  Vdown = EC + n * (dEcP_dn + ddeltaEc_dn) - (zeta+1.0) * ddeltaEc_dzeta;

  if (isnan(Vup))
    Vup = 0.0;
  if (isnan(Vdown))
    Vdown = 0.0;
}

#ifdef NOUNDERSCORE 
#define FORT(name) name
#else
#define FORT(name) name_
#endif 

extern "C" void FORT(exccor)(double &n, double &zeta, double &exc, double &vxc, 
			     double &vpol, int &type, int &Macdonald_Vosko);

void FortranExCorr(scalar  nup, scalar  ndown,
		   scalar &Vup, scalar &Vdown)
{
  double n = nup + ndown;
  double zeta = (nup-ndown)/n;
  
  int type = 4;
  int Macdonald_Vosko=/*0*/1;
  
  double exc, vxc, vpol;
  FORT(exccor)(n, zeta, exc, vxc, vpol, type, Macdonald_Vosko);

  //  fprintf (stderr, "Fortran Exc = %12.8f\n", exc);

  Vup = vxc + vpol;
  Vdown = vxc - vpol;
}
  


double FortranXCE (scalar nup, scalar ndown)
{
  double n = nup + ndown;
  double zeta = (nup-ndown)/n;
  
  int type = 4;
  int Macdonald_Vosko=/*0*/1;
  
  double exc, vxc, vpol;
  FORT(exccor)(n, zeta, exc, vxc, vpol, type, Macdonald_Vosko);
  return (exc);
}
  
  



//  scalar
//  LDA::ExchangePot(scalar r, const Array<RadialWF,1> &WFs)
//  {

//    return (0.0);
//  }


//  scalar 
//  LDA::CorrelationPot(scalar r, const Array<RadialWF,1> &WFs)
//  {

//    return (0.0);

//  }
