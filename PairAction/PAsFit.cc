#include "PAFit.h"
#include "../Splines/BicubicSpline.h"
#include "../SpecialFunctions/LegendrePoly.h"

const double Rho0Min = 1.0e-4;

/// The following routines are used only if we are creating fits, not
/// using them.
#ifdef MAKE_FIT
void PAsFitClass::ReadParams(IOSectionClass &inSection)
{
  assert(inSection.OpenSection ("qGrid"));
  qgrid = ReadGrid (inSection);
  inSection.CloseSection();
  assert(inSection.OpenSection ("yGrid"));
  ygrid = ReadGrid (inSection);
  inSection.CloseSection();
  assert(inSection.ReadVar("Order", Order));
  Coefs.resize(Order+1);
  Pn.resize(Order+1);
  GridsAreMine = true;
  UsePBC = inSection.ReadVar ("Box", Box);
}

void PAsFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
{
  outSection.WriteVar ("Type", "sfit");
  outSection.NewSection("qGrid");
  qgrid->Write(outSection);
  outSection.CloseSection();
  outSection.NewSection("yGrid");
  ygrid->Write(outSection);
  outSection.CloseSection();
  outSection.WriteVar("Order", Order);
}


class sFitIntegrand
{
private:
  double q, z;
  double r, rp;
  double sMin, sMax;
public:
  Rho &rho;
  bool IsdU;
  int k;

  Array<double,1> Ul, dUl;
  
  void Setqz(double newq, double newz,
	     double smin, double smax)
  {
    q = newq; z = newz;
    sMin = smin; sMax = smax;
    
    r = q + 0.5*z;
    rp = q - 0.5*z;
    rho.U_lArray(r,rp, Ul, dUl);
  }

  inline double operator()(double x)
  {
    double s=sMin+0.5*(sMax-sMin)*(x+1.0);
    double costheta;
    
    if ((r*rp)==0.0)
      costheta = 0.0;
    else
      costheta = (r*r + rp*rp - s*s)/(2.0*r*rp);
    
    costheta = min(1.0, costheta);
    costheta = max(-1.0, costheta);
    double U, dU;
    rho.UdU(r,rp,costheta, Ul, dUl, U, dU);
    if (isnan(U)) {
      cerr << "r = " << r << " rp = " << rp << endl;
      cerr << "costheta = " << costheta << endl;
      cerr << "x = " << x << endl;
    }
    // U = 0.0;
    //cerr << "r = " << r << " rp = " << rp << endl;
    //cerr << "smin = " << sMin << " smax = " << sMax << endl;
    //cerr << "x = " << x << endl;
    //cerr << "costheta = " << costheta << endl;
    if (!IsdU)
      return (0.5*(2.0*k+1.0)*U*LegendrePoly(k,x));
    else
      return (0.5*(2.0*k+1.0)*dU*LegendrePoly(k,x));
  }
  
  sFitIntegrand(Rho &myRho) : rho(myRho)
  {
    Ul.resize(rho.U_ls.size());
    dUl.resize(rho.U_ls.size());
  };
};


/// y = |z/z_max|
/// z_max = 
void PAsFitClass::AddFit (Rho &rho)
{
  const double Tolerance = 1.0e-7;
  NumBetas++;
  Usplines.resizeAndPreserve(NumBetas);
  dUsplines.resizeAndPreserve(NumBetas);
  sMax.resizeAndPreserve(NumBetas);

  double lambda = Particle1.lambda + Particle2.lambda;
  double beta = rho.Beta();
  sMax(NumBetas-1) = sqrt(-4.0*lambda*beta*log(Rho0Min));

  int numq = qgrid->NumPoints;
  int numy = ygrid->NumPoints;
  Array<double,3> Umat(numq, numy,Order+1), dUmat(numq, numy, Order+1);
  sFitIntegrand integrand(rho);
  for (int qi=0; qi<numq; qi++) {
    cerr << "qi = " << qi+1 << " of " << numq << ".\n";
    double q = (*qgrid)(qi);
    double zmax = 0.999999*min(2.0*q, sMax(NumBetas-1));
    for (int yi=0; yi<numy; yi++) {
      //cerr << "  yi = " << yi+1 << " of " << numy << ".\n";
      double z = zmax*(*ygrid)(yi);
      double smin = z;
      double smax = min (2.0*q, sMax(NumBetas-1));
      //cerr << "smin = " << smin << " smax = " << smax << endl;
      integrand.Setqz(q,z,smin,smax);
      for(int k=0; k<=Order; k++) {
	//cerr << " k = " << k << endl;
	integrand.k=k;
	integrand.IsdU=false;
	GKIntegration<sFitIntegrand,GK15> Uintegrator(integrand);
	// Accept a relative OR absolute toleraance of 1.0e-7
	Umat(qi, yi, k) = Uintegrator.Integrate(-1.0, 1.0, 
						beta*Tolerance, Tolerance, 
						false);

	integrand.IsdU = true;
	GKIntegration<sFitIntegrand,GK15> dUintegrator(integrand);
	dUmat(qi, yi, k) = dUintegrator.Integrate(-1.0, 1.0, Tolerance,
						  Tolerance, false);
      }
    }
  }
  // Initialize the bicubic splines
  Usplines(NumBetas-1).Init(qgrid,ygrid,Umat);
  dUsplines(NumBetas-1).Init(qgrid,ygrid,dUmat);
}



void PAsFitClass::Error(Rho &rho, double &Uerror, double &dUerror)
{
  int level = (int)floor(log(rho.Beta()/SmallestBeta)/log(2.0)+ 0.5);

  double U2err = 0.0;
  double dU2err = 0.0;
  double weight;
  FILE *Uxdat = fopen ("Ux.dat", "w");
  FILE *Ufdat = fopen ("Uf.dat", "w");
  FILE *dUxdat = fopen ("dUx.dat", "w");
  FILE *dUfdat = fopen ("dUf.dat", "w");
  FILE *tdat = fopen ("t.dat", "w");
  FILE *costhetadat = fopen ("costheta.dat", "w");
  FILE *ydat = fopen ("y.dat", "w");
  LinearGrid qgrid2(qgrid->Start, qgrid->End, 20);
  for (int qi=0; qi<qgrid2.NumPoints; qi++) {
    double q = qgrid2(qi);
    // HACK
    //q = 1.0;
    double zmax = 0.9999*min(2.0*q,sMax(level));
    LinearGrid zgrid(0.0, zmax, 20);
    for (int zi=0; zi<zgrid.NumPoints; zi++) {
      double z = zgrid(zi);
      double y = z/zmax;
      double smax = 0.9999*min(2.0*q,sMax(level));
      LinearGrid sgrid(z, smax, 20);
      for (int si=0; si<sgrid.NumPoints; si++) {
	double s = sgrid(si);
	//cerr << "q = " << q << " z = " << z << " s = " << s << endl;
	double t = s/smax;
	double w = exp(-s*s/(4.0*rho.lambda*rho.Beta()));
	double Uex, dUex, Ufit, dUfit;
	double r, rp, costheta;
	r  = q+0.5*z;
	rp = q-0.5*z;
	if (q == 0.0)
	  costheta = 1.0;
	else
	  costheta = (r*r + rp*rp - s*s)/(2.0*r*rp); 
	costheta = min(costheta,1.0);
	costheta = max(costheta,-1.0);


	rho.UdU(r, rp, costheta, Uex, dUex);
	if (!isnan(Uex) && !isnan(dUex)) {
	  Ufit = U(q, z, s*s, level);
	  dUfit = dU(q, z, s*s, level);
	  U2err += w*(Uex-Ufit)*(Uex-Ufit);
	  dU2err += w*(dUex-dUfit)*(dUex-dUfit);
	  weight += w;
	}
	fprintf (Uxdat, "%1.16e ", Uex);
	fprintf (Ufdat, "%1.16e ", Ufit);
	fprintf (dUxdat, "%1.16e ", dUex);
	fprintf (dUfdat, "%1.16e ", dUfit);
	fprintf (tdat, "%1.16e ", t);
	fprintf (costhetadat, "%1.16e ", costheta);
	fprintf (ydat, "%1.16e ", y);
      }
      fprintf (Uxdat, "\n");
      fprintf (Ufdat, "\n");
      fprintf (dUxdat, "\n");
      fprintf (dUfdat, "\n");
      fprintf (tdat, "\n");
      fprintf (ydat, "\n");
      fprintf (costhetadat, "\n");
    }
  }
  fclose (Uxdat); fclose(Ufdat); fclose(tdat); fclose(ydat); 
  fclose(costhetadat);
  Uerror = sqrt(U2err/weight);
  dUerror = sqrt(dU2err/weight);
}


void PAsFitClass::WriteFits (IOSectionClass &outSection)
{
  Array<double,3> Umat(qgrid->NumPoints, ygrid->NumPoints, Order+1); 
  Array<double,3> dUmat(qgrid->NumPoints, ygrid->NumPoints,Order+1); 
  double beta = SmallestBeta;
  for (int bi=0; bi<NumBetas; bi++) {
    cerr << "Writing Beta "<< bi+1 << " of " << NumBetas << ":\n";
    outSection.NewSection("Fit");
    outSection.WriteVar ("beta", beta);
    double smax = sMax(bi);
    outSection.WriteVar ("sMax", smax);
    for (int qi=0; qi<qgrid->NumPoints; qi++)
      for (int yi=0; yi<ygrid->NumPoints; yi++) 
	for (int j=0; j<=Order; j++) {
	  Umat(qi,yi,j)  =  Usplines(bi)(qi,yi,j);
	  dUmat(qi,yi,j) = dUsplines(bi)(qi,yi,j);
	}
    outSection.WriteVar ("Umat", Umat);
    outSection.WriteVar ("dUmat", dUmat);
    outSection.CloseSection();
    beta *= 2.0;
  }
}

#endif


double PAsFitClass::U(double q, double z, double s2, int level)
{
  z = fabs(z);
  double qmax = qgrid->End*1.000001;
  double zmax = 1.000001*min (2.0*q, sMax(level));
  double smax = 1.000001*min (2.0*q, sMax(level));
  double smin = z;
  double s=sqrt(s2);
  double x;

  if ((q<=qmax)&&(z<=zmax)&&(s<=smax)) {
    if (q == 0) {
      Usplines(level)(0.0,0.0,Coefs);
      x = 0.0;
    }
    else {
      double y = z/zmax;
      Usplines(level)(q,y, Coefs);
      x = (s-smin)/(smax-smin)*2.0 - 1.0;
    }

    LegendrePoly(x, Pn);
    // Now do summation
    double sum=0.0;
    for (int k=0; k<=Order; k++)
      sum += Pn(k)*Coefs(k);
    return (sum);
  }
  else {
    double beta = SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    double r = q+0.5*z;
    double rp = q-0.5*z;
    return (0.5*beta*(Potential->V(r)+Potential->V(rp)));
  }
}

double PAsFitClass::dU(double q, double z, double s2, int level)
{
  z = fabs(z);
  double qmax = qgrid->End*1.000001;
  double zmax = 1.000001*min (2.0*q, sMax(level));
  double smax = 1.000001*min (2.0*q, sMax(level));
  double smin = z;
  double s=sqrt(s2);
  double x;

  if ((q<=qmax)&&(z<=zmax)&&(s<=smax)) {
    if (q == 0) {
      dUsplines(level)(0.0,0.0,Coefs);
      x = 0.0;
    }
    else {
      double y = z/zmax;
      dUsplines(level)(q,y, Coefs);
      x = (s-smin)/(smax-smin)*2.0 - 1.0;
    }

    LegendrePoly(x, Pn);
    // Now do summation
    double sum=0.0;
    for (int k=0; k<=Order; k++)
      sum += Pn(k)*Coefs(k);
    return (sum);
  }
  else {
    double beta = SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    double r = q+0.5*z;
    double rp = q-0.5*z;
    return (0.5*(Potential->V(r)+Potential->V(rp)));
  }
}



bool PAsFitClass::Read (IOSectionClass &in,
			double smallestBeta, int numBetas)
{
  NumBetas = numBetas;
  SmallestBeta = smallestBeta;
  // Resize
  Usplines.resize(NumBetas);
  dUsplines.resize(NumBetas);
  Array<double,3> temp;

  // Read Particles;
  assert(in.OpenSection("Particle1"));
  Particle1.Read(in);
  in.CloseSection();
  assert(in.OpenSection("Particle2"));
  Particle2.Read(in);
  in.CloseSection();
  lambda = Particle1.lambda + Particle2.lambda;

  // Read Potential;
  assert(in.OpenSection("Potential"));
  Potential = ReadPH(in);
  in.CloseSection();

  // Read the fits
  assert(in.OpenSection("Fits"));
  // Read the qgrid
  assert(in.OpenSection("qGrid"));
  qgrid = ReadGrid (in);
  in.CloseSection();
  assert(in.OpenSection("yGrid"));
  ygrid = ReadGrid (in);
  in.CloseSection();
  GridsAreMine=true;
  // Read Order
  assert(in.ReadVar("Order", Order));
  Coefs.resize(Order+1);
  Pn.resize(Order+1);
  sMax.resize(NumBetas);

  double desiredBeta = smallestBeta;
  for (int betaIndex=0; betaIndex<NumBetas; betaIndex++) {
    int i=0;
    int numTemps = in.CountSections("Fit");
    bool found=false;
    while ((i<numTemps) && !found) {
      double beta;
      assert (in.OpenSection("Fit", i));
      assert (in.ReadVar("beta", beta));
      if ((fabs(beta-desiredBeta)/desiredBeta) < 1.0e-12)
	found = true;
      else {
	in.CloseSection();
	i++;
      }
    }
    if (!found) {
    cerr << "Couldn't find beta = " << desiredBeta 
	 << " in fit file.  Exitting\n";
    exit(1);
    }
    // Now read the fit coefficents
    double smax;
    assert(in.ReadVar("sMax", smax));
    sMax(betaIndex) = smax;
    assert(in.ReadVar("Umat", temp));
    Usplines(betaIndex).Init(qgrid,ygrid,temp);
    assert(in.ReadVar("dUmat", temp));
    dUsplines(betaIndex).Init(qgrid,ygrid,temp);
    in.CloseSection(); // "Fit"
    desiredBeta *= 2.0;
  }
  in.CloseSection(); // "Fits"
  return true;
}



