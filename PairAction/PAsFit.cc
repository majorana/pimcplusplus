#include "PAFit.h"
#include "../Splines/BicubicSpline.h"

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
  GridsAreMine = true;
  UsePBC = inSection.ReadVar ("Box", Box);
}

void PAsFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
{
  outSection.WriteVar ("Type", "sfit");
  outSection.NewSection("qGrid");
  qgrid->Write(outSection);
  outSection.CloseSection();
  outSection.NewSection("ygrid");
  ygrid->Write(outSection);
  outSection.CloseSection();
  outSection.WriteVar("Order", Order);
}


/// y = |z/z_max|
/// z_max = 
void PAsFitClass::AddFit (Rho &rho)
{
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
  for (int qi=0; qi<numq; qi++) {
    double q = (*qgrid)(qi);
    for (int yi=0; yi<numy; yi++) {
            
      double r, rp,costheta;
      double U, dU;
      rho.UdU(r,rp,costheta, U, dU);
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
  FILE *sdat = fopen ("s.dat", "w");
  FILE *costhetadat = fopen ("costheta.dat", "w");
  FILE *qdat = fopen ("q.dat", "w");
  LinearGrid qgrid2(qgrid->Start, qgrid->End, 315);
  for (int qi=0; qi<qgrid2.NumPoints; qi++) {
    double q = qgrid2(qi);
    LinearGrid sgrid(0.0, 2.0*q, 350);
    for (int si=0; si<sgrid.NumPoints; si++) {
      double s = sgrid(si);
      double w = exp(-s*s/(4.0*rho.lambda*rho.Beta()));
      double Uex, dUex, Ufit, dUfit;
      double costheta;
      if (q == 0.0)
	costheta = 1.0;
      else
	costheta = 1.0 - s*s/(2.0*q*q); 
      rho.UdU_Coulomb(q, q, costheta, Uex, dUex);
      Ufit = U(q, 0.0, s*s, level);
      dUfit = dU(q, 0.0, s*s, level);
      U2err += w*(Uex-Ufit)*(Uex-Ufit);
      dU2err += w*(dUex-dUfit)*(dUex-dUfit);
      weight += w;
      fprintf (Uxdat, "%1.16e ", Uex);
      fprintf (Ufdat, "%1.16e ", Ufit);
      fprintf (dUxdat, "%1.16e ", dUex);
      fprintf (dUfdat, "%1.16e ", dUfit);
      fprintf (sdat, "%1.16e ", s);
      fprintf (costhetadat, "%1.16e ", costheta);
      fprintf (qdat, "%1.16e ", q);
    }
    fprintf (Uxdat, "\n");
    fprintf (Ufdat, "\n");
    fprintf (dUxdat, "\n");
    fprintf (dUfdat, "\n");
    fprintf (sdat, "\n");
    fprintf (qdat, "\n");
    fprintf (costhetadat, "\n");
  }
  fclose (Uxdat); fclose(Ufdat); fclose(sdat); fclose(qdat); 
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
    outSection.NewSection("Fit");
    outSection.WriteVar ("beta", beta);
    for (int qi=0; qi<qgrid->NumPoints; qi++)
      for (int yi=0; yi<ygrid->NumPoints; yi++) 
	for (int j=0; j<=Order; j++) {
	  Umat(qi,yi,j) = Usplines(bi)(qi,yi,j);
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
  if (q <= (qgrid->End*1.0000001)) {
    if (q == 0) 
      Usplines(level)(0.0,0.0,Coefs);
    else {
      double t = sqrt(s2)/(2.0*q);
      Usplines(level)(q,t, Coefs);
    }
    // Now do summation
  }
  else {
    // Coulomb action is independent of z
    double beta = SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    return (beta*Potential->V(q));
  }
}

double PAsFitClass::dU(double q, double z, double s2, int level)
{
  if (q <= (qgrid->End*1.0000001)) {
    if (q == 0)
      dUsplines(level)(0.0,0.0,Coefs);
    else {
      double t = sqrt(s2)/(2.0*q);
      dUsplines(level)(q,t,Coefs);
    }
    // Now do summation
  }
  else {
    // Coulomb action is independent of z
    return (Potential->V(q));
  }
}




bool PAsFitClass::Read (IOSectionClass &in,
				double smallestBeta, int NumBetas)
{
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



