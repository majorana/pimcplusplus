#include "PAFit.h"
#include "../Splines/BicubicSpline.h"

/// The following routines are used only if we are creating fits, not
/// using them.
#ifdef MAKE_FIT
void PAcoulombBCFitClass::ReadParams(IOSectionClass &inSection)
{
  assert(inSection.OpenSection ("qGrid"));
  qgrid = ReadGrid (inSection);
  inSection.CloseSection();
  assert(inSection.OpenSection ("tGrid"));
  tgrid = ReadGrid (inSection);
  inSection.CloseSection();
  GridsAreMine = true;
  UsePBC = inSection.ReadVar ("Box", Box);
}

void PAcoulombBCFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
{
  outSection.WriteVar ("Type", "coulombBCfit");
  outSection.NewSection("qGrid");
  qgrid->Write(outSection);
  outSection.CloseSection();
  outSection.NewSection("tGrid");
  tgrid->Write(outSection);
  outSection.CloseSection();
}



void PAcoulombBCFitClass::AddFit (Rho &rho)
{
  NumBetas++;
  Usplines.resizeAndPreserve(NumBetas);
  dUsplines.resizeAndPreserve(NumBetas);

  int numq = qgrid->NumPoints;
  int numt = tgrid->NumPoints;
  Array<double,2> Umat(numq, numt), dUmat(numq, numt);
  for (int qi=0; qi<numq; qi++) {
    double q = (*qgrid)(qi);
    for (int ti=0; ti<numt; ti++) {
      double t = (*tgrid)(ti);
      double s = 2.0*q*t;

      double costheta;
      if (q == 0.0)
	costheta = 1.0;
      else
	costheta = 1.0 - s*s/(2.0*q*q);
      costheta = min(1.0, costheta);
      costheta = max(-1.0, costheta);
      
      double U, dU;
      rho.UdU_Coulomb(q,q,costheta, U, dU);
      Umat(qi,ti) = U;
      dUmat(qi,ti) = dU;
    }
  }
  // Initialize the bicubic splines
  Usplines(NumBetas-1).Init(qgrid,tgrid,Umat);
  dUsplines(NumBetas-1).Init(qgrid,tgrid,dUmat);
}



void PAcoulombBCFitClass::Error(Rho &rho, double &Uerror, double &dUerror)
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


void PAcoulombBCFitClass::WriteFits (IOSectionClass &outSection)
{
  Array<double,2> Umat(qgrid->NumPoints, tgrid->NumPoints); 
  Array<double,2> dUmat(qgrid->NumPoints, tgrid->NumPoints); 
  double beta = SmallestBeta;
  for (int bi=0; bi<NumBetas; bi++) {
    outSection.NewSection("Fit");
    outSection.WriteVar ("beta", beta);
    for (int qi=0; qi<qgrid->NumPoints; qi++)
      for (int ti=0; ti<tgrid->NumPoints; ti++) {
	Umat(qi,ti) = Usplines(bi)(qi,ti);
	dUmat(qi,ti) = dUsplines(bi)(qi,ti);
      }
    outSection.WriteVar ("Umat", Umat);
    outSection.WriteVar ("dUmat", dUmat);
    outSection.CloseSection();
    beta *= 2.0;
  }
}

#endif


double PAcoulombBCFitClass::U(double q, double z, double s2, int level)
{
  if (q <= (qgrid->End*1.0000001)) {
    if (q == 0)
      return (Usplines(level)(0.0,0.0));
    else {
      double t = sqrt(s2)/(2.0*q);
      return (Usplines(level)(q,t));
    }
  }
  else {
    // Coulomb action is independent of z
    double beta = SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    return (beta*Pot->V(q));
  }
}

double PAcoulombBCFitClass::dU(double q, double z, double s2, int level)
{
  if (q <= (qgrid->End*1.0000001)) {
    if (q == 0)
      return (dUsplines(level)(0.0,0.0));
    else {
      double t = sqrt(s2)/(2.0*q);
      return (dUsplines(level)(q,t));
    }
  }
  else {
    // Coulomb action is independent of z
    return (Pot->V(q));
  }
}




bool PAcoulombBCFitClass::Read (IOSectionClass &in,
				double smallestBeta, int NumBetas)
{
  SmallestBeta = smallestBeta;
  // Resize
  Usplines.resize(NumBetas);
  dUsplines.resize(NumBetas);
  Array<double,2> temp;

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
  Pot = ReadPotential(in);
  in.CloseSection();

  // Read the fits
  assert(in.OpenSection("Fits"));
  // Read the qgrid
  assert(in.OpenSection("qGrid"));
  qgrid = ReadGrid (in);
  in.CloseSection();
  assert(in.OpenSection("tGrid"));
  tgrid = ReadGrid (in);
  in.CloseSection();
  GridsAreMine=true;
  // Read Order

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
    Usplines(betaIndex).Init(qgrid,tgrid,temp);
    assert(in.ReadVar("dUmat", temp));
    dUsplines(betaIndex).Init(qgrid,tgrid,temp);
    in.CloseSection(); // "Fit"
    desiredBeta *= 2.0;
  }
  in.CloseSection(); // "Fits"
  return true;
}


bool PAcoulombBCFitClass::IsLongRange()
{
  return true;
}


double PAcoulombBCFitClass::Vlong(double q, int level)
{
  if (q <= 0.0)
    return 2.0/sqrt(M_PI)*Z1Z2*alpha;
  else 
    return Z1Z2/q*erf(alpha*q);
}

double PAcoulombBCFitClass::Vlong_k(double boxVol, double k, int level)
{
  if (k <= 0.0)
    k = 1.0e-30;
  return 4.0*M_PI/(boxVol*k*k)*exp(-k*k/(4.0*alpha*alpha));
}

double PAcoulombBCFitClass::dVlong(double q, int level)
{
  return 0.0;
}


void PAcoulombBCFitClass::DoBreakup(const dVec &box, const Array<dVec,1> &kVecs)
{
  // Calculate the cutoff parameter
  double minL, boxVol;
  boxVol = minL = box[0];
  for (int i=1; i<NDIM; i++) {
    minL = min(minL, box[i]);
    boxVol *= box[i];
  }
  alpha = 7.0/minL;

  // First, subtract off long-ranged part from the bicubic spline to
  // make them short-ranged.
  for (int level=0; level<NumBetas; level++) {
    double tau = SmallestBeta * pow(2.0, level);
    for (int qi=0; qi<Usplines(level).Nx; qi++) {
      double q = (*Usplines(level).Xgrid)(qi);
      double v = Vlong(q, level);
      for (int ti=0; ti<Usplines(level).Ny; ti++) {
	Usplines(level)(qi,ti) -= tau*v;
	dUsplines(level)(qi,ti) -= v;
      }
    }
  }


  // Now, calculate the k-space parts
  Ulong_k.resize(NumBetas, kVecs.size());
  dUlong_k.resize(NumBetas, kVecs.size());
  U_RPA_long_k.resize(NumBetas, kVecs.size());
  dU_RPA_long_k.resize(NumBetas, kVecs.size());
  for (int level=0; level<NumBetas; level++) {
    double tau = SmallestBeta * pow(2.0, level);
    Ulong_0(level)  = tau*Vlong(0.0, level);
    dUlong_0(level) = Vlong(0.0, level);
    for (int ki=0; ki<kVecs.size(); ki++) {
      double k = sqrt(dot(kVecs(ki), kVecs(ki)));
      Ulong_k = tau*Vlong_k(boxVol, k, level);
      dUlong_k = Vlong_k(boxVol, k, level);
    }
  }
}

