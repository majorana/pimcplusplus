#include "PAFit.h"

const double Rho0Min = 1.0e-3;

/// The following routines are used only if we are creating fits, not
/// using them.
#ifdef MAKE_FIT
void PAtricubicFitClass::ReadParams(IOSectionClass &inSection)
{
  assert(inSection.OpenSection ("qGrid"));
  qgrid = ReadGrid (inSection);
  inSection.CloseSection();
  assert(inSection.OpenSection ("yGrid"));
  ygrid = ReadGrid (inSection);
  inSection.CloseSection();
  assert(inSection.OpenSection ("tGrid"));
  tgrid = ReadGrid (inSection);
  inSection.CloseSection();
  assert(fabs(ygrid->End-1.0) < 1.0e-10);
  assert(fabs(tgrid->End-1.0) < 1.0e-10);
  GridsAreMine = true;
  UsePBC = inSection.ReadVar ("Box", Box);
}

void PAtricubicFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
{
  if (Comm.MyProc() == 0) {
    outSection.WriteVar ("Type", "tricubicfit");
    outSection.NewSection("qGrid");
    qgrid->Write(outSection);
    outSection.CloseSection();
    outSection.NewSection("yGrid");
    ygrid->Write(outSection);
    outSection.CloseSection();
    outSection.NewSection("tGrid");
    tgrid->Write(outSection);
    outSection.CloseSection();
  }
}


/// y = |z/z_max|
/// z_max = 
void PAtricubicFitClass::AddFit (Rho &rho)
{
  NumBetas++;
  Usplines.resizeAndPreserve(NumBetas);
  dUsplines.resizeAndPreserve(NumBetas);
  sMax.resizeAndPreserve(NumBetas);
  sMaxInv.resizeAndPreserve(NumBetas);

  double lambda = Particle1.lambda + Particle2.lambda;
  double beta = rho.Beta();
  sMax(NumBetas-1) = sqrt(-4.0*lambda*beta*log(Rho0Min));
  sMaxInv(NumBetas-1) = 1.0/sMax(NumBetas-1);

  int numq = qgrid->NumPoints;
  int numy = ygrid->NumPoints;
  int numt = tgrid->NumPoints;
  Array<double,3> Umat(numq,numy,numt), dUmat(numq,numy,numt);
  Array<double,1> Ul, dUl;


  /*  double q=4.0;
  for (double z=0.0; z<1.0; z+=0.0005) {
    double r = q+0.5*z;
    double rp = q-0.5*z;
    double x = rho.Transform.r2x(r);
    double xp = rho.Transform.r2x(rp);
    rho.U_lArray(r,rp,Ul,dUl);
    double U, dU;
    rho.UdU(r,rp,1.0, Ul, dUl, U, dU);
    fprintf (stderr, "%1.5e %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e\n", 
	     z, U, r, rp, x, xp,Ul(0),Ul(1));
	     }*/

  for (int qi=0; qi<numq; qi++) {
    double q = (*qgrid)(qi);
    double smax = min(2.0*q, sMax(NumBetas-1));
    for (int yi=0; yi<numy; yi++) {
      double y = (*ygrid)(yi);
      double z = y*smax;
      double r = q+0.5*z; 
      double rp= q-0.5*z;
      if (r == 0.0)
	r = 1.0e-8;
      if (rp == 0.0)
	rp = 1.0e-8;
      
      rho.U_lArray(r,rp,Ul,dUl);
      for (int ti=0; ti<numt; ti++) {
	double t = (*tgrid)(ti);
	double s = z + (smax-z)*t;
	double costheta;
	if ((r*rp)==0.0)
	  costheta = 1.0;
	else
	  costheta = (r*r + rp*rp - s*s)/(2.0*r*rp);
	costheta = min(1.0, costheta);
	costheta = max(-1.0, costheta);	
	double U, dU;
	rho.UdU(r,rp,costheta, Ul, dUl, U, dU);
	if (isnan(U))
	  fprintf (stderr, "NAN in U at (qi,yi,ti) = (%d,%d,%d)\n",
		   qi, yi, ti);
	if (isnan(dU))
	  fprintf (stderr, "NAN in dU at (qi,yi,ti) = (%d,%d,%d)\n",
		   qi, yi, ti);
	Umat(qi,yi,ti) = U;
	dUmat(qi,yi,ti) = dU;
      }
    }
  } 
  Usplines(NumBetas-1).Init(qgrid,ygrid,tgrid,Umat);
  dUsplines(NumBetas-1).Init(qgrid,ygrid,tgrid,dUmat);
}



void PAtricubicFitClass::Error(Rho &rho, double &Uerror, double &dUerror)
{
  int level = (int)floor(log(rho.Beta()/SmallestBeta)/log(2.0)+ 0.5);

  double U2err = 0.0;
  double dU2err = 0.0;
  double weight = 0.0;
  FILE *Uxdat = fopen ("Ux.dat", "w");
  FILE *Ufdat = fopen ("Uf.dat", "w");
  FILE *dUxdat = fopen ("dUx.dat", "w");
  FILE *dUfdat = fopen ("dUf.dat", "w");
  FILE *tdat = fopen ("t.dat", "w");
  FILE *costhetadat = fopen ("costheta.dat", "w");
  FILE *ydat = fopen ("y.dat", "w");
  LinearGrid qgrid2(qgrid->Start, 0.999*qgrid->End, 20);
  Array<double,1> Ul, dUl;
  for (int qi=0; qi<qgrid2.NumPoints; qi++) {
    double q = qgrid2(qi);
    // HACK
    //q = 1.0;
    double zmax = 0.9999999*min(2.0*q,sMax(level));
    LinearGrid zgrid(0.0, zmax, 20);
    for (int zi=0; zi<zgrid.NumPoints; zi++) {
      double z = zgrid(zi);
      double y = z/zmax;
      double smax = 0.9999999*min(2.0*q,sMax(level));
      //cerr << "smin = "  << z << " smax = " << smax << endl;
      LinearGrid sgrid(z, smax, 100);
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
	
	//cerr << "costheta = " << costheta << endl;

	costheta = min(costheta,1.0);
	costheta = max(costheta,-1.0);

	rho.U_lArray(r,rp,Ul, dUl);
	rho.UdU(r, rp, costheta, Ul, dUl, Uex, dUex);
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


void PAtricubicFitClass::WriteFits (IOSectionClass &outSection)
{
  if (Comm.MyProc() == 0) {
    Array<double,3> Umat(qgrid->NumPoints,ygrid->NumPoints,tgrid->NumPoints); 
    Array<double,3> dUmat(qgrid->NumPoints,ygrid->NumPoints,tgrid->NumPoints); 
    double beta = SmallestBeta;
    for (int bi=0; bi<NumBetas; bi++) {
      cerr << "Writing Beta "<< bi+1 << " of " << NumBetas << ":\n";
      outSection.NewSection("Fit");
      outSection.WriteVar ("beta", beta);
      double smax = sMax(bi);
      outSection.WriteVar ("sMax", smax);
      for (int qi=0; qi<qgrid->NumPoints; qi++)
	for (int yi=0; yi<ygrid->NumPoints; yi++) 
	  for (int ti=0; ti<tgrid->NumPoints; ti++) {
	    Umat(qi,yi,ti)  =  Usplines(bi)(qi,yi,ti);
	    dUmat(qi,yi,ti) = dUsplines(bi)(qi,yi,ti);
	  }
      outSection.WriteVar ("Umat", Umat);
      outSection.WriteVar ("dUmat", dUmat);
      outSection.CloseSection();
      beta *= 2.0;
    }
  }
}

#endif


double PAtricubicFitClass::U(double q, double z, double s2, int level)
{
  z = fabs(z);
  double qmax = qgrid->End;
  double smax = min (2.0*q, sMax(level));
  double zmax = smax;
  double smin = z;
  double s=sqrt(s2);
  double x;

  if ((q<=qmax)&&(z<=zmax)&&(s<=smax)) {
    if (q == 0.0) 
      return(Usplines(level)(0.0,0.0,0.0));
    else {
      double y = z/zmax;
      double t = (s-z)/(smax-z);
      return(Usplines(level)(q,y,t));
    }
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



double PAtricubicFitClass::dU(double q, double z, double s2, int level)
{
  z = fabs(z);
  double qmax = qgrid->End;
  double smax = min (2.0*q, sMax(level));
  double zmax = smax;
  double smin = z;
  double s=sqrt(s2);
  double x;

  if ((q<=qmax)&&(z<=zmax)&&(s<=smax)) {
    if (q == 0.0) 
      return(dUsplines(level)(0.0,0.0,0.0));
    else {
      double y = z/zmax;
      double t = (s-z)/(smax-z);
      return(dUsplines(level)(q,y,t));
    }
  }
  else {
    double r = q+0.5*z;
    double rp = q-0.5*z;
    return (0.5*(Potential->V(r)+Potential->V(rp)));
  }
}





bool PAtricubicFitClass::Read (IOSectionClass &in,
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
  assert(in.OpenSection("tGrid"));
  tgrid = ReadGrid (in);
  in.CloseSection();
  GridsAreMine=true;
  sMax.resize(NumBetas);
  sMaxInv.resize(NumBetas);

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
    sMaxInv(betaIndex) = 1.0/smax;
    assert(in.ReadVar("Umat", temp));
    Usplines(betaIndex).Init(qgrid,ygrid,tgrid,temp);
    assert(in.ReadVar("dUmat", temp));
    dUsplines(betaIndex).Init(qgrid,ygrid,tgrid,temp);
    in.CloseSection(); // "Fit"
    desiredBeta *= 2.0;
  }
  in.CloseSection(); // "Fits"
  return true;
}



