#include "PAFit.h"
#include "../SpecialFunctions/HermitePoly.h"
#include "../Fitting/Fitting.h"

/// The following routines are used only if we are creating fits, not
/// using them.
#ifdef MAKE_FIT
void PAcoulombFitClass::ReadParams(IOSectionClass &inSection)
{
  assert(inSection.OpenSection ("qGrid"));
  qgrid = ReadGrid (inSection);
  GridIsMine = true;
  assert(inSection.ReadVar ("Order", Order));
  Ucoefs.resize(Order+1);
  UsePBC = inSection.ReadVar ("Box", Box);
}

void PAcoulombFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
{
  outSection.WriteVar ("Type", "coulombfit");
  outSection.WriteVar ("Order", Order);
  outSection.NewSection("qGrid");
  qgrid->Write(outSection);
  outSection.CloseSection();
}



void PAcoulombFitClass::AddFit (Rho &rho)
{
  NumBetas++;
  Uj.resizeAndPreserve(NumBetas);
  dUj.resizeAndPreserve(NumBetas);

  FILE *Udebug = fopen ("U.dat", "w");
  FILE *sdebug = fopen ("s.dat", "w");

  double lambda = rho.lambda;
  double beta = rho.Beta();
  const int N = 400;
  int M = Order+1;
  Array<double,2>  UCoefs(qgrid->NumPoints, M);
  Array<double,2> dUCoefs(qgrid->NumPoints, M);

  Array<double,1> Uexact(N), dUexact(N), sigma(N);
  Array<double,2> Basis(N,M-1);
  Array<double,1> Hn(2*M), errors(M), Ui(M), dUi(M);
  for (int qi=0; qi<qgrid->NumPoints; qi++) {
    double q = (*qgrid)(qi);
    LinearGrid sgrid(0.0, 2.0*q, N);
    double U0, dU0;
    rho.UdU_Coulomb (q, q, 1.0, U0, dU0);
    for (int si=0; si<sgrid.NumPoints; si++) {
      double s = sgrid(si);
      double costheta;
      if (q == 0.0)
	costheta = 1.0;
      else
	costheta = 1.0 - s*s/(2.0*q*q);
      costheta = min(1.0, costheta);
      costheta = max(-1.0, costheta);
      // Compute exact U and dU
      double Uex, dUex;
      rho.UdU_Coulomb (q, q, costheta, Uex, dUex);
      Uexact(si) = Uex-U0;
      fprintf (Udebug, "%1.16e ", Uex-U0);
      fprintf (sdebug, "%1.16e ", s);
      fflush (Udebug);
      fflush (sdebug);
      dUexact(si) = dUex-dU0;
      if (isnan(Uexact(si))) {
	cerr << "NAN is Uexact for q = " << q << " s = " << s << endl;
	abort();
      }
      if (isnan(dUexact(si))){
	cerr << "NAN in dUexact for q = " << q << "s = " << s << endl;
	abort();
      }
      // Compute weight for point
      double rho = exp(-s*s/(4.0*lambda*beta));
      if (rho > 1e-8)
	sigma(si) = 1.0/sqrt(rho);
      else
	sigma(si) = 1.0/sqrt(1e-8);
      // Compute basis functions
      double t = s / sqrt(4.0*lambda*beta);
      for (int i=1; i<M; i++) 
	Basis(si,i-1) = pow(s, 2*i);
    }
    // Now do fits
    LinFitSVD( Uexact, sigma, Basis,  Ui, errors, 1.0e-12);
    LinFitSVD(dUexact, sigma, Basis, dUi, errors, 1.0e-12);
    UCoefs(qi,0) = U0;
    dUCoefs(qi,0) = dU0;
    UCoefs(qi,Range(1,M-1)) = Ui;
    dUCoefs(qi,Range(1,M-1)) = dUi;
    fprintf (Udebug, "\n");
    fflush (Udebug);
    fprintf (sdebug, "\n");
    fflush (sdebug);
  } 
    
  fclose (Udebug);
  fclose (sdebug);

  // Initialize splines
  Uj(NumBetas-1).Init(qgrid, UCoefs);
  dUj(NumBetas-1).Init(qgrid, dUCoefs);  
}

void PAcoulombFitClass::Error(Rho &rho, double &Uerror, double &dUerror)
{
  int level = (int)floor(log(rho.Beta()/SmallestBeta)/log(2.0)+ 0.5);

  double U2err = 0.0;
  double dU2err = 0.0;
  double weight;
  FILE *Uxdat = fopen ("Ux.dat", "w");
  FILE *Ufdat = fopen ("Uf.dat", "w");
  FILE *sdat = fopen ("s.dat", "w");
  FILE *qdat = fopen ("q.dat", "w");
  for (int qi=0; qi<qgrid->NumPoints; qi++) {
    double q = (*qgrid)(qi);
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
      fprintf (sdat, "%1.16e ", s);
      fprintf (qdat, "%1.16e ", q);
    }
    fprintf (Uxdat, "\n");
    fprintf (Ufdat, "\n");
    fprintf (sdat, "\n");
    fprintf (qdat, "\n");
  }
  fclose (Uxdat); fclose(Ufdat); fclose(sdat); fclose(qdat);
  Uerror = sqrt(U2err/weight);
  dUerror = sqrt(dU2err/weight);
}


void PAcoulombFitClass::WriteFits (IOSectionClass &outSection)
{
  Array<double,2> UCoefs(qgrid->NumPoints, Order+1); 
  Array<double,2> dUCoefs(qgrid->NumPoints, Order+1); 
  double beta = SmallestBeta;
  for (int i=0; i<NumBetas; i++) {
    outSection.NewSection("Fit");
    outSection.WriteVar ("beta", beta);
    for (int j=0; j<qgrid->NumPoints; j++)
      for (int k=0; k<(Order+1); k++) {
	UCoefs(j,k) = Uj(i)(j,k);
	dUCoefs(j,k) = dUj(i)(j,k);
      }
    outSection.WriteVar ("Ucoefs", UCoefs);
    outSection.WriteVar ("dUcoefs", dUCoefs);
    outSection.CloseSection();
    beta *= 2.0;
  }
}

#endif

double PAcoulombFitClass::U(double q, double z, double s2, int level)
{
  if (q < qgrid->End) {
    Uj(level)(q, Ucoefs);
    double s2j = 1.0;
    double Usum = 0.0;
    for (int j=0; j<=Order; j++){
      Usum += Ucoefs(j)*s2j;
      s2j *= s2;
    }
    return (Usum);
  }
  else {
    double beta = SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    // Coulomb action is independent of z
    return (beta*Potential->V(q));
  }
}

double PAcoulombFitClass::dU(double q, double z, double s2, int level)
{
  if (q < qgrid->End) {
    dUj(level)(q, Ucoefs);
    double s2j = 1.0;
    double dUsum = 0.0;
    for (int j=0; j<=Order; j++){
      dUsum += Ucoefs(j)*s2j;
      s2j *= s2;
    }
    return (dUsum);
  }
  else
    return 0.0;
}




bool PAcoulombFitClass::Read (IOSectionClass &in,
			      double smallestBeta, int NumBetas)
{
  // Resize
  Uj.resize(NumBetas);
  dUj.resize(NumBetas);
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
  Potential = ReadPH(in);
  in.CloseSection();

  // Read the fits
  assert(in.OpenSection("Fits"));
  // Read the qgrid
  assert(in.OpenSection("qGrid"));
  qgrid = ReadGrid (in);
  in.CloseSection();
  GridIsMine=true;
  // Read Order
  assert (in.ReadVar ("Order", Order));
  Ucoefs.resize(Order+1);
  dUcoefs.resize(Order+1);

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
    assert(in.ReadVar("Ucoefs", temp));
    Uj(betaIndex).Init(qgrid,temp);
    assert(in.ReadVar("dUcoefs", temp));
    dUj(betaIndex).Init(qgrid,temp);
    in.CloseSection(); // "Fit"
    desiredBeta *= 2.0;
  }
  in.CloseSection(); // "Fits"
  return true;
}



// void Rho::WriteFit (string fileName)
// {
//   FILE *fout = fopen (fileName.c_str(), "w");
//   FILE *fpout = fopen ("var.dat", "w");
//   assert (fout != NULL);
//   int lmax = U_ls.size()-1;
//   LinearGrid qgrid (1e-4, 7.5, 101);
  
//   // HACK
//   double myBeta = 0.125;
//   //double myBeta = Beta();

//   for (int qindex=0; qindex<qgrid.NumPoints; qindex++)
//     {
//       cerr << " q = " << qindex << endl;
//       double q = qgrid(qindex);
//       double smax = min (2.0*q, sqrt(-4.0*lambda*myBeta*log(1.0e-4)));
//       LinearGrid zgrid(-smax, smax, 101);
//       LinearGrid sgrid(0.0, smax, 101);
//       for (int zindex=0; zindex<zgrid.NumPoints; zindex++)
// 	{
// 	  double z = zgrid(zindex);
// 	  double r = q + 0.5*z;
// 	  double rp = q - 0.5*z;
// 	  if (r <= 0.0) r = 1.0e-5;
// 	  if (rp <= 0.0) rp = 1.0e-5;
// 	  double x = Transform.r2x(r);
// 	  double xp = Transform.r2x(rp);
// 	  Array<double,1> Ulvec(lmax+1), dUlvec(lmax+1);
// 	  // First, calculate the diagonal parts
// 	  U_lArray(r,r, Ulvec, dUlvec);
// 	  double Ur, Urp, dUr, dUrp;
// 	  UdU(r,r,1.0,Ulvec, dUlvec, Ur, dUr);
// 	  U_lArray(rp,rp, Ulvec, dUlvec);
// 	  UdU(r,r,1.0,Ulvec, dUlvec, Ur, dUr);
// // 	  for (int l=0; l<=lmax; l++)
// // 	    Ularray(l) = U_ls(l).U(x,x);
// // 	  double Udiag_r = U(r,r,1.0, Ularray);

// // 	  for (int l=0; l<=lmax; l++)
// // 	    Ularray(l) = U_ls(l).U(xp,xp);
// // 	  double Udiag_rp = U(rp,rp,1.0, Ularray);

// 	  // Now calculate off-diagonal elements
// // 	  for (int l=0; l<=lmax; l++)
// // 	    Ularray(l) = U_ls(l).U(x,xp);
// 	  U_lArray(r,rp, Ulvec, dUlvec);
// 	  //LinearGrid sgrid(fabs(r-rp)+1e-9, r+rp-1e-9, 101);
// 	  for (int sindex=0; sindex<sgrid.NumPoints; sindex++)
// 	    {
// 	      double s = sgrid(sindex);	      
// 	      double costheta;
// 	      if ((r<=0.0) || (rp<=0.0))
// 		costheta = 1.0;
// 	      else
// 		costheta = (r*r + rp*rp - s*s)/(2.0*r*rp);
// 	      //if (costheta > 1.0)  costheta = 1.0;
// 	      //if (costheta < -1.0) costheta = -1.0;
// 	      double theta = acos (costheta);
// 	      double Uval, dUval, Ufitval, dUfitval, Ufit2val, dUfit2val;
// 	      if (isnan(theta))
// 		{
// // 		  cerr << "NAN in theta.\n"
// // 		       << "costheta = " << costheta << endl
// // 		       << "r = " << r << " rp = " << rp << endl 
// // 		       << "q = " << q << " z = " << z << " s = " << s << endl;
// 		  Uval = NAN;	   dUval = NAN;
// 		  Ufitval = NAN;   dUfitval = NAN;
// 		  Ufit2val = NAN;  dUfit2val = NAN;
// 		}
// 	      else
// 		{
// 		  if ((r < grid->End) && (rp < grid->End)) {
// 		    UdU(r,rp,costheta, Ulvec, dUlvec, Uval, dUval);
// 		    //Uval = U(r,rp,costheta, Ularray);
// 		    //Uval -= 0.5*(Udiag_r + Udiag_rp);
// 		    Ufitval = Ufit(q, s, z);
// 		    //Ufitval -= 0.5*(Udiag_r + Udiag_rp);
// 		    Ufit2val = Ufit2(q,s,z);
// 		    dUfit2val = Ufit2.dU(q,s,z);
// 		    //Ufit2val -= 0.5*(Udiag_r + Udiag_rp);
// 		  }
// 		  else {
// 		    cerr << "r or rp outside grid.\n";
// 		    Uval = 0.0;
// 		  }
// 		}
// 	      fprintf (fout, "%1.14e %1.14e %1.14e %1.14e %1.14e %1.14e\n", 
// 		       Uval, Ufitval, Ufit2val, dUval, dUfitval, dUfit2val);
// 	      fprintf (fpout, "%1.6e %1.6e %1.6e %1.8e\n", q, z, s,
// 		       exp(-s*s/(2.0*beta)));
// 	    }
// 	}
//     }
//   fclose (fout);
//   fclose(fpout);
// }
