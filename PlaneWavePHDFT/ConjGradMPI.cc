#include "ConjGrad2.h"
#include "../MatrixOps/MatrixOps.h"
#include "../MPI/Communication.h"

void ConjGradMPI::Setup()
{
  int numBands = Bands.rows();
  int N = H.GVecs.size();

  // Figure out which bands I'm responsible for
  int numProcs = Communicator.NumProcs();
  int myProc = Communicator.MyProc();
  int band = 0;
  for (int proc=0; proc<numProcs; proc++) {
    int procBands = numBands/numProcs + ((numBands % numProcs)>proc);
    if (proc==myProc) {
      MyFirstBand = band;
      MyLastBand = band + procBands-1;
      band += procBands;
    }
  }

  // Now allocate memory for everything
  Bands.resize (numBands, N);
  Energies.resize(numBands);
  cnext.resize (N);
  Hc.resize (N);
  Phi.resize (N);
  Phip.resize(N);
  Phipp.resize(N);
  Xi.resize(N);
  Eta.resize(N);
  Bands = 0.0;
  Vec3 kBox = H.GVecs.GetkBox();
  double maxk = max(kBox[0], max(kBox[1], kBox[2]));
  double maxE = 1.0*maxk*maxk;
  int numNonZero;
  for (int band=0; band<numBands; band++) {
    numNonZero = 0;
    c.reference(Bands(band,Range::all()));
    c = 0.0;
    for (int i=0; i<N; i++) {
      double dist2 = dot(H.GVecs(i), H.GVecs(i));
      c(i) = exp(-2.0*dist2);
    }
    Normalize(c);
  }
  InitBands();
  IsSetup = true;
}

void ConjGradMPI::InitBands()
{
  int numBands = Bands.rows();
  //  int numVecs = 6 * numBands;
  int numVecs = 4 * numBands;
  assert (numVecs <= H.GVecs.size());
  Array<complex<double>,2> Hmat(numVecs, numVecs);
  Array<complex<double>,2> EigVecs(numBands, numVecs);
  Array<double,1> EigVals (numBands);

  H.Vion->Vmatrix(Hmat);
  cerr << "kPoint = " << H.kPoint << endl;
  for (int i=0; i<numVecs; i++) {
    Vec3 Gpk = H.GVecs(i) + H.kPoint;
    Hmat(i,i) += 0.5 * dot (Gpk,Gpk);
  }
  SymmEigenPairs (Hmat, numBands, EigVals, EigVecs);
  
  for (int i=0; i<numBands; i++)
    cerr << "Mini energy(" << i << ") = " << 27.211383*EigVals(i) << endl;

  // Now put results in Bands
  Bands = 0.0;
  for (int band=0; band<numBands; band++) 
    for (int i=0; i<numVecs; i++)
      Bands(band, i) = EigVecs(band, i);
}

double ConjGradMPI::CalcPhiSD()
{
  H.Apply(c, Hc);
  E0 = realconjdot (c, Hc);
  cerr << "E = " << E0 << endl;
  Phip = E0*c - Hc;
  double nrm = norm (Phip);
  Normalize (Phip);  
  return nrm;
}

void ConjGradMPI::Precondition()
{
  double Tinv = 1.0/T;
  for (int i=0; i<c.size(); i++) {
    double x = 0.5*dot(H.GVecs(i), H.GVecs(i))*Tinv;
    double num = 27.0 + 18.0*x +12.0*x*x + 8.0*x*x*x;
    double denom = num + 16.0*x*x*x*x;
    Eta(i) = (num/denom)* Xi(i);
  }
}


// Returns the norm of the residual Hc - Ec
double ConjGradMPI::CalcPhiCG()
{
  Hc = 0.0;
  H.Kinetic.Apply (c, Hc);
  T = realconjdot (c, Hc);
  H.Vion->Apply (c, Hc);
  E0 = realconjdot (c, Hc);
  Xi = E0*c - Hc;
  double residualNorm = norm (Xi);
  /// Orthonalize to other bands here
  zVec &Xip = Xi;
  Orthogonalize2 (Bands, Xip, CurrentBand);

  Precondition();

  // Now, orthogonalize to all bands
  // rename for clarity
  zVec &Etap = Eta;
  Orthogonalize2 (Bands, Etap, CurrentBand);
  // Etap = Eta - conjdot (c, Eta)*c;

  // Compute conjugate direction
  complex<double> etaxi = conjdot(Etap, Xip);
  complex<double> gamma; 
  if (EtaXiLast != complex<double>(0.0, 0.0)) 
    gamma = etaxi/EtaXiLast;
  else
    gamma = 0.0;
  EtaXiLast = etaxi;
  
  Phi = Etap + gamma * Phi;
  
  // Orthogonalize to present band
  complex<double> cPhi = conjdot(c,Phi);
  Phip = Phi - cPhi * c;
  Normalize (Phip);
  Energies(CurrentBand) = E0;
  return residualNorm;
}


void ConjGradMPI::Solve(int band)
{
  CurrentBand = band;
  if (!IsSetup)
    Setup();
  c.reference (Bands(band,Range::all()));
  EtaXiLast = 0.0;
  Orthogonalize2 (Bands, c, band-1);
  Normalize(c);

  double Elast = 1.0e100;
  double residualNorm = 1.0;
  //  while (fabs (Elast - Energies(band)) > Tolerance) {
  int iter=0;
  while ((residualNorm > 1.0e-8) && (iter < 50)) {
    Elast = Energies(band);
    // First, calculate conjugate gradient direction
    residualNorm = CalcPhiCG();
    
    // Now, pick optimal theta for 
    double dE_dtheta = 2.0*realconjdot(Phip, Hc);

    H.Apply (Phip, Hc);
    double d2E_dtheta2 = 2.0*(realconjdot(Phip, Hc) - E0);
    double thetaMin = 0.5*atan(-dE_dtheta/(0.5*d2E_dtheta2));

    double costhetaMin, sinthetaMin;
    sincos(thetaMin, &sinthetaMin, &costhetaMin);
    c = costhetaMin*c + sinthetaMin*Phip;
    iter++;
  }
  if (residualNorm > 1.0e-8)
    cerr << "Warning:  conjugate gradient residual norm = " 
	 << residualNorm << endl;
  cerr << "# of iterations = " << iter << endl;
}

void
ConjGradMPI::PrintOverlaps()
{
  cerr << "Overlaps = \n";
  zVec x, y;
  for (int i=0; i<Bands.rows(); i++) {
    x.reference (Bands(i, Range::all()));
    for (int j=0; j<Bands.rows(); j++) {
      y.reference (Bands(j,Range::all()));
      double s = realconjdot(x,y);
      fprintf (stderr, "%12.6e ", s);
    }
    fprintf (stderr, "\n");
  }
}


void
ConjGradMPI::CollectBands()
{
  Communicator.AllGatherRows(Bands);
}
