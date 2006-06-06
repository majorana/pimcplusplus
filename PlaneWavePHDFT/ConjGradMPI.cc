#include "ConjGradMPI.h"
#include "../MatrixOps/MatrixOps.h"
#include "../MPI/Communication.h"

void ConjGradMPI::Setup()
{
  int numBands = Bands.rows();
  int N      = H.GVecs.size();

  // Figure out which bands I'm responsible for
  int numProcs = BandComm.NumProcs();
  int myProc = BandComm.MyProc();
  int band = 0;
  for (int proc=0; proc<numProcs; proc++) {
    int procBands = numBands/numProcs + ((numBands % numProcs)>proc);
    if (proc==myProc) {
      MyFirstBand = band;
      MyLastBand = band + procBands-1;
    }
    band += procBands;
  }

  // Now allocate memory for everything
  Bands.resize (numBands, N);
  lastPhis.resize(numBands, N);
  Energies.resize(numBands);
  Residuals.resize(numBands);
  EtaXiLast.resize(numBands);
  EtaXiLast = complex<double>(0.0, 0.0);
  cnext.resize (N);
  Hc.resize (N);
  Phi.resize (N);
  Phip.resize(N);
  Phipp.resize(N);
  Xi.resize(N);
  Eta.resize(N);
  T.resize(N);
    
  Bands = 0.0;
  Vec3 kBox = H.GVecs.GetkBox();
  double maxk = max(kBox[0], max(kBox[1], kBox[2]));
  double maxE = 1.0*maxk*maxk;
  int numNonZero;
  InitBands();
  IsSetup = true;
}

void ConjGradMPI::InitBands()
{
  int numBands = Bands.rows();
  int numVecs = 8 * numBands;
  assert (numVecs <= H.GVecs.size());
  Array<complex<double>,2> Hmat(numVecs, numVecs);
  Array<complex<double>,2> EigVecs(numBands, numVecs);
  Array<double,1> EigVals (numBands);

  H.Vion->Vmatrix(Hmat);
  perr << "kPoint = " << H.kPoint << endl;
  for (int i=0; i<numVecs; i++) {
    Vec3 Gpk = H.GVecs(i) + H.kPoint;
    Hmat(i,i) += 0.5 * dot (Gpk,Gpk);
  }
  SymmEigenPairs (Hmat, numBands, EigVals, EigVecs);

  if (BandComm.MyProc() == 0)
    for (int i=0; i<numBands; i++)
      perr << "Mini energy(" << i << ") = " << 27.211383*EigVals(i) << endl;

  // Now put results in Bands
  for (int band=0; band<numBands; band++)
    for (int vec=0; vec<Bands.cols(); vec++)
      Bands(band,vec) = 1.0e-6*(drand48()-0.5);
  for (int band=0; band<numBands; band++) 
    for (int i=0; i<numVecs; i++)
      Bands(band, i) = EigVecs(band, i);
  GramSchmidt(Bands);
  // We have to broadcast in order to make sure everyone is starting
  // with the exact same bands.
  BandComm.Broadcast(0,Bands);
}

double ConjGradMPI::CalcPhiSD()
{
  H.Apply(c, Hc);
  E0 = realconjdot (c, Hc);
  perr << "E = " << E0 << endl;
  Phip = E0*c - Hc;
  double nrm = norm (Phip);
  Normalize (Phip);  
  return nrm;
}

void ConjGradMPI::Precondition()
{
  double Tinv = 1.0/T(CurrentBand);
  Vec3 k = H.kPoint;
  for (int i=0; i<c.size(); i++) {
    double x = 0.5*dot(H.GVecs(i)+k, H.GVecs(i)+k)*Tinv;
    double num = 27.0 + 18.0*x +12.0*x*x + 8.0*x*x*x;
    double denom = num + 16.0*x*x*x*x;
    Eta(i) = (num/denom)* Xi(i);
  }
}

// Returns the norm of the residual Hc - Ec
double 
ConjGradMPI::CalcPhiCG()
{
  //  Hc.reference(Hcs(CurrentBand, Range::all()));
  Hc = 0.0;
  H.Kinetic.Apply (c, Hc);
  T = realconjdot (c, Hc);
  if (VHXC.size() != 0)
    H.Vion->Apply (c, Hc, VHXC);
  else
    H.Vion->Apply (c, Hc);
  E0 = realconjdot (c, Hc);
  // Steepest descent direction: (5.10)
  Xi = E0*c - Hc;
  double residualNorm = norm (Xi);
  /// Orthonalize to other bands lower than me here (5.12)
  zVec &Xip = Xi;
  OrthogLower (Bands, Xip, CurrentBand);

  Precondition();

  // Now, orthogonalize to all bands, (including present band);
  // rename for clarity (5.18)
  zVec &Etap = Eta;
  OrthogLower(Bands, Etap, CurrentBand+1);

  // Compute conjugate direction (5.20)
  complex<double> etaxi = conjdot(Etap, Xip);
  complex<double> gamma; 
  if (EtaXiLast(CurrentBand) != complex<double>(0.0, 0.0)) 
    gamma = etaxi/EtaXiLast(CurrentBand);
  else
    gamma = 0.0;

  EtaXiLast(CurrentBand) = etaxi;
  
  // (5.19)
  Phi = Etap + gamma * lastPhis(CurrentBand,Range::all());
  lastPhis(CurrentBand,Range::all()) = Phi;

  // Orthogonalize to all lower bands band: (5.21)
  Phip = Phi;
  //  OrthogExcluding(Bands, Phip, -1);
  OrthogLower(Bands, Phip, CurrentBand+1);
  // (5.22)
  Normalize (Phip);

  Energies(CurrentBand) = E0;
  return residualNorm;
}

inline double max (const Array<double,1> &v)
{
  double mval = v(0);
  for (int i=1; i<v.size(); i++)
    if (v(i) > mval)
      mval = v(i);
  return mval;
}

double
ConjGradMPI::Iterate()
{
  if (!IsSetup)
    Setup();
  for (int band=MyFirstBand; band<=MyLastBand; band++) {
    CurrentBand = band;
    c.reference     (Bands(band, Range::all()));
    Residuals(band) = CalcPhiCG();

    // Now, pick optimal theta for 
    double dE_dtheta = 2.0*realconjdot(Phip, Hc);
    if (VHXC.size() != 0)
      H.Apply (Phip, Hc, VHXC);
    else
      H.Apply (Phip, Hc);
    double d2E_dtheta2 = 2.0*(realconjdot(Phip, Hc) - E0);
    double thetaMin = 0.5*atan(-dE_dtheta/(0.5*d2E_dtheta2));

    double costhetaMin, sinthetaMin;
    sincos(thetaMin, &sinthetaMin, &costhetaMin);
    c = costhetaMin*c + sinthetaMin*Phip;
  }
  CollectBands();
  GramSchmidt(Bands);
  // CheckOverlaps();
  return max(Residuals);
}



void ConjGradMPI::Solve()
{
  /// Hamiltonian has changed, so we can't assume we're conjugate to
  /// previous directions.
  EtaXiLast = complex<double>(0.0, 0.0);
  int iter = 0;
  double residual = 1.0;
  while ((iter < 100) && (residual > Tolerance)) {
    cerr << "iter = " << iter << endl;
    cerr << "Energies = " << Energies << endl;
    residual = Iterate();
    iter++;
  }

  if (residual > Tolerance)
    perr << "Warning:  conjugate gradient residual norm = " 
	 << residual << endl;
  perr << "# of iterations = " << iter << endl;
}

void
ConjGradMPI::PrintOverlaps()
{
  perr << "Overlaps = \n";
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
ConjGradMPI::CheckOverlaps()
{
  zVec x, y;
  for (int i=0; i<Bands.rows(); i++) {
    x.reference (Bands(i, Range::all()));
    for (int j=0; j<Bands.rows(); j++) {
      y.reference (Bands(j,Range::all()));
      double s = realconjdot(x,y);
      if (i==j)
	assert (fabs(s-1.0)<1.0e-12);
      else
	if (fabs(s) > 1.0e-12) {
	  perr << "Overlap(" << i << "," << j << ") = " << s << endl;
	  abort();
	}
     }
   }
}

void
ConjGradMPI::CollectBands()
{
  BandComm.AllGatherRows(Bands);
  BandComm.AllGatherVec(Residuals);
  BandComm.AllGatherVec(Energies);
}


