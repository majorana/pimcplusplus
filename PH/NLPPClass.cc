#include "NLPPClass.h"
#include "../Integration/GKIntegration.h"
#include "../MatrixOps/MatrixOps.h"

bool
NLPPClass::IsNonlocal()
{
  return true;
}

void
NLPPClass::Read(IOSectionClass &in)
{
  assert(in.OpenSection("pseudo"));
  assert(in.ReadVar("AtomicNumber", AtomicNumber));
  assert(in.ReadVar("LocalChannel", lLocal));
  assert(in.ReadVar("ValenceCharge", Zion));
  
  int numChannels = in.CountSections("lChannel");
  vector<Array<double,1> > vl(numChannels), ul(numChannels);
  vector<double> rc(numChannels);
  assert (in.OpenSection("PotentialGrid"));
  Grid *grid = ReadGrid (in);
  in.CloseSection(); // "PotentialGrid"
  for (int i=0; i<numChannels; i++) {
    assert (in.OpenSection("lChannel", i));
    int l;
    assert (in.ReadVar("l", l));
    if (l > numChannels) {
      cerr << "Skipped channels in NLPPClass read.\n";
      abort();
    }
    assert (in.ReadVar("Vl", vl[l]));
    assert (in.ReadVar("ul", ul[l]));
    assert (in.ReadVar("Cutoff", rc[l]));
    in.CloseSection (); // "lChannel"
  }
  
  Vl.resize(numChannels);
  for (int l=0; l<numChannels; l++) {
    Vl[l].l = l;
    Vl[l].V.Init (grid, vl[l]);
    Vl[l].u.Init (grid, ul[l]);
    Array<double,1> deltav(vl[l].size());
    deltav = vl[l] - vl[lLocal];
    Vl[l].DeltaV.Init (grid, deltav);
    Vl[l].rc = rc[l];
    Vl[l].R0 = 1.75*Vl[l].rc;
  }

  in.CloseSection(); // "pseudo"
}


void
NLPPClass::Write(IOSectionClass &out)
{
  cerr << "NLPPClass::Write not implemented.\n";
  abort();
}


void
NLPPClass::SetupProjectors(double G_max, double G_FFT)
{
  for (int i=0; i<Vl.size(); i++) 
    if (Vl[i].l != lLocal) 
      Vl[i].SetupProjector(G_max, G_FFT);

  FILE *fout = fopen ("zeta_r.dat", "w");
  for (double r=0.0; r<50.0; r+=0.001)
    fprintf (fout, "%12.16e %12.16e %12.16e\n", r, 
	     Vl[0].zeta_r(r), Vl[1].zeta_r(r));
  fclose (fout);
  
  fout = fopen ("zeta_q.dat", "w");
  for (double q=0.0; q<G_FFT; q+=0.001)
    fprintf (fout, "%12.16e %12.16e %12.16e\n", q, 
	     Vl[0].zeta_q(q), Vl[1].zeta_q(q));
  fclose (fout);

  fout = fopen ("chi_q.dat", "w");
  for (double q=0.0; q<G_FFT; q+=0.001)
    fprintf (fout, "%12.16e %12.16e %12.16e\n", q, 
	     Vl[0].chi_q(q), Vl[1].chi_q(q));
  fclose (fout);

  fout = fopen ("chi_r.dat", "w");
  for (double r=0.0; r<50.0; r+=0.001)
    fprintf (fout, "%12.16e %12.16e %12.16e\n", r, 
	     Vl[0].chi_r(r), Vl[1].chi_r(r));
  fclose (fout);

}


double
ChannelPotential::A(double q, double qp)
{
  // HACK to avoid using L'Hospital's rule
  if (q == qp)
    q = qp+1.0e-6;

  if (l == 0) 
    return -1.0/(q*q - qp*qp) * 
      q * qp *(q*cos(q*R0)*sin(qp*R0) - qp*cos(qp*R0)*sin(q*R0));
  else if (l == 1)
    return q*qp/(q*q - qp*qp) *
      (q*cos(qp*R0)*sin(q*R0) - qp*cos(q*R0)*sin(qp*R0))
      - 1.0/R0 * sin(q*R0)*sin(qp*R0);
  else if (l == 2)
    return sqrt(-1.0);
  else 
    return sqrt(-1.0);
}

// double
// ChannelPotential::A(double q, double qp, double R0)
// {
//   double sum = 0.0;
//   double delta = 0.01;
//   for (double r=0.0; r<R0; r+= delta)
//     sum += jl(l,q*r)*jl(l,qp*r)*r*r;
//   sum *= delta*q*q*qp*qp;
//   return sum;
// }



void
ChannelPotential::SetupProjector (double G_max, double G_FFT)
{
  Grid &grid = *u.grid;
  // First, compute zeta_r, normalization, and E_KB
  Job = NORM;
  GKIntegration<ChannelPotential> integrator(*this);
  double norm = integrator.Integrate(0.0, grid.End, 1.0e-12);
  ProjectorNorm = 1.0/sqrt(norm);
  Job = EKB;
  double E_KB = norm/integrator.Integrate(0.0, grid.End, 1.0e-12);
  cerr << "l = " << l << "  Norm is " << norm 
       << "  E_KB is " << E_KB << "  R0 = " << R0 << endl;
  
  // Compute zeta(r)
  Array<double,1> zeta(grid.NumPoints);
  zeta(0) = ProjectorNorm * DeltaV(0)*u(1.0e-8)*1.0e8;
  for (int i=1; i<grid.NumPoints; i++) 
    zeta(i) = ProjectorNorm * DeltaV(i)*u(i)/grid(i);
  zeta_r.Init (&grid, zeta);

  // Compute zeta(q)
  Job = ZETA_Q;
  qGrid.Init (0.0, G_FFT, 1000);
  zeta.resize(qGrid.NumPoints);
  for (int i=0; i<qGrid.NumPoints; i++) {
    qCurrent = qGrid(i);
    zeta(i) = integrator.Integrate(0.0, grid.End, 1.0e-12);
  }
  zeta_q.Init (&qGrid, zeta);
  chi_q.Init (&qGrid, zeta);

  double gamma = G_FFT - G_max;
  // Zero out zeta_q above gamma;
  for (int i=0; i<qGrid.NumPoints; i++)
    chi_q(i) = (qGrid(i) >= gamma) ? 0.0 : zeta_q(i);

  // Now for the magic:  We adjust chi_q between G_max and gamma so
  // that the real-space oscillations outside R0 are damped out
  // See King-Smith et al, PRB 44 13063
  // Find index of gamma
  int gammaIndex = qGrid.ReverseMap(gamma);
  int G_maxIndex  = qGrid.ReverseMap(G_max)+1;
  double delta = qGrid(1)-qGrid(0);
  int nb = gammaIndex - G_maxIndex;
  Array<double,1> b(nb), x(nb);
  // First, create the b vector
  b = 0.0;
  for (int i=0; i<nb; i++) {
    double q = qGrid(G_maxIndex+i);
    for (int j=0; j<G_maxIndex; j++) {
      double qp = qGrid (j);
      b(i) -= delta * A(q, qp)*zeta_q(j);
    }
  }
  // Now, create the M matrix
  Array<double,2> M(nb, nb);
  for (int i=0; i<nb; i++) {
    double q = qGrid(G_maxIndex+i);
    for (int j=0; j<nb; j++) {
      double qp = qGrid(G_maxIndex+j);
      M(i,j) = delta*A(q, qp);
    }
    M(i,i) -= 0.5*M_PI*q*q;
  }
//   if (l==0) {
//     FILE *fout = fopen ("M.dat", "w");
//     for (int i=0; i<nb; i++) {
//       for (int j=0; j<nb; j++) 
// 	fprintf (fout, "%24.16e ", M(i,j));
//       fprintf (fout, "\n");
//     }
//     fclose (fout);
//   }
  // Now solve Mx = b
  Array<int,1> perm;
  double sign;
  LUdecomp (M, perm, sign);
  LUsolve (M, perm, b);
  x = b;
  for (int i=0; i<nb; i++) 
    chi_q(G_maxIndex+i) = x(i);
  chi_q(G_maxIndex+nb) = 0.0;
  
  // Now transform back to real-space, computing chi(r)
  Job = CHI_R;
  Array<double,1> chi(grid.NumPoints);
  for (int i=0; i<grid.NumPoints; i++) {
    rCurrent = grid(i);
    chi(i) = integrator.Integrate(0.0, qGrid.End, 1.0e-8);
  }
  chi_r.Init (&grid, chi);

  // Finally, check to see if chi_r is small outside R0
  Job = CHECK_CHI_R;
  double norm2 = integrator.Integrate(0.0, qGrid.End, 1.0e-8);
  double error = integrator.Integrate( R0, qGrid.End, 1.0e-8);
  if (error > 1.0e-10)
    cerr << "Fractional error in real-space projection = "
	 << (error / norm2) << endl;
}


double
NLPPClass::V(double r)
{
  return Vl[lLocal].V(r);
}

double
NLPPClass::dVdr(double r)
{
  return Vl[lLocal].V.Deriv(r);
}

double
NLPPClass::d2Vdr2(double r)
{
  return Vl[lLocal].V.Deriv2(r);
}
