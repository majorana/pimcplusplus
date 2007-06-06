#include "NonlocalClass.h"
#include "../PathDataClass.h"

NonlocalClass::NonlocalClass (PathDataClass &pathData) :
  ActionBaseClass (pathData)
{
  // Nothing else for now

}

inline double 
LegendrePll (int l, double x) {
  if (l==0)
    return 1.0;
  else {
    double sqt = std::sqrt(1.0-x)*std::sqrt(1.0+x);
    double val = 1.0;
    double dblfact = 1.0;
    for (int i=1; i<=l; i++) {
      val *= -sqt;
      val *= dblfact;
      dblfact += 2.0;
    }
    return val;
  }
}

inline double 
LegendrePlm (int l, int m, double x) {
  if (m < 0) {
    m = std::abs (m);
    double posval = LegendrePlm (l, m, x);
    double sign = (m%2==0) ? 1.0 : -1.0;
    double mfact = 1.0;
    for (int i=2; i<=(l-m); i++)
      mfact *= static_cast<double>(i);
    double pfact = 1.0;
    for (int i=2; i<=(l+m); i++)
      pfact *= static_cast<double>(i);
    return posval * sign*mfact/pfact;
  }
  // Now we can assume that m is not negative
  double pmm = LegendrePll (m, x);
  double pmp1m = x*(2*m+1)*pmm;
  
  if (m == l) 
    return pmm;
  else if (l==(m+1))
    return pmp1m;
  else { // Use recursive formula
    double Plm2m = pmm;
    double Plm1m = pmp1m;
    double Pl;
    for (int i=m+2; i<=l;  i++) {
      Pl = (1.0/static_cast<double>(i-m)) * (x*(2*i-1)*Plm1m - (i+m-1)*Plm2m);
      Plm2m = Plm1m;
      Plm1m = Pl;
    }
    return Pl;
  }
}

inline std::complex<double> 
Ylm(int l, int m, const TinyVector<double,3>& r)
{
  double costheta = r[0];
  double phi   = std::atan2(r[2],r[1]);
  int lmm = l - m;
  int lpm = l + m;
  double mfact = 1.0;
  double pfact = 1.0;
  for (int i=lmm; i>0; i--)
    mfact *= static_cast<double>(i);
  for (int i=lpm; i>0; i--)
    pfact *= static_cast<double>(i);
  double prefactor = std::sqrt (static_cast<double>(2*l+1)*mfact/(4.0*M_PI*pfact));
  prefactor *= LegendrePlm (l, m, costheta);
  return std::complex<double>(prefactor*std::cos(m*phi), prefactor*std::sin(m*phi));
}


void
NonlocalClass::SetQuadratureRule (int rule)
{
  int nk;
  double w;
  typedef enum {SINGLE, TETRA, OCTA, ICOSA} SymmType;
  SymmType symmetry;
  int lexact;
  double A, B, C, D;
  A = B = C = D = 0.0;
  vector<Vec3> points;
  vector<double> weights;

  switch (rule) {
  case 1:
    nk = 1;   symmetry = SINGLE; lexact = 0;
    A = 1.0;
    break;
  case 2:
    nk = 4;   symmetry = TETRA;  lexact = 2;
    A=0.25;
    break;
  case 3:
    nk = 6;   symmetry = OCTA;   lexact = 3;
    A=1.0/6.0;
    break;
  case 4:
    nk = 12;  symmetry = ICOSA;  lexact = 5;
    A = 1.0/12.0;
    B = 1.0/12.0;
    break;
  case 5:
    nk = 18;  symmetry = OCTA;   lexact = 5;
    A = 1.0/30.0; 
    B = 1.0/15.0;
    break;
  case 6:
    nk = 26;  symmetry = OCTA;   lexact = 7;
    A = 1.0  / 21.0;
    B = 4.0  / 105.0;
    C = 27.0 / 840.0;
    break;
  case 7:
    nk = 50;  symmetry = OCTA;   lexact = 11;
    A = 4.0/315.0;
    B = 64.0/2835.0;
    C = 27.0/1280.0;
    D = 14641.0/725760.0;
    break;
  default:
    cerr << "Unrecognized spherical quadrature rule " << rule << ".";
    abort();
  }
  
  // First, build a_i, b_i, and c_i points
  vector<Vec3> a, b, c, d;
  double p = 1.0/std::sqrt(2.0);
  double q = 1.0/std::sqrt(3.0);
  double r = 1.0/std::sqrt(11.0);
  double s = 3.0/std::sqrt(11.0);
  
  if (symmetry == SINGLE) {
    a.push_back (Vec3(1.0, 0.0, 0.0));
  }
  else if (symmetry == TETRA) {
    a.push_back(Vec3( q, q, q));
    a.push_back(Vec3( q,-q,-q));
    a.push_back(Vec3(-q, q,-q));
    a.push_back(Vec3(-q,-q, q));
  }
  else if (symmetry == OCTA) {
    a.push_back(Vec3( 1.0, 0.0, 0.0));
    a.push_back(Vec3(-1.0, 0.0, 0.0));
    a.push_back(Vec3( 0.0, 1.0, 0.0));
    a.push_back(Vec3( 0.0,-1.0, 0.0));
    a.push_back(Vec3( 0.0, 0.0, 1.0));
    a.push_back(Vec3( 0.0, 0.0,-1.0));
    
    b.push_back(Vec3(   p,   p, 0.0));
    b.push_back(Vec3(   p,  -p, 0.0));
    b.push_back(Vec3(  -p,   p, 0.0));
    b.push_back(Vec3(  -p,  -p, 0.0));
    b.push_back(Vec3(   p, 0.0,   p));
    b.push_back(Vec3(   p, 0.0,  -p));
    b.push_back(Vec3(  -p, 0.0,   p));
    b.push_back(Vec3(  -p, 0.0,  -p));
    b.push_back(Vec3( 0.0,   p,   p));
    b.push_back(Vec3( 0.0,   p,  -p));
    b.push_back(Vec3( 0.0,  -p,   p));
    b.push_back(Vec3( 0.0,  -p,  -p));
    
    c.push_back(Vec3(   q,   q,   q));
    c.push_back(Vec3(   q,   q,  -q));
    c.push_back(Vec3(   q,  -q,   q));
    c.push_back(Vec3(   q,  -q,  -q));
    c.push_back(Vec3(  -q,   q,   q));
    c.push_back(Vec3(  -q,   q,  -q));
    c.push_back(Vec3(  -q,  -q,   q));
    c.push_back(Vec3(  -q,  -q,  -q));
    
    d.push_back(Vec3(   r,   r,   s));
    d.push_back(Vec3(   r,   r,  -s));
    d.push_back(Vec3(   r,  -r,   s));
    d.push_back(Vec3(   r,  -r,  -s));
    d.push_back(Vec3(  -r,   r,   s));
    d.push_back(Vec3(  -r,   r,  -s));
    d.push_back(Vec3(  -r,  -r,   s));
    d.push_back(Vec3(  -r,  -r,  -s));
    
    d.push_back(Vec3(   r,   s,   r));
    d.push_back(Vec3(   r,   s,  -r));
    d.push_back(Vec3(   r,  -s,   r));
    d.push_back(Vec3(   r,  -s,  -r));
    d.push_back(Vec3(  -r,   s,   r));
    d.push_back(Vec3(  -r,   s,  -r));
    d.push_back(Vec3(  -r,  -s,   r));
    d.push_back(Vec3(  -r,  -s,  -r));
    
    d.push_back(Vec3(   s,   r,   r));
    d.push_back(Vec3(   s,   r,  -r));
    d.push_back(Vec3(   s,  -r,   r));
    d.push_back(Vec3(   s,  -r,  -r));
    d.push_back(Vec3(  -s,   r,   r));
    d.push_back(Vec3(  -s,   r,  -r));
    d.push_back(Vec3(  -s,  -r,   r));
    d.push_back(Vec3(  -s,  -r,  -r));
  }
  else if (symmetry == ICOSA) {
    double t, p;  // theta and phi
    // a points
    t = 0.0;  p=0.0;
    a.push_back(Vec3(std::cos(t),std::sin(t)*std::cos(p),std::sin(t)*std::sin(p)));
    t = M_PI;  p=0.0;
    a.push_back(Vec3 (std::cos(t),std::sin(t)*std::cos(p),std::sin(t)*std::sin(p)));
    // b points
    for (int k=0; k<5; k++) {
      t = std::atan(2.0);          p = (double)(2*k+0)*M_PI/5.0;
      b.push_back(Vec3 (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
      t = M_PI-std::atan(2.0);     p = (double)(2*k+1)*M_PI/5.0;
      b.push_back(Vec3 (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
    }
    // c points
    double t1 = std::acos ((2.0+std::sqrt(5.0)) / std::sqrt(15.0+6.0*std::sqrt(5.0)));
    double t2 = std::acos (      1.0            / std::sqrt(15.0+6.0*std::sqrt(5.0)));
    for (int k=0; k<5; k++) {
      t = t1; p = (double)(2*k+1)*M_PI/5.0;
      c.push_back(Vec3 (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
      t = t2; p = (double)(2*k+1)*M_PI/5.0;
      c.push_back(Vec3 (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
      t = M_PI - t1; p = (double)(2*k+0)*M_PI/5.0;
      c.push_back(Vec3 (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
      t = M_PI - t2; p = (double)(2*k+0)*M_PI/5.0;
      c.push_back(Vec3 (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
    }
  }
  // Now, construct rule
  if (std::fabs(A) > 1.0e-10) 
    for (int i=0; i<a.size(); i++) {
      points.push_back(a[i]);
      weights.push_back(A);
    }
  if (std::fabs(B) > 1.0e-10) 
    for (int i=0; i<b.size(); i++) {
      points.push_back(b[i]);
      weights.push_back(B);
    }
  if (std::fabs(C) > 1.0e-10) 
    for (int i=0; i<c.size(); i++) {
      points.push_back(c[i]);
      weights.push_back(C);
    }
  if (std::fabs(D) > 1.0e-10) 
    for (int i=0; i<d.size(); i++) {
      points.push_back(d[i]);
      weights.push_back(D);
    }
  
  QuadPoints.resize(points.size());
  QuadWeights.resize(points.size());
  for (int i=0; i<points.size(); i++) {
    QuadPoints(i)  = points[i];
    QuadWeights(i) = weights[i];
  }
  
  // Allocate storage for wave function ratios
  //  pp_nonloc->resize_warrays(nk,NumNonLocal,Lmax);
  
  // Finally, check the rule for correctness
  assert (QuadPoints.size()  == nk);
  assert (QuadWeights.size() == nk);
  double wSum = 0.0;
  for (int k=0; k < nk; k++) {
    Vec3 r = QuadPoints(k);
    double nrm = dot(r,r);
    assert (std::fabs(nrm-1.0) < 1.0e-14);
    wSum += QuadWeights(k);
  }
  assert (std::fabs(wSum - 1.0) < 1.0e-14);
  // Check the quadrature rule
  CheckQuadratureRule(lexact);
}


void
NonlocalClass::CheckQuadratureRule(int lexact)
{
  Array<Vec3,1> &grid = QuadPoints;
  Array<double,1> &w = QuadWeights;
  for (int l1=0; l1<=lexact; l1++) 
    for (int l2=0; l2 <= (lexact-l1); l2++) 
      for (int m1=-l1; m1<=l1; m1++)
	for (int m2=-l2; m2<=l2; m2++) {
	  complex<double> sum(0.0, 0.0);
	  for (int k=0; k<grid.size(); k++) {
	    complex<double> v1 = Ylm(l1, m1, grid(k));
	    complex<double> v2 = Ylm(l2, m2, grid(k));
	    sum += 4.0*M_PI*w(k) * conj(v1)*v2;
	  }
	  double re = real (sum);
	  double im = imag (sum);
	  if ((l1==l2) && (m1==m2)) 
	    re -= 1.0;
	  if ((std::fabs(im) > 1.0e-14) || (std::fabs(re) > 1.0e-14)) {
	    cerr << "Broken spherical quadrature for " << grid.size() << "-point rule.\n" << endl;
	    abort();
	  }
	  // 	    fprintf (stderr, "(l1,m1,l2m,m2) = (%2d,%2d,%2d,%2d)  sum = (%20.16f %20.16f)\n",
	  // 	     l1, m1, l2, m2, real(sum), imag(sum));
	  
	}
}

void
NonlocalClass::Setup (FixedPhaseClass *fixedPhase)
{
  FixedPhase = fixedPhase;
  // Find potential
  NLPP = dynamic_cast<NLPPClass*>
    (PathData.Actions.PairMatrix(FixedPhase->IonSpeciesNum,
				 FixedPhase->UpSpeciesNum)->Pot);
  if (NLPP == NULL) {
    cerr << "The potential passed to Nonlocal action is not an NLPP.\n";
    abort();
  }

  // Compute DeltaV splines for nonlocal channels


  // Setup quadrature rule -- first, check all the rules
  for (int i=1; i<=7; i++)
    SetQuadratureRule (i);
  // Set it to a reasonable default
  SetQuadratureRule (4);
  ScaledPoints.resize(QuadPoints.size());
  WFratios.resize(QuadPoints.size());
  Legendre.resize(NLPP->NumChannels());
  DeltaV.resize(NLPP->NumChannels());
  SpeciesClass &up   = PathData.Path.Species(FixedPhase->UpSpeciesNum);
  SpeciesClass &down = PathData.Path.Species(FixedPhase->DownSpeciesNum);
  Electrons.resize (up.NumParticles + down.NumParticles);
  for (int i=0; i<up.NumParticles; i++)
    Electrons(i) = i+up.FirstPtcl;
  for (int i=0; i<down.NumParticles; i++)
    Electrons(i+up.NumParticles) = i + down.FirstPtcl;
}

void
NonlocalClass::Read (IOSectionClass &in)
{
  

}

void
NonlocalClass::NearestIon (int slice, int ptcl, Vec3 &ionpos, double &dist)
{
  PathClass &Path = PathData.Path;
  Vec3 r = Path (slice, ptcl);
  SpeciesClass &ions = Path.Species(FixedPhase->IonSpeciesNum);
  int first = ions.FirstPtcl;  int last = ions.LastPtcl;
  Vec3 disp;

  ionpos = Path(slice,first);
  disp = r - ionpos;
  Path.PutInBox(disp);
  dist = dot (disp, disp);
  for (int ion=first+1; ion <= last; ion++) {
    disp = r - Path(slice, ion);
    Path.PutInBox(disp);
    double tryDist = dot (disp, disp);
    if (tryDist < dist) {
      dist = tryDist;
      ionpos = Path(slice, ion);
    }
  }
  dist = sqrt (dist);
}

void
NonlocalClass::ScaleQuadPoints (Vec3 ionpos, double dist)
{
  // Setup rotation matrix

  // Scale and translate points;
  for (int i=0; i<ScaledPoints.size(); i++) {
    ScaledPoints(i) = ionpos + dist * QuadPoints(i);
    PathData.Path.PutInBox(ScaledPoints(i));
  }
}

double
NonlocalClass::SingleAction(int slice1, int slice2,
			    const Array<int,1> &activeParticles, int level)
{
  PathClass &Path = PathData.Path;
  int skip = 1<<level;
  double U = 0.0;
  // Find maximum cutoff radius
  double max_rc = 0.0;
  for (int l=0; l<NLPP->NumChannels(); l++)
    if (l != NLPP->LocalChannel())
      if (NLPP->Getrc(l) > max_rc)
	max_rc = NLPP->Getrc(l);
  
  double levelTau = ldexp (Path.tau, level);

  // Outer loop overs slices
  for (int slice=slice1; slice <= slice2; slice++) {
    double prefactor = (slice==slice1 || slice==slice2) ? 0.5 : 1.0;
    prefactor *= levelTau/(4.0*M_PI);
    for (int ptclIndex=0; ptclIndex < activeParticles.size(); ptclIndex++) {
      int ptcl = activeParticles(ptclIndex);
      Vec3 ionpos;
      double dist;
      NearestIon (slice, ptcl, ionpos, dist);
      if (dist < max_rc) {
	ScaleQuadPoints (ionpos, dist);
	double distInv = 1.0/dist;
	// Compute the wave function ratios
	FixedPhase->CalcWFratios (slice, ptcl, ScaledPoints, WFratios);
	// Now loop over l-channels
	for (int l=0; l<NLPP->NumChannels(); l++)
	  DeltaV(l) = NLPP->GetDeltaV(l,dist);
	DeltaV(NLPP->LocalChannel()) = 0.0;
	for (int pi=0; pi<ScaledPoints.size(); pi++) {
	  double costheta = distInv*distInv*dot(ScaledPoints(pi), ionpos);
	  double P_l, P_lm1, P_lp1;
	  P_l = 1.0;
	  P_lm1 = 0.0;
	  for (int l=0; l<NLPP->NumChannels(); l++) {
	    double dl = (double)l;
	    P_lp1 = ((2.0*dl+1.0)*costheta*P_l - dl*P_lm1)/(1.0+dl);
	    U += prefactor * DeltaV(l) * P_l * real(WFratios(pi)) * QuadWeights(pi);
	    P_lm1 = P_l;
	    P_l = P_lp1;
	  }
	}
      }
    }
  }
  //  cerr << "U = " << U << endl;
  return U;
}


double 
NonlocalClass::d_dBeta (int slice1, int slice2, int level)
{
  double levelTau = ldexp (Path.tau, level);
  double U = SingleAction (slice1, slice2, Electrons, level);
  return U / levelTau;
}

void 
NonlocalClass::GradAction (int slice1, int slice2, const Array<int,1> &ptcls,
			   int level, Array<dVec,1> &gradVec)
{


}

string 
NonlocalClass::GetName()
{
  return "Nonlocal";
}
