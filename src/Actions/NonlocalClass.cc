#include "NonlocalClass.h"

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
  QuadPoints.clear();
  QuadWeights.clear();

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
      QuadPoints.push_back(a[i]);
      QuadWeights.push_back(A);
    }
  if (std::fabs(B) > 1.0e-10) 
    for (int i=0; i<b.size(); i++) {
      QuadPoints.push_back(b[i]);
      QuadWeights.push_back(B);
    }
  if (std::fabs(C) > 1.0e-10) 
    for (int i=0; i<c.size(); i++) {
      QuadPoints.push_back(c[i]);
      QuadWeights.push_back(C);
    }
  if (std::fabs(D) > 1.0e-10) 
    for (int i=0; i<d.size(); i++) {
      QuadPoints.push_back(d[i]);
      QuadWeights.push_back(D);
    }
  
  // Allocate storage for wave function ratios
  //  pp_nonloc->resize_warrays(nk,NumNonLocal,Lmax);
  
  // Finally, check the rule for correctness
  assert (QuadPoints.size()  == nk);
  assert (QuadWeights.size() == nk);
  double wSum = 0.0;
  for (int k=0; k < nk; k++) {
    Vec3 r = QuadPoints[k];
    double nrm = dot(r,r);
    assert (std::fabs(nrm-1.0) < 1.0e-14);
    wSum += QuadWeights[k];
  }
  assert (std::fabs(wSum - 1.0) < 1.0e-14);
  // Check the quadrature rule
  CheckQuadratureRule(lexact);
}


void
NonlocalClass::CheckQuadratureRule(int lexact)
{
  vector<Vec3> &grid = QuadPoints;
  vector<double> &w = QuadWeights;
  for (int l1=0; l1<=lexact; l1++) 
    for (int l2=0; l2 <= (lexact-l1); l2++) 
      for (int m1=-l1; m1<=l1; m1++)
	for (int m2=-l2; m2<=l2; m2++) {
	  complex<double> sum(0.0, 0.0);
	  for (int k=0; k<grid.size(); k++) {
	    complex<double> v1 = Ylm(l1, m1, grid[k]);
	    complex<double> v2 = Ylm(l2, m2, grid[k]);
	    sum += 4.0*M_PI*w[k] * conj(v1)*v2;
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
  // Compute DeltaV splines for nonlocal channels

  // Setup quadrature rule -- first, check all the rules
  for (int i=1; i<=7; i++)
    SetQuadratureRule (i);
  // Set it to a reasonable default
  SetQuadratureRule (4);
}

void
NonlocalClass::Read (IOSectionClass &in)
{
  

}

double
NonlocalClass::SingleAction(int slice1, int slice2,
			    const Array<int,1> &activeParticles, int level)
{
  return 0.0;
}


double 
NonlocalClass::d_dBeta (int slice1, int slice2, int level)
{
  return 0.0;
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
