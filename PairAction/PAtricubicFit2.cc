#include "PAFit.h"
#include "../Fitting/Fitting.h"
#include <gsl/gsl_sf.h>

const double Rho0Min  = 1.0e-4;

/// The following routines are used only if we are creating fits, not
/// using them.
#ifdef MAKE_FIT
void 
PAtricubicFit2Class::ReadParams(IOSectionClass &inSection)
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
}

void 
PAtricubicFit2Class::WriteBetaIndependentInfo (IOSectionClass &outSection)
{
  if (Comm.MyProc() == 0) {
    outSection.WriteVar ("Type", "tricubicfit2");
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


class SCintegrand2
{
  Vec2 r, rp, rpp;
  Rho &rho;
  LinearGrid qgrid;
  CubicSpline Udiag, dUdiag;
public:
  bool IsdU;
  double operator()(double x)
  {
    rpp = r + x*(rp-r);
    double dist = sqrt(rpp[0]*rpp[0]+rpp[1]*rpp[1]);
    //    return (PH->V(dist));
    if ((dist >= rho.grid->Start) && (dist <= rho.grid->End)) 
      return (IsdU ? dUdiag(dist) : Udiag(dist));
    else {
      double V = rho.Pot->V(dist);
      return (IsdU ? V : rho.Beta()*V);
    }
  }
  void Set (double rmag, double rpmag, double costheta)
  {
    r[0] = rmag;
    r[1] = 0.0;
    rp[0] = rpmag * costheta;
    rp[1] = rpmag * sqrt(1.0-costheta*costheta);
  }
  SCintegrand2(Rho &rho_) : rho(rho_)
  {
    int N=1000;
    double qmin = rho.grid->Start;
    double qmax = rho.grid->End;
    Array<double,1> Ud(N), dUd(N);
    qgrid.Init(qmin, qmax, N);
    double sigma = sqrt (2.0*rho.lambda*rho.Beta());
    for (int i=0; i<N; i++) {
      double q = qgrid(i);
      /// HACK HACK HACK
      /// We correct for poor extrapolation of the tail for
      /// PH's at small beta.  The 22 is empirically determined.
      if ((q > 3.5) && (q/sigma > 22.0)) {
	double V = rho.Pot->V(q);
	Ud(i)  = rho.Beta()*V;
	dUd(i) = V;
      }
      else {
	double U, dU;
	// Calculate diagonal action
	rho.UdU(q, q, 1.0, U, dU);
	Ud(i) = U;
	dUd(i) = dU;
      }
    }
    Udiag.Init(&qgrid, Ud);
    dUdiag.Init(&qgrid, dUd);
  }
};

class USemiclassical2
{
private:

  SCintegrand2 SC;
  double beta;
public:

  double U(double r, double rp, double costheta)
  {
    SC.Set(r, rp, costheta);
    SC.IsdU = false;
    GKIntegration<SCintegrand2,GK15> Integrator(SC);
    double Uavg = Integrator.Integrate(0.0, 1.0, 1.0e-7,
				       1.0e-7, false);
    return (Uavg);
  }
  double dU(double r, double rp, double costheta)
  {
    SC.Set(r, rp, costheta);
    SC.IsdU = true;
    GKIntegration<SCintegrand2,GK15> Integrator(SC);
    double dUavg = Integrator.Integrate(0.0, 1.0, 1.0e-7,
					1.0e-7, false);
    return (dUavg);
  }
  USemiclassical2 (Rho &rho, double beta_) : SC(rho)
  {
    beta = beta_;
  }
};


void PAtricubicFit2Class::DoFit (Rho &rho)
{
  Usplines.resize(1);
  dUsplines.resize(1);
  sMax.resize(1);   sMaxInv.resize(1);

  double lambda = Particle1.lambda + Particle2.lambda;
  SmallestBeta = rho.Beta();
  double beta = SmallestBeta;
  sMax(0)  = sqrt(-4.0*lambda*beta*log(Rho0Min));
  sMaxInv(0) = 1.0/sMax(0);

  cerr << "sMax = "  << sMax(0) << endl;

  int numq = qgrid->NumPoints;
  int numy = ygrid->NumPoints;
  int numt = tgrid->NumPoints;
  Array<double,3> Umat(numq,numy,numt), dUmat(numq,numy,numt);
  Array<double,1> Ul, dUl;
  USemiclassical2 Usemi(rho, beta);

  double qmax = max (3.5, 22.0*sqrt(2.0*beta*lambda));
  cerr << "qmax = " << qmax << endl;

  for (int qi=0; qi<numq; qi++) {
    double q = (*qgrid)(qi);
    cerr << "qi = " << qi << " of " << numq << " q = " << q << endl;
    for (int yi=0; yi<numy; yi++) {
      double y = (*ygrid)(yi);
      double z, s;
      yt2zs (q, y, 0.0, 0, z, s);
      double r = q+0.5*z; 
      double rp= q-0.5*z;
      if (r == 0.0)
	r = 1.0e-8;
      if (rp == 0.0)
	rp = 1.0e-8;
      
      rho.U_lArray(r,rp,Ul,dUl);
      
      for (int ti=0; ti<numt; ti++) {
	double t = (*tgrid)(ti);
	yt2zs (q, y, t, 0, z, s);
	double costheta;
	if ((r*rp)==0.0)	
	  costheta = 1.0;
	else
	  costheta = (r*r + rp*rp - s*s)/(2.0*r*rp);
	costheta = min(1.0, costheta);
	costheta = max(-1.0, costheta);	
	double U, dU;
	// HACK HACK HACK
	if (q > qmax) {
	  U  = Usemi.U (r,rp,costheta);
	  dU = Usemi.dU(r,rp,costheta);
	}
	else {
	  rho.UdU(r,rp,costheta, Ul, dUl, U, dU);
	  if (isnan(U)) {
	    //U = Usemi.U(r,rp,costheta);
	    fprintf (stderr, "NAN in U at (qi,yi,ti) = (%d,%d,%d)\n",
		     qi, yi, ti); 
	    fprintf (stderr, "(q, z, s) = (%1.5e %1.5e %1.5e)\n", 
		     q, z, s);
	    fprintf (stderr, "exp(-s^2/4lb) = %1.4e\n", 
		     exp(-s*s/(4.0*lambda*beta)));
	  }
	  if (isnan(dU)) {
	    //dU = Usemi.dU(r,rp,costheta);
	    fprintf (stderr, "NAN in dU at (qi,yi,ti) = (%d,%d,%d)\n",
		     qi, yi, ti);
	    fprintf (stderr, "(q, z, s) = (%1.5e %1.5e %1.5e)\n", 
		     q, z, s);
	  }
	}
	Umat(qi,yi,ti) = U;
	dUmat(qi,yi,ti) = dU;
      }
    } 
  }
  Usplines(0).Init(qgrid,ygrid,tgrid,Umat);
  dUsplines(0).Init(qgrid,ygrid,tgrid,dUmat);
}



void 
PAtricubicFit2Class::Error(Rho &rho, double &Uerror, double &dUerror)
{
  int level = (int)floor(log(rho.Beta()/SmallestBeta)/log(2.0)+ 0.5);

  double lambda = Particle1.lambda + Particle2.lambda;
  double beta = rho.Beta();

  double U2err = 0.0;
  double dU2err = 0.0;
  double weight = 0.0;
  FILE *Uxdat = fopen ("Ux.dat", "w");
  FILE *Ufdat = fopen ("Uf.dat", "w");
  FILE *dUxdat = fopen ("dUx.dat", "w");
  FILE *dUfdat = fopen ("dUf.dat", "w");
  FILE *tdat = fopen ("t.dat", "w");
  FILE *sdat = fopen ("s.dat", "w");
  FILE *costhetadat = fopen ("costheta.dat", "w");
  FILE *ydat = fopen ("y.dat", "w");
  LinearGrid qgrid2(qgrid->Start, 0.999*qgrid->End, 20);
  Array<double,1> Ul, dUl;
  for (int qi=0; qi<qgrid2.NumPoints; qi++) {
    double q = qgrid2(qi);
    double zmax = min(sMax(0),2.0*q);
    LinearGrid zgrid(0.0, zmax, 21);
    for (int zi=0; zi<zgrid.NumPoints; zi++) {
      double z = zgrid(zi);
      LinearGrid sgrid(z, zmax, 21);
      for (int si=0; si<sgrid.NumPoints; si++) {
	double s = sgrid(si);
	double y, t;
	zs2yt(q,z,s,0,y,t);
	double w = exp(-s*s/(4.0*rho.lambda*rho.Beta()));
	double Uex, dUex, Ufit, dUfit;
	double r, rp, costheta;
	r  = q+0.5*z;
	rp = q-0.5*z;
	if ((r*rp)==0.0)
	  costheta = 1.0;
	else
	  costheta = (r*r + rp*rp - s*s)/(2.0*r*rp); 
	
	costheta = min(costheta,1.0);
	costheta = max(costheta,-1.0);

	rho.U_lArray(r,rp,Ul, dUl);
	rho.UdU(r, rp, costheta, Ul, dUl, Uex, dUex);

	Ufit = U(q, z, s*s, level);
	dUfit = dU(q, z, s*s, level);
	if (!isnan(Uex) && !isnan(dUex)) {
	  if (s <= zmax) {
	    U2err += w*(Uex-Ufit)*(Uex-Ufit);
	    dU2err += w*(dUex-dUfit)*(dUex-dUfit);
	    weight += w;
	  }
	}
	fprintf (Uxdat, "%1.16e ", Uex);
	fprintf (Ufdat, "%1.16e ", Ufit);
	fprintf (dUxdat, "%1.16e ", dUex);
	fprintf (dUfdat, "%1.16e ", dUfit);
	fprintf (tdat, "%1.16e ", t);
	fprintf (sdat, "%1.16e ", s);
	fprintf (costhetadat, "%1.16e ", costheta);
	fprintf (ydat, "%1.16e ", y);
      }
      fprintf (Uxdat, "\n");
      fprintf (Ufdat, "\n");
      fprintf (dUxdat, "\n");
      fprintf (dUfdat, "\n");
      fprintf (tdat, "\n");
      fprintf (sdat, "\n");
      fprintf (ydat, "\n");
      fprintf (costhetadat, "\n");
    }
  }
  fclose (Uxdat); fclose(Ufdat); 
  fclose (dUxdat); fclose(dUfdat);fclose(tdat); fclose(ydat); fclose(sdat);
  fclose(costhetadat);
  Uerror = sqrt(U2err/weight);
  dUerror = sqrt(dU2err/weight);
}


void 
PAtricubicFit2Class::WriteFit(IOSectionClass &outSection)
{
  if (Comm.MyProc() == 0) {
    Array<double,3> Umat(qgrid->NumPoints,ygrid->NumPoints,tgrid->NumPoints); 
    Array<double,3> dUmat(qgrid->NumPoints,ygrid->NumPoints,tgrid->NumPoints); 
    double beta = SmallestBeta;
    outSection.NewSection("Fit");
    outSection.WriteVar ("beta", beta);
    for (int qi=0; qi<qgrid->NumPoints; qi++)
      for (int yi=0; yi<ygrid->NumPoints; yi++) 
	for (int ti=0; ti<tgrid->NumPoints; ti++) {
	  Umat(qi,yi,ti)  =  Usplines(0)(qi,yi,ti);
	  dUmat(qi,yi,ti) = dUsplines(0)(qi,yi,ti);
	}
    outSection.WriteVar ("Umat", Umat);
    outSection.WriteVar ("dUmat", dUmat);
    outSection.WriteVar ("sMax", sMax(0));
    outSection.CloseSection();
  }
}

#endif


// double PAtricubicFit2Class::U(double q, double z, double s2, int level)
// {
//   z = fabs(z);
//   double qmax = qgrid->End;
//   double smax = min (2.0*q, sMax(level));
//   double zmax = smax;
//   double smin = z;
//   double s=sqrt(s2);
//   double x;

//   if ((q<=qmax)&&(z<=zmax)&&(s<=smax)) {
//     if (q == 0.0) 
//       return(Usplines(level)(0.0,0.0,0.0));
//     else {
//       double y = z/zmax;
//       double t = (s-z)/(smax-z);
//       return(Usplines(level)(q,y,t));
//     }
//   }
//   else {
//     double beta = SmallestBeta;
//     for (int i=0; i<level; i++)
//       beta *= 2.0;
//     double r = q+0.5*z;
//     double rp = q-0.5*z;
//     return (0.5*beta*(Pot->V(r)+Pot->V(rp)));
//   }
// }

double 
PAtricubicFit2Class::U(double q, double z, double s2, int level)
{
  z = fabs(z);
  double y, t;
  double s = sqrt(s2);
  zs2yt (q, z, s, level, y, t);

  if (q<=qgrid->End) {
    if (q == 0.0) 
      return(Usplines(level)(0.0,0.0,0.0));
    else {
      if (t < 1.0) {
	double spline = Usplines(level)(q,y,t);
	if (isnan (spline)) {
	  cerr << "q = " << q << " s = " << s << " z = " << z << endl;
	  cerr << "NAN in spline!!!!!!!!\n";
	}
	return(spline);
      }
      else if (y < 1.0) {
	double val   = Usplines(level)(q, y, 1.0);
	double deriv = Usplines(level).d_dz(q,y,1.0);
	return (val + deriv*(t-1.0));
      }
      else {
	double val   = Usplines (level)(q, 1.0, 1.0);
	double deriv = Usplines (level).d_dy(q,1.0,1.0);
	return (val + deriv*(y-1.0));
      }
    }
  }
  else {
    double beta = SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    double r = q+0.5*z;
    double rp = q-0.5*z;
    double prim  = 0.5*beta*(Pot->V(r)+Pot->V(rp));
    if (isnan(prim))
      cerr << "NAN in prim!!!!!\n";
    return (prim);
  }
}

double PAtricubicFit2Class::dU(double q, double z, double s2, int level)
{
  z = fabs(z);
  double y, t;
  double s = sqrt(s2);
  zs2yt (q, z, s, level, y, t);

  if (q<=qgrid->End) {
    if (q == 0.0) 
      return(dUsplines(level)(0.0,0.0,0.0));
    else {
      if (t < 1.0) {
	double spline = dUsplines(level)(q,y,t);
	return(spline);
      }
      else if (y < 1.0) {
	double val   = dUsplines(level)(q, y, 1.0);
	double deriv = dUsplines(level).d_dz(q,y,1.0);
	return (val + deriv*(t-1.0));
      }
      else {
	double val   = dUsplines (level)(q, 1.0, 1.0);
	double deriv = dUsplines (level).d_dy(q,1.0,1.0);
	return (val + deriv*(y-1.0));
      }
    }
  }
  else {
    double r = q+0.5*z;
    double rp = q-0.5*z;
    return (0.5*(Pot->V(r)+Pot->V(rp)));
  }
}



// double PAtricubicFit2Class::dU(double q, double z, double s2, int level)
// {
//   z = fabs(z);
//   double qmax = qgrid->End;
//   double smax = min (2.0*q, sMax(level));
//   double zmax = smax;
//   double smin = z;
//   double s=sqrt(s2);
//   double x;

//   if ((q<=qmax)&&(z<=zmax)&&(s<=smax)) {
//     if (q == 0.0) 
//       return(dUsplines(level)(0.0,0.0,0.0));
//     else {
//       double y = z/zmax;
//       double t = (s-z)/(smax-z);
//       return(dUsplines(level)(q,y,t));
//     }
//   }
//   else {
//     double r = q+0.5*z;
//     double rp = q-0.5*z;
//     return (0.5*(Pot->V(r)+Pot->V(rp)));
//   }
// }





bool PAtricubicFit2Class::Read (IOSectionClass &in,
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
  Pot = ReadPotential(in);
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
    assert(in.ReadVar("Umat", temp));
    Usplines(betaIndex).Init(qgrid,ygrid,tgrid,temp);
    assert(in.ReadVar("dUmat", temp));
    dUsplines(betaIndex).Init(qgrid,ygrid,tgrid,temp);
    // Read sMax
    assert(in.ReadVar("sMax", sMax(betaIndex)));
    sMaxInv(betaIndex) = 1.0/sMax(betaIndex);
    in.CloseSection(); // "Fit"
    desiredBeta *= 2.0;
  }
  in.CloseSection(); // "Fits"
  return true;
}



bool PAtricubicFit2Class::IsLongRange()
{
  return true;
}


// double PAtricubicFit2Class::Vlong_k(double boxVol, double k, int level)
// {
//   if (k <= 0.0)
//     k = 1.0e-30;
//   return 4.0*M_PI/(boxVol*k*k)*exp(-k*k/(4.0*alpha*alpha));
// }

// double PAtricubicFit2Class::dVlong(double q, int level)
// {
//   return 0.0;
// }


// void PAtricubicFit2Class::DoBreakup(const dVec &box, const Array<dVec,1> &kVecs)
// {
//   // Calculate the cutoff parameter
//   double minL, boxVol;
//   boxVol = minL = box[0];
//   for (int i=1; i<NDIM; i++) {
//     minL = min(minL, box[i]);
//     boxVol *= box[i];
//   }
//   alpha = 7.0/minL;

//   // First, subtract off long-ranged part from the bicubic spline to
//   // make them short-ranged.
//   for (int level=0; level<NumBetas; level++) {
//     double tau = SmallestBeta * pow(2.0, level);
//     for (int qi=0; qi<Usplines(level).Nx; qi++) {
//       double q = (*Usplines(level).Xgrid)(qi);
//       double v = Vlong(q, level);
//       for (int ti=0; ti<Usplines(level).Ny; ti++) {
// 	Usplines(level)(qi,ti) -= tau*v;
// 	dUsplines(level)(qi,ti) -= v;
//       }
//     }
//   }


//   // Now, calculate the k-space parts
//   Ulong_k.resize(NumBetas, kVecs.size());
//   dUlong_k.resize(NumBetas, kVecs.size());
//   U_RPA_long_k.resize(NumBetas, kVecs.size());
//   dU_RPA_long_k.resize(NumBetas, kVecs.size());
//   for (int level=0; level<NumBetas; level++) {
//     double tau = SmallestBeta * pow(2.0, level);
//     Ulong_0(level)  = tau*Vlong(0.0, level);
//     dUlong_0(level) = Vlong(0.0, level);
//     for (int ki=0; ki<kVecs.size(); ki++) {
//       double k = sqrt(dot(kVecs(ki), kVecs(ki)));
//       Ulong_k = tau*Vlong_k(boxVol, k, level);
//       dUlong_k = Vlong_k(boxVol, k, level);
//     }
//   }
// }

double PAtricubicFit2Class::Udiag (double q, int level)
{
  if (q <= (qgrid->End*1.0000001)) 
    return Usplines(level)(q, 0.0, 0.0);
  else {
    double beta = SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    return (beta*Pot->V(q));
  }
}

double PAtricubicFit2Class::Udiag_p (double q, int level)
{
  if (q <= (qgrid->End*1.0000001)) 
    return Usplines(level).d_dx(q, 0.0, 0.0);
  else {
    // Coulomb action is independent of z
    double beta = SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    return (beta*Pot->dVdr(q));
  }
}

double PAtricubicFit2Class::Udiag_pp (double q, int level)
{
  if (q <= (qgrid->End*1.0000001)) 
    return Usplines(level).d2_dx2(q, 0.0, 0.0);
  else {
    // Coulomb action is independent of z
    double beta = SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    return (beta*Pot->d2Vdr2(q));
  }
}



double PAtricubicFit2Class::dUdiag (double q, int level)
{
  if (q <= (qgrid->End*1.0000001)) 
    return dUsplines(level)(q, 0.0, 0.0);
  else // Coulomb action is independent of z
    return (Pot->V(q));
}

double PAtricubicFit2Class::dUdiag_p (double q, int level)
{
  if (q <= (qgrid->End*1.0000001)) 
    return dUsplines(level).d_dx(q, 0.0, 0.0);
  else // Coulomb action is independent of z
    return (Pot->dVdr(q));
}

double PAtricubicFit2Class::dUdiag_pp (double q, int level)
{
  if (q <= (qgrid->End*1.0000001)) 
    return dUsplines(level).d2_dx2(q, 0.0, 0.0);
  else  // Coulomb action is independent of z
    return (Pot->d2Vdr2(q));
}


double PAtricubicFit2Class::Xk_U (double k, int level) 
{
  double C0 = -4.0*M_PI/(k*k) * cos(k*rcut);
  double C1 =  4.0*M_PI/k * (gsl_sf_Si(k*rcut) - 0.5*M_PI);
  double C2 =  4.0*M_PI/k * (k*gsl_sf_Ci(k*rcut) - sin(k*rcut)/rcut);

  return (Ucoefs(0,level)*C0 + Ucoefs(1,level)*C1 + Ucoefs(2,level)*C2);
}

double PAtricubicFit2Class::Xk_dU (double k, int level) 
{
  double C0 = -4.0*M_PI/(k*k) * cos(k*rcut);
  double C1 =  4.0*M_PI/k * (gsl_sf_Si(k*rcut) - 0.5*M_PI);
  double C2 =  4.0*M_PI/k * (k*gsl_sf_Ci(k*rcut) - sin(k*rcut)/rcut);

  return (dUcoefs(0,level)*C0 + dUcoefs(1,level)*C1 + dUcoefs(2,level)*C2);
}

double PAtricubicFit2Class::Xk_V (double k)
{
  double C0 = -4.0*M_PI/(k*k) * cos(k*rcut);

  return (Z1Z2*C0); 

}

/// HACK HACK HACK HACK -- this we need to fix!
double PAtricubicFit2Class::Vk (double k)
{
  double C0 = 4.0*M_PI/(k*k);
  return (Z1Z2*C0); 
}


void PAtricubicFit2Class::Setrc (double rc)
{
  rcut = rc;
  const int numPoints=1000;

  int numLevels = Usplines.size();
  
  Array<double,1> coefs, error;
  Array<double,1> y(numPoints), sigma(numPoints);
  Array<double,2> F(numPoints,3);
  double deltar = (qgrid->End - rcut)/(numPoints-1);
  
  // Setup basis functions:
  for (int i=0; i<numPoints; i++) {
    double r = rcut + i*deltar;
    F(i,0) = 1.0/r;
    F(i,1) = 1.0/(r*r);
    F(i,2) = 1.0/(r*r*r);
    sigma(i) = 1.0;
  }

  // Do U fits;
  Ucoefs.resize(3, numLevels);
  for (int level=0; level<numLevels; level++) {
    for (int i=0; i<numPoints; i++) {
      double r = rcut + i*deltar;
      y(i) = Udiag (r, level);
    }
    LinFitSVD (y, sigma, F, coefs, error, 1.0e-15);
    for (int i=0; i<coefs.size(); i++)
      Ucoefs(i, level) = coefs(i);
    //    cerr << "Ucoefs = " << coefs << endl;
  }

  // Do dU fits;
  dUcoefs.resize(3, numLevels);
  for (int level=0; level<numLevels; level++) {
    for (int i=0; i<numPoints; i++) {
      double r = rcut + i*deltar;
      y(i) = dUdiag (r, level);
    }
    LinFitSVD (y, sigma, F, coefs, error, 1.0e-15);
    for (int i=0; i<coefs.size(); i++)
      dUcoefs(i, level) = coefs(i);
  }

  // Do V fits;
  Vcoefs.resize(3);
  for (int i=0; i<numPoints; i++) {
    double r = rcut + i*deltar;
    y(i) = V(r);
  }
  LinFitSVD (y, sigma, F, coefs, error, 1.0e-15);
  for (int i=0; i<coefs.size(); i++)
    Vcoefs(i) = coefs(i);
  Vcoefs = 0.0;
  Vcoefs(0) = Z1Z2; Vcoefs(1) = 0.0; Vcoefs(2) = 0.0;
}
