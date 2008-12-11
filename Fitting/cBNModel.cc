#include "cBNModel.h"

double
cBNModel::FindV(double P, double T)
{
  double Vlo = 10.0;
  double Vhi = 150.0;
  double au2GPa = 29421.01;

  while ((Vhi - Vlo) > 1.0e-10) {
    double Vtry = 0.5*(Vlo + Vhi);
    double Ptry = Static.Pressure(Vtry) +
      au2GPa * Phonon.P(Vtry,T);
    if (Ptry > P)
      Vlo = Vtry;
    else
      Vhi = Vtry;
  }
  return 0.5*(Vlo+Vhi);
}



void
cBNModel::SetStatic (string fname)
{
  TinyVector<double,4> params(80.0, 0.51, 0.0155, 3.7), errors;  
  Static.SetParams (params);
  NonlinearFitClass<4,VinetEOSClass> fitter(Static);
  Array<double,1> V, E, sigma;
  
  IOSectionClass in;
  assert(in.OpenFile (fname));
  assert(in.ReadVar("V0",  params[0]));
  assert(in.ReadVar("E0",  params[1]));
  assert(in.ReadVar("B0",  params[2]));
  assert(in.ReadVar("B0p", params[3]));
  params[2] /= 29421.01;
	
  assert(in.ReadVar("V", V));
  assert(in.ReadVar("E", E));
  assert(in.ReadVar("Sigma", sigma));
  cerr << "V = " << V << endl;
  cerr << "E = " << E << endl;
  cerr << "Sigma = " << sigma << endl;

  fprintf (stdout, "V0  = %12.8f\n", params[0]);
  fprintf (stdout, "Ec  = %12.8f\n", params[1]);
  fprintf (stdout, "B0  = %12.8f GPa\n", params[2] * 29421.01);
  Static.SetParams (params);
  fprintf (stdout, "B0p = %12.8f\n", params[3]);
  fitter.Fit (V, E, sigma, params);

  Array<double,2> &covar = fitter.GetCovariance();
  errors[0] = sqrt (covar(0,0));
  errors[1] = sqrt (covar(1,1));
  errors[2] = sqrt (covar(2,2));
  errors[3] = sqrt (covar(3,3));
  fprintf (stdout, "V0  = %12.8f +/- %12.8f\n", params[0], errors[0]);
  fprintf (stdout, "E0  = %12.8f +/- %12.8f\n", params[1], errors[1]);
  fprintf (stdout, "B0  = %12.8f +/- %12.8f\n", 
	   params[2]* 29421.01, errors[2]* 29421.01);
  Static.SetParams(params);
  fprintf (stdout, "B0p = %12.8f +/- %12.8f\n", params[3], errors[3]);
  fprintf (stdout, "Delta = %1.5f\n", 0.0);
  double V0 = params[0];
  double latConst = cbrt(4.0*V0);
  double latConstError = cbrt(4.0)/(3.0*cbrt(V0)*cbrt(V0)) * errors[0];
  fprintf (stdout, "Lattice const. = %1.6f +/- %1.6f bohr radii\n",
	   latConst, latConstError);
  fprintf (stdout, "               = %1.6f +/- %1.6f angstrom\n",  
	   0.52917721*latConst, 0.52917721*latConstError);
}


void
cBNModel::SetPhonon (string fname)
{
  FILE *fin = fopen (fname.c_str(), "r");
  assert (fin != NULL);
  std::vector<double> F_T, U_T;
  std::vector<std::vector<double> > F_VT, U_VT;
  std::vector<double> V, T;
  bool done = false;

  double lastv = 0.0;
  bool first = true;

  double v,t,f,u;
  while (!done) {
    int retval = fscanf (fin, "%lf %lf %lf %lf", &v, &t, &f, &u);
    if (retval == 4) {
      if (v != lastv) {
	if (first) {
	  first = false;
	}
	else {
	  F_VT.push_back (F_T);
	  U_VT.push_back (U_T);
	  V.push_back(lastv);
	  F_T.clear();
	  U_T.clear();
	  T.clear();
	}
      }
      T.push_back (t);
      F_T.push_back (f);
      U_T.push_back (u);
      lastv = v;
    }
    else
      done = true;
  }
  F_VT.push_back(F_T);
  U_VT.push_back(U_T);
  fclose (fin);
  
  int numV = V.size();
  int numT = T.size();
  fprintf (stderr, "numT = %d\n", numT);
  double startV = V[0];
  double endV = V[V.size()-1];
  LinearGrid *Vgrid = new LinearGrid (startV, endV, numV);
  double startT = T[0];
  double endT   = T[T.size()-1];
  LinearGrid *Tgrid = new LinearGrid (startT, endT, numT);

  Array<double,2> F (numV, numT);
  Array<double,2> U (numV, numT);
  for (int iV=0; iV<numV; iV++) {
    //fprintf (stderr, "V = %1.2f\n", (*Vgrid)(iV));
    for (int iT=0; iT<numT; iT++) {
      //    fprintf (stderr, "T = %1.2f\n", (*Tgrid)(iV));
      F(iV, iT) = F_VT[iV][iT];
      U(iV, iT) = U_VT[iV][iT];
    }
  }
  Phonon.SetF(Vgrid, Tgrid, F);
  Phonon.SetU(Vgrid, Tgrid, U);

  fin = fopen (fname.c_str(), "r");
  lastv = 0;
  done=false;
  double F0 = 0.0;
  double U0 = 0.0;
  FILE *out = fopen ("ModelF2.dat", "w");
  while (!done) {
    double v,t,f,u;
    int retval = fscanf (fin, "%lf %lf %lf %lf", &v, &t, &f, &u);
    if (retval != 4)
      done = true;
    else {
      if (v != lastv)
	F0 = f; U0 = u;
      fprintf (out, "%1.9f  %5.1f   %12.14f  %12.14f\n",
	       v, t, Phonon.F(v,t), Phonon.U(v,t));
      lastv = v;
    }
  }
  fclose(out);
  fclose (fin);

  out = fopen ("TestF.dat", "w");
  for (double V=40.0; V<95.0; V+=0.1) {
    for (double T=0.0; T<2500.0; T+= 3.32342341)
      fprintf (out, "%1.12e ", Phonon.F(V,T));
    fprintf (out, "\n");
  }
  fclose(out);
}


main(int argc, char **argv)
{
  cBNModel model;
  if (argc > 2) {
    model.SetStatic (argv[1]);
    model.SetPhonon (argv[2]);
    FILE *fout = fopen ("300K_Model.dat", "w");
    for (double V=40.0; V<=95.0; V+=0.01)
      fprintf (fout, "%5.2f %10.6f\n", V, model.P(V, 300.0));
    fclose(fout);
  }
}
