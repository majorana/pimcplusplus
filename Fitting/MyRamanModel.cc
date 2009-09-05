#include "GoncharovRamanModel.h"
#include "MyRamanModel.h"
#include "NonlinearFitTemp.h"
#include "../Random/Random.h"
#include <iostream>
#include <fstream>
#include <strstream>
#include <cstdio>
#include <vector>


void FitRamanError(string fname)
{
  ifstream fin;
  fin.open(fname.c_str());
  vector<double> Pvec, Tvec, nuvec, sigmavec;


  while (!fin.eof()) {
    char line[500];
    fin.getline(line,500);
    if (line[0] == '#') continue;

    istrstream sin(line);
    double V, T, P, nu, sigma;
    sin >> V >> T >> P >> nu >> sigma;
    Pvec.push_back(P);
    Tvec.push_back(T);
    nuvec.push_back(nu);
    sigmavec.push_back(sigma);
  }
  int N = Pvec.size();
  Array<double,1> nu(N), sigma(N);
  Array<TinyVector<double,2>,1> PT(N);
  for (int i=0; i<N; i++) {
    nu(i) = nuvec[i];
    sigma(i) = sigmavec[i];
    PT(i) = TinyVector<double,2>(Pvec[i],Tvec[i]);
  }

  const int numSamples = 100;
  Array<double,1> nuSum(1001), nu2Sum(1001);
  nuSum = 0.0;
  nu2Sum = 0.0;

  GoncharovRamanModel GoncharovModel, DatchiModel;
  MyRamanModel model;
  TinyVector<double,9> Gparams (3.48, -1.75e-4, -8.55e-8,
				325.6, -1.11e-2, -8.9e-6,
				1060.2, -0.012, -1.57e-5);
  TinyVector<double,9> Dparams (1.9597, 0.0, 0.0,
				452.63, -0.035640, 0.0,
				1058.3, -0.0100, -1.42e-5);
  GoncharovModel.SetParams(Gparams);
  DatchiModel.SetParams(Dparams);


  CommunicatorClass comm;
  RandomClass rand(comm);
  rand.Init();
  for (int sample=0; sample < numSamples; sample++) {
    Array<double,1> nuPert(nu.size());
    for (int i=0; i<nu.size(); i++)
      nuPert(i) = nu(i) + rand.LocalGaussian(sigma(i));

    TinyVector<double,7> myparams (2.995, 1e-3, 1e-4,
				   325.6, -1.11e-4, -8.9e-8, 
				   1062.3, -100.0, 4000.0);
    
    model.SetParams(myparams);
    
    NonlinearFitClass<7, MyRamanModel, TinyVector<double,2> > fitter(model);
    TinyVector<double,9> Goncharov_params;
    TinyVector<double,7> params;
    
    model.SetParams(myparams);
    params = myparams;
    cerr << "Before fitter.Fit.\n";
    fitter.Fit (PT, nuPert, sigma, params);
    cerr << "After fitter.Fit.\n";
    model.SetParams(params);
    for (int iP=0; iP<=1000; iP++) {
      double P = (double)iP;
      TinyVector<double,2> PT(P, 300.0);
      double my_nu = model (PT);
      nuSum(iP) += my_nu;
      nu2Sum(iP) += my_nu * my_nu;
    }
  }

  // Write out results
  FILE *fout = fopen ("AllRamanModels.dat", "w");
  fprintf (fout, "# Pressure   nu_Datchi   nu_Goncharov  nu_me  sigma\n");
  for (int iP=0; iP<=1000; iP++) {
    double P = (double)iP;
    TinyVector<double,2> PT(P, 300.0);
    double my_nu = nuSum(iP) / (double)numSamples;
    double sigma = std::sqrt(nu2Sum(iP) / (double)numSamples - my_nu * my_nu);
    PT[0] = P;
    double Goncharov_nu  = GoncharovModel(PT);
    double Datchi_nu  = DatchiModel(PT);
    fprintf (fout, "%10.4f  %12.4f %12.4f %12.4f %12.14f\n",
	     P, Datchi_nu, Goncharov_nu, my_nu, sigma);
  }
  fclose (fout); 
}


void FitRaman(string fname)
{
  ifstream fin;
  fin.open(fname.c_str());
  vector<double> Pvec, Tvec, nuvec, sigmavec;


  while (!fin.eof()) {
    char line[500];
    fin.getline(line,500);
    if (line[0] == '#') continue;

    istrstream sin(line);
    double V, T, P, nu, sigma;
    sin >> V >> T >> P >> nu >> sigma;
    const double Timport = 1000.0;
    if ((T < 5000.0 && P > -10.0)) {
      Pvec.push_back(P);
      Tvec.push_back(T);
      nuvec.push_back(nu);
      sigmavec.push_back(sigma);
    }
  }
  int N = Pvec.size();
  Array<double,1> nu(N), sigma(N);
  Array<TinyVector<double,2>,1> PT(N);
  for (int i=0; i<N; i++) {
    nu(i) = nuvec[i];
    sigma(i) = sigmavec[i];
    PT(i) = TinyVector<double,2>(Pvec[i],Tvec[i]);
  }

  GoncharovRamanModel GoncharovModel, DatchiModel;
  MyRamanModel model;
  TinyVector<double,9> Gparams (3.48, -1.75e-4, -8.55e-8,
				325.6, -1.11e-2, -8.9e-6,
				1060.2, -0.012, -1.57e-5);
  TinyVector<double,9> Dparams (1.9597, 0.0, 0.0,
				452.63, -0.035640, 0.0,
				1058.3, -0.0100, -1.42e-5);
  TinyVector<double,7> myparams (2.874,
				 395, 200.0, 200.0,   //2771.1,  96.26,
				 1073.58, -154.58, 1400.0);

  GoncharovModel.SetParams(Gparams);

  FILE *fcheck = fopen ("GonFreqs.dat", "w");
  for (double P=0.0; P<1000.0; P++) {
    double nu = GoncharovModel(TinyVector<double,2>(P,300.0));
    fprintf (fcheck, "%6.2f %12.8f\n", P, nu);
  }
  fclose(fcheck);
    
  NonlinearFitClass<7, MyRamanModel, TinyVector<double,2> > fitter(model);
  TinyVector<double,9>  Goncharov_params;
  TinyVector<double,7> params;

  GoncharovModel.SetParams(Gparams);
  DatchiModel.SetParams(Dparams);
  model.SetParams(myparams);
  params = myparams;
  fitter.Fit (PT, nu, sigma, params);
  model.SetParams(params);
  fprintf (stderr, "b0     = %12.5e\n", params[0]);
  //  fprintf (stderr, "b1     = %12.5e\n", params[1]);
  // fprintf (stderr, "b2     = %12.5e\n", params[2]);

  fprintf (stderr, "R0     = %12.5e\n", params[1]);
  fprintf (stderr, "R1     = %12.5e\n", params[2]);
  fprintf (stderr, "R2     = %12.5e\n", params[3]);

  fprintf (stderr, "nu0    = %12.5e\n", params[4]);
  fprintf (stderr, "nu1    = %12.5e\n", params[5]);
  fprintf (stderr, "nu2    = %12.5e\n", params[6]);

  FILE *fout = fopen ("AllRamanModelsNew.dat", "w");
  
  for (double P=0; P<=1000.0; P+=1.0) {
    TinyVector<double,2> PT(P, 300.0);
    double my_nu = model (PT);
    PT[0] = P;
    double Goncharov_nu  = GoncharovModel(PT);
    double Datchi_nu  = DatchiModel(PT);
    fprintf (fout, "%10.4f  %12.4f %12.4f %12.4f\n",
	     P, Datchi_nu, Goncharov_nu, my_nu);
  }
  fclose (fout); 


  cerr << "params = " << params << endl;
}



main(int argc, char **argv)
{
  if (argc > 1) {
    //FitRamanError(argv[1]);
    FitRaman(argv[1]);
  }
}
