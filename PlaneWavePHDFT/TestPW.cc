#include "ConjGrad.h"

main()
{
  Vec3 box(10.0, 11.0, 12.0);
  IOSectionClass in;
  in.OpenFile ("NaUnscreenedPH_Feb18_05.h5");
    //in.OpenFile ("NaPH_US_March1_05b.h5");
  Potential *ph = ReadPotential(in);
  in.CloseFile();
  Hamiltonian H(box, 4.0, 1.0, *ph);
  H.PHFFT.Setup();
  int numBands = 6;
  ConjGrad CG(H,numBands);
  clock_t start, end;
  start = clock();
  for (int band=0; band<numBands; band++) 
    for (int j=0; j<25; j++) {
      CG.Iterate(band);
      fprintf (stderr, "Energy(%d) = %14.10f\n", band, CG.Energies(band));
    }
//   for (int band=0; band<numBands; band++)
//     cerr << "Energies = ";
//   for (int band=0; band<numBands; band++)
//     fprintf (stderr, "%10.6f ", CG.Energies(band));
//   cerr << endl;
//   }
  end = clock();
  CG.PrintOverlaps();

  FILE *fout = fopen ("HRho_k10.0.dat", "w");
  Vec3 r(0.0, 0.0, 0.0);
  for (double x=-0.5*box[0]; x<=0.5*box[0]; x+=0.01) {
    r[0] = x;
    fprintf (fout, "%1.12e ", x);
    for (int band=0; band<numBands; band++) {
      complex<double> psi(0.0, 0.0);
      for (int i=0; i<H.GVecs.size(); i++) {
	double phase = dot (H.GVecs(i), r);
	double s,c;
	sincos (phase,&s, &c);
	complex<double> z(c, s);
	psi += z * CG.Bands(band,i);
      }
    
      psi /= sqrt(H.GVecs.GetBoxVol());
      fprintf (fout, "%1.12e ", real(conj(psi)*psi));
    }
    fprintf (fout, "\n");
  }
  fclose (fout);
  

  fprintf (stderr, "Time = %1.3f\n", 
	   (double)(end-start)/(double)CLOCKS_PER_SEC);

}
