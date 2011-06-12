#include <stdio.h>

main()
{

  for (int ix=0; ix<8; ix++) {
    double x = 0.125*ix + 0.0625;
    char fname[100];
    snprintf (fname, 100, "twists512_%d.in", ix);
    FILE *fout = fopen (fname, "w");
    
    fprintf (fout, "  string Type = \"FIXEDPHASE\";\n");
    fprintf (fout, "  string IonSpecies  = \"Na\";\n");
    fprintf (fout, "  string UpSpecies   = \"eup\";\n");
    fprintf (fout, "  string DownSpecies = \"edown\";\n");
    fprintf (fout, "  double kCut = 6.0;\n");
    fprintf (fout, "  bool UseMDExtrap = true;\n");
    fprintf (fout, "  bool UseLDA = true;\n");
    fprintf (fout, "  bool UseSubspaceRotation = true;\n");
    fprintf (fout, "  int SmearOrder = 2;\n");
    fprintf (fout, "  double SmearWidth = 0.015; \n\n");

    fprintf (fout, "  Array<double,2> TwistAngles(64,3) = [ ");
    for (int iy=0; iy<8; iy++) {
      double y = 0.125*iy + 0.0625;
      for (int iz=0; iz<8; iz++) {
	double z = 0.125*iz + 0.0625;
	if ((iy != 0) || (iz != 0))
	  fprintf (fout, "                                        ");
	fprintf (fout, "%1.4f, %1.4f, %1.4f", x, y, z);
	if (!((iy==7)&&(iz==7)))
	  fprintf (fout, ",\n");
      }
    }
    fprintf (fout, " ];\n");
    fclose (fout);
  }
}
      
