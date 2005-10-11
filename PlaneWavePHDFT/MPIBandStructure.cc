#include "MPIBandStructure.h"

void
MPIBandStructureClass::Read(IOSectionClass &in)
{
  Communicator.SetWorld();
  // Read the pseudoHamiltonian
  assert (in.OpenSection ("Potential"));
  PH = ReadPotential (in);
  in.CloseSection();

  // Read the box
  Array<double,1> box;
  assert (in.ReadVar ("Box", box));
  Box[0]=box(0);   Box[1]=box(1);    Box[2]=box(2);
  
  // Read ion positions
  Array<double,2> tmpIons;
  assert (in.ReadVar ("Positions", tmpIons));
  Rions.resize(tmpIons.rows());
  for (int i=0; i<Rions.size(); i++)
    for (int j=0; j<3; j++)
      Rions(i)[j] = tmpIons(i,j);

  // Read the number of bands
  assert (in.ReadVar ("NumBands", NumBands));
	  
  // Read the k-point list
  Array<double,2> tmpkPoints;
  assert (in.ReadVar("kPoints", tmpkPoints));
  kPoints.resize(tmpkPoints.rows());
  for (int i=0; i<kPoints.size(); i++)
    for (int j=0; j<3; j++)
      kPoints(i)[j] = tmpkPoints(i,j);

  // Read the number of interpolation points per kPoint
  assert (in.ReadVar ("InterpPoints", InterpPoints));

  // Read the k-cutoff
  assert(in.ReadVar("kCut", kCut));
  System = new MPISystemClass(NumBands, Communicator);

  // Read the output file name
  assert (in.ReadVar ("OutFilename", OutFilename));

  // Setup the plane-wave system
  System->Setup (Box, kPoints(0), kCut, *PH, true);
  System->SetIons(Rions);


};


void
MPIBandStructureClass::CalcBands()
{
  FILE *fout = fopen (OutFilename.c_str(), "w");
  assert (fout != NULL);
 
  for (int ki=0; ki<kPoints.size()-1; ki++) 
    for (int i=0; i<InterpPoints; i++) {
      double alpha = (double)i/(double)InterpPoints;
      Vec3 k = (1.0-alpha)*kPoints(ki) + alpha*kPoints(ki+1);
      fprintf (fout, "%1.6e %1.6e %1.6e ", k[0], k[1], k[2]);
      System->Setk(k);
      System->DiagonalizeH();
      for (int band=0; band<NumBands; band++)
	fprintf (fout, "%1.16e ", System->GetEnergy(band));
      fprintf (fout, "\n");
      fflush (fout);
    }
  fclose(fout);
}


main(int argc, char **argv)
{
  COMM::Init(argc, argv);
  if (argc < 2) {
    cout << "Usage:\n";
    cout << "BandStructure myfile.in\n"; 
  }
  else {
    IOSectionClass in;
    assert (in.OpenFile(argv[1]));
    MPIBandStructureClass BandStructure;
    BandStructure.Read(in);
    in.CloseFile();
    BandStructure.CalcBands();
  }
}
