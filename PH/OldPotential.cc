#include "PH.h"

//////////////////////////////////////////////////////////////////
//                FullCorePotential Routines                    //
//////////////////////////////////////////////////////////////////

void
FullCorePotential::Read(char *FName)
{
  strncpy (FileName, FName, 1000);
  FILE *fin;
  if ((fin = fopen(FileName, "r")) == NULL)
    {
      cerr << "Cannot open " << FileName << " for reading. Exitting.";
      exit(1);
    }

  fscanf (fin, " %lf ", &Z);
  Grid *TempGrid;
  TempGrid = ReadGrid(fin);
  if (GridInitialized)
    delete(V.grid);
  GridInitialized = 1;
  scalar dummy;
  Array<scalar, 1> PotVals(TempGrid->NumPoints);
  for (int i=0; i<TempGrid->NumPoints; i++)
    fscanf(fin, " %lf %lf ", &PotVals(i), &dummy);
  V.Init(TempGrid, PotVals, 5.0e30, 5.0e30);
  fclose(fin);
}
