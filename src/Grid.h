#ifndef GRID_H
#define GRID_H

#include "Blitz.h"
#include "InputFile.h"

//Ken's Grid Class

/// The different types of grids that we currently allow
typedef enum {NONE, LINEAR, OPTIMAL, LOG} GridType;


/// Parent class for all grids
class Grid
{
 protected:
  /// Contains the grid points 
  Array<scalar,1> grid;
 public:
  /// First and last grid points
  scalar Start, End;

  /// Number of points in the grid
  int NumPoints;

  /// The i'th point in the grid
  inline scalar operator()(int i) const
    {
      return (grid(i));
    }


  /// Returns the type of the grid (i.e. linear, optimal, etc)
  virtual GridType Type()
    { return (NONE); }

  /// Reads the parameters in the file.
  virtual void ReadParams (FILE *fin)
  {
    cerr << "Should never get here.\n";
    exit(1);
  }
   
  ///Write the parameters in the file
  virtual void WriteParams (FILE *fout)
  {
    cerr << "Should never get here.\n";
    exit(1);
  }

 
  /// Writes the grid paramaters in inputfile format
  virtual void WriteInput (FILE *fout)
  {
    cerr << "Should never get here.\n";
    exit(1);
  }

  ///Returns the index of the nearest point below r. 
  virtual int ReverseMap (scalar r)
  {
    cerr << "Should never get here.\n";
    exit(1);
    return (0);
  }



};


/// Linear Grid inherets from Grid.  
class LinearGrid : public Grid
{
 private:
  /// The value between successive grid points.
  scalar delta;
 public:
  /// Returns the type of the grid (in this case LINEAR)
  GridType Type()
    { return (LINEAR); }

  /// Reads the paramaters from the file.
  void ReadParams (FILE *fin)
    {
      fscanf (fin, " %lf  %lf %d ", &Start, &End, &NumPoints);
      grid.resize(NumPoints);
      delta = (End-Start)/(scalar)(NumPoints-1);
      for (int i=0; i<NumPoints; i++)
	grid(i) = Start + (scalar)i*delta;
    }
  ///Writes the paramaters to the file.
  void WriteParams (FILE *fout)
    {
      fprintf (fout, "%1.16e %1.16e %d\n", Start, End, NumPoints);
    }

  /// Writes the grid paramaters in inputfile format  
  void WriteInput (FILE *fout)
  {
    fprintf (fout, "  {\n");
    fprintf (fout, "    Type = \"LINEAR\";\n");
    fprintf (fout, "    Start = %1.16e;\n", Start);
    fprintf (fout, "    End = %1.16e;\n", End);
    fprintf (fout, "    NumPoints = %d;\n", NumPoints);
    fprintf (fout, "  }\n");
  }

  /// Returns the index of the nearest point below r. 
  int ReverseMap(scalar r)
    {
      return ((int)floor((r-Start)/delta));
    }

  /// Initializes the linear grid.
  inline void Init(scalar start, scalar end, int numpoints)
  {
      Start=start; End=end; NumPoints=numpoints;
      grid.resize(NumPoints);
      delta = (End-Start)/(scalar)(NumPoints-1);
      for (int i=0; i<NumPoints; i++)
	grid(i) = Start + (scalar)i*delta;
  }

  /// Useless constructor
  LinearGrid ()
  {
    // Do nothing
  }

  /// Constructor that sets the number of points, start and end point
  /// of the original grid 
  LinearGrid (scalar start, scalar end, int numpoints)
    {
      Init (start, end, numpoints);
    }
};



class OptimalGrid : public Grid
{
 private:
  scalar a, b;
 public:

  GridType Type()
    { return (OPTIMAL); }

  void ReadParams (FILE *fin)
    {
      fscanf (fin, " %lf  %lf %d ", &a, &b, &NumPoints);

      Start = a * (exp(b) - 1.0);
      End   = a * (exp(b*NumPoints) - 1.0);
      
      grid.resize(NumPoints);
      
      for (int i=0; i<NumPoints; i++)
	grid(i) = a*(exp(b*(i+1))-1.0);
    }

  void WriteParams (FILE *fout)
    {
      fprintf (fout, "%1.16e %1.16e %d\n", a, b, NumPoints);
    }


  void WriteInput (FILE *fout)
  {
    fprintf (fout, "  {\n");
    fprintf (fout, "    Type = \"OPTIMAL\";\n");
    fprintf (fout, "    a = %1.16e;\n", a);
    fprintf (fout, "    b = %1.16e;\n", b);
    fprintf (fout, "    NumPoints = %d;\n", NumPoints);
    fprintf (fout, "  }\n");
  }

  int ReverseMap(scalar r)
    {
      if ((r/a) < 1e-6)
	return ((int)floor(r/(a*b)+0.5)-1);
      else
	return((int)floor(log(r/a + 1.0)/b + 0.5) -1);
    }

  double Geta() const
    { return (a); }

  double Getb() const
    { return (b); }

  OptimalGrid ()
  {
    // Do nothing
  }


  OptimalGrid (int numpoints, scalar rmax, scalar bval)
  {
    NumPoints = numpoints;
    b = bval;
    End = rmax;
    a = End / (exp(b*(scalar)NumPoints) - 1.0);  
    Start = a * (exp(b) - 1.0);
    grid.resize(NumPoints);
      
    for (int i=0; i<NumPoints; i++)
      grid(i) = a*(exp(b*(i+1))-1.0);
  }

  OptimalGrid (scalar aval, scalar bval, int numpoints)
  {
    a = aval;
    b = bval;
    NumPoints = numpoints;
    Start = a * (exp(b) - 1.0);
    End   = a * (exp(b*NumPoints) - 1.0);
      
    grid.resize(NumPoints);
      
    for (int i=0; i<NumPoints; i++)
      grid(i) = a*(exp(b*(i+1))-1.0);
  }

  OptimalGrid (scalar Z, scalar rmax)
    {
      a = 4.34e-6/Z;
      //a = 4.0e-2;
      b = 0.002304;
      //b = 0.004;

      NumPoints = (int)ceil(log(rmax/a+1.0)/b);
      cerr << "NumPoints = " << NumPoints << "\n";
      b = log(rmax/a+1.0)/(scalar)NumPoints;
      Start = a * (exp(b) - 1.0);
      End = rmax;
      //End   = a * (exp(b*NumPoints) - 1.0);
      
      grid.resize(NumPoints);
      
      for (int i=0; i<NumPoints; i++)
	{
	  grid(i) = a*(exp(b*(i+1))-1.0);
	  //fprintf (stdout, "%1.12e\n", grid(i));
	}
    }
};





class LogGrid : public Grid
{
 public:
  scalar Z, r0, Spacing;

  GridType Type()
    { return (LOG); }

  void ReadParams (FILE *fin)
    {
      fscanf (fin, " %d %lf %lf %lf ", &NumPoints, &Z, &r0, &Spacing);

      Start = r0 / Z;
      End = r0/Z * pow(Spacing, (scalar) (NumPoints-1));
      
      grid.resize(NumPoints);
      
      for (int i=0; i<NumPoints; i++)
	grid(i) = r0/Z * pow (Spacing, (scalar) i);
    }

  void WriteParams (FILE *fout)
    {
      fprintf (fout, " %d %1.16e %1.16e %d\n", NumPoints, Z, r0, Spacing);
    }

  void WriteInput (FILE *fout)
  {
    fprintf (fout, "  {\n");
    fprintf (fout, "    Type = \"LOG\";\n");
    fprintf (fout, "    r0 = %1.16e;\n", r0/Z);
    fprintf (fout, "    Spacing = %1.16e;\n", Spacing);
    fprintf (fout, "    NumPoints = %d;\n", NumPoints);
    fprintf (fout, "  }\n");
  }

  int ReverseMap(scalar r)
    {
      return ((int)(floor(log(Z*r/r0)/log(Spacing))));
    }

  LogGrid ()
  {
    // Do nothing
  }

  LogGrid (scalar R0, scalar spacing, int numpoints)
  {
    NumPoints = numpoints;
    Z = 1.0; r0 = R0; Spacing = spacing;
    Start = r0;
    End = r0 * pow(Spacing, (scalar) NumPoints-1);
    grid.resize (NumPoints);
    
    for (int i=0; i<NumPoints; i++)
      grid(i) = r0 * pow(Spacing, (scalar) i);
  }


  LogGrid (int numpoints, scalar z, scalar R0, scalar spacing)
    {
      NumPoints = numpoints;
      Z = z; r0 = R0; Spacing = spacing;

      Start = r0 / Z;
      End = r0/Z * pow(Spacing, (scalar) (NumPoints-1));
      
      grid.resize(NumPoints);
      
      for (int i=0; i<NumPoints; i++)
	grid(i) = r0/Z * pow (Spacing, (scalar) i);
    }
};

Grid *ReadGrid (InputBuffer &SectionBuf);

//  class FlexGrid
//  {
//  private:
//    int Initialized;
//  public:
//    Grid *grid;
  
//    inline void Read(FILE *fin)
//    {
//      char TypeString[500];
//      fgets(TypeString, 500, fin);
//      TypeString[strlen(TypeString)-1] = '\0';

//      if (!strcmp(TypeString, "LINEAR"))
//        grid = new LinearGrid;
//      else if (!strcmp(TypeString, "OPTIMAL"))
//        grid = new OptimalGrid;
//      else
//        {
//  	cerr << "Error: Unrecognized grid type " << TypeString << "\n";
//  	exit (1);
//        }
//      grid.ReadParams(fin);
//      Initialized = 1;
//    }
//    inline scalar operator()(int i)
//    {
//      if (!Initialized)
//        {
//  	cerr << "Grid not initilized in FlexGrid::operator().\n";
//  	exit (1);
//        }
//      return ((*grid)(i));
//    }
//    inline int ReverseMap(scalar r)
//    {
//      if (!Initialized)
//        {
//  	cerr << "Grid not initilized in FlexGrid::operator().\n";
//  	exit (1);
//        }
//      return (grid->ReverseMap(r));
//    }

//    inline FlexGrid()
//    {
//      Initialized = 0;
//    }
//    inline FlexGrid~()
//    {
//      if (Initialized)
//        delete grid;
//    }
//  };


inline Grid *ReadGrid(FILE *fin)
{
  Grid *grid;
  char TypeString[500];
  fgets(TypeString, 500, fin);
  TypeString[strlen(TypeString)-1] = '\0';
  
  if (!strcmp(TypeString, "LINEAR"))
    grid = new LinearGrid;
  else if (!strcmp(TypeString, "OPTIMAL"))
    grid = new OptimalGrid;
  else
    {
      cerr << "Error: Unrecognized grid type " << TypeString << "\n";
      exit (1);
    }
  grid->ReadParams(fin); 
  return (grid);
}

Grid *ReadGrid (InputBuffer &SectionBuf);

#endif
