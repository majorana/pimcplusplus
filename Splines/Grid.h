#ifndef GRID_H
#define GRID_H

#include "../Blitz.h"
#include "../IO/InputFile.h"
#include "../IO/InputOutput.h"

//Ken's Grid Class

/// The different types of grids that we currently allow
typedef enum {NONE, LINEAR, OPTIMAL, OPTIMAL2, LOG, CLUSTER} GridType;


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
  inline double* data()
  {
    return grid.data();
  }
  inline Array<double,1>& Points()
  {
    return grid;
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

  virtual void Write (IOSectionClass &out) = 0;
  virtual void Read  (IOSectionClass &inSection) = 0;
};


/// Linear Grid inherets from Grid.  
class LinearGrid : public Grid
{
 private:
  /// The value between successive grid points.
  scalar delta, deltainv;
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
      return ((int)floor((r-Start)*deltainv));
    }

  /// Initializes the linear grid.
  inline void Init(scalar start, scalar end, int numpoints)
  {
      Start=start; End=end; NumPoints=numpoints;
      grid.resize(NumPoints);
      delta = (End-Start)/(scalar)(NumPoints-1);
      deltainv = 1.0/delta;
      for (int i=0; i<NumPoints; i++)
	grid(i) = Start + (scalar)i*delta;
  }

  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Points", grid); 
    outSection.WriteVar ("Type", "Linear");
    outSection.WriteVar ("Start", Start);
    outSection.WriteVar ("End", End);
    outSection.WriteVar ("NumPoints", NumPoints);
  }
  void Read (IOSectionClass &inSection)
  {
    assert(inSection.ReadVar("Start", Start));
    assert(inSection.ReadVar("End", End));
    assert(inSection.ReadVar("NumPoints", NumPoints));
    Init (Start, End, NumPoints);
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


/// The OptimalGrid class stores a grid which has linear spacing at
/// the origin and exponential spacing further out.  It has the
/// analytic form \f[r_k = a\left(e^{kb}-1\right)\f].
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

  /// Returns a parameter
  double Geta() const
    { return (a); }

  /// Returns b parameter
  double Getb() const
    { return (b); }

  OptimalGrid ()
  {
    // Do nothing
  }

  /// This form of the constructor takes the number of points, the
  /// maximum radius and the value of b.
  void Init (int numpoints, scalar rmax, scalar bval)
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

  OptimalGrid (int numPoints, scalar rmax, scalar bval)
  { Init (numPoints, rmax, bval); }

  void Init (scalar aval, scalar bval, int numPoints)
  {
    a = aval;
    b = bval;
    NumPoints = numPoints;
    Start = a * (exp(b) - 1.0);
    End   = a * (exp(b*NumPoints) - 1.0);
      
    grid.resize(NumPoints);
      
    for (int i=0; i<NumPoints; i++)
      grid(i) = a*(exp(b*(i+1))-1.0);
  }


  /// This form of the constructor takes a, b, and the number of points.
  OptimalGrid (scalar aval, scalar bval, int numPoints)
  { 
    Init (aval, bval, numPoints);
  }

  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Points", grid); 
    outSection.WriteVar ("Type", "Optimal");
    outSection.WriteVar ("a", a);
    outSection.WriteVar ("b", b);
    outSection.WriteVar ("NumPoints", NumPoints);
  }

  void Read (IOSectionClass &inSection)
  {
    double aval, bval;
    int numPoints;
    if (inSection.ReadVar("a", aval)) {
      assert(inSection.ReadVar("b", bval));
      assert(inSection.ReadVar("NumPoints", numPoints));
      Init (aval,bval,numPoints);
    }
    else {
      double Z, rmax;
      assert(inSection.ReadVar("Z", Z));
      assert(inSection.ReadVar("rmax", rmax));
      Init (Z, rmax);
    }
  }
  
  void InitRatio (double end, double ratio, int numpoints)
  {
    End = end; 
    NumPoints = numpoints;
    
    b = log(ratio)/(double)(numpoints-2);
    a = end/(exp(b*(double)(numpoints-1)) - 1);
      
    grid.resize(NumPoints);
      
    for (int i=0; i<NumPoints; i++)
      grid(i) = a*(exp(b*i)-1.0);
  }



  /// This form of the constructor takes a nuclear charge and a
  /// maxmimum radius and chooses an appropriate number of points for
  /// that atom.
  void Init (scalar Z, scalar rmax)
    {
      a = 4.34e-6/Z;
      //a = 4.0e-2;
      b = 0.002304;
      //b = 0.004;

      NumPoints = (int)ceil(log(rmax/a+1.0)/b);
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
  OptimalGrid (scalar Z, scalar rmax)
  { Init (Z, rmax); }

};


/// The OptimalGrid class stores a grid which has linear spacing at
/// the origin and exponential spacing further out.  It has the
/// analytic form \f[r_k = a\left(e^{kb}-1\right)\f].
class OptimalGrid2 : public Grid
{
 private:
  double a, b, c;
  double Ratio;
 public:

  GridType Type()
    { return (OPTIMAL2); }

  void ReadParams (FILE *fin)
    {    }

  void WriteParams (FILE *fout)
    {    }


  void WriteInput (FILE *fout)
  {  }

  int ReverseMap(scalar r)
    {
      r -= c;
      if ((r/a) < 1e-6)
	return ((int)floor(r/(a*b)+0.5));
      else
	return((int)floor(log(r/a + 1.0)/b + 0.5));
    }

  /// Returns a parameter
  double Geta() const
    { return (a); }

  /// Returns b parameter
  double Getb() const
    { return (b); }

  OptimalGrid2 ()
  {
    // Do nothing
  }

  /// This form of the constructor takes the number of points, the
  /// maximum radius and the value of b.
  OptimalGrid2 (int numpoints, scalar rmax, scalar bval)
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

  void Init (double start, double end, double ratio, int numpoints)
  {
    Start = start;
    End = end; 
    Ratio = ratio;
    NumPoints = numpoints;
    
    b = log(ratio)/(double)(numpoints-2);
    c = Start;
    a = (end - c)/(exp(b*(double)(numpoints-1)) - 1);
      
    grid.resize(NumPoints);
      
    for (int i=0; i<NumPoints; i++)
      grid(i) = c + a*(exp(b*i)-1.0);
  }


  /// This form of the constructor takes a, b, and the number of points.
  OptimalGrid2 (double start, double end, double ratio, int numpoints)
  { 
    Init (start, end, ratio, numpoints);
  }

  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Points", grid); 
    outSection.WriteVar ("Type", "Optimal2");
    outSection.WriteVar ("Start", Start);
    outSection.WriteVar ("End", End);
    outSection.WriteVar ("Ratio", Ratio);
    outSection.WriteVar ("NumPoints", NumPoints);
  }

  void Read (IOSectionClass &inSection)
  {
    double start, end, ratio;
    int numPoints;
    assert(inSection.ReadVar("Start", start));
    assert(inSection.ReadVar("End", end));
    assert(inSection.ReadVar("Ratio", ratio));
    assert(inSection.ReadVar("NumPoints", numPoints));
    Init (start,end,ratio,numPoints);
  }
};


/// LogGrid is a function whose gridpoints increase exponentially with
/// the index.  That is, it has the analytic form
/// \f[ r_k = \frac{r_0}{Z} \Delta^k.\f]  It is appropriate for
/// functions which change rapidly near the origin but vary smoothly
/// further out.
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

  void Init (scalar R0, scalar spacing, int numpoints)
  {
    NumPoints = numpoints;
    Z = 1.0; r0 = R0; Spacing = spacing;
    Start = r0;
    End = r0 * pow(Spacing, (scalar) NumPoints-1);
    grid.resize (NumPoints);
    
    for (int i=0; i<NumPoints; i++)
      grid(i) = r0 * pow(Spacing, (scalar) i);
  }



  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Points", grid); 
    outSection.WriteVar ("Type", "Log");
    outSection.WriteVar ("r0", r0);
    outSection.WriteVar ("Spacing", Spacing);
  }

  void Read (IOSectionClass &inSection)
  {
    double tempr0, tempSpacing;
    int  tempNumPoints;
    assert (inSection.ReadVar("r0", tempr0));
    assert (inSection.ReadVar("Spacing", tempSpacing));
    assert (inSection.ReadVar("NumPoints", tempNumPoints));
    Init (tempr0, tempSpacing, tempNumPoints);
  }

  LogGrid (scalar R0, scalar spacing, int numpoints)
  {
    Init  (R0, spacing, numpoints);
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




/// ClusterGrid is a function whose gridpoints are clustered tightly
/// around the origin.
class ClusterGrid : public Grid
{
private:
  double x0, dri, rr;

 public:
  double Start, End, Cluster;

  GridType Type()
    { return (CLUSTER); }

  void ReadParams (FILE *fin)
    {    }

  void WriteParams (FILE *fout)
    {    }

  void WriteInput (FILE *fout)
  {  }

  int ReverseMap(scalar r)
    {
      return ((int)floor (dri/(r-rr) -1.0 + x0));
    }

  void Init (double start, double end, double cluster, int numpoints)
  {
    Start = start; End = end; Cluster = cluster;
    NumPoints = numpoints;
    
    x0 = (NumPoints - Cluster)/(1.0-Cluster);
    dri = -(End-Start)*((double) NumPoints-x0)*(1.0-x0) / 
      ((double)NumPoints-1.0);
    rr = Start - dri/(1.0-x0);
    grid.resize(NumPoints);
    for (int i=0; i<NumPoints; i++)
      grid(i) = rr+dri/(i+1.0-x0);
  }



  void Write (IOSectionClass &outSection)
  {
    outSection.WriteVar ("Points", grid); 
    outSection.WriteVar ("Type", "Cluster");
    outSection.WriteVar ("Start", Start);
    outSection.WriteVar ("End", End);
    outSection.WriteVar ("Cluster", Cluster);
    outSection.WriteVar ("NumPoints", NumPoints);
  }

  void Read (IOSectionClass &inSection)
  {
    double start, end, cluster;
    int numpoints;
    assert (inSection.ReadVar("Start", start));
    assert (inSection.ReadVar("End", end));
    assert (inSection.ReadVar("Cluster", cluster));
    assert (inSection.ReadVar("NumPoints", numpoints));
    Init (start, end, cluster, numpoints);
  }

  ClusterGrid (double start, double end, double cluster, int numpoints)
  {
    Init  (start, end, cluster, numpoints);
  }
  ClusterGrid () 
  { /* Do nothing */ }
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

inline Grid* ReadGrid (IOSectionClass &inSection)
{
  string Type;
  assert (inSection.ReadVar ("Type", Type));

  Grid *newGrid;
  if (Type == "Linear")
    newGrid = new LinearGrid;
  else if (Type == "Optimal")
    newGrid = new OptimalGrid;
  else if (Type == "Optimal2")
    newGrid = new OptimalGrid2;
  else if (Type == "Log")
    newGrid = new LogGrid;
  else if (Type == "Cluster")
    newGrid = new ClusterGrid;
  else
    {
      cerr << "Unrecognized Grid type " << Type << "\n";
      exit(1);
    }
  newGrid->Read(inSection);
  return (newGrid);
}




#endif
