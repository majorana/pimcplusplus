#ifndef GRID_H
#define GRID_H

#include "Blitz.h"
#include "InputFile.h"

typedef enum {NONE, LINEAR, OPTIMAL, LOG} GridType;

class Grid
{
 protected:
  Array<scalar,1> grid;
 public:
  scalar Start, End;
  int NumPoints;

  inline scalar operator()(int i) const
    {
      return (grid(i));
    }

  virtual GridType Type()
    { return (NONE); }

  virtual void ReadParams (FILE *fin)
  {
    cerr << "Should never get here.\n";
    exit(1);
  }
  
  virtual void WriteParams (FILE *fout)
  {
    cerr << "Should never get here.\n";
    exit(1);
  }

  virtual void WriteInput (FILE *fout)
  {
    cerr << "Should never get here.\n";
    exit(1);
  }

  virtual int ReverseMap (scalar r)
  {
    cerr << "Should never get here.\n";
    exit(1);
    return (0);
  }
  inline scalar Interp(const Array<scalar,1> &Vals, scalar r, int order)
    {
      if (/*(r<Start) ||*/ (r>(End*1.0000000001)))
	{
	  cerr << "r = " << r << "outside grid range in Interp.\n";
	  exit(1);
	}
      int NumVals = order+1;
      int StartIndex, EndIndex;
      if (NumPoints&1)
	StartIndex = ReverseMap(r) - NumVals/2;
      else
	StartIndex = ReverseMap(r) - NumVals/2 + 1;
      EndIndex = StartIndex + order;
      if (StartIndex < 0)
	{
	  EndIndex -= StartIndex;
	  StartIndex = 0;
	}
      if (EndIndex >= NumPoints)
	{
	  StartIndex -= (EndIndex-(NumPoints-1));
	  EndIndex = NumPoints-1;
	}

      scalar delta =  grid(StartIndex) - grid(EndIndex);
      for (int i=StartIndex; i<=EndIndex; i++)
	if (fabs((grid(i)-r)/delta) < 1.0e-7)
	  return (Vals(i));

      // Now do a standard, nth order interpolation
      scalar Sum=0.0;

      for (int i=StartIndex; i<=EndIndex; i++)
	{
	  scalar Numer, Denom;
	  Numer = Denom = 1.0;
	  for (int j=StartIndex; j<=EndIndex; j++)
	    if (i!=j) 
	      {
		Numer *= (r-grid(j));
		Denom *= (grid(i) - grid(j));
	      }
	  Sum += Numer/Denom * Vals(i);
	}
      return(Sum);
    }  

  inline scalar Deriv (const Array<scalar,1> &Vals, scalar r, int order)
    {
      if (/*(r<Start) ||*/ (r>(End*1.000001)))
	{
	  cerr << "r = " << r << "outside grid range in Interp.\n";
	  // exit(1);
	}

      
      int index = ReverseMap(r);
      //cerr << "r = " << r << " grid(index) = " << grid(index) << "  index = " << index << "\n";
      if (index < 0)
	index = 0;
      if (index > (NumPoints-2))
	index = NumPoints-2;
      
      return ((Vals(index+1)-Vals(index)) / (grid(index+1)-grid(index)));
    

      //  int NumVals = order+1;
//        int StartIndex, EndIndex;
//        if (NumPoints&1)
//  	StartIndex = ReverseMap(r) - NumVals/2;
//        else
//  	StartIndex = ReverseMap(r) - NumVals/2 + 1;
//        EndIndex = StartIndex + order;
//        if (StartIndex < 0)
//  	{
//  	  EndIndex -= StartIndex;
//  	  StartIndex = 0;
//  	}
//        if (EndIndex >= NumPoints)
//  	{
//  	  StartIndex -= (EndIndex-(NumPoints-1));
//  	  EndIndex = NumPoints-1;
//  	}

//        scalar delta =  grid(StartIndex) - grid(EndIndex);

      /*for (int m=StartIndex; m<=EndIndex; m++)
	      if (fabs((grid(m)-r)/delta) < 1.0e-3)
	  {
	    //cerr << m << "\n";
	    scalar Sum=0.0;
	    for (int i=StartIndex; i<=EndIndex; i++)
	      {
		if (i == m)
		  {
		    scalar Numer, Denom;
		    Numer = 0.0;
		    Denom = 1.0;
		    for (int j=StartIndex; j<=EndIndex; j++)
		      {
			scalar NumerFact = 1.0;
			for (int k=StartIndex; k<=EndIndex; k++)
			  if (k!=i)
			    if (k!=j)
			      NumerFact *= (r-grid(j));
			if (i!=j) 
			  {
			    Numer += NumerFact;
			    Denom *= (grid(i) - grid(j));
			  }
		      }
		    Sum += Numer/Denom * Vals(i);
		  }
		else
		  {
		    scalar Numer = 1.0;
		    scalar Denom = 1.0;
		    for (int j=StartIndex; j<=EndIndex; j++)
		      {
			if (i!=j)
			  {
			    if (j != m)
			      Numer *= -(r - grid(j));
			    Denom *= (grid(j) - grid(i));
			  }
		      }
		    Sum += (Numer/Denom) * Vals(i);
		  } 
	      }
	    return (Sum);
	    }*/
		  
            

      // Now do a standard, nth order interpolation

      //  scalar Sum=0.0;
//        for (int i=StartIndex; i<=EndIndex; i++)
//  	{
//  	  scalar Numer, Denom;
//  	  Numer = 0.0;
//  	  Denom = 1.0;
//  	  for (int j=StartIndex; j<=EndIndex; j++)
//  	    {
//  	      scalar NumerFact = 1.0;
//  	      for (int k=StartIndex; k<=EndIndex; k++)
//  		if (k!=i)
//  		  if (k!=j)
//  		    NumerFact *= (r-grid(j));
//  	      if (i!=j) 
//  		{
//  		  Numer += NumerFact;
//  		  Denom *= -(grid(i) - grid(j));
//  		}
//  	    }
//  	  Sum += Numer/Denom * Vals(i);
//  	}
//        return(Sum);
    }
};

class LinearGrid : public Grid
{
 private:
  scalar delta;
 public:

  GridType Type()
    { return (LINEAR); }

  void ReadParams (FILE *fin)
    {
      fscanf (fin, " %lf  %lf %d ", &Start, &End, &NumPoints);
      grid.resize(NumPoints);
      delta = (End-Start)/(scalar)(NumPoints-1);
      for (int i=0; i<NumPoints; i++)
	grid(i) = Start + (scalar)i*delta;
    }

  void WriteParams (FILE *fout)
    {
      fprintf (fout, "%1.16e %1.16e %d\n", Start, End, NumPoints);
    }

  void WriteInput (FILE *fout)
  {
    fprintf (fout, "  {\n");
    fprintf (fout, "    Type = \"LINEAR\";\n");
    fprintf (fout, "    Start = %1.16e;\n", Start);
    fprintf (fout, "    End = %1.16e;\n", End);
    fprintf (fout, "    NumPoints = %d;\n", NumPoints);
    fprintf (fout, "  }\n");
  }

  int ReverseMap(scalar r)
    {
      return ((int)floor((r-Start)/delta));
    }

  inline void Init(scalar start, scalar end, int numpoints)
  {
      Start=start; End=end; NumPoints=numpoints;
      grid.resize(NumPoints);
      delta = (End-Start)/(scalar)(NumPoints-1);
      for (int i=0; i<NumPoints; i++)
	grid(i) = Start + (scalar)i*delta;
  }

  LinearGrid ()
  {
    // Do nothing
  }
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
