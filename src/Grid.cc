#include "Grid.h"


Grid *ReadGrid (InputBuffer &SectionBuf)
{
  Grid *grid;
  char TypeString[50];
  int Success;
  
  if (!SectionBuf.ReadVar ("Type", TypeString, 50))
    Abort("Could not find \"Type\" variable in section \"Grid\"");
  
  SectionBuf.Rewind();
  if (!strcmp(TypeString, "OPTIMAL"))
    {
      scalar a, b, Z, rmax;
      int NumPoints;
      if (SectionBuf.ReadVar ("a", a))
	{
	  SectionBuf.Rewind();
	  if (SectionBuf.ReadVar("b", b))
	    {
	      SectionBuf.Rewind();
	      if (SectionBuf.ReadVar ("NumPoints", NumPoints))
		  grid = new OptimalGrid(a, b, NumPoints);
	      else
		return NULL;
	    }
	  else
	    return NULL;
	}
      else
	{
	  SectionBuf.Rewind();
	  if (SectionBuf.ReadVar ("Z", Z))
	    {
	      SectionBuf.Rewind();
	      if (SectionBuf.ReadVar ("rmax", rmax))
		grid = new OptimalGrid (Z, rmax);
	      else
		return NULL;
	    }
	  else
	    {
	      SectionBuf.Rewind();
	      if (SectionBuf.ReadVar ("NumPoints", NumPoints))
		{
		  SectionBuf.Rewind();
		  if (SectionBuf.ReadVar ("rmax", rmax))
		    {
		      SectionBuf.Rewind();
		      if (SectionBuf.ReadVar("b", b))
			grid = new OptimalGrid (NumPoints, rmax, b);
		      else
			return NULL;
		    }
		  else
		    return NULL;
		}
	      else return NULL;
	    }
	}
    }
  else if (!strcmp(TypeString, "LINEAR"))
    {
      scalar Start, End;
      int NumPoints;
      if (SectionBuf.ReadVar ("Start", Start))
	{
	  SectionBuf.Rewind();
	  if (SectionBuf.ReadVar ("End", End))
	    {
	      if (SectionBuf.ReadVar ("NumPoints", NumPoints))
		grid = new LinearGrid (Start, End, NumPoints);
	      else
		return NULL;
	    }
	  else
	    return NULL;
	}
      else
	return NULL;
    }
  else if (!strcmp(TypeString, "LOG"))
    {
      scalar r0, Spacing;
      int NumPoints;
      if (SectionBuf.ReadVar ("r0", r0))
	{
	  SectionBuf.Rewind();
	  if (SectionBuf.ReadVar("Spacing", Spacing))
	    {
	      SectionBuf.Rewind();
	      if (SectionBuf.ReadVar("NumPoints", NumPoints))
		grid = new LogGrid (r0, Spacing, NumPoints);
	      else return NULL;
	    }
	  else return NULL;
	}
      else return NULL;
     
    }
  else
    Abort ("Unknown Type in section \"Grid\".");
  return (grid);
}
