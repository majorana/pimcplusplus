#include "InputOutput.h"
#include <stdio.h>
#include <stdlib.h>

//////////////////////////////////////////////////////////////
// Read a file into the buffer, ignoring C++ style comments.//
//////////////////////////////////////////////////////////////
int
InputSectionClass::Open(string FileName)
{
  int InQuotes=0;

  FILE *fin;
  if ((fin = fopen (FileName.c_str(), "r")) == NULL)
    return (0);

  int count = 0;

  char ch;
  do
    {
      ch = fgetc(fin);
      if (ch != EOF)
	{
	  if (InQuotes)
	    {
	      if (ch == '\"')
		InQuotes = 0;
	      
	      count++;
	    }
	  else if (ch == '\"')
	    {
	      InQuotes = 1;
	      count++;
	    }
	  else if (ch == '/')
	    {
	      char ch2 = fgetc(fin);
	      if (ch2 == '/')
		{
		  char ch3;
		  do
		    {
		      ch3 = fgetc(fin);
		    } while ((ch3 != EOF) && (ch3 != '\n'));
		}

	      else if (ch2 == '*')
		{
		  int found = 0;
		  char ch3;
		  ch3 = fgetc(fin);
		  while (!found)
		    {
		      if (ch3 == '*')
			{
			  ch3 = fgetc(fin);
			  if (ch3 == '/')
			    found = 1;
			}
		      else
			ch3 = fgetc(fin);
		    }
		}
	      else
		{
		  count ++;
		  if (ch2 != EOF)
		    count++;

		}
	    }
	  else if (ch == '#')
	    {
	      char str[500];
	      fgets(str, 500, fin);
	      int len = strlen (str);
	      if (!strncmp(str, "include", 7))
	      {
		int i = 7;
		while ((i<len) && ((str[i]==' ') || (str[i]=='\t')))
		  i++;
		if (str[i] != '\"')
		  {
		    cerr << "Badly formed #include directive " 
			 << "in InputBuffer::Read().\n";
		    exit(1);
		  }
		else
		  {
		    i++;
		    char IncludeName[500];
		    int j=0;
		    while ((str[i] != '\"') && (i < len))
		      {
			IncludeName[j] = str[i];
			i++;
			j++;
		      }
		    IncludeName[j] = '\0';
		    if (str[i] != '\"')
		      {
			cerr << "Badly formed #include directive " 
			     << "in InputBuffer::Read().\n";
			exit(1);
		      }
		    InputBuffer IncludeBuf;
		    if (!IncludeBuf.Read(IncludeName))
		      {
			cerr << "Error reading included file "
			     << IncludeName << " in InputBuffer::Read().\n";
			exit(1);
		      }
		    //IncludeBuf.Write(stderr);
		    count += IncludeBuf.Size();
		  }
	      }
	    }
	  else
	    count++;
	}
    } while (ch != EOF);
			 
  // Now resize the buffer
  Resize (count);
  
  // Now actually read it in;
  fseek(fin, (long)0, SEEK_SET);
  count = 0;
  
  do
    {
      ch = fgetc(fin);
      if (ch != EOF)
	{
	  if (InQuotes)
	    {
	      if (ch == '\"')
		InQuotes = 0;
	      
	      buffer(count) = ch;
	      count++;
	    }
	  else if (ch == '\"')
	    {
	      InQuotes = 1;
	      buffer(count) = ch;
	      count++;
	    }
	  else if (ch == '/')
	    {
	      char ch2 = fgetc(fin);
	      if (ch2 == '/')
		{
		  char ch3;
		  do
		    {
		      ch3 = fgetc(fin);
		    } while ((ch3 != EOF) && (ch3 != '\n'));
		}
	      else if (ch2 == '*')
		{
		  int found = 0;
		  char ch3;
		  ch3 = fgetc(fin);
		  while (!found)
		    {
		      if (ch3 == '*')
			{
			  ch3 = fgetc(fin);
			  if (ch3 == '/')
			    found = 1;
			}
		      else
			ch3 = fgetc(fin);
		    }
		}
	      else
		{
		  buffer(count) = ch;
		  count ++;
		  if (ch2 != EOF)
		    {
		      buffer(count) = ch2;
		      count++;
		    }
		}
	    }
	  else if (ch == '#')
	    {
	      char str[500];
	      fgets(str, 500, fin);
	      int len = strlen (str);
	      if (!strncmp(str, "include", 7))
	      {
		int i = 7;
		while ((i<len) && ((str[i]==' ') || (str[i]=='\t')))
		  i++;
		if (str[i] != '\"')
		  {
		    cerr << "Badly formed #include directive " 
			 << "in InputBuffer::Read().\n";
		    exit(1);
		  }
		else
		  {
		    i++;
		    char IncludeName[500];
		    int j=0;
		    while ((str[i] != '\"') && (i < len))
		      {
			IncludeName[j] = str[i];
			i++;
			j++;
		      }
		    IncludeName[j] = '\0';
		    if (str[i] != '\"')
		      {
			cerr << "Badly formed #include directive " 
			     << "in InputBuffer::Read().\n";
			exit(1);
		      }
		    InputBuffer IncludeBuf;
		    IncludeBuf.Read(IncludeName);
		    for (int k=0; k<IncludeBuf.Size(); k++)
		      {
			buffer(count) = IncludeBuf(k);
			count++;
		      }
		  }
	      }

	    }
	  else
	    {
	      buffer(count) = ch;
	      count++;
	    }
	}
    } while (ch != EOF);
  fclose (fin);
  return (1);
}		      
