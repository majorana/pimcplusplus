#ifndef DavidPAClass_H
#define DavidPAClass_H

#include "PAFitBase.h"
#include "../Blitz.h"
#include <fstream>
#include <iostream>
//#include "../../PathDataClass.h"


///This is the pair action class. It uses the following formula in
///order to calculate the pair action
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
class DavidPAClass : public PairActionFitClass 
{

 private:
  
  ///Holds the Ukj coefficients for a given q
  Array<double,1> TempukjArray;
  Array<double,1> TempdukjArray;

  //  DistanceTableClass *DistanceTable;
  /// Skips to the next string in the file whose substring matches skipToString
  inline string SkipTo(ifstream &infile, string skipToString);
  /// Reads a Fortran 3 tensor
  inline void ReadFORTRAN3Tensor(ifstream &infile, Array<double,3> &tempUkj);
 public:
  Array<double,1> Potential; 
  string type1,type2;
  inline bool Read(IOSectionClass &IOSection,double desiredTau, int numLevels);
  inline void Print();
  double DesiredTau;
  int TauPos;
  int NumLevels;
  int NumTau;
  /// This stores the coefficients of the expansion specified above.
  /// The array index is over the levels.  You call it with the q
  /// value and a temporary array to get all of the values in that
  /// column. 
  Array<MultiCubicSpline,1> ukj; ///<(level )
  ///Same as ukj but stores the beta derivatives.
  Array<MultiCubicSpline,1> dukj; ///<(level )
  /////  MultiCubicSpline Pot;
  /// Calculate the U(s,q,z) value when given s,q,z and the level 
  void calcUsqz(double s,double q,double z,int level,
		       double &U, double &dU, double &V);
  /// This is the order of the fit to use. 
  int n;
  /// This is the temperature 
  double tau;
  /// Function to read David's squarer file input.
  void ReadDavidSquarerFile(string DMFile);
  double U (double q, double z, double s2, int level);
  double dU(double q, double z, double s2, int level);
  double V(double r);
  bool IsLongRange(); 
  void DoBreakup (const dVec &box, const Array<dVec,1> &kVecs);
#ifdef MAKE_FIT
  void ReadParams  (IOSectionClass &inSection);
  void WriteBetaIndependentInfo (IOSectionClass &outSection);
  /// Returns weighter RMS error
  void Error (Rho &rho, double &Uerror, double &dUerror);
  void AddFit (Rho &rho);
  void WriteFits(IOSectionClass &outSection);
#endif
};



inline bool DavidPAClass::Read(IOSectionClass &in,double x, int y)
{
  string fileName;
  DesiredTau=x;
  NumLevels=y;
  // Read Particles;
  assert(in.OpenSection("Fits"));
  assert(in.OpenSection("Particle1"));
  Particle1.Read(in);
  in.CloseSection();
  assert(in.OpenSection("Particle2"));
  Particle2.Read(in);
  in.CloseSection();
  lambda = Particle1.lambda + Particle2.lambda;
  

  //  assert(in.ReadVar("tau",DesiredTau));
  //  assert(in.ReadVar("MaxLevels",NumLevels));
  assert(in.ReadVar("Daviddmfile",fileName));
  ReadDavidSquarerFile(fileName.c_str());
  //  assert(in.ReadVar("type1",type1));
  //  assert(in.ReadVar("type2",type2));
  in.CloseSection();
  return true;
}
     

inline string DavidPAClass::SkipTo(ifstream &infile,string skipToString)
{
  string lineString;
  do{
    getline(infile,lineString);

  } while (lineString.find(skipToString,0)==string::npos && !infile.eof());

  if (infile.eof()){
    cerr<<"ERROR!!!  No "<<skipToString<<" found in Davids squarer file\n";
  }
  return lineString;

}


inline void DavidPAClass::ReadFORTRAN3Tensor(ifstream &infile,Array<double,3> &tempUkj)
{

   
  for (int counterTau=0;counterTau<tempUkj.extent(2);counterTau++){
    for (int counterUkj=0;counterUkj<tempUkj.extent(1);counterUkj++){
      for (int counterGridPoints=0;counterGridPoints<tempUkj.extent(0);counterGridPoints++){

	infile>>tempUkj(counterGridPoints,counterUkj,counterTau);
	//	tempUkj(counterGridPoints,counterUkj,counterTau)=0;
      }
    }
  }

}


inline bool IsDigit(char c)
{
  return (((c>='0') && (c<='9')) || (c == '-'));
}

inline bool IsNumberChar(char c)
{
  bool IsNumber = false;
  
  IsNumber = IsNumber || ((c>='0')&&(c<='9'));
  IsNumber = IsNumber || (c=='-');
  IsNumber = IsNumber || (c=='.');
  return IsNumber;
}



inline int GetNextInt(string &s)
{
  int i=0;
  int size = s.size();
  while ((i<size) && (!IsDigit(s[i])))
    i++;
  s = s.substr (i, size-i);

  istringstream sstream(s);

  int num;
  sstream >> num;
  int pos = sstream.tellg();
  if (pos == -1)
    s = "";
  else
    s = s.substr(pos, s.size()-pos);

  return (num);
}


inline string GetNextWord(string &s)
{
  int counter=0;
  while (s[counter] == ' '){
    counter++;
  }
  int startSpot=counter;
  while (s[counter] != ' '){
    counter++;
  }
  string toReturn=s.substr(startSpot,counter-startSpot);
  s=s.substr(counter,s.size()-counter);
  return toReturn;
}

inline double GetNextDouble(string &s)
{
  int i=0;
  int size = s.size();
  while ((i<size) && (!IsNumberChar(s[i])))
    i++;
  s = s.substr (i, size-i);

  istringstream sstream(s);

  double num;
  sstream >> num;
  int pos = sstream.tellg();
  if (pos == -1)
    s = "";
  else
    s = s.substr(pos, s.size()-pos);

  return (num);
}

inline void DavidPAClass::Print()
{ cerr<<"hi\n";}
            



#endif
