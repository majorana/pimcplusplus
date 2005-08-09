#ifndef PATH_CLASS_H
#define PATH_CLASS_H

#include <Common/IO/InputOutput.h>
#include "MirroredClass.h"
#include "SpeciesClass.h"
#include <Common/Random/Random.h>
#include <Common/MPI/Communication.h>

class ActionsClass;

///The number of time slices is the number of slices on this processor.
///In all cases this processor shares a time slice with the processor 
///ahead of it and behind it. The convention for the shared slices
///is that the processor owns its first but not its last slice.
class PathClass
{
private:  
  /// Path stores the position of all the particles at all time
  /// slices.  The order for access is timeslice, particle
  Mirrored2DClass<dVec> Path;
  /// This stores the path at the reference slice.
  /// Stores what species a particle belongs to
  Array<int,1> SpeciesNumber;
  Array<SpeciesClass *,1> SpeciesArray;
  int MyNumSlices;

  ////////////////////////
  /// Class References ///
  ////////////////////////
  ActionsClass &Actions;

  /////////////////////
  /// Misc. Helpers ///
  /////////////////////
  void LeviFlight (Array<dVec,1> &vec, double lambda);
  void ReadOld(string fileName, bool replicate);
  ////////////////////////////////
  /// Boundary conditions stuff //
  ////////////////////////////////
  dVec IsPeriodic;

  dVec Box, BoxInv;
  dVec kBox; //kBox(i)=2*Pi/Box(i)

  ///////////////////////////////////////////////
  /// k-space stuff for long-range potentials ///
  ///////////////////////////////////////////////
 public:
  /// This is the maximum number of k vectors in each direction
  TinyVector<int,NDIM> MaxkIndex;
  /// Stores the radius of the k-space sphere we sum over
  double kCutoff;
  /// Allocates and sets up the k-vectors.
  void SetupkVecs3D();

  void SetupkVecs2D();
 private:
  ///These two preceding function/array only need to be used in order to
  ///make the readings of David's file work correctly (sorry,
  ///Ken...I'll find out how to make it more transparent in a bit)
  void SortRhoK();

  /// Stores indices into C array for fast computation of rho_k
  Array<TinyVector<int,NDIM>,1> kIndices;
  /// This stores e^{i\vb_i \cdot r_i^\alpha}
  TinyVector<Array<complex<double>,1>,NDIM> C;

  /// Thus function finds the two closest particles in species
  /// speciesNum and swaps their positions.  This is used to flip the
  /// sign of the nodal determinant.
  void SwapClosest(int speciesNum);
public:
  Array<int,1> MagKint;
  Array<double,1> MagK;

  CommunicatorClass &Communicator;

  Mirrored1DClass<dVec> RefPath;
  void BroadcastRefPath();
  /// True if we need k-space sums for long range potentials.
  bool LongRange;
  ///True if we are doing long range in David's way
  bool DavidLongRange;
  /// Stores the actual k vectors
  Array<dVec,1> kVecs;
  /// This holds the density of particles in k-space.  Indexed by
  /// (slice, species, k-vector).  Defined as
  /// \rho^\alpha_k = \sum_i e^{i\mathbf{k}\cdot\mathbf{r}_i^\alpha}
  /// Stores the kvectors needed for the reciporical space sum.
  /// Stores only half the vectors because of k/-k symmetry.
  Mirrored3DClass< complex<double> > Rho_k;
  void CalcRho_ks_Slow(int slice, int species);  
  void CalcRho_ks_Fast(int slice, int species);  

private:
  //  int RefSliceCheck;
  void ShiftPathData(int sliceToShift);
  void ShiftRho_kData(int sliceToShift);
public:
  /// Stores the position of the reference slice w.r.t. time slice 0
  /// on processor 0.
  int RefSlice;

  Mirrored1DClass<int> Permutation;
  /// This function accumulates the total permutation vector
  /// from all of the processors individual permutation vector.
  /// Only processor 0 gets the result
  void TotalPermutation (Array<int,1> &permVec);
  RandomClass &Random;
  int TotalNumSlices;
  double tau; //we need to set this still
  /// A scratch array to hold a boolean indicating whether we've
  /// looped over this particle yet
  Array<bool,1> DoPtcl;

  //jgadd
  int numMol;  // we want to store the number of molecules; it's loaded in the .cc file
  Array<int,1> MolRef;

  inline void  SetBox (dVec box);
  inline const dVec& GetBox();
  inline const dVec& GetkBox();
  inline double GetVol();
  inline double Getkc();
  inline void  SetPeriodic(TinyVector<bool,NDIM> period);
  inline dVec  GetPeriodic() const { return IsPeriodic; }

  //////////////////////////////////
  /// TimeSlice parallelism stuff //
  //////////////////////////////////

  /// Returns the range of time slices that a processor holds.
  inline void SliceRange (int proc, int &slice1, int &slice2);
  /// Returns which processor owns the given slice
  inline int SliceOwner (int slice);


  /////////////////////////////////
  /// Displacements / Distances ///
  /////////////////////////////////
  inline void DistDisp (int slice, int ptcl1, int ptcl2,
			double &dist, dVec &disp);
  inline void DistDisp (int sliceA, int sliceB, int ptcl1, int ptcl2,
			double &distA, double &distB,dVec &dispA, dVec &dispB);
  void RefDistDisp (int slice, int refPtcl, int ptcl,
		    double &dist, dVec &disp);
  //  inline double Distance (int slice, int ptcl1, int ptcl2);Not used?
  inline dVec Velocity (int sliceA, int sliceB, int ptcl);
  inline void PutInBox (dVec &v);


  //////////////////////////
  /// Data manipulations ///
  //////////////////////////
  inline const dVec& operator() (int slice, int ptcl) const;
  inline dVec& operator() (int slice, int ptcl);
  inline void SetPos (int slice, int ptcl, const dVec& r);
  inline int NumParticles();
  inline int NumTimeSlices();
  inline int GetRefSlice() const;

  void MoveJoin(int oldJoin, int newJoin);      
  void ShiftData(int sliceToShift);
  void AcceptCopy(int startTimeSlice,int endTimeSlice, 
		  const Array <int,1> &activeParticle);
  void RejectCopy(int startTimeSlice,int endTimeSlice, 
		  const Array <int,1> &activeParticle );


  /////////////////////////////
  /// Species Manipulations ///
  /////////////////////////////
  inline int ParticleSpeciesNum(int ptcl);
  inline SpeciesClass& ParticleSpecies(int ptcl);
  inline SpeciesClass& Species(int speciesNum);
  inline int SpeciesNum (string name);
  inline void AddSpecies (SpeciesClass *newSpecies);
  inline int NumSpecies();


  //////////////////////////
  /// IO and allocations ///
  //////////////////////////
  void Read(IOSectionClass &inSection);
  void Allocate();
  void InitPaths(IOSectionClass &inSection);

  /// This class will create a new brownian random walk for
  /// species(speciesNum).  If the species is fermion, it will do its
  /// best to construct a node-avoiding walk that is reasonable.
  void NodeAvoidingLeviFlight (int speciesNum, Array<dVec,1> &initialPoints);

  void SetupClones();
  inline PathClass(CommunicatorClass &communicator,
		   RandomClass &random,
		   ActionsClass &actions);
  friend void SetupPathNaCl(PathClass &path);
  friend void SetupPathZincBlend(PathClass &path);
  friend void SetupPathSimpleCubic(PathClass &path);

  //////////////////////////
  /// Fermions           ///
  //////////////////////////
  inline bool HasFermions(const Array<int,1>& activeParticles);


  //////////////////////////
  /// Open Loops         ///
  //////////////////////////
  MirroredClass<int> OpenPtcl;
  MirroredClass<int> OpenLink;
  bool OpenPaths;
  int OpenSpeciesNum;
  void InitOpenPaths();
  void DistanceToTail();
  MirroredClass<int> Weight;
};

inline bool PathClass::HasFermions(const Array<int,1>& activeParticles)
{
  bool HasFermion=false;
  for (int ptclIndex=0;ptclIndex<activeParticles.size();ptclIndex++){
    int ptcl=activeParticles(ptclIndex);
    HasFermion = HasFermion || ParticleSpecies(ptcl).GetParticleType();
  }
  return HasFermion;
}





inline int PathClass::SpeciesNum (string name)
{
  int i=0;
  while ((i<SpeciesArray.size()) && (SpeciesArray(i)->Name != name))
    i++;
  if (i == SpeciesArray.size())
    return -1;
  else
    return i;
}


/// Return what species type a particle belongs to;
inline int 
PathClass::ParticleSpeciesNum(int ptcl)
{
  return (SpeciesNumber(ptcl));
}
/// Return a reference to the species that particle ptcl belongs to
inline SpeciesClass& 
PathClass::ParticleSpecies(int ptcl) 
{
  return *(SpeciesArray(SpeciesNumber(ptcl)));
}
/// Return a species references
inline SpeciesClass& 
PathClass::Species(int speciesNum)
{
  return (*(SpeciesArray(speciesNum)));
}


/// Returns the number of particle Species
inline int 
PathClass::NumSpecies() 
{ 
  return SpeciesArray.size();
}

inline int 
PathClass::NumParticles() 
{ 
  return Path.cols()-OpenPaths;
}

///The number of time slices is the number of slices on this processor.
///In all cases this processor shares a time slice with the processor 
///ahead of it and behind it. The convention for the shared slices
///is that the processor owns its first but not its last slice.
inline int 
PathClass::NumTimeSlices() 
{ 
  return Path.rows();
}

/// Returns the position of the reference slice w.r.t. slice 0 on this
/// processor 0.
inline int 
PathClass::GetRefSlice() const
{
  return RefSlice;
}


/// Returns the position of particle ptcl at time slice timeSlice
inline const dVec& 
PathClass::operator() (int slice, int ptcl) const
{ 
  return Path(slice, ptcl); 
}

/// Returns the position of particle ptcl at time slice timeSlice
inline dVec& 
PathClass::operator() (int slice, int ptcl)
{ 
  return Path(slice, ptcl); 
}

inline void 
PathClass::SetPos(int slice, int ptcl, const dVec &r)
{
  (*this)(slice, ptcl) = r;
}

inline void 
PathClass::AddSpecies (SpeciesClass *newSpecies)
{
  int numSpecies = SpeciesArray.size();
  /// Add an element for the new species
  SpeciesArray.resizeAndPreserve(numSpecies+1);
  SpeciesArray(numSpecies) = newSpecies;
}

inline 
PathClass::PathClass (CommunicatorClass &communicator,
			     RandomClass &random,
			     ActionsClass &actions) : 
  Communicator(communicator), Random(random), Actions(actions)
{
  //      NumSpecies = 0;
  TotalNumSlices=0;

  //  OpenPaths=true;
  OpenPaths=false; //turns off open loops (Should be read at some poitn)
  Weight=1;
}


// <<<<<<< .mine
// =======
// ///Must start diagonal or this is goign to break because of where
// ///you are starting the openlink at 0
// inline void 
// PathClass::InitOpenPaths()
// {
//   cerr<<"Starting to initialize"<<endl;
//   if (OpenPaths){
//     cerr<<"openpaths"<<endl;
//     SetMode(OLDMODE);
//     OpenLink=NumTimeSlices()-1;
//     cerr<<"Set up the open link"<<endl;
//     OpenPtcl=Species(OpenSpeciesNum).FirstPtcl;
//     cerr<<"set up the open particle"<<endl;
//     SetMode(NEWMODE);
//     OpenLink=NumTimeSlices()-1;
//     OpenPtcl=0;
//     cerr<<"Preparing for moving things in path around"<<endl;
//     Path[OLDMODE]((int)OpenLink,NumParticles())=Path[OLDMODE]((int)OpenLink,(int)OpenPtcl);
//     cerr<<"Moved the first thing"<<endl;
//     Path[NEWMODE]((int)OpenLink,NumParticles())=Path[NEWMODE]((int)OpenLink,(int)OpenPtcl);
//     cerr<<"Moved the second thing"<<endl;
//     Path[OLDMODE](0,NumParticles())=Path[OLDMODE](0,(int)OpenPtcl);
//     cerr<<"Moved the first thing"<<endl;
//     Path[NEWMODE](0,NumParticles())=Path[NEWMODE](0,(int)OpenPtcl);
// >>>>>>> .r523

////////////////////////////////
/// Boundary conditions stuff //
////////////////////////////////

inline void 
PathClass::SetBox (dVec box)
{
  Box = box;
  for (int i=0; i<NDIM; i++){
    BoxInv(i) = 1.0/box(i);
    kBox(i)=2*M_PI*BoxInv(i);
  }
  
}

inline const dVec& 
PathClass::GetBox()
{
  return Box;
}

inline const dVec& 
PathClass::GetkBox()
{
  return kBox;
}

inline double 
PathClass::GetVol()
{
  double  vol=1.0;
  for (int i=0;i<NDIM;i++){
    vol*=Box(i);
  }
  return vol;
}

inline double 
PathClass::Getkc()
{
  if (LongRange)
    return kCutoff;
  else
    return 0.0;
}


inline void 
PathClass::SetPeriodic(TinyVector<bool,NDIM> period)
{
  for (int i=0; i<NDIM; i++)
    IsPeriodic(i) = period(i) ? 1.0 : 0.0;
}

inline void 
PathClass::DistDisp (int slice, int ptcl1, int ptcl2,
			       double &dist, dVec &disp)
{
  disp = Path(slice, ptcl2) -Path(slice, ptcl1);
  
  for (int i=0; i<NDIM; i++) {
    double n = -floor(disp(i)*BoxInv(i)+0.5);
    disp(i) += n*IsPeriodic(i)*Box(i);
  }
  dist = sqrt(dot(disp,disp));

#ifdef DEBUG
  dVec DBdisp = Path(slice, ptcl2) -Path(slice, ptcl1);
  for (int i=0; i<NDIM; i++) {
    while (DBdisp(i) > 0.5*Box(i))
      DBdisp(i) -= Box(i);
    while (DBdisp(i) < -0.5*Box(i)) 
      DBdisp(i) += Box(i);
    if (fabs(DBdisp(i)-disp(i)) > 1.0e-12){ 
      cerr<<DBdisp(i)<<" "<<disp(i)<<endl;
    }
    //    assert (fabs(DBdisp(i)-disp(i)) < 1.0e-12);
  }
#endif
}


inline void 
PathClass::DistDisp (int sliceA, int sliceB, int ptcl1, int ptcl2,
		     double &distA, double &distB, 
		     dVec &dispA, dVec &dispB)
{  
  //  //  bool changePtcl1=(OpenPaths && sliceB==(int)OpenLink && ptcl1==(int)OpenPtcl);
  //  //  ptcl1=ptcl1*!changePtcl1+NumParticles()*changePtcl1;
  //  //  bool changePtcl2=(OpenPaths && sliceB==(int)OpenLink && ptcl2==(int)OpenPtcl);
  //  //  ptcl2=ptcl2*!changePtcl2+NumParticles()*changePtcl2;
  dispA = Path(sliceA, ptcl2) - Path(sliceA,ptcl1);
  if (OpenPaths && sliceB==(int)OpenLink && ptcl1==(int)OpenPtcl){
    dispB=Path(sliceB,ptcl2)-Path(sliceB,NumParticles());
  }
  else if (OpenPaths && sliceB==(int)OpenLink && ptcl2==(int)OpenPtcl){
    dispB=Path(sliceB,NumParticles())-Path(sliceB,ptcl1);
  }
  else{
    dispB = Path(sliceB, ptcl2) - Path(sliceB,ptcl1);
  }
  //  dVec tempDispB;
  //  dVec tempDispBN;
  //  dVec dispBNew;
//   cerr << "A1 = " << Path(sliceA,ptcl1) << endl;
//   cerr << "A2 = " << Path(sliceA,ptcl2) << endl;
//   cerr << "B1 = " << Path(sliceB,ptcl1) << endl;
//   cerr << "B2 = " << Path(sliceB,ptcl2) << endl;
  int m;
  for (int i=0; i<NDIM; i++) {
    double n = -floor(dispA(i)*BoxInv(i)+0.5);
    dispA(i) += n*IsPeriodic(i)*Box(i);
//     double m = -floor(dispB(i)*BoxInv(i)+0.5);
//     dispB(i) += m*IsPeriodic(i)*Box(i);
    double mNew=-floor((dispA(i)-dispB(i))*BoxInv(i)+0.5);
    dispB(i)-= mNew*IsPeriodic(i)*Box(i);
//     // HACK HACK HACK
//     m=0;
//     tempDispB(i) = dispB(i)+m*IsPeriodic(i)*Box(i);
//     tempDispBN(i)=dispB(i)-m*IsPeriodic(i)*Box(i);
//     while (fabs(dispA(i)-tempDispB(i))>Box(i)/2.0 &&
// 	   fabs(dispA(i)-tempDispBN(i))>Box(i)/2.0){
//       m=m+1;
//       tempDispB(i) = dispB(i)+m*IsPeriodic(i)*Box(i);
//       tempDispBN(i)=dispB(i)-m*IsPeriodic(i)*Box(i);
//     }
//     if (fabs(dispA(i)-tempDispB(i))<=Box(i)/2.0)
//       dispB(i)=tempDispB(i);
//     else if (fabs(dispA(i)-tempDispBN(i))<=Box(i)/2.0)
//       dispB(i)=tempDispBN(i);
//     else cerr<<"ERROR! ERROR! ERROR!"<<endl;
//     if (fabs(dispBNew(i)-dispB(i))>1e-12){
//       cerr<<"dispBNew and dispB are not the same!\n";
//     }
//     //    double m = -floor(dispB(i)*BoxInv(i)+0.5);
//     //    dispB(i) += m*IsPeriodic(i)*Box(i);
//     //    cerr << "n = " << n << endl;
  }
//   cerr << "dispA = " << dispA << endl;
//   cerr << "dispB = " << dispB << endl;
  distA = sqrt(dot(dispA,dispA));
  distB = sqrt(dot(dispB,dispB));

#ifdef GARBAGEDEBUG
  dVec DBdispA = Path(sliceA, ptcl2) -Path(sliceA, ptcl1);
  dVec DBdispB = Path(sliceB, ptcl2) -Path(sliceB, ptcl1);
  for (int i=0; i<NDIM; i++) {
    while (DBdispA(i) > 0.5*Box(i)) 
      DBdispA(i) -= Box(i);
    while (DBdispA(i) < -0.5*Box(i)) 
      DBdispA(i) += Box(i);
    while ((DBdispB(i)-DBdispA(i)) > 0.5*Box(i))
      DBdispB -= Box(i);
    while ((DBdispB(i)-DBdispA(i)) < -0.5*Box(i))
      DBdispB += Box(i);
  }
//   cerr << "DBdispA = " << DBdispA << endl;
//   cerr << "DBdispB = " << DBdispB << endl;
  for (int i=0; i<NDIM; i++) {
    assert (fabs(DBdispA(i)-dispA(i)) < 1.0e-12);
    assert (fabs(DBdispB(i)-dispB(i)) < 1.0e-12);
  }
#endif
}

inline dVec 
PathClass::Velocity (int sliceA, int sliceB, int ptcl)
{

  //  bool changePtcl=(OpenPaths && sliceB==(int)OpenLink && ptcl==(int)OpenPtcl);
  //  ptcl=ptcl*!changePtcl+NumParticles()*changePtcl;
  dVec vel;
  if (OpenPaths && sliceB==(int)OpenLink && ptcl==(int)OpenPtcl){
    vel=Path(sliceB,NumParticles())-Path(sliceA,ptcl);
  }
  else{
    vel = Path(sliceB, ptcl) - Path(sliceA,ptcl);
  }
  PutInBox(vel);
  

/* #ifdef DEBUG */
/*  dVec DBvel = Path(sliceB,ptcl) - Path(sliceA,ptcl); */
/*    for (int dim=0; dim<NDIM; dim++) */
/*      { */
/*        while (DBvel[dim] > (0.5*Box[dim])) */
/*   	DBvel[dim] -= Box[dim]; */
/*        while (DBvel[dim] < -(0.5*Box[dim])) */
/* 	 DBvel[dim] += Box[dim]; */
/*        if (fabs(DBvel[dim]-vel[dim])<1e-12){ */
/* 	 cerr<<DBvel[dim]<<" "<<vel[dim]<<endl; */
/*        } */
/*        assert(fabs(DBvel[dim]-vel[dim])<1e-12); */
/*      } */

/* #endif */

  return vel;
}

inline void 
PathClass::PutInBox (dVec &v)
{
#ifdef DEBUG
  dVec Dv=v;
#endif
  for (int i=0; i<NDIM; i++) {
    double n = -floor(v(i)*BoxInv(i)+0.5);
    v(i) += n*IsPeriodic(i)*Box(i);
  }
/* #ifdef DEBUG */
/*   for (int dim=0; dim<NDIM; dim++){ */
/*     while (Dv[dim] > (0.5*Box[dim])) */
/*       Dv[dim] -= Box[dim]; */
/*     while (Dv[dim] < (-(0.5*Box[dim]))) */
/*       Dv[dim] += Box[dim]; */
/*     if (fabs(Dv[dim]-v[dim])>1e-12){ */
/* 	 cerr<<Dv[dim]<<" "<<v[dim]<<" "<<Box(dim)<<" "<<dim<<BoxInv(dim) */
/* 	     <<endl; */
/*     } */
/*     assert(fabs(Dv[dim]-v[dim])<1e-12); */
/*   } */
/* #endif  */
  
}


inline void 
PathClass::SliceRange(int proc, int &start, int &end)
{
  end = 0;
  int nProcs = Communicator.NumProcs();
  for (int i=0; i<=proc; i++) {
    start = end;
    int numSlices = TotalNumSlices/nProcs+(i<(TotalNumSlices % nProcs));
    end = start + numSlices;
  }
}

inline int 
PathClass::SliceOwner(int slice)
{
  int proc = 0;
  int nProcs = Communicator.NumProcs();
  for (int i=0; i<Communicator.NumProcs(); i++) {
    int slice1, slice2;
    SliceRange (i, slice1, slice2);
    if ((slice1 <= slice) && (slice2 >= slice))
      proc = i;
  }
  return proc;
}


#endif
