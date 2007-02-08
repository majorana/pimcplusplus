/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "MoleculeInteractionsClass.h"
#include "../PathDataClass.h"
#include "../Moves/MoveUtils.h"

MoleculeInteractionsClass::MoleculeInteractionsClass (PathDataClass &pathData) :
  ActionBaseClass (pathData)
{
  elementary_charge = 1.602*pow(10.0,-19);
  N_Avogadro = 6.022*pow(10.0,23.0);
  kcal_to_joule = 4184;
  epsilon_not = 8.85*pow(10.0,-12);
  angstrom_to_m = pow(10.0,10);
  SI = 1/(4*M_PI*epsilon_not);
  k_B = 1.3807*pow(10.0,-23);
  erg_to_eV = 1.6*pow(10.0,12);
  joule_to_eV = pow(10.0,19)/1.602;
  bohr_per_angstrom = 1.890359;

	// radial cutoffs for ST2 modulation function
  // these are in angstrom
	RL = 2.0160;
	RU = 3.1287;

  // these are initialized in Read() now
	// SPC/F2 intramolecular potential parameters
	//rho = 2.361;
	//D = 0.708;
	//alpha = 108.0*M_PI/180.0;
	//R_OH_0 = 1.0;
	//R_HH_0 = 2*R_OH_0*sin(alpha/2);
	//b = 1.803;
	//c = -1.469;
	//d = 0.776;
	//// conversion factor for mdyn*angstrom^-1 --> kcal*mol^-1*angstrom^-2
	//Dyn2kcal = 143.929;

	ReadComplete = false;
}

string MoleculeInteractionsClass::GetName(){
	return("MoleculeInteractionsClass");
}

double MoleculeInteractionsClass::SingleAction (int startSlice, int endSlice, 
	       const Array<int,1> &activeParticles, int level){

	//cerr << "MoleculeInteractions::Action__________________ for slices " << startSlice << " to " << endSlice;// << endl;
	assert(ReadComplete);
  if(startSlice == 0 && endSlice == 0){
    startSlice -= 1;
    endSlice += 1;
  }
	bool IsAction = true;
	double TotalU = ComputeEnergy(startSlice+1, endSlice-1, activeParticles, level, TruncateAction, IsAction);
	//cerr << " RETURNING " << TotalU*PathData.Path.tau << endl;
  return(TotalU*PathData.Path.tau);
}

double MoleculeInteractionsClass::d_dBeta (int startSlice, int endSlice,  int level)
{
	//cerr << "MoleculeInteractions::d_dBeta__________________ for slices " << startSlice << " to " << endSlice;// << endl;
	assert(ReadComplete);
  Array<int,1> activeParticles(PathData.Path.NumParticles());
  for (int i=0;i<PathData.Path.NumParticles();i++)
    activeParticles(i)=i;
	
	bool IsAction = false;
	double TotalU = ComputeEnergy(startSlice, endSlice-1, activeParticles, level, TruncateEnergy, IsAction);
	//cerr << " RETURNING " << TotalU << endl;
  return TotalU;
}

double MoleculeInteractionsClass::ComputeEnergy(int startSlice, int endSlice, 
	       const Array<int,1> &activeParticles, int level, bool with_truncations, bool isAction){

  //cerr << isAction << " MoleculeInteraction computing SPC energy over slices " << startSlice << " to " << endSlice << endl;
	Updated = false;
  for (int counter=0; counter<Path.DoPtcl.size(); counter++)
    Path.DoPtcl(counter)=true;

  double TotalU = 0.0;
	double TotalLJ = 0.0;
	double TotalCharge = 0.0;
	double TotalSpring = 0.0;
	double TotalKinetic = 0.0;
	double TotalHarmonic = 0.0;
  int numChangedPtcls = activeParticles.size();
  int skip = 1<<level;

	//cerr << "Hello.  This is MoleculeInteractionsClass::ComputeEnergy.  ActivePtcls are " << activeParticles << endl;
	//cerr << "I'm doing Intramolecular " << IntraMolecular << " and with_truncations " << with_truncations << " and with modulation " << withS << endl;
	//cerr << "CUTOFF is " << CUTOFF << endl;
  for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){
    int ptcl1 = activeParticles(ptcl1Index);
    Path.DoPtcl(ptcl1) = false;
		//cerr << "DoPtcl bools are after killing " << ptcl1 << Path.DoPtcl << endl;
    int species1=Path.ParticleSpeciesNum(ptcl1);

		//	Calculate Lennard-Jones interaction
    if (Interacting(species1,0)){
			//cerr << "LJ for ptcl " << ptcl1 << " of species " << species1 << endl;
  		for (int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++){
    		int species2=Path.ParticleSpeciesNum(ptcl2);
				//cerr << "Considering LJ between " << ptcl1 << " and " << ptcl2 << " for which DoPtcl is " << Path.DoPtcl(ptcl2) << "...";
    		if (Interacting(species2,0) && Path.DoPtcl(ptcl2)){
					//cerr << " Going to compute between " << ptcl1 << " and " << ptcl2;
	  			for (int slice=startSlice; slice<=endSlice; slice+=skip){
	    			dVec r;
	    			double rmag;
	    			PathData.Path.DistDisp(slice, ptcl1, ptcl2, rmag, r);
						//cerr << "... rmag is " << rmag;
						//if(ptcl1==Path.MolRef(ptcl1) && ptcl2==Path.MolRef(ptcl2) && !Updated(ptcl1,ptcl2)){
						//	Updated(ptcl1,ptcl2) = Updated(ptcl2,ptcl1) = true;
						//	COMTable(ptcl1,ptcl2) = COMTable(ptcl2,ptcl1) = rmag;
						//	COMVecs(ptcl1,ptcl2) = r;
						//	COMVecs(ptcl2,ptcl1) = -1*r;
						//}
						//	disregard interactions outside spherical cutoff
            if (rmag <= CUTOFF){
							//cerr << "<"<<CUTOFF;
              double rinv = 1.0/rmag;
              double sigR = PathData.Species(species1).Sigma*rinv;
              double sigR6 = sigR*sigR*sigR*sigR*sigR*sigR;
							double offset = 0.0;
							if(with_truncations){
  							double sigma_over_cutoff = PathData.Species(species1).Sigma/CUTOFF;
  							offset = pow(sigma_over_cutoff,12) - pow(sigma_over_cutoff,6);
							}
							//cerr << "; using offset " << offset << endl;
	      			double lj = conversion*4*PathData.Species(species1).Epsilon*(sigR6*(sigR6-1) - offset); // this is in kcal/mol 
							TotalLJ += lj;
	      			TotalU += lj;
							//cerr  << TotalU << " added " << lj << " from LJ interaction between " << ptcl1 << " and " << ptcl2 << endl;
            }
						//cerr << "... outside" << endl;
	  			}
				}
				else{
					//cerr << "Skipped" << endl;
				}
      }
    }
		if(Interacting(species1,1)){
      /// calculating coulomb interactions

      /// hack: we're going to sum over the first set of images too
  		for (int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++){
    		int species2=Path.ParticleSpeciesNum(ptcl2);
    		if (Interacting(species2,1) && Path.DoPtcl(ptcl2)){
					//	don't compute intramolecular interactions
					//  unless told to
					if(IntraMolecular || Path.MolRef(ptcl1)!=Path.MolRef(ptcl2)){
	  				for (int slice=startSlice;slice<=endSlice;slice+=skip){
	    				double Ormag;
	    				dVec Or;
	    				PathData.Path.DistDisp(slice,Path.MolRef(ptcl1),Path.MolRef(ptcl2),Ormag,Or);
              dVec L = Path.GetBox();
						  for (int x=-1; x<=1; x++) {
						    for (int y=-1; y<=1; y++) {
						      for (int z=-1; z<=1; z++) {
                    Or(0) += x*L(0);
                    Or(1) += y*L(1);
                    Or(2) += z*L(2);
                    Ormag = Mag(Or);
							      // implement spherical cutoff
            	      //double Ormag = COMSeparation(slice,ptcl1,ptcl2);
            	      if (Ormag <= CUTOFF){
                      dVec r = Or - (Path(slice,ptcl1) - Path(slice,Path.MolRef(ptcl1)))
                        + (Path(slice,ptcl2) - Path(slice,Path.MolRef(ptcl2)));
                      double rmag = Mag(r);
							      	double truncate = 0.0;
                    	double ptclCutoff = 1.0;
							      	if(with_truncations){
							      		truncate = 1.0;
							      		//ptclCutoff = CalcCutoff(ptcl1,ptcl2,slice,CUTOFF);
							      		ptclCutoff = Mag(r - Or + Scale(Or,CUTOFF));
							      	}
							      	double modulation=1.0;
							      	if(withS)
							      		modulation = S(Ormag);
                    	double coulomb = conversion*prefactor*PathData.Species(species1).pseudoCharge
							      									*PathData.Species(species2).pseudoCharge*modulation
                                      *(1.0/rmag - truncate/ptclCutoff);
							      	TotalCharge += coulomb;
	      			      	TotalU += coulomb;
							      	//cerr  << TotalU << " added " << coulomb << " from charge-charge interaction between " << ptcl1 << " and " << ptcl2 << endl;
  	                }
                  }
                }
						  }
  	        }
					}
  	    }
  	  }
  	}
		if(Interacting(species1,2) && (ptcl1 == PathData.Path.MolRef(ptcl1))){
			// quadratic intramolecular potential
			// SPC/F2; see Lobaugh and Voth, JCP 106, 2400 (1997)
			double spring = 0.0;
	  	for (int slice=startSlice;slice<=endSlice;slice+=skip){
				vector<int> activeP(0);
				activeP.push_back(ptcl1);
				//cerr << " INTRA added " << ptcl1;
  			for (int index2=1; index2<PathData.Path.MolMembers(ptcl1).size(); index2++){
					int ptcl2 = PathData.Path.MolMembers(ptcl1)(index2);
    			int species2=Path.ParticleSpeciesNum(ptcl2);
    			if (Interacting(species2,2) && Path.DoPtcl(ptcl2)){
						activeP.push_back(ptcl2);
						//cerr << " INTRA added " << ptcl2;
					}
				}
        assert(activeP.size() == 3);
				dVec r;
				double ROH1, ROH2, RHH;
				PathData.Path.DistDisp(slice,activeP[0],activeP[1],ROH1,r);
				PathData.Path.DistDisp(slice,activeP[0],activeP[2],ROH2,r);
				PathData.Path.DistDisp(slice,activeP[1],activeP[2],RHH,r);

				double term1 = rho*rho*D*((ROH1 - R_OH_0)*(ROH1 - R_OH_0) + (ROH2 - R_OH_0)*(ROH2 - R_OH_0));
				double term2 = 0.5*b*(RHH - R_HH_0)*(RHH - R_HH_0);
				double term3 = c*(ROH1 + ROH2 - 2*R_OH_0)*(RHH - R_HH_0);
				double term4 = d*(ROH1 - R_OH_0)*(ROH2 - R_OH_0);
				//cerr << "slice " << slice << " INTRA: ROH1 " << ROH1 << " ROH2 " << ROH2 << " RHH " << RHH << " term1 " << term1 << " term2 " << term2 << " term3 " << term3 << " term4 " << term4 << endl;
				spring = conversion*Dyn2kcal*(term1 + term2 + term3 + term4);
				outfile << spring << " " << term1 << " " << term2 << " " << term3 << " " << term4 << endl;
				TotalSpring += spring;
				TotalU += spring;
			}
		}
		if(Interacting(species1,3)){
			// compute kinetic action: imaginary time spring term
			// Based on KineticClass::SingleAction, execpt that displacements are WRT the molecule COM

			//cerr << "species " << species1 << ", ptcl " << ptcl1 << " computing kinetic...";
		  double TotalK = 0.0;
		  int skip = 1<<level;
		  double levelTau = Path.tau* (1<<level);
		  double lambda = lambdas(species1);
		  if (lambda != 0.0){
		    double FourLambdaTauInv=1.0/(4.0*lambda*levelTau*levelTau);
        if(isAction){
		      for (int slice=startSlice-1; slice<= endSlice;slice+=skip) {
		        dVec vel;
				  	//vel = COMVelocity(slice, slice+skip, ptcl1);
				  	vel = PathData.Path.Velocity(slice, slice+skip, ptcl1);
			
		        double GaussProd = 1.0;
		        for (int dim=0; dim<NDIM; dim++) {
				  		double GaussSum=0.0;
				  		for (int image=-NumImages; image<=NumImages; image++) {
				  		  double dist = vel[dim]+(double)image*Path.GetBox()[dim];
				  		  GaussSum += exp(-dist*dist*FourLambdaTauInv);
				  		}
				  		GaussProd *= GaussSum;
		        }
				  	TotalK -= log(GaussProd);
		      }
        } else {
		      for (int slice=startSlice; slice <= endSlice;slice+=skip) {
		        dVec vel;
				  	//vel = COMVelocity(slice, slice+skip, ptcl1);
				  	vel = PathData.Path.Velocity(slice, slice+skip, ptcl1);
			
		        double GaussProd = 1.0;
		        for (int dim=0; dim<NDIM; dim++) {
				  		double GaussSum=0.0;
				  		for (int image=-NumImages; image<=NumImages; image++) {
				  		  double dist = vel[dim]+(double)image*Path.GetBox()[dim];
				  		  GaussSum += exp(-dist*dist*FourLambdaTauInv);
				  		}
				  		GaussProd *= GaussSum;
		        }
				  	TotalK -= log(GaussProd);
		      }
        }
		  }
			TotalKinetic += TotalK;
			TotalU += TotalK;
			//cerr << TotalKinetic << endl;
		}
		if(Interacting(species1,4)){
			// harmonic intermolecular potential
			// hard-wired parameters for ST2 water dimer
			// could be generalized
			double harmonic = 0.0;
			double omega = 26; // ps^-1
			// for Rossky ST2 in kcal/mol*s^2 m^-2
			//double m_H2O = 0.043265; // ???? this seems off by a factor of 2
      //double m_H2O = 0.02163; // kcal/mol*ps^2 angstrom^-2
      double m_H2O = 0.0060537; // kcal/mol*ps^2 bohr^-2
			// Lobaugh & Voth in amu
			//double m_H2O = 9.0;
			//double R0 = 2.85; // angstrom
			double R0 = 5.39; // bohr
  		for (int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++){
    		int species2=Path.ParticleSpeciesNum(ptcl2);
    		if (Interacting(species2,4) && Path.DoPtcl(ptcl2)){
	  			for (int slice=startSlice;slice<=endSlice;slice+=skip){
  					dVec COMr;
  					double COMrmag;
  					PathData.Path.DistDisp(slice, ptcl1, ptcl2, COMrmag, COMr);
						//cerr << "Harmonic: COMrmag is " << COMrmag << " between " << ptcl1 << " and " << ptcl2 << " and R0 is " << R0 << endl;
						harmonic = conversion*0.5*m_H2O*omega*omega*(COMrmag - R0)*(COMrmag - R0);
						TotalHarmonic += harmonic;
						TotalU += harmonic;
					}
				}
			}
		}
	}
  //cerr << " returning energy " << TotalU << " and action " << TotalU*PathData.Path.tau << endl;
	//cerr << "Empirical potential contributions: " << TotalU << " " << TotalLJ << " " << TotalCharge << " " << TotalSpring << " " << TotalHarmonic << " " << TotalKinetic;// << endl;
	// write additional output file

	//if(!isAction && special){
	//	outfile << TotalU << " " << TotalLJ << " " << TotalCharge << " " << TotalSpring << " " << TotalHarmonic << " " << TotalKinetic << endl;
	//}
  return (TotalU);
}


double MoleculeInteractionsClass::CalcCutoff(int ptcl1, int ptcl2, int slice, double Rcmag){
  // get oxygen particle ids
  int Optcl1 = Path.MolRef(ptcl1);
  int Optcl2 = Path.MolRef(ptcl2);
  // get vectors of oxygens
  dVec O1 = Path(slice,Optcl1);
  dVec O2 = Path(slice,Optcl2);
  // get vector between oxygens
  dVec Roo;
  double Ormag;
	//if(!Updated(Optcl1,Optcl2)){
  	Path.DistDisp(slice, Optcl1, Optcl2, Ormag, Roo);
	//	Updated(Optcl1,Optcl2) = Updated(Optcl2,Optcl1) = true;
	//	COMTable(Optcl1,Optcl2) = COMTable(Optcl2,Optcl1) = Ormag;
	//	COMVecs(Optcl1,Optcl2) = Roo;
	//	COMVecs(Optcl2,Optcl1) = -1*Roo;
	//} else {
	//	Ormag = COMTable(Optcl1,Optcl2);
	//	Roo = COMVecs(Optcl1,Optcl2);
	//}
  dVec Rc = Scale(Roo,Rcmag);
  // get constituent coordinates WRT oxygen COM
  dVec P1 = Path(slice,ptcl1);
  dVec P2 = Path(slice,ptcl2);
  P1 -= O1;
  P2 -= O2;
  // solve for vector between constituents (if molecule 2 is at cutoff radius)
  dVec R12 = Rc + P2 - P1;
  // obtain cutoff magnitude
  double r12 = Mag(R12);
  return r12;
}

double MoleculeInteractionsClass::COMSeparation (int slice,int ptcl1,int ptcl2)
{
  // get oxygen particle ids
  int Optcl1 = Path.MolRef(ptcl1);
  int Optcl2 = Path.MolRef(ptcl2);
  dVec Or;
  double Ormag;

	//if(!Updated(Optcl1,Optcl2)){
  	PathData.Path.DistDisp(slice, Optcl1, Optcl2, Ormag, Or);
//		Updated(Optcl1,Optcl2) = Updated(Optcl2,Optcl1) = true; COMTable(Optcl1,Optcl2) = COMTable(Optcl2,Optcl1) = Ormag;
//		COMVecs(Optcl1,Optcl2) = Or;
//		COMVecs(Optcl2,Optcl1) = -1*Or;
//	} else {
//		Ormag = COMTable(Optcl1,Optcl2);
//	}

  return Ormag;
}

double MoleculeInteractionsClass::S(double r)
{
  double mod;
  if(r<RL)
    mod = 0.0;
  else if(r>RU)
    mod = 1.0;
  else{
    double diff1 = r - RL;
    double diff2 = 3*RU - RL - 2*r;
    double diff3 = RU - RL;
    mod = diff1*diff1*diff2/(diff3*diff3*diff3);
  }
  return mod;
}

dVec MoleculeInteractionsClass::COMVelocity(int sliceA, int sliceB, int ptcl){
	dVec p1 = Path(sliceB, ptcl);
	dVec p2 = Path(sliceA, ptcl);
	if(ptcl != PathData.Path.MolRef(ptcl)){
		int COMptcl = PathData.Path.MolRef(ptcl);
 		p1 -= Path(sliceB, COMptcl);
 		p2 -= Path(sliceA, COMptcl);
	}
 	dVec vel = p1 - p2;
  PathData.Path.PutInBox(vel);
	return vel;
}

void MoleculeInteractionsClass::Read (IOSectionClass &in)
{
	if(!ReadComplete){
		cerr << "In MoleculeInteractionsClass::Read" << endl;
		// DEFAULTS
		prefactor = SI*angstrom_to_m*elementary_charge*
								elementary_charge*N_Avogadro/kcal_to_joule;

    // SPC/F2 intramolecular potential parameters
    // Lobaugh and Voth, JCP 106 2400 (1997)
    // these are the stated parameters from the paper in units of angstrom
	  rho = 2.361;
	  D = 0.708;
	  alpha = 108.0*M_PI/180.0;
	  R_OH_0 = 1.0;
	  R_HH_0 = 2*R_OH_0*sin(alpha/2);
	  b = 1.803;
	  c = -1.469;
	  d = 0.776;
	  // conversion factor for mdyn*angstrom^-1 --> kcal*mol^-1*angstrom^-2
	  Dyn2kcal = 143.929;

    string units = "angstrom";
    in.ReadVar("Units",units);
    if(units == "bohr"){
      prefactor *= bohr_per_angstrom;

      // these are converted into bohr from angstrom
	    rho = 1.24897;
	    D = 1.338;
	    //R_OH_0 = 1.8904;
      // optimized value
	    R_OH_0 = 1.8696;
	    //R_HH_0 = 2*R_OH_0*sin(alpha/2);
      // optimized value
	    R_HH_0 = 2.98;
	    b = 0.9538;
	    c = -0.7771;
	    d = 0.4105;
	    // conversion factor for mdyn*bohr^-1 --> kcal*mol^-1*bohr^-2
	    Dyn2kcal = 76.138;
    }

		CUTOFF = Path.GetBox()(0)/2;
		IntraMolecular = false;
		withS = true;
		TruncateAction = true;
		TruncateEnergy = false;

    conversion = 1.0;
		if(in.ReadVar("Prefactor",conversion)){
      cerr << "Setting conversion factor to " << conversion << ". Default units are kcal/mol with length in " << units << endl;
    }
		in.ReadVar("Cutoff",CUTOFF);
		in.ReadVar("Modulated",withS);
		in.ReadVar("Intramolecular",IntraMolecular);
		in.ReadVar("TruncateAction",TruncateAction);
		in.ReadVar("TruncateEnergy",TruncateEnergy);

		Interacting.resize(PathData.NumSpecies(),5);
		Interacting = false;
		LJSpecies.resize(0);
		ChargeSpecies.resize(0);
		SpringSpecies.resize(0);
		KineticSpecies.resize(0);
		in.ReadVar("LJSpecies",LJSpecies);
		in.ReadVar("ChargeSpecies",ChargeSpecies);
		in.ReadVar("IntraMolecularSpecies",SpringSpecies);
		in.ReadVar("KineticActionSpecies",KineticSpecies);
		in.ReadVar("QuadraticSpecies",QuadSpecies);

		cerr << "Read LJSpec " << LJSpecies << " and chargeSpec " << ChargeSpecies << " etc..." << endl;
		if(KineticSpecies.size() > 0){
			assert(in.ReadVar("Lambdas",lambdas));
			cerr << "Read lambda for each species: " << lambdas << endl;
		}

		for(int s=0; s<LJSpecies.size(); s++){
			Interacting(Path.SpeciesNum(LJSpecies(s)), 0) = true;
			cerr << "Setting " << LJSpecies(s) << " to interact via LJ!!" << endl;
		}
		for(int s=0; s<ChargeSpecies.size(); s++)
			Interacting(Path.SpeciesNum(ChargeSpecies(s)), 1) = true;
		for(int s=0; s<SpringSpecies.size(); s++)
			Interacting(Path.SpeciesNum(SpringSpecies(s)), 2) = true;
		for(int s=0; s<KineticSpecies.size(); s++)
			Interacting(Path.SpeciesNum(KineticSpecies(s)), 3) = true;
		for(int s=0; s<QuadSpecies.size(); s++)
			Interacting(Path.SpeciesNum(QuadSpecies(s)), 4) = true;
		cerr << "MoleculeInteractions::Read I have loaded Interacting table " << Interacting << endl;

		// make sure something is filled
		bool empty = true;
		for(int i=0; i<Interacting.rows(); i++){
			for(int j=0; j<Interacting.cols(); j++){
				if(Interacting(i,j) != 0){
					empty = (empty && false);
				}
			}
		}
		assert(!empty);
				
	
		Updated.resize(PathData.Path.numMol);
		COMTable.resize(PathData.Path.numMol);
		COMVecs.resize(PathData.Path.numMol);

		special = false;
		in.ReadVar("ExtraOutput",special);
		if(special){
			string filename;
			assert(in.ReadVar("File",filename));
			outfile.open(filename.c_str());
			//outfile << "# Total LJ Coulomb Intramolecular Quadratic Kinetic" << endl;
			outfile << "# Intramolecular Int1 Int2 Int3 Int4" << endl;
		}
			
		
		ReadComplete = true;
	}
}

void MoleculeInteractionsClass::SetNumImages (int num)
{
	NumImages = num;
}
