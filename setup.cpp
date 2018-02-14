/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * setup.cpp
 * Created by ljohnson on 5/31/09.
 *	Setup functions for Robinson MC code. These are simply intended to clean up the formerly-bloated MCmain.cpp
 *	Anything that the user doesn't edit should be here unless it would require passing a lot of variables around.
 *	Anything the user edits should remain in MCmain.cpp
 *
 * updated by Andreas Tillack on Jan 25, 2011
 * - multiple molecule type support
 */

#include "setup.h"

inline double EPScondition(double EPS_CUT, double C, double r, double delta, unsigned int n)
{
	return EPS_CUT-(C/qpwr(r-delta,n)-C/qpwr(r,n));
}

// numerically solves EPS=C*(1.0/(r-delta)^n-1.0/r^n) for r
inline double optcut(double EPS_CUT, double C, double delta, unsigned int n)
{
	if(C>EPS_CUT){
		double min=delta+EPS;
		double max=1.0/EPS;
		double result=(min+max)/2.0;
		while(fabs(max-min)*2.0>EPS){
			if((EPScondition(EPS_CUT,C,min,delta,n)>0.0) || (EPScondition(EPS_CUT,C,max,delta,n)<0.0)){ // safety check
				cout << "ERROR: Can not find zero crossing.\n";
				return 0.0;
			}
			if(EPScondition(EPS_CUT,C,result,delta,n)<0.0) min=result; else max=result;
			result=(min+max)/2.0;
		}
		return result;
	} else return 0.0;
}

inline double conditionLJ(double sigma, double r, unsigned int n, unsigned int m)
{
	return n/sigma*qpwr(sigma/r,n+1)-m/sigma*qpwr(sigma/r,m+1);
}

inline double LJmin(double sigma, unsigned int n, unsigned int m)
{
	double min=sigma+EPS;
	double max=1.0/EPS;
	double result=(min+max)/2.0;
	while(fabs(max-min)*2.0>EPS){
		if((conditionLJ(sigma,min,n,m)>0.0) || (conditionLJ(sigma,max,n,m)<0.0)){ // safety check
			cout << "ERROR: Can not find zero crossing.\n";
			return sigma;
		}
		if(conditionLJ(sigma,result,n,m)<0.0) min=result; else max=result;
		result=(min+max)/2.0;
	}
	return result;
}

inline double EPSconditionLJ(double EPS_CUT, double C, double sigma, double r, double delta, unsigned int n, unsigned int m)
{
	return EPS_CUT+C*((qpwr(sigma/(r-delta),m)-qpwr(sigma/(r-delta),n))-(qpwr(sigma/r,m)-qpwr(sigma/r,n)));
}

// numerically solves EPS=C*[((sigma/(r-delta))^m-(sigma/(r-delta))^n)-((sigma/r)^m-(sigma/r)^n] for r
inline double optcutLJ(double EPS_CUT, double C, double sigma, double delta, unsigned int n, unsigned int m)
{
	if(C>EPS_CUT){
		double min=2.0*LJmin(sigma,n,m)-sigma+delta+EPS;
		double max=min+1.0/EPS;
		double result=(min+max)/2.0;
		while(fabs(max-min)*2.0>EPS){
			if((EPSconditionLJ(EPS_CUT,C,sigma,min,delta,n,m)>0.0) || (EPSconditionLJ(EPS_CUT,C,sigma,max,delta,n,m)<0.0)){ // safety check
				cout << "ERROR: Can not find zero crossing.\n";
				return 0.0;
			}
			if(EPSconditionLJ(EPS_CUT,C,sigma,result,delta,n,m)<0.0) min=result; else max=result;
			result=(min+max)/2.0;
		}
		return result;
	} else return 0.0;
}

/// List of physical constants. Since kB isn't likely to change in the near future, the user doesn't need to edit these.
void phys_configuration(Config_Data* configuration)
{
	//Load physical constants
	configuration->trackchange = 0; //counter for changing E and T. Should always be 0 to start.
	configuration->changesteps = 0; //steps where E and T are changed. Must be same length as arrays below.
	configuration->kT = kB*configuration->T;
	configuration->Ekin=0.0;
	for(unsigned int i=0; i<configuration->num_element_types; i++) configuration->Ekin += configuration->element_types[i]->number*configuration->element_types[i]->dof/2.0*configuration->kT; //elements
	for(unsigned int i=0; i<configuration->num_groups; i++) configuration->Ekin += configuration->groups[i]->number*configuration->groups[i]->Type->dof/2.0*configuration->kT; //groups
	configuration->beta = 1.0/configuration->kT;
}

/// Load vdW parameters into configuration
void setup_vdw(Config_Data* configuration)
{
	//Set up vdw cutoff
	double bltemp;
	if (configuration->latticetype < 3)
	{
		bltemp = configuration->boxlength[0];
		if (configuration->boxlength[1] < bltemp) bltemp = configuration->boxlength[1];
		if (configuration->boxlength[2] < bltemp) bltemp = configuration->boxlength[2];
		bltemp /= 2.0;
	}
	else if (configuration->latticetype > 2)
	{
		bltemp = configuration->spherecylr;
		if (configuration->latticetype == 4)
		{
			if (configuration->boxlength[2]/2.0 < bltemp) bltemp = configuration->boxlength[2]/2.0;
		}
	}
	if (configuration->ljcut > bltemp)
	{
		configuration->ljcut = bltemp;
		configuration->ljcut2 = bltemp*bltemp;
		printf("Warning: Lennard-Jones cutoff is too large. Resizing to half of the box length (%f)\n", bltemp);
	}
	
	//Calculate parameters for LJ correction.
	unsigned int i,j;
	double m = (double)(configuration->LJexp[0])/2.0;
	double n = (double)(configuration->LJexp[1])/2.0;
	double cfshift = pow(n/m,1/(n-m)); //shift used with solvent parameters
	double cr = n/m;
	double cr_scale = pow(cr,cr/(cr-1.0))/(cr-1.0); //scaling coefficient (4 for 12-6 potential)
	configuration->fshift = cfshift;
	configuration->r = cr_scale;
	// configuration->r is used to scale Lennard-Jones interaction only, hence can be used for energy scaling
	configuration->r *= configuration->energy_scale[2];
	
	//Precalculate vdw epsilon and ro values, as well as average distances for if Touch is not used
	//Will need to be replaced with an arbitrary matrix if more than two types of particles are allowed
	// -- done :-) AT Jan 25, 2011 late at night/early morning
	configuration->pre_sigma2 = new double[configuration->num_element_types*configuration->num_element_types];
	configuration->pre_eps = new double[configuration->num_element_types*configuration->num_element_types+2*configuration->num_element_types]; // reserve space for eps values of LJ wall
	configuration->pre_touch = new bool[configuration->num_element_types*configuration->num_element_types];
	for(i=0; i<configuration->num_element_types; i++){ // not the most efficient code, but since it only gets called once ...
		for(j=0; j<configuration->num_element_types; j++){
			unsigned int ikk=i*configuration->num_element_types+j;
			// Mixing rule for LJ radius: 0 .. arithmetic (sigma_ij = r_i+r_j), 1 .. geometric (sigma_ij = 2*sqrt(r_i*r_j)), 2 .. sixth power (sigma_ij = 2*((r_i^6+r_j^6)/2)^(1/6) (default: 0)
			switch(configuration->LJ_radius_mixing){
				default:
				case 0:
					configuration->pre_sigma2[ikk]=qqpwr(configuration->element_types[i]->radius+configuration->element_types[j]->radius,2);
#if DEBUG_LEVEL>3
					cout << i << ", " << j << " -> " << ikk << ": (" << configuration->element_types[i]->radius << " + " << configuration->element_types[j]->radius << ")^2 = " << configuration->pre_sigma2[ikk] << "\n";
#endif
					break;
				case 1:
					configuration->pre_sigma2[ikk]=4.0*configuration->element_types[i]->radius*configuration->element_types[j]->radius;
					break;
				case 2:
					configuration->pre_sigma2[ikk]=4.0*cbrt((qqpwr(configuration->element_types[i]->radius,6)+qqpwr(configuration->element_types[j]->radius,6))/2.0);
					break;
			}
			configuration->pre_eps[ikk]=determine_epsilon(configuration->element_types[i],configuration->element_types[j]->rT)*determine_epsilon(configuration->element_types[j],configuration->element_types[i]->rT);
			if(configuration->LJ_epsilon_sixthpower) configuration->pre_eps[ikk]*=2*qqpwr(configuration->element_types[i]->rT*configuration->element_types[j]->rT,3)/(qqpwr(configuration->element_types[i]->rT,6)+qqpwr(configuration->element_types[j]->rT,6));
			// Precalculate whether or not touch is used for a given pair of particle types
			bool sphere_i=((fabs(configuration->element_types[i]->saxes.vec[0]-configuration->element_types[i]->saxes.vec[1])<EPS) && (fabs(configuration->element_types[i]->saxes.vec[1]-configuration->element_types[i]->saxes.vec[2])<EPS));
			bool sphere_j=((fabs(configuration->element_types[j]->saxes.vec[0]-configuration->element_types[j]->saxes.vec[1])<EPS) && (fabs(configuration->element_types[j]->saxes.vec[1]-configuration->element_types[j]->saxes.vec[2])<EPS));
			configuration->pre_touch[ikk]=!(sphere_i && sphere_j);
		}
	}
	double maxq=0.0;
	double maxdipole=0.0;
	double maxeps=0.0;
	double maxsigma=0.0;
	for(i=0; i<configuration->num_element_types; i++){
		if(configuration->element_types[i]->number>0){
			if(configuration->element_types[i]->Vvdw>maxeps) maxeps=configuration->element_types[i]->Vvdw;
			if(maxd(configuration->element_types[i]->saxes.vec,3)>maxsigma) maxsigma=maxd(configuration->element_types[i]->saxes.vec,3);
			if(fabs(fabs(configuration->element_types[i]->dipole)-maxdipole)>EPS) maxdipole=fabs(configuration->element_types[i]->dipole);
			for(j=0; j<configuration->element_types[i]->nr_charges; j++){
				if(fabs(configuration->element_types[i]->q[j])>maxq) maxq=fabs(configuration->element_types[i]->q[j]);
			}
		}
		configuration->pre_eps[configuration->num_element_types*configuration->num_element_types+2*i]=sqrt(configuration->element_types[i]->Vvdw*configuration->LJwall_a);
		configuration->pre_eps[configuration->num_element_types*configuration->num_element_types+2*i+1]=sqrt(configuration->LJwall_b/configuration->LJwall_a);
	}
	for(i=0; i<configuration->num_groups; i++){
		if(configuration->groups[i]->number>0){
//			double overlap_lambda=0.5*configuration->groups[i]->Type->adjust_overlap;
			for(j=0; j<configuration->groups[i]->nr_elements; j++){
				Element_Type* et=configuration->group_elements[configuration->groups[i]->elements[j]].MyType;
/*				if(configuration->groups[i]->levelofdetail && (overlap_lambda>0.0)){ // take care of overlaps for LOD groups
					unsigned int type=configuration->group_elements[configuration->groups[i]->elements[j]].mytype;
					et->Vvdw*=1.0-overlap_lambda*configuration->groups[i]->Type->overlaps[j];
					for(unsigned int k=0; k<configuration->num_element_types; k++){
						configuration->pre_eps[type*configuration->num_element_types+k]*=sqrt(1.0-overlap_lambda*configuration->groups[i]->Type->overlaps[j]);
						configuration->pre_eps[k*configuration->num_element_types+type]*=sqrt(1.0-overlap_lambda*configuration->groups[i]->Type->overlaps[j]);
					}
				}*/
				if(et->Vvdw>maxeps) maxeps=et->Vvdw;
				if(maxd(et->saxes.vec,3)>maxsigma) maxsigma=maxd(et->saxes.vec,3);
				if(fabs(et->dipole)>maxdipole) maxdipole=fabs(et->dipole);
				for(unsigned int k=0; k<et->nr_charges; k++) if(fabs(et->q[k])>maxq) maxq=fabs(et->q[k]);
			}
		}
	}
	maxsigma*=2.0;
	configuration->opt_cutoffs=Vec4(optcut(configuration->ES_cutoff_precision,maxq*maxq/configuration->n2,configuration->maxtrans/2.0,1),
					optcut(configuration->ES_cutoff_precision,maxq*maxdipole/configuration->n2,configuration->maxtrans/2.0,2),
					optcut(configuration->ES_cutoff_precision,maxdipole*maxdipole/configuration->n2,configuration->maxtrans/2.0,3),
					optcutLJ(configuration->LJ_cutoff_precision,configuration->r*maxeps,maxsigma,configuration->maxtrans/2.0,configuration->LJexp[0],configuration->LJexp[1]));
	if(configuration->cut) cout << "Optimal cutoff distances (charge-charge, charge-dipole, dipole-dipole, LJ): " << configuration->opt_cutoffs.V4Str(',') << "\n";
	bool adjust_optcutoffs=false;
	if(configuration->opt_cutoffs.vec[0]>configuration->escut){
		configuration->opt_cutoffs.vec[0]=configuration->escut;
		adjust_optcutoffs=true;
	}
	if(configuration->opt_cutoffs.vec[1]>configuration->escut){
		configuration->opt_cutoffs.vec[1]=configuration->escut;
		adjust_optcutoffs=true;
	}
	if(configuration->opt_cutoffs.vec[2]>configuration->escut){
		configuration->opt_cutoffs.vec[2]=configuration->escut;
		adjust_optcutoffs=true;
	}
	if(configuration->opt_cutoffs.vec[3]>configuration->ljcut){
		configuration->opt_cutoffs.vec[3]=configuration->ljcut;
		adjust_optcutoffs=true;
	}
	if(adjust_optcutoffs && configuration->cut){
		cout << "\tAdjusted to user-specified cutoff distances: " << configuration->opt_cutoffs.V4Str(',') << "\n";
	}
	if(configuration->cut){
		if(configuration->opt_cutoffs.vec[3]>EPS) configuration->ljcut=configuration->opt_cutoffs.vec[3];
		configuration->ljcut2=configuration->ljcut*configuration->ljcut;
		if((maxq>EPS) && (maxd(configuration->opt_cutoffs.vec,3)>EPS)){
			configuration->escut=maxd(configuration->opt_cutoffs.vec,3);
			cout << "\tFree charge carriers are present\n\t\tDetermining electrostatics cutoff from maximum of charge-charge, charge-dipole, and dipole-dipole cutoff: " << configuration->escut << "\n";
		} else{
			if(configuration->opt_cutoffs.vec[2]>EPS){
				configuration->escut=configuration->opt_cutoffs.vec[2];
				cout << "\tDetermine electrostatics cutoff from dipole-dipole cutoff: " << configuration->escut << "\n";
			}
		}
		cout << "-> Cutoff distances used (ES, LJ): " << configuration->escut << ", " << configuration->ljcut << "\n";
		configuration->escut2=configuration->escut*configuration->escut;
	}
	
	double nnscale = 1.25; //Scaling factor for outer edge of nearest-neighbor shell.
	
	//Set up Touch if it is used
	if(configuration->vdwtype>=4){
		//Warn about performance - I couldn't resist but build on Lewis' spelling of "performance" ;-) AT
		printf("Warning: Touch is enabled — 'Preformance will be slower', says Ford Prefect. Marvin agrees.\n");
		
		configuration->touchtrunc = 1; //Uses Touch only at short range for efficiency
		
		//Calculate Touch truncation distance
		if(configuration->touchtrunc){
			//Determine maximum and minimum dimensions of two oids
			bool adjust_cut=false;
			double miT, maT;
			double trs = 1.5; //Touch ratio scaling factor. Adjust only if preformance is poor.
			configuration->tcut2 = new double[configuration->num_element_types*configuration->num_element_types];
			int ij;
			for(i=0; i<configuration->num_element_types; i++){
				for(j=0; j<configuration->num_element_types; j++){
					ij=i*configuration->num_element_types+j;
					miT = min(mind(configuration->element_types[i]->saxes.vec,3),mind(configuration->element_types[j]->saxes.vec,3));
					maT = max(maxd(configuration->element_types[i]->saxes.vec,3),maxd(configuration->element_types[j]->saxes.vec,3));
					//Set Touch cutoff based on eccentricity of ellipses. Cutoff will get larger as ellipses become more eccentric.
					configuration->tcut2[ij] = 2.0*trs*maT*maT/miT;
					configuration->tcut2[ij] *= configuration->tcut2[ij];
					if(configuration->tcut2[ij] > configuration->ljcut2){
						configuration->tcut2[ij] = configuration->ljcut2;
						adjust_cut=true;
#if DEBUG_LEVEL>2
						if(i<j) cout << configuration->element_types[i]->name << " <-> " << configuration->element_types[j]->name << " touch radius is larger than LJ cutoff (" << sqrt(configuration->tcut2[ij]) << ")\n";
#endif
					}
				}
			}
			if(adjust_cut){
				cout << "Warning: Touch truncation radius is larger than LJ cutoff for";
#if DEBUG_LEVEL>2
				cout << " element interactions listed above.";
#else
				cout << " some elements.";
#endif
				cout << " Using Touch for these LJ interactions.\n";
			}
		}
	}
	//Define nearest neighbor shell maximum
	double maxshell=configuration->element_types[0]->saxes.vec[0];
	for(i=1; i<configuration->num_element_types; i++){
		maxshell=max(maxshell,maxd(configuration->element_types[i]->saxes.vec,3));
	}
	
	configuration->nndmax = 2.0*nnscale*maxshell;
	
	//Check for invalid solvent parameters
	if((configuration->vdwtype==2) && (fabs(configuration->Solvent[0])<EPS)){
		printf("Solvent[0] cannot be 0 if implicit solvent is used.\n");
		exit(1);
	} else{
		if((configuration->vdwtype==1) || (configuration->vdwtype==4) || (configuration->vdwtype==6)){
			double ljratio = (double)(configuration->LJexp[1])/(double)(configuration->LJexp[0]);
			if(fabs(ljratio-2)>EPS){
				printf("Simple LJ potential requires repulsive exponent to be twice attractive component (e.g. 6 and 12)\n");
				exit(1);
			}
		}
	}
}

/// Configure electrostatics (including reaction field)
void setup_electrostatics(Config_Data* configuration)
{
	//Set up electrostatic cutoff
	double bltemp;
	if (configuration->latticetype < 3)
	{
		bltemp = configuration->boxlength[0];
		if (configuration->boxlength[1] < bltemp) bltemp = configuration->boxlength[1];
		if (configuration->boxlength[2] < bltemp) bltemp = configuration->boxlength[2];
		bltemp /= 2.0;
	}
	else if (configuration->latticetype > 2)
	{
		bltemp = configuration->spherecylr;
		if (configuration->latticetype == 4)
		{
			if (configuration->boxlength[2] < bltemp/2.0) bltemp = configuration->boxlength[2]/2.0;
		}
	}
	if (configuration->escut > bltemp)
	{
		configuration->escut = bltemp;
		configuration->escut2 = bltemp*bltemp;
		cout << "WARNING: Electrostatic cutoff is too large. Resizing to half of the box length (" << bltemp << ")\n";
	}
	//Check geometry of simulation cell
	if((fabs(configuration->epsilon-configuration->n2)>EPS) || configuration->dyneps){
		if(configuration->latticetype==3){
			configuration->rfcut = configuration->spherecylr;
		} else{
			if(configuration->cut){
				configuration->rfcut = configuration->escut;
			} else{
				cout << "Reaction field requires a spherical cutoff or simulation cell!\n";
				exit(1);
			}
		}
	} else configuration->rfcut = 1.0; //this can be any nonzero number; it just gets multiplied by 0 - LEJ 12/22
	
	if(!configuration->dyneps && configuration->selfrf){
		cout << "WARNING: Self reaction field is on, but reaction field off. Switching self reaction field off.\n";
		configuration->selfrf=false;
	}
	if(configuration->selfrf) cout << "Self reaction field is on.\n";
	
	//Check polarizable medium dielectric constant n2
	if(fabs(configuration->n2)>EPS){
		configuration->in2 = 1.00/configuration->n2;
	} else{
		cout << "Cannot define cavity dielectric as zero!\n";
		exit(1);
	}
	// since 1/n^2 (in2 is only used for that, not for reaction field n^2) is used to scale electrostatic interaction, energy scaling can be applied here too
	configuration->in2 *= configuration->energy_scale[1];
	configuration->mucorr = 1.00;
	if(fabs(configuration->n2-1.0)>EPS){
		printf("Polarizable medium active.\n");
		if(configuration->polarize){
			cout << "Onsager polarizability active.\n";
			configuration->mucorr = (configuration->n2 + 2.0)/3.0;
			configuration->in2 *= qqpwr(configuration->mucorr,2);
		}
	}
	
	//Check poling field if present
	if(configuration->Efield.V3Norm()>EPS){
		if (configuration->realfield){
			printf("Poling field defined at constant voltage.\n");
		} else printf("Poling field defined as field inside cavity.\n");
	}
}

/// Check that the box is valid and calculate its volume
void check_box(Config_Data* configuration)
{
	//Check to make sure PBCs are being used properly
	if (configuration->latticetype == 3 && (configuration->PBCs[0] == 1 || configuration->PBCs[1] == 1 || configuration->PBCs[2] == 1)){
		printf("PBCs are not yet supported for spherical simulation cells.\n");
		exit(1);
	} else if (configuration->latticetype == 4 && (configuration->PBCs[0] == 1 || configuration->PBCs[1] == 1)){
			printf("PBCs are currently only supported on the z-axis for cylindrical simulation cells.\n");
			exit(1);
		} else if (configuration->PBCs[0] == 0 && configuration->PBCs[1] == 0 && configuration->PBCs[2] == 0 && configuration->latticetype < 3){
			printf("Warning: PBCs are off! Simulating a finite cube of liquid.\n");
		}
}

/// Check to make sure number of particles, simulation length, etc, are valid
void final_validation(Config_Data* configuration)
{
	//Make sure the simulation is worth running/runnable
	bool rot_trans_test=false;
	bool average_test=false;
	for(unsigned int i=0; i<configuration->num_element_types; i++){
		if(configuration->element_types[i]->still==0) rot_trans_test=1;
		if(configuration->element_types[i]->calculate_order) average_test=true;
		if ((configuration->element_types[i]->saxes.vec[0] < 0.0) || (configuration->element_types[i]->saxes.vec[1] < 0.0) || (configuration->element_types[i]->saxes.vec[2] < 0.0)){
			printf("Invalid dimension for %s-type oids. All dimensions must be positive.\n",configuration->element_types[i]->name.c_str());
			exit(1);
		}
		if(configuration->vdwtype<4){
			if(((fabs(configuration->element_types[i]->saxes.vec[0]-configuration->element_types[i]->saxes.vec[1])>EPS) || (fabs(configuration->element_types[i]->saxes.vec[0]-configuration->element_types[i]->saxes.vec[2])>EPS) || (fabs(configuration->element_types[i]->saxes.vec[1]-configuration->element_types[i]->saxes.vec[2])>EPS)) && (configuration->element_types[i]->number>0)){
				printf("%s-type oids are non-spherical and require Touch (vdw type >= 4)\n",configuration->element_types[i]->name.c_str());
				exit(1);
			}
		}
	}
	if(configuration->n_oids == 0){
		printf("The simulation doesn't contain any particles!\n");
		exit(1);
	}
	else if (rot_trans_test == 0)
	{
		printf("Running a simulation where no particles can rotate or translate is rather pointless.\n");
		exit(1);
	}
	else if (average_test == 0)
	{
		printf("Rather odd, but OK: No elements specified for order calculation.\n");
	}
	else if (configuration->anneal == 1 && ((configuration->steps - configuration->laststep) < configuration->randsteps*2))
	{
		printf("Temperature is not stable during averaging period.\n");
	}
	// Now is a good time to set the maximum cutoff
	adjust_cut2plus(configuration);
}

/// A couple of basic checks for the molecules in the box — can probably be merged with other functions.
void check_oids(MC_Elements &BeOids)
{
	//Determine the number of particles actually averaged over
	int to_average = 0;
	
	for(unsigned int i = 0; i < BeOids.configuration->n_oids; i++){
		if(BeOids[i]->MyType->calculate_order) to_average++;
	}
	
	BeOids.configuration->n_avg = to_average;
}

inline string Int2CenterStr(int i, unsigned int width)
{
	string result="";
	string s=int2str(i);
	unsigned int slen=s.length();
	for(unsigned int k=0; k<(width-slen)/2; k++){
		result+=" ";
		s+=" ";
	}
	
	if(result.length()+s.length()!=width) result+=" "+s; else result+=s;
	return result;
}

/*!
 * Summary file generator
 * Generate machine-readable summary file. This no longer requires a batch job. LEJ 02/17/2010
 */
void open_summary(Config_Data* configuration)
{
	time_t rawtime;
	struct tm * start_time;
	time (&rawtime);
	start_time = localtime (&rawtime);
	
	//Calculate corrected dipole moments, density, and reduced dipole density y
	double rho = (double)(configuration->N)/configuration->V;
	
	//Generate summary file header with molecule and simulation info needed for post-processing
	string summfn = configuration->fileout + "_summary_" + int2str(configuration->runid) + ".dat";
	FILE* summ = fopen(summfn.c_str(), "wa");
	fprintf(summ, "Summary file opened at %s \n", asctime(start_time));
	fprintf(summ, "#SIMULATION SUMMARY\n\n");
	double muz_total=0;
	fprintf(summ, "#Present element type properties:\n#name\tid #\tsemiaxes [Ang]\tLJ epsilon [perg]\tN\tcharges 4xn-tupel\tdipole vector\n");
	for(unsigned int i=0; i<configuration->num_element_types; i++){
		if(configuration->element_types[i]->number>0){
			string charges="{";
			if(configuration->element_types[i]->nr_charges>0){
				for(unsigned int q=0; q<configuration->element_types[i]->nr_charges; q++){
					if(q>0) charges+="|";
					charges+=double2str(configuration->element_types[i]->q[q]/e_in_esu)+"; "+configuration->element_types[i]->q_pos[q].V3Str(',');
				}
				charges+="}";
			} else charges="-";
			fprintf(summ, "%s	%i	%i	(%s)	%f	%s	(%s)\n", configuration->element_types[i]->name.c_str(), i+1, configuration->element_types[i]->number, (configuration->element_types[i]->saxes.V3Str(',')).c_str(), configuration->element_types[i]->Vvdw, charges.c_str(), (configuration->element_types[i]->initial_dipole.V3Str(',')).c_str());
			muz_total+=(double)configuration->element_types[i]->number*configuration->element_types[i]->dipole*configuration->element_types[i]->dipole;
		}
	}
	fprintf(summ, "\n#Present group properties\n#name\tid #\tsemiaxes [Ang]\tLJ epsilon [perg]\tN\tcharges 4xn-tupel\tdipole vector\n");
	for(unsigned int i=0; i<configuration->num_groups; i++){
		Element_Group* group=configuration->groups[i];
		if(group->number>0){
			fprintf(summ, "#Group: %s\n",group->Type->name.c_str());
			for(unsigned int j=0; j<group->nr_elements; j++){
				Element* element=&configuration->group_elements[group->elements[j]];
				bool not_done=true;
				for(unsigned int d=0; d<j; d++){
					if(element->mytype==configuration->group_elements[group->elements[d]].mytype){
						not_done=false;
						break;
					}
				}
				if(not_done){
					string charges="{";
					if(element->MyType->nr_charges>0){
						for(unsigned int q=0; q<element->MyType->nr_charges; q++){
							if(q>0) charges+="|";
							charges+=double2str(element->MyType->q[q]/e_in_esu)+"; "+element->MyType->q_pos[q].V3Str(',');
						}
						charges+="}";
					} else charges="-";
					fprintf(summ, "%s	%i	%i	(%s)	%f	%s	(%s)\n", element->MyType->name.c_str(), element->mytype+1, group->number, (element->MyType->saxes.V3Str(',')).c_str(), element->MyType->Vvdw, charges.c_str(), (element->MyType->initial_dipole.V3Str(',')).c_str());
					muz_total+=(double)element->MyType->number*element->MyType->dipole*element->MyType->dipole;
				}
			}
		}
	}
	fprintf(summ, "\n#Present element upper LJ epsilon interaction matrix in perg\n");
	string nrs="";
	bool first=true;
	unsigned int ii=0;
	unsigned int jj=0;
	for(unsigned int i=0; i<configuration->num_element_types; i++){
		if(configuration->element_types[i]->number>0){
			unsigned int i_nr=0;
			for(unsigned int j=0; j<configuration->num_element_types; j++){
				unsigned int ikk=i*configuration->num_element_types+j;
				if(configuration->element_types[j]->number>0){
					if(i_nr==0){
						fprintf(summ, "%i", i+1);
						ii++;
						jj=1;
					}
					if(first) nrs+="\t"+Int2CenterStr(j+1,8);
					if(jj<ii) fprintf(summ, "\t    -   "); else fprintf(summ, "\t%f", configuration->pre_eps[ikk]);
					jj++;
					i_nr++;
				}
			}
			for(unsigned int g=0; g<configuration->num_groups; g++){
				Element_Group* group=configuration->groups[g];
				if(group->number>0){
					for(unsigned int j=0; j<group->nr_elements; j++){
						Element* element=&configuration->group_elements[group->elements[j]];
						bool not_done=true;
						for(unsigned int d=0; d<j; d++){
							if(element->mytype==configuration->group_elements[group->elements[d]].mytype){
								not_done=false;
								break;
							}
						}
						if(not_done){
							unsigned int ikk=i*configuration->num_element_types+element->mytype;
							if(i_nr==0){
								fprintf(summ, "%i", i+1);
								ii++;
								jj=1;
							}
							if(first) nrs+="\t"+Int2CenterStr(element->mytype+1,8);
							if(jj<ii) fprintf(summ, "\t    -   "); else fprintf(summ, "\t%f", configuration->pre_eps[ikk]);
							jj++;
							i_nr++;
						}
					}
				}
			}
			if(i_nr>0) fprintf(summ, "\n");
			first=false;
		}
	}
	for(unsigned int g=0; g<configuration->num_groups; g++){
		Element_Group* groupi=configuration->groups[g];
		if(groupi->number>0){
			for(unsigned int gi=0; gi<groupi->nr_elements; gi++){
				Element* elementi=&configuration->group_elements[groupi->elements[gi]];
				bool not_done=true;
				for(unsigned int d=0; d<gi; d++){
					if(elementi->mytype==configuration->group_elements[groupi->elements[d]].mytype){
						not_done=false;
						break;
					}
				}
				if(not_done){
					unsigned int i=elementi->mytype;
					unsigned int i_nr=0;
					for(unsigned int j=0; j<configuration->num_element_types; j++){
						unsigned int ikk=i*configuration->num_element_types+j;
						if(configuration->element_types[j]->number>0){
							if(i_nr==0){
								fprintf(summ, "%i", i+1);
								ii++;
								jj=1;
							}
							if(first) nrs+="\t"+Int2CenterStr(j+1,8);
							if(jj<ii) fprintf(summ, "\t    -   "); else fprintf(summ, "\t%f", configuration->pre_eps[ikk]);
							jj++;
							i_nr++;
						}
					}
					for(unsigned int g=0; g<configuration->num_groups; g++){
						Element_Group* group=configuration->groups[g];
						if(group->number>0){
							for(unsigned int j=0; j<group->nr_elements; j++){
								Element* element=&configuration->group_elements[group->elements[j]];
								not_done=true;
								for(unsigned int d=0; d<j; d++){
									if(element->mytype==configuration->group_elements[group->elements[d]].mytype){
										not_done=false;
										break;
									}
								}
								if(not_done){
									unsigned int ikk=i*configuration->num_element_types+element->mytype;
									if(i_nr==0){
										fprintf(summ, "%i", i+1);
										ii++;
										jj=1;
									}
									if(first) nrs+="\t"+Int2CenterStr(element->mytype+1,8);
									if(jj<ii) fprintf(summ, "\t    -   "); else fprintf(summ, "\t%f", configuration->pre_eps[ikk]);
									jj++;
									i_nr++;
								}
							}
						}
					}
					if(i_nr>0) fprintf(summ, "\n");
					first=false;
				}
			}
		}
	}
	if(nrs!="") nrs+="\n";
	nrs+="\n";
	fprintf(summ, "%s",nrs.c_str());
	double y = 4*pi*muz_total/(9.0*configuration->V*configuration->kT);
	fprintf(summ, "RHO	%f\n", rho);
	fprintf(summ, "T	%f\n", configuration->T);
	fprintf(summ, "N2	%f\n", configuration->n2);
	fprintf(summ, "Y	%f\n", y);
	fclose(summ);
}

