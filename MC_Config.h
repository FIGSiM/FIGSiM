/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * updated by LEJ, 01/15/2011
 *
 * updated by AT Jan 28, 2011
 *	- added ability to read values from files (format "value name = value")
 *	- added comments for autocreation of documentation through Doxygen
*/

#ifndef INCLUDED_MC_CONFIG
#define INCLUDED_MC_CONFIG

#ifdef OPENBABEL
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#endif

#ifndef INCLUDED_CONFIG
#include "Config.h"
#endif
#ifndef INCLUDED_CONFIGREADER
#include "ConfigReader.h"
#endif

const unsigned nr_bondorder_types=5;
const string bondorder_types[nr_bondorder_types]={"am","ar","du","un","nc"};

typedef struct _Element Element; // pulling a rabbit from the cylinder ;-)
typedef struct _Element_Type Element_Type;
typedef struct _Interaction Interaction;
typedef struct _Interaction_Potential Interaction_Potential;
typedef struct _Element_Group Element_Group;
typedef struct _Level_of_Detail Level_of_Detail;

typedef double Vstore[NUM_V_STORE];

typedef struct _Traj_EP{
	unsigned int element_type;
	int group_type;
	double speed;
	Vec3 position;
	Vec3 velocity;
	Vec4 rotation_vector;
} Traj_EP;

// struct and class definitions
typedef struct _Config_Data{
	string configfile;	///< configuration filename
	string configdir;	///< configuration file directory
	string trajfile;	///< trajectory filename
	
	bool fit2lod;		///< flag to indicate that configuration is being called through <fit2lod>
	bool use_trajectory;	///< there is a trajectory file to use
	__int32_t placement_rngseed;	///< random number generator seed used for element creation
	bool restart_calculations; ///< restart calculations from scratch at resume
	bool restart_volume;	///< restart volume shrinking process upon resume
	unsigned int last_step;	///< last step stored in trajectory file
	int resume_step;	///< what step to resume from (-1 means last)
	string trajectoryfile;	///< filename of trajectory file
	string trajectorydir;	///< trajectory file directory
	unsigned int trajectorynr; ///< subnumber in case trajectory file grows too big
	
	//RNG memory. Must not be constant! (RUNTIME)
#ifdef USE_CMWC4096
	__uint32_t rngseed;
#else
	__int32_t rngseed;
	__int32_t* idum;
#endif
	
	//Batch job data
	int runid;	///< ID number of current simulation
	string fileout; ///< Root of name for output file
	string UseMap; ///< Root of name for output file
	
	// Element type related constants
	unsigned int num_element_types;	///< number of element types
	Element_Type** element_types;	///< holds the element types
	
	// Group related constants
	unsigned int num_groups;	///< number of groups
	Element_Group** groups;	///< groups defined in the configuration
	Element* group_elements;	///< holds element-archetypes which are associated to the groups
	
	Vec3* group_centers;	///< holds the center of individual group archetype's bounding sphere
	double* group_radii;	///< holds the radius of individual group archetype's bounding sphere
	double* max_group_dipoles;	///< maximum dipole per defined group (if group_dipoles is true overall dipole for whole group, if not overall combined dipole from individual elements)
	
	// Interaction potential related constants
	unsigned int num_potentials;	///< number of interaction potentials
	Interaction_Potential** potentials;	///< array of pointers -- holds interaction potentials
	
	// Level-of-Detail related constants
	unsigned int num_levelofdetail;	///< number of level of detail sections
	Level_of_Detail** lods;
	
	double Mass;		///< precalculated overall system mass (if masses defined for elements)
	double max_dipole;	///< precalculated maximum dipole moment based on configuration
	double dipoleQ;	///< precalculated maximum Q for dipole moments of individual element types
	double quadrupole_charge_factor;	///< adjust quadrupol charges by this factor > 1 in order to space quadrupole charges closer
	double group_max_centerdist;	///< maximum distance from center of all groups in simulation box (used to determine neighborhood radius)
	double max_neighbor_trans;	///< maximum translation size in neighborhood
	double rT;	///< global test sphere size for LOD ellipsoid calculations (default: 1.65)
	
	double LJwall_a;
	double LJwall_b;
	double LJwall_width;
	double LJwall_xm;
	double LJwall_x;
	double LJwall_epsilon;
	double LJwall_mirror_factor;
	bool LJwall_fixed;
	bool LJwall_calc;
	
	//precalculated touch, rho, and eps parameters for LJ interactions of all molecule type
	//is created dynamically once we know how many there are
	double* pre_eps;	///< Lennard-Jones epsilons nxn array for all molecule types
	double* pre_sigma2;	///< Lennard-Jones sigma squared
	bool* pre_touch;	///< Control flags for whether or not to use touch
	
	//Boolean flags for configuring simulation type
	bool use_gpu;	///< Enables usage of GPU (in case OpenCL version is run, no effect on MiniCL thread fallback)
	bool NpT;	///< Enables NpT calculation mode
	bool adjust_internal;	///< adjust internal energies (VLJ and VES) until system reaches target volume
	bool user_rT;	///< if true indicates that the user set the rT value
	bool test;	///< Enables test mode if true
	bool cut;	///< Enables potential truncation if true
	bool cut_internal;	///< Enables potential truncation (VLJ and VES) for group internal calculations if true
	bool ions_present;	///< internal variable set true when ions are present in simulation
	bool PBCs[3];	///< Enables x,y,z periodic boundary conditions
	bool dogr;	///< Enables calculation of correlation functions if true
	bool dyneps;	///< Enables SCRF if true
	bool selfrf;	///< Enables reaction field contribution from originating element/group (if group_dipole is set)
	bool chkcoord;	///< Dumps intermediate coordinate files if true
	bool touchtrunc;	///< Uses an additional, smaller truncation radius for steric anisotropy if true 
	bool doxm;	///< Calculates extended order parameters (wrt x,y,z) if true
	bool offctrmu;	///< Molecular dipoles are allowed to be off-center if true (SLOW!)
	bool realfield;	///< Poling field defined at constant potential (Maxwell field) if true
	bool anneal;	///< Simulated annealing used if true
	bool polarize;	///< Dipole moments enhanced by Lorentz cavity factor if true
	bool auto_volume;	///< if true automatically extend volume if running out of space
	bool smooth_animation;	///< wether animation interpolates smoothly (true) or not (false) -- for large rotational changes smooth animation can lead to wierd behavior (oids moving unphysically)
	bool time;	///< define time using boltzmann-maxwell distribution and use in MC algorithm
	bool randsteps_nointernal;	///< if true, groups are not randomized internally (aka stiff groups) during randomization
	bool randsteps_novolume;	///< if true, no volume adjustments from initial placement volume are undertaken during randomization (neither in NpT nor constant V - for constant V, the volume drop happens after randomization in half the randomization steps)
	bool dyneps_varm;	///< if true, var(M) instead of M^2 is used in reaction field epsilon calculations when no poling field is defined (default: false)
	bool transition;	///< if true, transition from large volume, no LJ dispersion to target volume, full LJ dispersion mimicking solvent evaporation
	bool NpT_lnV_change;	///< if true, NpT volume move is V' = exp( lnV +/- dV_max*[-1,1] )
	bool NpT_adjust_EScut;	///< if true, adjust electrostatics cutoff distance with NpT volume step
	bool NpT_relative_coords;	///< if true (default), adjust element coordinates to be at the same location relative to the box (scale coordinates with box), otherwise apply boundary conditions
	bool LJ_interaction_area;	///< if true (default is false to stay consistent with old code), adjust ellipsoid Lennard-Jones epsilon based on cross-sectional area of interellipsoid distance vector normal plane with ellipsoids
	bool LJ_interaction_area_fit;	///< if true (default is false), adjust IA values to (on average) match the textured values
	bool LJ_adjust_width;		///< if true (default is false) average width for each fit ellipsoid is used
	bool LJ_epsilon_sixthpower;	///< if true, use sixth power LJ epsilon mixing rule (e_ij = (2*sqrt(e_i*e_j)*sigma_i^3*sigma_j^3)/(sigma_i^6+sigma_j^6), otherwise use geometric rule e_ij = sqrt(e_i*e_j) (default: 0)
	bool collective_moves;		///< if true (default is false) use collective moves for bound things
	bool VmuE_squared;		///< if true, calculate charge/dipole-Efield interaction as -(mu*E)^2/(|mu||E|) (default: false, VmuE = -mu*E)
	
	//Integer parameters - none of these should be adjusted during run
	unsigned int vdwtype;	///< Sets LJ interaction type
	unsigned int loader;	///< Sets loading algorithm used for building the box (default: 2 (random))
	int LJexp[2];	///< Sets Lennard-Jones exponents
	double LJlambda;	///< scaling factor for dispersive (attractive) Lennard-Jones term (default: 1.0)
	unsigned int latticetype;	///< Sets lattice type for initial box/sphere (default: 1 (simple cubic))
	
	unsigned int n_oids;	///< Total number of elements in system
	unsigned int n_groups;	///< Total number of groups in system
	unsigned int n_group_oids;	///< Total number of elements in all groups
	unsigned int N;		///< N is number of individual elements and groups (use for things like NkT ...)
	unsigned int N_NpT;	///< N for NpT move (is either N or N+1)
	unsigned int max_neighbors;	///< maximum number of neighbors possible
	unsigned int NpT_move_nr;	///< attemp NpT move every NpT_move_nr moves (defaults to n_oids)
	unsigned int add_slowly;	///< if >0 add individual elements or groups one by one every add_slowly cycles starting from the center (default: 0, add everything at beginning)
	unsigned int LJ_radius_mixing;	///< Mixing rule for LJ radius: 0 .. arithmetic (sigma_ij = r_i+r_j), 1 .. geometric (sigma_ij = 2*sqrt(r_i*r_j)), 2 .. sixth power (sigma_ij = 2*((r_i^6+r_j^6)/2)^(1/6) (default: 0)
	double Ekin;		///< average E_kin = dof/2*N*kT (is in pico-ergs=10^(-19) J)
	
	// Parameters relevant for group/element placement
	double rmax_g;	///< maximum bounding sphere radius of all groups (as defined in configuration file)
	double a_group;	///< fcc lattice constant (size in Angström)
	unsigned int nr_fcc_cells;	///< how many fcc unit cells in total in system volume
	unsigned int fcc[3];	///< number of fcc cells fitting in user defined box in each dimension
	
	double c_group[3];	///< radius of space in between fcc lattice elements
	
	double rmax_e;	///< maximum bounding sphere radius of all element types
	double a_element;	///< sc lattice constant (size in Angström)
	unsigned int nr_sc_cells;	///< how many sc unit cells in total in system volume (if only elements and no groups are specified)
	unsigned int sc[3];	///< number of simple cubic cells in each direction inside empty space in fcc unit cell
	unsigned int elements_per_fcc;	///< elements fitting inside empty space per fcc unit cell
	
	double scale[3];	///< scaling factors of unit cells to fit in overall system volume
	double energy_scale[4];	///< scaling factor for individual potential energies (0..VmuE, 1..VES, 2..VLJ, 3..VG)
	
	unsigned int n_avg;	///< Total number of molecules being averaged (SETUP)
	unsigned int xmfreq;	///< Frequency of extended means calculations
	unsigned int grinc;	///< Number of grid elements used for discretizing correlation functions
	unsigned int grfreq;	///< Frequency for correlation function calculations
	unsigned int dyneps_average_cycles;	///< Number of cycles over which to average M, M^2, or calculate var(M) over for reaction field epsilon calculations (see also dyneps_varm), defaults to zero
	unsigned int correlation_length;	///< correlation length to use block average (should be obtained from proper correlation function)
	unsigned int running_average;	///< Number of blocks to use +/- for running averages
	unsigned int steps;	///< Total steps in a simulation
	unsigned int randsteps;	///< Number of steps used for randomization at beginning of simulations
	unsigned int laststep;	///< Number of steps used for averaging at end of simulation
	unsigned int stepsize_average;	///< Number of steps over which stepsizes are averaged (running average) to determine optimum stepsizes per element/group - No average if set to 0 (default)
	unsigned int transition_start;	///< Number of steps after which transition (currently: volume down to targetV, LJ lambda up to 1.0) starts
	unsigned int transition_cycles;	///< Width (in cycles) of transition
	double group_rand_frac;	///< fraction of the total number of elements in a group being randomized after rot/trans of whole group
	
	//Parameters that may be adjusted during run (most currently aren't adjustable, though)
	double epsilon;	///< Initial reaction field dielectric constant (dynamic value is NOT stored here)
	double boxlength[3];	///< Length of the x,y,z dimensions of the box
	double inv_boxlength[3];	///< Inverse length of the x,y,z dimensions of the box
	double spherecylr;	///< Radius if a spherical or cylinderical box is used
	double V;	///< adjusted volume to comfortably fit elements
	double pext;	///< external pressure (use in NpT mode)
	double targetV;	///< specified volume of box
	double transition_start_volume;	///< start volume for transition
	double transition_target_volume;	///< target volume for transition
	double transition_delta_n2;	///< adjusts 1/n^2 for electrostatics calculations from 1/(n^2+dn^2) to 1/n^2
	double transition_LJ_start;	///< LJ attractive fraction starting point (default: 0)
	double transition_LJ_end;	///< LJ attractive fraction end point (default: 1)
	double transition_LJ_delta;	///< LJ attractive fraction difference (calculated from transition_LJ_start and transition_LJ_end)
	double rmax;	///< Cutoff distance for correlation functions
	double refdr;	///< Size of each grid element for correlation functions
	double maxtrans;	///< Maximum distance for a translation move
	double maxrot;	///< Maximum angle for a rotational move
	double maxdV;	///< Maximum change in volume for NpT move
	double n2;	///< Square of refractive index, used for implicit polarizability
	double in2;	///< Inverse of square of refractive index (SETUP)
	
	double mucorr;	///< Lorentz field fractor for scaling dipole moments (SETUP)
	double Rand_Time; ///< Time the randomization step took (MC_Elements::evolve_system)
	double Runtime;
	bool storetime;
	bool output_code;	///< Store simulation source code alongside trajectory file (default: true)
	
	//dynamic E and T adjustment
	unsigned int trackchange;	///< Counter for temperature changes, currently always 0
	unsigned int changesteps;	///< Steps at which T is changed, not currently used
	Vec3 Efield;	///< Initial poling field strength. How this is defined depends on cnsts->realfield
	bool noEfield;	///< true if no field specified
	bool rand_Efield;	///< if true, external field interaction is calculated during randomization
	unsigned int rotate_Efield_steps;	///< number of steps after which to rotate system so total dipole moment direction matches electric field (do nothing if zero)
	double T;	///< Initial absolute temperature
	double kT;	///< self-explanatory
	double beta; ///< as in 1/kT
	
	//Parameters that should not be adjusted during run
	double Solvent[3];	///< Width, amplitude, and rmind for implicit solvent correction
	double stilleuler[3];	///< Euler angles for non-rotating molecules
	double stillrot[3][3];	///< Rotation matrix for non-rotating molecules
	double nndist;	///< Nearest-neighbor distance for molecules (density^(-1/3))
	Vec4 opt_cutoffs;	///< Optimal cutoff distances (when maximum potential difference is on order of EPS) for (q-q, q-mu, mu-mu, LJ) interactions
	double ljcut;	///< Lennard-Jones cutoff distance
	double escut;	///< Electrostatics cutoff distance
	double ljcut2;	///< Lennard-Jones cutoff distance squared
	double escut2;	///< Electrostatics cutoff distance squared
	double cut2plus;	///< Largest cutoff radius squared plus maximum displacement squared (for neighbor determination)
	double ES_cutoff_precision;	///< Another way of determining cutoff distances: energy difference cutoff for likely move (maxtrans/2)
	double LJ_cutoff_precision;	///< Another way of determining cutoff distances: energy difference cutoff for likely move (maxtrans/2)
	double rfcut;	///< Reaction field cavity radius (generally equal to cnsts->escut) (SETUP)
	double* tcut2;	///< Touch cutoff square array (SETUP)
	double nndmax;	///< Approximate radius of nearest-neighbor shell	(SETUP)
	
	double r;
	//LDBR parameters. Can these be constant? (not currently used)
	double fshift;
	double Vvdwboundary;
	double uBoundary;
	double rboundary;
	double uSigma;
} Config_Data;

/*!
This class is used to store the system configuration as well as any value needed by multiple
classes or functions. All variables are assigned values through GetFromFile which reads and parses
a textfile containing the setup of the system.
*/
class MC_Config : public ConfigReader
{
public:
	MC_Config();
	virtual ~MC_Config();
	
	void GetFromFile(const char* filename);
	bool GetAllElementProperties(Traj_EP** elements, unsigned int* steps, double* V, double* time, unsigned int &nr_steps){ return GetAllElementProperties(elements,steps,V,time,nr_steps,0); };
	bool GetAllElementProperties(Traj_EP** elements, unsigned int* steps, double* V, double* time, unsigned int &nr_steps, unsigned int start_step){ return GetAllElementProperties(elements,steps,V,NULL,time,nr_steps,start_step); };
	bool GetAllElementProperties(Traj_EP** elements, unsigned int* steps, double* V, double* Xm, double* time, unsigned int &nr_steps, unsigned int start_step);
	bool GetElementProperties(const unsigned step, Traj_EP* elements);
	bool GetStatistics(const unsigned int uptostep, Vstore* potentials, Vec3* dipoles, Vec3* cosmeans, double* Vs, double* msmoved, unsigned int* accepted, unsigned int* tries){ return GetStatistics(0,uptostep,potentials,dipoles,cosmeans,Vs,msmoved,accepted,tries); };
	bool GetStatistics(const unsigned int initialstep, const unsigned int uptostep, Vstore* potentials, Vec3* dipoles, Vec3* cosmeans, double* Vs, double* msmoved, unsigned int* accepted, unsigned int* tries);
	bool VLJatomisticZero(Vec3 &center, unsigned int* combine_elements, unsigned int count, Vec3 direction, Vec3 &r, double rT);
	bool dVLJatomisticZero(Vec3 &center, unsigned int* combine_elements, unsigned int count, Vec3 direction, Vec3 &r, double rT);
	
	Config_Data parameters;
	unsigned int max_nr_components, max_levelofdetail, max_LOD_elements;
protected:
	int nr_group_elements; ///< running total of group elements created
	
	virtual char* IncludeFile(readfile* file, char* include_type, unsigned int &pos, unsigned int &include_position, string name, string &subname, bool combine_sections, string* inc_files);
	
	void LoadSimParams(char* conf);
	void LoadElementTypes(char* conf, int nr){ LoadElementTypes(conf,nr,true); }; // load charges by default
	void LoadElementTypes(char* conf, int nr, bool load_charges);
	bool LoadGroup(char* conf, int nr, bool rerun); ///< loads a group of elements, returns true if succesful (if not succesful, the group consists of other groups and needs to be created later)
	void LoadPotential(char* conf, int nr, string specifier);
	void CombineElements(Element_Group* group, Element_Group* newgroup, unsigned int* combine_elements, unsigned int count, string specifier, double origV, double origEps, double rT);
	void SetEllipsoidLinks(Element_Group* group, Element_Group* newgroup, int* combine_elements, unsigned int* counters);
	void SetEllipsoidCharges(Element_Group* group, Element_Group* newgroup);
	void LoadLevelofDetail(char* conf, Element_Group* group, int nr);
	
	void SetSystemProperties();
	void SetSystemVolume();
	void AutoPotentials(Element_Group* group);
	void AssignFixedElements(Element_Group* group);
	
	char* Mol2Convert(char* content, char* conf, string groupname);
	char* OpenBabelConvert(char* content, char* conf, string groupname, string extension);
	char* PDBConvert(char* content, char* conf, string groupname);
	void GetSpecialDistances(double* special, int nr);
	
	double GroupLJdisp(Element_Group* group, unsigned int* deleted, unsigned int nr_deleted, bool* connection_site, unsigned int start, unsigned int divider, double theta, Vec3 axis, Vec3 location);
	Element_Type* NewType(Element_Type* source);
	void CopyElement(Element* source, Element* dest, int offset){ CopyElement(source,dest,offset,true); }
	void CopyElement(Element* source, Element* dest, int offset, bool copyproperties);
	void AddGroup2Group(Element_Group* group, Element_Group* add_group, unsigned int basegroup_start);
};


inline void update_volume(Config_Data* configuration, double newV)
{
	if(fabs(newV-configuration->V)>EPS){ // adjust parameters only when needed
		double change=cbrt(newV/configuration->V); // how much each dimension needs to change
		switch(configuration->latticetype){
			case 1:
			case 2: // Rectangular box
				configuration->boxlength[0]*=change;
				configuration->inv_boxlength[0]=1.0/configuration->boxlength[0];
				configuration->boxlength[1]*=change;
				configuration->inv_boxlength[1]=1.0/configuration->boxlength[1];
				configuration->boxlength[2]*=change;
				configuration->inv_boxlength[2]=1.0/configuration->boxlength[2];
				break;
			case 3: // Spherical volume
				configuration->spherecylr*=change;
				break;
			case 4: // Cylindrical volume
				configuration->spherecylr*=change;
				configuration->boxlength[2]*=change;
				configuration->inv_boxlength[2]=1.0/configuration->boxlength[2];
				break;
			default: // Default condition
				cout << "Cannot change volume. Lattice type " << configuration->latticetype << " not recognized. Exiting.\n";
				exit(1);
		}
		configuration->nndist*=change;
		configuration->V=newV;
	}
}

inline void delete_from_array(unsigned int** array, unsigned int &nr_elements, unsigned int which)
{
	if((which<nr_elements) && (nr_elements>0)){
		// shift element array one down (onto element number to be delete) and resize array
		nr_elements--;
		for(unsigned int i=which; i<nr_elements; i++) (*array)[i]=(*array)[i+1];
		*array=(unsigned int*)realloc(*array,nr_elements*sizeof(unsigned int));
		if(!(*array)){ // should not happen (shrinking should always work) but we'd be in big trouble if it did, so ...
			cout << "ERROR: The universe just folded onto itself while losing matter/antimatter.\n";
			exit(42);
		}
	} else{
		cout << "ERROR: Trying to delete entry " << which << " from an array with " << nr_elements << " entries.\n";
		exit(4);
	}
}

inline void delete_from_array(int** array, unsigned int &nr_elements, unsigned int which)
{
	if((which<nr_elements) && (nr_elements>0)){
		// shift element array one down (onto element number to be delete) and resize array
		nr_elements--;
		for(unsigned int i=which; i<nr_elements; i++) (*array)[i]=(*array)[i+1];
		*array=(int*)realloc(*array,nr_elements*sizeof(int));
		if(!(*array)){ // should not happen (shrinking should always work) but we'd be in big trouble if it did, so ...
			cout << "ERROR: The universe just folded onto itself while losing matter/antimatter.\n";
			exit(42);
		}
	} else{
		cout << "ERROR: Trying to delete entry " << which << " from an array with " << nr_elements << " entries.\n";
		exit(4);
	}
}

inline void delete_from_array(double** array, unsigned int &nr_elements, unsigned int which)
{
	if((which<nr_elements) && (nr_elements>0)){
		// shift element array one down (onto element number to be delete) and resize array
		nr_elements--;
		for(unsigned int i=which; i<nr_elements; i++) (*array)[i]=(*array)[i+1];
		*array=(double*)realloc(*array,nr_elements*sizeof(double));
		if(!(*array)){ // should not happen (shrinking should always work) but we'd be in big trouble if it did, so ...
			cout << "ERROR: The universe just folded onto itself while losing matter/antimatter.\n";
			exit(42);
		}
	} else{
		cout << "ERROR: Trying to delete entry " << which << " from an array with " << nr_elements << " entries.\n";
		exit(4);
	}
}

MC_Config* GetConfig(); ///< returns pointer to global_Config* object defined for program-wide use

const unsigned short arbitrary_element_names_nr = 16;
const string arbitrary_element_names[] = {"arthur dent", "ford prefect", "trillian", "zaphod beeblebrox", "marvin", "slartibartfast", "deep thought", "42", "zarniwoop", "fenchurch", "mark ii", "random frequent flyer dent", "tricia mcmillian", "agrajag", "emily saunders", "colin"};

const unsigned short arbitrary_group_names_nr = 10;
const string arbitrary_group_names[] = {"Pink Floyd", "Led Zeppelin", "Queen", "The Beatles", "The Rolling Stones", "The Who", "Radiohead", "Air", "Morcheeba", "The Doors"};

#endif

