/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * updated by AT March 30, 2011
 */

#ifndef INCLUDED_MC_ELEMENTS
#define INCLUDED_MC_ELEMENTS

#include <signal.h>
#include <omp.h>

#ifndef INCLUDED_CONFIG
#include "Config.h"
#endif

#ifndef INCLUDED_MC_CONFIG
#include "MC_Config.h"
#endif

using namespace std;

typedef struct _Element Element; ///< magic happens here: forward declaration of Element structure (a.k.a. pulling a rabbit out of the hat)

typedef struct _special_distance{
	unsigned int distance; // 4 bytes
	float factor; // float is good enough for prefactor and 4-byte length allows for better array accesibility
} special_distance;

typedef struct _Element_Group_Type{
	string name; ///< group name
	Vec3 group_dipole_color; ///< color of dipole(s) of this particular group type
	
	unsigned int dof; ///< degrees of freedom particular group (defaults to non-linear 6 - vibrations are taken care of by potential energy)
	unsigned int rand_elements; ///< how many elements are randomized during each internal configuration randomization
	unsigned int nr_movable; ///< total number of elements allowed to be moved inside group
	double nr_rand_per_cycle; ///< how many times per MC cycle on average the group's elements are randomized (defaults to sqrt(# of movable group elements))
	double inv_sqrt_nrpc; ///< 1.0/nr_rand_per_cycle
	
	bool group_dipole; ///< if true add charges and dipole moments of all group elements in order to get group dipole moment (default: false)
	bool rand_independent; ///< if true group elements are randomized independent of each other during group randomization (default: false)
	bool allow_bond_bend; ///< allow bond bending even if no bend potential is assigned
	bool still; ///< Whether or not the group can move
	bool rot_notrans; ///< Whether or not the group can rotate but not translate
	bool add_start_LOD; ///< parameter for super groups which controls the import of groups when the base name is used as either a group's start LOD level (if greater than 0) or as the base level, default is 1
	bool show_just_dipoles; ///< display only dipole moments in x3d animation (either uses component dipoles which have "component_show_dipole.<nr> = 1", or falls back to use whole group's dipole)
	bool is_virtual; ///< indicates if group is a "virtual" (a collection of groups and elements where no internal randomization shall take part)
	
	double mass;
	double density;
	double packing_density;
	double Volume;
	double charge;
	double mincharge;
	double maxcharge;
	double ES_scale;
	
	bool calculate_order;
	
	bool visual_PBCs;
	bool label_elements;
	double transparency;
	
	double* overlaps;
	double adjust_overlap; // constant to adjust half overlap used
	
	bool* elements_in_ring;
	unsigned int nr_rings;
	unsigned int** rings;
	
	bool* fixed_elements; ///< if true, particular element (relative to group start) is in a bond fixed with another element
	
	unsigned int nr_connection_sites;
	unsigned int* connection_sites;
	
	unsigned int bond_range; ///< number of bonds counted from a given element after which nonbonding interactions are switched on again
	special_distance range_factors[SPECIAL_DISTANCES]; ///< array containing Lennard-Jones and Electrostatic potential prefactors for special distances
	__uint32_t* range_field; ///< 4-bit bitfield holding information if elements are within bond_range
	
	Level_of_Detail* LOD;
} Element_Group_Type;

/// This struct is a container representing a group of elements
typedef struct _Element_Group{
	unsigned int number; ///< group number (to be able to distinguish from each other) -- doubles as number to be created in configuration
	unsigned int type; ///< group type is number of group in configuration to access resources like center, etc. ...
	unsigned int levelofdetail;
	
	double VES;
	double VLJ;
	double VG;
	double inv_nr_elements; ///< inverse number of elements in group
	double rT; ///< test sphere radius
	
	Vec3 center; ///< group center (geometric) -- is updated each time group is randomized
	
	unsigned int nr_neighbors;
	unsigned int old_nr_neighbors;
	unsigned int nr_neighbors_above;
	bool fixed;	///< group is there but can not move or rotate
	unsigned int nr_elements; ///< number of elements in group
	unsigned int* elements; ///< actual numbers (in class array) of elements in group
	
	unsigned int nr_potentials; ///< number of interaction potentials in this group
	Interaction_Potential** potentials; ///< array of pointers to all potentials in this group
	
	Element_Group_Type* Type;
} Element_Group;

typedef struct _Level_of_Detail{
	unsigned int levels; ///< nr of levels detail (0 downto this number of levels)
	unsigned int start_level; ///< simulation uses this level to begin with (default: 0 - group as specified)
	unsigned int end_level; ///< last level of detail simulation uses (default: number of levels specified)
	bool gyration_sphere; ///< if enabled, gyration radius is used to construct spherical object rather than ellipsoidal one (default: false)
	bool symmetry_center; ///< if enabled, use most symmetric center with respect to fully atomistic potential (default: true)
	bool symmetry_axes; ///< if enabled, calculate ellipsoid axes as average 90 degree cones with respect to fully atomistic potential boundaries NOTE: If both symmetry_axes and match_volume are false, use gyration tensor results as is (default: true)
	bool match_volume; ///< if enabled, Lennard-Jones volume originally excluded will be preserved (default: true)
	bool match_epsilon; ///< if enabled, Lennard-Jones epsilon is calculated from average epsilon around LOD ellipsoids based on fully atomistic model
	bool* keep_original_potentials; ///< controls if original potentials are kept for LOD models
	bool* use_epsilon_texture;
	
	bool* visual_original; ///< if true, base model is visualized along with LOD model
	double* transparency; ///< transparency setting (0..1) for each LOD level
	double* inside_transparency; ///< transparency value (0..1) for fully-atomistic "inside" of LOD 3D plots (default: 0.0 - no transparency)
	double scale_original_vdw; ///< scale factor for base model ellipsoids
	bool show_internal_bonds; ///< wether to show bonds of the underlying fully-atomistic model
	bool internal_charge_colors; ///< if true, base model element colors are derived from charges
	
	double** reduce_electrostatics; ///< charges leftover after charge reduction (currently -1 (keep original charges), 0 (no monopole), 1 (one charge, one dipole), <-1 (shrink charge distances, preserve dipole))
	double* distances; ///< distance at which each level activates (not used currently)
	
	Element_Group** groups; ///< associated group per level
	
	unsigned int order_dipole_component; ///< which component's dipole to use for order calculations (default: 0 = whole group)
	unsigned int nr_components; ///< number of components defined
	unsigned int* nr_elements; ///< number of elements per component
	unsigned int** component_elements; ///< 2D array containing element numbers per component (first index component nr, second element number in component)
	int* element_in_component; ///< which component number is a particular base group element in (-1 means it is in no component)
	
	bool* component_show_dipole; ///< visualize dipole moments of components (not implemented yet)
	Vec3* component_color; ///< color of component in RGB
	double* component_transparency; ///< component transparencies
	string* component_texture; ///< textures for components
	bool* component_texture_is_dipole; ///< visualize dipole moments of components (not implemented yet)
	
	unsigned int* nr_ellipsoids; ///< total number of ellipsoids per LOD level
	unsigned int** ellipsoid_counts; ///< stores information for array "element_groups" below, first index is LOD level (LOD 1 = 0), second index is respective ellipsoid nr (value is second index of "element_groups" yielding how many base elements are in that particular ellipsoid)
	int** element_groups; ///< 2D array holding information which underlying model elements (start count at *1*) belong to each ellipsoid, first index is LOD level (LOD 1 = 0), second index is relative to count entry (is int flex tupel: count_1, nr 1, ..., nr count_1; count_2 , nr 1 , ... nr count_2 ; ...)
	int** element_in_ellipsoid; ///< maps underlying model element number to specific LOD model ellipsoid (counting starts at *0*)
	
	double** volumes; ///< volumes of each LOD ellipsoid (first index LOD level (LOD 1 = 0), second index ellipsoid number)
	double** epsilons; ///< Lennard-Jones epsilons for each LOD ellipsoid (first index LOD level (LOD 1 = 0), second index ellipsoid number)
	
	Vec3* placement_deltas; ///< placement difference to base model for each LOD level
} Level_of_Detail;

/// This struct holds information for interaction potentials
typedef struct _Interaction_Potential{
	string* name;	///< points to name of this interaction potential
	unsigned short type; ///< type of potential (e.g. stretch, bend, rotate, 3-body interaction, etc.)
	unsigned short n; ///< how many elements does this potential connect (depends on type)
	unsigned short nr_parameters; ///< how many parameters do we have
	
	unsigned int* partners; ///< numbers of elements (in group array) *or* element types involved (depends on particular potential and if potential is group-local or global)
	unsigned int* links; ///< nxn matrix containing link nr in partners (links[i*n+k] is link nr on partner i pointing to k)
	double* parameters; ///< points to array of potential parameters (size depends on particular potential)
} Interaction_Potential;

/// This struct holds information for an individual interaction (currently used for bonds)
typedef struct _Interaction{
	bool fixed; ///< if true then bond is fixed (no rotating, bending, or stretching allowed)
	bool allow_bond_stretch; ///< is true if bond stretch potential has been assigned
	bool allow_bond_bend; ///< is true if bond bend potential has been assigned
	double bond_length; ///< fixed bond length (if allow_bond_stretch is false)
	double bond_order; ///< if positive it's the bond_order, if negative it's the index for the bondorder_types defined in MC_Config.h
	
	unsigned int nr_potentials; ///< number of interaction potentials acting on this interaction
	Interaction_Potential** potentials; ///< array of pointers to potentials associated with this interaction (e.g. bond-stretch, bend, rotate, 3-body interaction, etc.)
	unsigned int partner; ///< what number is the element the interaction is with (in group)
	unsigned int back_link; ///< what interaction number on partner element belongs to the back link
	
	Vec3* location; ///< NULL if initial_location is (0,0,0) -- location of interaction attachment with respect to Element center (gets updated by rot/trans functions)
	Vec3* normal; ///< NULL if initial_normal is (0,0,0) -- normal direction of interaction starting from attachement location (auto-updated by rot/trans functions)
	Vec3* tangent; ///< NULL if initial_tangent is (0,0,0) -- tangent of interaction starting from attachement location (auto-updated by rot/trans functions)
	
	Vec3* initial_location; ///< location of interaction from center
	Vec3* initial_normal; ///< normal of interaction from location
	Vec3* initial_tangent; ///< tangent of interaction from location
	
	double potential_map[3]; ///< maps all potentials to stretch, bend, and torsion to bias group rebuilding
} Interaction;

/// This struct holds information for each individual molecule type
typedef struct _Element_Type{
	unsigned int number; ///< how many to create individually (does not contain how many are used in a group)
	unsigned int dof; ///< degrees of freedom particular type (defaults to 5)
	int archetype; ///< archetype (element type number)
	
	Vec3 saxes; ///< Semiaxes of ellipse
	Vec3 saxes2; ///< Semiaxes of ellipse squared
	Vec3 invsaxes2; ///< Inverse of semiaxes of ellipse squared
	double invsaxesvolume; ///< inverse volume of box spanned by semiaxes
	Vec3 mu_pos; ///< Dipole position inside ellipse
	
	double dipole; ///< Initial dipole moment
	Vec3 initial_dipole; ///< Initial dipole moment vector
	
	double Vvdw; ///< LJ epsilon for ellipse
	double sqrtVvdw; ////< square root of LJ epsilon for ellipse
	double inv_avg_area; ////< inverse average oid area (based on average semiaxes)
	double avg_width; ///< average width (in terms of semi-axis) of LJ potential (0 means don't use)
	double rT; ///< test sphere radius based on either physical size *or* LJ potential width
	double radius; ///< average physical size in the far-range
	double sphericity; ///< sphericity
	double IA_average; ///< average of interaction area code
	double charge; ///< overall residual charge of element
	double* Vvdw_coefficients;
	double* IA_coefficients;
	double mass; ///< mass of element type in amu -- used for center of mass and rot/trans bias of groups (more mass => less rot/trans)
	Mat33 MOI; ///< moment of inertia tensor
	
	bool hasmu; ///< Whether or not the ellipse has a dipole in it
	bool still; ///< Whether or not the Element can move
	bool rot_notrans; ///< Whether or not the Element can rotate but not translate
	bool calculate_order; ///< Whether or not an Element is included in averages for order parameter calculations
	
// Visualizer properties
	Vec3 color; ///< color of element in RGB space
	double transparency; ///< transparency of element (0 (opaque) - 1 (completely transparent))
	bool label_elements;
	bool texture_is_dipole;
	string texture; ///< texture file name (leave empty for no texture)
	
	string name; ///< name of this element type
	unsigned int thistype; ///< Type number
	unsigned int nr_Vvdw_coefficients;
	
	unsigned int nr_charges; ///< number of charges inside ellipse
	Vec3* q_pos; ///< Charge position(s) inside ellipse
	double* q; ///< Charge magnitude(s) in esu for ellipse
	
	double* eps_texture; ///< sqrt of epsilon texture for all directions (will usually be NULL)
	double* r_diff; ///< difference from ellipsoid shape in r-direction (NULL if unused)
	double* width_diff; ///< difference from ellipsoid LJ well depth width in r-direction (NULL if unused)
} Element_Type;

/// This struct holds information for one element in the MC simulation
typedef struct _Element{
	Vec3 dipole; ///< Dipole moment (gets updated each time ellipse changes rotational orientation)
	Mat33 rot; ///< Rotation matrix to rotate the Element from [0 0 1] to its lab frame position
	Vec3 center; ///< Cordinates for the center of the Element
	
	Vec3 ds;
	Vec3 current_s;
	Vec3 dtheta;
	double Ekin;
	double max_step;
	double VE;
	double delta; ///< Scale factor for ellipsoid semi-axes
	double gamma_quarter; ///< one quarter of this ellipsoid's gamma value (used for VdW potential of type (a/r)^(12+gamma)-(a/r)^6)
	
	unsigned int nr_trials;
	unsigned int nr_moves;
	unsigned int nr_neighbors;
	unsigned int old_nr_neighbors;
	unsigned int nr_neighbors_above;
	bool fixed;	///< element is just there but can not move or rotate
	
	unsigned int mytype; ///< Type of an individual Element
	unsigned int nr_interactions; ///< nr of bonds on this Element
	
	Element_Type* MyType; ///< points to the respective element type
	Vec3* q_pos; ///< Up to date charge position(s) inside ellipse
	Interaction* interactions; ///< interactions (currently, read: bonds) this particular Element has with other Elements
	Element_Group* group; ///< points to the group the Element belongs to
} Element;

typedef struct _group_save{
	Vec3 dipole;
	Mat33 rot;
	Vec3 center;
	double gamma_quarter;
	Vec3* q_pos;
	Vec3* interaction_location;
	Vec3* interaction_normal;
	Vec3* interaction_tangent;
} group_save;

typedef struct _group_rebuild{
	unsigned int* element_storage;
	unsigned int* visited; // stores number (plus one) of previous element before getting to current element
	unsigned int* to_link; // link number of previous element leading to current element
	unsigned int* randomized; // elements which have been randomized before current element (group relative)
	bool* en_route;
	Vec3* delta_trans; // translation of previous elements
	Mat33* delta_rot; // rotation matrix difference of previous elements
	double* theta2;
} group_rebuild;

static bool exitnow=false;

inline void RotateElement(Element* element, Vec3 &center, Mat33 &delta_rot)
{
	element->center=center+delta_rot*(element->center-center);
	element->rot=delta_rot*element->rot;
	element->dipole=element->rot*element->MyType->initial_dipole;
	// charges
	for(unsigned int j=0; j<element->MyType->nr_charges; j++) element->q_pos[j]=element->rot*element->MyType->q_pos[j];
	// links
	for(unsigned int j=0; j<element->nr_interactions; j++){
		if(element->interactions[j].location) *element->interactions[j].location=element->rot*(*element->interactions[j].initial_location);
		if(element->interactions[j].normal) *element->interactions[j].normal=element->rot*(*element->interactions[j].initial_normal);
		if(element->interactions[j].tangent) *element->interactions[j].tangent=element->rot*(*element->interactions[j].initial_tangent);
	}
}

inline void RotateElement(Element* element, Mat33 &delta_rot)
{
	element->rot=delta_rot*element->rot;
	element->dipole=element->rot*element->MyType->initial_dipole;
	// charges
	for(unsigned int j=0; j<element->MyType->nr_charges; j++) element->q_pos[j]=element->rot*element->MyType->q_pos[j];
	// links
	for(unsigned int j=0; j<element->nr_interactions; j++){
		if(element->interactions[j].location) *element->interactions[j].location=element->rot*(*element->interactions[j].initial_location);
		if(element->interactions[j].normal) *element->interactions[j].normal=element->rot*(*element->interactions[j].initial_normal);
		if(element->interactions[j].tangent) *element->interactions[j].tangent=element->rot*(*element->interactions[j].initial_tangent);
	}
}

/*!
MC_Elements class - stores MC_Elements parameters, constructors, and functions for rotating and translating elements. Fundamental
class for systems being simulated; the entire ensemble is stored in an array of Element
\warning For execution speed, no safety checks are performed when accessing i-th element, so be careful
*/
class MC_Elements{
public:
	Config_Data* configuration; ///< points to configuration
	
	MC_Elements(Config_Data* conf); ///< constructor
	~MC_Elements(); ///< destructor
	
	unsigned int AddElement(const unsigned int type); ///< Add an Element to the ensemble
	Element_Group* AddGroup(const unsigned int type); ///< Add an element group to the ensemble
	
	double Rotate(const unsigned int idx); ///< Randomly rotate a MC_Element with the constraints of configuration->maxrot
	double Rotate(Element_Group* group, Vec3 &center); ///< Randomly rotate group around center with the constraints of configuration->maxrot/number of group members
	double Translate(const unsigned int idx); ///< Randomly translate a MC_Element in a cubic box within the constraints of configuration->maxtrans
	double Translate_Sphere(const unsigned int idx); ///< Same, but spherical box
	double Translate_Cylinder(const unsigned int idx); ///< Same, but cylindrical box
	
	void Rotate(const unsigned int idx, double roll, double yaw, double pitch); ///< Rotate element around its center by (r)oll, (y)aw, (p)itch angles (think of a plane moving in z-direction, x points up, y to the left)
	void Translate(const unsigned int idx, Vec3 delta); ///< Translate an element by delta
	void Translate(Element_Group* group, Vec3 delta); ///< Translate group by delta
	void Rotate(Element_Group* group, Vec3 &center, double roll, double yaw, double pitch); ///< Rotate group around center by (r)oll, (y)aw, (p)itch angles (think of a plane moving in z-direction, x points up, y to the left)
	inline void undo_PBCs(Vec3 &ref, Vec3 &rvec);
	inline void apply_PBCs(Vec3 &rvec); ///< apply periodic boundary conditions to vector
	/// Get i-th Element
	inline Element* operator[](const unsigned int i){ return &Elements[i]; };
	
	void evolve_system(){ evolve_system(false); }
	void benchmark_system(){ evolve_system(true); }
private:
	void evolve_system(bool benchmark);
	static void exit_gracefully(int signum){
		time_t rawtime;
		time (&rawtime);
		char* timestring = asctime(localtime(&rawtime));
		timestring[strlen(timestring)-1]=0;
		signal(signum,MC_Elements::exit_gracefully); // reset signal
		cout << "\nExiting gracefully (received signal " << signum << " on " << timestring << ").\n\n";
		exitnow=true;
	}
	static void exit_not(int signum){
		time_t rawtime;
		time (&rawtime);
		char* timestring = asctime(localtime(&rawtime));
		timestring[strlen(timestring)-1]=0;
		signal(signum,MC_Elements::exit_not); // reset signal
		cout << "\nReceived signal " << signum << " on " << timestring << ". Is this important?\n -> Probably not -- crossing fingers and continuing.\n";
	}
	static void continue_execute(int signum){
		signal(signum,MC_Elements::continue_execute); // reset signal
		cout << "\n... continuing:\n";
	}
	double (MC_Elements::*LJpair)(const unsigned int kk, const unsigned int i); /// Lennard-Jones pair function to use
	Element* Elements;
	Element_Group** Groups; // contains pointers to all allocated groups
	double dt;
	double smallest_volume;
	unsigned int max_nr_neighbors;
	bool neighbor_switch;
	unsigned int* neighbor_store;
	double* pair_energy_store;
	unsigned int* energy_index_store;
	unsigned int* neighbors;
	double* pair_energies;
	unsigned int* energy_indices;
	bool* bound_neighbors; // list of bound elements for collective moves
	unsigned int N, nr_elements, nr_group_elements, elements_allocated, nr_groups, max_group_nrelements, max_charges, nr_individual_elements, nr_bound_neighbors, nr_N_to_kk;
	group_rebuild rebuild_storage;
	group_save** group_storage; // array of arrays - one for each group type
	
	cl_context Context;
	cl_command_queue Queue;
	cl_device_id deviceID;
	cl_program program;
	cl_kernel distance_kernel, lj_kernel, es_kernel, metropolis_kernel;
	
	Vec3* cosmeans; // Stores order parameters
	Vec3* Ms; // Stores total dipole moment of system
	double* Vs; // Stores volume of system
	Vstore* totalV; // Stores potential energy components for systems
	
	double* trans_store;
	double* rot_store;
	unsigned int* running_average;
	
	double msmoved[2];
	unsigned int accept[3]; // Move acceptance counters
	__uint64_t nr_tries_initial, nr_tries_equi, nr_tries;
	unsigned int NpT_accept[3]; // Move acceptance counters
	__uint64_t NpT_tries_initial, NpT_tries_equi, NpT_tries;
	double simulation_time, dt_acc;
	
	cl_double rf_correction;
	cl_double rf_neutralization;
	Vec3* Rmus;
	double* Rdist2;
	double* Rdist2LJ;
	double* Rmirror;
	bool* LJwallMirror;
	double LJlambda; // current value of LJ lambda (dispersive term scaling factor)
	double selfLJlambda; // current value of LJ lambda (dispersive term scaling factor) for internal energies
	Vec3 Evec;
	double epsRF, V_diff;
	bool changevolume;
	unsigned int traj_laststep, step_pos;
	bool traj_created;
	
	void CreateGroupStorage();
	void FreeGroupStorage();
	
	bool setEandT(unsigned int step);
	inline double calc_lambda();
	inline double calc_lambda(unsigned int step);
	inline void update_eps(Vec3 &M, unsigned int step);
	inline void update_rf_correction();
	inline void get_means(Vec3 &means, Vec3 &M);
	void store_trajectory(unsigned int step);
	
	inline void initialize_neighbors();
	inline void extend_neighbors();
	inline void find_neighbors();
	inline bool get_distances(const unsigned int i); ///< fill distance array with distances from i-th Element's center (periodic boundary conditions apply)
	inline bool get_distances_above_index(const unsigned int i); ///< fill distance array with distances from i-th Element's center and above (periodic boundary conditions apply)
	inline Interaction* get_link(const unsigned int i, const unsigned int k);
	inline Vec3 dipole_distance(const unsigned int i, const unsigned int k);
	
	bool recalcall(double *totalVs, bool skip_charges);
	double calc_charge_interaction();
	
	inline double Vwall(const unsigned int i, bool skip_charges, bool &success);
	
	inline double VmuE(const unsigned int i);
	
	inline double phi_q_mu(const unsigned int k, Vec3 &thermu, double &RFp, bool updateRF);
	inline double phi_q(const unsigned int k, Vec3 &thermu, double &RFp, bool updateRF);
	inline Vec3 E_mu_q(const unsigned int k, Vec3 &thermu, Vec3 &RF, bool updateRF);
	inline Vec3 E_mu(const unsigned int k, Vec3 &thermu, Vec3 &RF, bool updateRF);
	
	inline double self_phi_q_mu(const unsigned int k, Vec3 &thermu, double &RFp, bool updateRF);
	inline double self_phi_q(const unsigned int k, Vec3 &thermu, double &RFp, bool updateRF);
	inline Vec3 self_E_mu_q(const unsigned int k, Vec3 &thermu, Vec3 &RF, bool updateRF);
	inline Vec3 self_E_mu(const unsigned int k, Vec3 &thermu, Vec3 &RF, bool updateRF);
	
	double VSS_notouch(const unsigned int kk, const unsigned int i);
	double simpleVLJ_notouch(const unsigned int kk, const unsigned int  i);
	double VLJ_notouch(const unsigned int kk, const unsigned int i);
	double simpleVLJ_touch(const unsigned int kk, const unsigned int i);
	double modulate_touch(const unsigned int kk, const unsigned int i);
	double VSE_touch(const unsigned int kk, const unsigned int i);
	double nmVLJ_notouch(const unsigned int kk, const unsigned int  i);
	inline double touch(const unsigned int kk, const unsigned int i);
	
	inline void Rotate(Element_Group* group, Vec3 &center, Mat33& rot);
	inline void Rotate(Element_Group* group, Mat33& rot);
	inline void Rotate(Element* element, Vec3 &center, Mat33& rot);
	inline Mat33 Rotate(Element* element, Vec3 &center, double &theta);
	inline double GroupShuffle(Element_Group* group, const unsigned int kk);
	inline double ShrinkVolume(unsigned int step);
	void find_bound_neighbors(unsigned int kk, unsigned int level);
	void find_bound_neighbors(unsigned int kk){ find_bound_neighbors(kk,0); }
	
	bool selfInteraction(Element_Group* group, double &VES, double* VLJ, bool skip_charges);
	bool OneStep(const unsigned int kk, double &Ves, double &VLJ, bool skip_charges);
	bool OneStepGroup(const unsigned int kk, double &Ves, double &VLJ, bool skip_charges);
	inline double NpT_step(double* Vstep, double* save_groupVG, unsigned int step, unsigned int eqsteps, bool &adjust_internal);
	double GroupPotentialEnergy(Element_Group* group);
};

inline void MC_Elements::Rotate(Element* element, Vec3 &center, Mat33& rot) ///< Rotate group around center with rotation matrix rot
{
	element->center=center+rot*(element->center-center);
	element->rot=rot*element->rot;
	
	element->dipole=element->rot*element->MyType->initial_dipole;
	// charges
	for(unsigned int j=0; j<element->MyType->nr_charges; j++) element->q_pos[j]=element->rot*element->MyType->q_pos[j];
	// links
	for(unsigned int j=0; j<element->nr_interactions; j++){
		if(element->interactions[j].location) *element->interactions[j].location=element->rot*(*element->interactions[j].initial_location);
		if(element->interactions[j].normal) *element->interactions[j].normal=element->rot*(*element->interactions[j].initial_normal);
		if(element->interactions[j].tangent) *element->interactions[j].tangent=element->rot*(*element->interactions[j].initial_tangent);
	}
}

inline void MC_Elements::Rotate(Element_Group* group, Vec3 &center, Mat33& rot) ///< Rotate group around center with rotation matrix rot
{
	// Apply to group elements
	Vec3 delta;
	Element* element;
	group->center.vec[0]=0.0; group->center.vec[1]=0.0; group->center.vec[2]=0.0;
	for(unsigned int i=0; i<group->nr_elements; i++){
		element=&Elements[group->elements[i]];
		// Location of element first
		delta=element->center-center;
		delta=rot*delta;
		element->center=center+delta;
		group->center+=element->center;
		// Rotate element
		RotateElement(element,rot);
	} // done
	group->center*=group->inv_nr_elements;
}

inline void MC_Elements::Rotate(Element_Group* group, Mat33& rot) ///< Rotate group with rotation matrix rot
{
	// Apply to group elements
	Vec3 delta;
	Element* element;
	for(unsigned int i=0; i<group->nr_elements; i++){
		element=&Elements[group->elements[i]];
		// Location of element first
		delta=element->center-group->center;
		delta=rot*delta;
		element->center=group->center+delta;
		// Rotate element
		RotateElement(element,rot);
	} // done
}

inline double IA0(Element_Type* type, Vec3 direction)
{
	double a=direction*direction;
	if(a>EPS*EPS){
		a = 1.0/a;
		double c=direction.vec[0]*direction.vec[0]*(type->invsaxes2.vec[1]*type->invsaxes2.vec[2])+direction.vec[1]*direction.vec[1]*(type->invsaxes2.vec[0]*type->invsaxes2.vec[2])+direction.vec[2]*direction.vec[2]*(type->invsaxes2.vec[0]*type->invsaxes2.vec[1]);
		// a*beta^2 + b*beta + c = 0
		// beta_+/- = -b/2a +/- sqrt(b^2 - 4*a*c)/2a
		// beta+ * beta_ = (-b + sqrt(b^2 - 4*a*c)) * (-b - sqrt(b^2 - 4*a*c))/(2a)^2 = (b^2 - (b^2 - 4ac))/4a^2 = 4ac/4a^2 = c/a
		return pi*type->inv_avg_area/sqrt(c*a);
	}
	return 0.0;
}

inline double IA(Element_Type* type, Vec3 direction)
{
	double a=direction*direction;
	if(a>EPS*EPS){
		a = 1.0/a;
		double c=direction.vec[0]*direction.vec[0]*(type->invsaxes2.vec[1]*type->invsaxes2.vec[2])+direction.vec[1]*direction.vec[1]*(type->invsaxes2.vec[0]*type->invsaxes2.vec[2])+direction.vec[2]*direction.vec[2]*(type->invsaxes2.vec[0]*type->invsaxes2.vec[1]);
		// a*beta^2 + b*beta + c = 0
		// beta_+/- = -b/2a +/- sqrt(b^2 - 4*a*c)/2a
		// beta+ * beta_ = (-b + sqrt(b^2 - 4*a*c)) * (-b - sqrt(b^2 - 4*a*c))/(2a)^2 = (b^2 - (b^2 - 4ac))/4a^2 = 4ac/4a^2 = c/a
		double IAval=pi*type->inv_avg_area/sqrt(c*a);
		if(type->IA_coefficients){
			if(IAval>type->IA_coefficients[0]){
				return type->IA_coefficients[1]+IAval*(type->IA_coefficients[2]+type->IA_coefficients[3]*IAval);
			} else{
				return type->IA_coefficients[5]+IAval*(type->IA_coefficients[6]+type->IA_coefficients[4]*IAval);
			}
		}
		return IAval;
	}
	return 0.0;
}

inline void adjust_cut2plus(Config_Data* configuration)
{
	if(configuration->cut){
		configuration->cut2plus=configuration->ljcut;
		if(configuration->escut2>configuration->ljcut2) configuration->cut2plus=configuration->rfcut;
		configuration->cut2plus+=configuration->max_neighbor_trans;
		configuration->cut2plus=qqpwr(configuration->cut2plus,2);
	}
}

inline double determine_epsilon(Element_Type* element_type, double rp)
{
	if(element_type->nr_Vvdw_coefficients){
		double epsilon=0.0;
		double xacc=1.0;
		for(unsigned int j=1; j<element_type->nr_Vvdw_coefficients; j++){
			epsilon+=element_type->Vvdw_coefficients[j]*xacc;
			xacc*=rp;
		}
		epsilon/=qqpwr(element_type->Vvdw_coefficients[0]+rp,element_type->nr_Vvdw_coefficients-2);
		return epsilon;
	} else return element_type->sqrtVvdw;
}

/// recursive algorithm to "walk" from element current to all linked ones and store shortest distances in storage (based on Dijkstra's algorithm)
inline void follow_links(int* storage, Element_Group* group, Element* elements, int current, int depth, int limit, unsigned int &statistics)
{
	statistics++;
	storage[current]=depth;
	if(depth<limit){ // stop walking when the distance limit has been reached
		int next;
		for(unsigned int i=0; i<elements[group->elements[current]].nr_interactions; i++){ // follow links in all directions
			next=elements[group->elements[current]].interactions[i].partner;
			if((storage[next]>depth+1) || (storage[next]<0)){ // only follow link if other end has not been visited or current path is shorter
				follow_links(storage,group,elements,next,depth+1,limit,statistics);
			}
		}
	}
}

/// non-recursive algorithm to "walk" from element current to all linked ones and store shortest distances in storage (based on Dijkstra's algorithm)
/// - this is the one currently used (slightly faster)
inline void follow_links2(int* element_storage, int* storage, Element_Group* group, Element* elements, int start, int limit, unsigned int &statistics, unsigned int* closed_nodes, unsigned int &cn_count)
{
	statistics++;
	int distance=0;
	storage[start]=distance;
	
	unsigned int next, newcount;
	unsigned int nr_elements=0;
	for(unsigned int i=0; i<elements[group->elements[start]].nr_interactions; i++){
		next=elements[group->elements[start]].interactions[i].partner;
		if(storage[next]<0){
			nr_elements++;
			element_storage[i<<1]=next;
		}
	}
	
	bool alternate=0;
	unsigned int current;
	cn_count=0;
	while((distance<limit) && (nr_elements>0)){
		distance++;
#if DEBUG_LEVEL>2
		cout << "distance: " << distance << "\n";
#endif
		newcount=0;
		for(unsigned int i=0; i<nr_elements; i++){
			current=element_storage[(i<<1)+alternate];
#if DEBUG_LEVEL>2
			cout << current+1 << "\n";
#endif
			storage[current]=distance;
			statistics++;
			for(unsigned int j=0; j<elements[group->elements[current]].nr_interactions; j++){
				next=elements[group->elements[current]].interactions[j].partner;
				if((storage[next]>distance) || (storage[next]<0)){ // only follow link if other end has not been visited or current path is shorter
					unsigned int id=(newcount<<1)+(!alternate);
					// find out if particular element is already in list of next elements ...
					bool found=false;
					for(unsigned int l=0; l<newcount; l++){ // if there are closed loops, no need to stay on multiple tracks at once
						if((int)next==element_storage[(l<<1)+(!alternate)]){
							found=true;
							break;
						}
					}
					if(!found){
						// find out if we've been there previously
						for(unsigned int l=0; l<nr_elements; l++){ // if there are closed loops, no need to stay on multiple tracks at once
							if((int)next==element_storage[(l<<1)+(alternate)]){
								found=true;
								break;
							}
						}
					}
					if(found){ // already been there (either shorter distance or same -- that's why path does not need to be retraced)
#if DEBUG_LEVEL>2
						cout << "closed loop between " << start+1 << " and " << next+1 << ".\n";
#endif
						if(closed_nodes){
#if DEBUG_LEVEL>2
							cout << "closed loop between " << start+1 << " and " << next+1 << ".\n";
#endif
							closed_nodes[cn_count]=next;
							cn_count++;
						}
					} else{
						element_storage[id]=next;
						newcount++;
					}
				}
			}
		}
		alternate=(!alternate); // flip-flop
		nr_elements=newcount;
#if DEBUG_LEVEL>2
		cout << "---\n";
#endif
	}
}

inline void DetermineGroupVolume(Element_Group* group, Element* elements)
{
	cout << "\t-> Determining excluded volume of <" << group->Type->name << "> using Monte-Carlo routine (may take a while)\n";
#if DEBUG_LEVEL>1
	double tstart = (double)clock();
#endif
	group->Type->overlaps = NULL;
	int* distance_matrix = NULL;
	if(group->Type->adjust_overlap>0.0){
		group->Type->overlaps=new double[group->nr_elements];
		int* element_storage = new int[group->nr_elements<<2];
		int* storage = new int[group->nr_elements+1];
		distance_matrix = new int[group->nr_elements*group->nr_elements];
		unsigned int cn_count;
		unsigned int steps=0;
		for(unsigned int i=0; i<group->nr_elements; i++){ // get distance from element i to j
			memset(storage,0xFF,group->nr_elements*sizeof(int)); // puts -1 in each field
			follow_links2(element_storage,storage,group,elements,i,group->nr_elements-1,steps,NULL,cn_count);
			for(unsigned int j=0; j<group->nr_elements; j++) distance_matrix[i*group->nr_elements+j]=storage[j];
		}
		delete[] element_storage;
		delete[] storage;
	}
	// Find bounding box first
	Vec3 box(0.0);
	Vec3 startpoint=elements[group->elements[0]].center;
	Vec3 box_corner=startpoint;
	for(unsigned int i=0; i<group->nr_elements; i++){
		if(group->Type->overlaps) group->Type->overlaps[i]=0.0;
		Element* curr=&elements[group->elements[i]];
		Vec3 dist=curr->center-startpoint;
		double maxradius=maxd(curr->MyType->saxes.vec,3);
		for(unsigned int j=0; j<3; j++){
			if(box_corner.vec[j]>dist.vec[j]-maxradius) box_corner.vec[j]=dist.vec[j]-maxradius;
			if(box.vec[j]<dist.vec[j]+maxradius) box.vec[j]=dist.vec[j]+maxradius;
		}
	}
	box-=box_corner;
	double boundingV=box.vec[0]*box.vec[1]*box.vec[2];
	unsigned int Nin=0;
	unsigned int N=0;
	unsigned int n_overlaps;
	double inside;
	unsigned int* already_in = new unsigned int[group->nr_elements];
#if DEBUG_LEVEL>1
	double var;
#endif
	do{
		Vec3 point=box_corner;
#ifdef USE_CMWC4096
		point.vec[0]+=ranQ()*box.vec[0];
		point.vec[1]+=ranQ()*box.vec[1];
		point.vec[2]+=ranQ()*box.vec[2];
#else
		point.vec[0]+=ran2(idum)*box.vec[0];
		point.vec[1]+=ran2(idum)*box.vec[1];
		point.vec[2]+=ran2(idum)*box.vec[2];
#endif
		point=point+startpoint;
		n_overlaps=0;
		for(unsigned int i=0; i<group->nr_elements; i++){
			Element* curr=&elements[group->elements[i]];
			if(Point_in_Ellipsoid(point,curr->MyType->saxes,curr->center,curr->rot)){
				if(n_overlaps==0){
					Nin++; // only single-count
					if(!group->Type->overlaps) break;
				}
				already_in[n_overlaps]=i;
				n_overlaps++;
			}
		}
		if(group->Type->overlaps){
			for(unsigned int i=0; i<n_overlaps; i++){
				for(unsigned int j=i+1; j<n_overlaps; j++){
					int dist=distance_matrix[already_in[i]*group->nr_elements+already_in[j]];
					if(dist>0) group->Type->overlaps[dist]+=1.0; else cout << "WARNING: Overlap between non-attached group elements " << already_in[i]+1 << " and " << already_in[i]+1 << ".\n";
				}
			}
		}
		N++;
		inside=(double)Nin/N;
	} while(((1.0-inside)/(N*inside)>1E-6) || (Nin<1000*group->nr_elements)); // first statement is testing for relative error of 1E-3 (0.1%)
	group->Type->Volume=inside*boundingV;
#if DEBUG_LEVEL>1
	var=boundingV*boundingV*inside*(1-inside)/N;
	double tend = (double)clock();
	cout << "\t\t-> Originally excluded volume is " << group->Type->Volume << " +/- " << sqrt(var) << " Angström³ (" << (unsigned int)N << " iterations, took " << (tend-tstart)/CLOCKS_PER_SEC << " s)\n";
#endif
	if(group->Type->overlaps){
#if DEBUG_LEVEL>1
		cout << "\t\t-> Overlap volume fractions for nearest neighbors: ";
#endif
		double* ovl_vol = new double[group->nr_elements];
		for(unsigned int i=0; i<group->nr_elements; i++) ovl_vol[i]=0.0;
		for(unsigned int i=0; i<group->nr_elements-1; i++){
			for(unsigned int j=i+1; j<group->nr_elements; j++){
				int dist=distance_matrix[i*group->nr_elements+j];
				if(dist>0) ovl_vol[dist]+=0.5*4.0*pi/3.0*(elements[group->elements[i]].MyType->saxes.vec[0]*elements[group->elements[i]].MyType->saxes.vec[1]*elements[group->elements[i]].MyType->saxes.vec[2]+elements[group->elements[j]].MyType->saxes.vec[0]*elements[group->elements[j]].MyType->saxes.vec[1]*elements[group->elements[j]].MyType->saxes.vec[2]);
			}
		}
		unsigned int dist_factors=0;
		for(unsigned int i=1; i<group->nr_elements; i++){
			group->Type->overlaps[i]/=(double)N;
			if(ovl_vol[i]>0.0) group->Type->overlaps[i]*=boundingV/ovl_vol[i];
			if(group->Type->overlaps[i]>EPS){
#if DEBUG_LEVEL>1
				if(i>1) cout << ", ";
				cout << group->Type->overlaps[i]*100 << "% (" << i << ".)";
#endif
				if(i>=group->Type->bond_range) dist_factors++;
			}
		}
#if DEBUG_LEVEL>1
		cout << "\n";
#endif
		if((group->Type->bond_range>0) && (group->nr_elements>1)){
			bool first=true;
			cout << "\t\t\t-> Adjusting distance dependent internal energy scaling factors: ";
			for(unsigned int i=1; i<group->nr_elements; i++){
				if(group->Type->overlaps[i]>EPS){
					if(i>=group->Type->bond_range){
						group->Type->range_factors[i].distance=i;
						group->Type->range_factors[i].factor=1.0-group->Type->overlaps[i];
						if(!first) cout << ", ";
						if(first) first=false;
						cout << group->Type->range_factors[i].factor*100 << "% (" << i << ".)";
					}
				}
			}
		}
		cout << "\n";
	}
}

inline void DetermineGroupRings(Element_Group* group, Element* elements)
{
#if DEBUG_LEVEL>0
	cout << "Determining closed loops for group <" << group->Type->name << ">\n";
#endif
	unsigned int steps=0;
	unsigned int rings=0;
	unsigned int** ring_storage=NULL; // array of ringe elements array (whose zero element is the number of elements in it)
	if(group->nr_elements>1){
		int* element_storage = new int[group->nr_elements<<2];
		int* storage = new int[group->nr_elements+1];
		int* back_storage = new int[group->nr_elements+1];
		unsigned int* closed_nodes = new unsigned int[group->nr_elements>>1];
		memset(storage,0xFF,group->nr_elements*sizeof(int)); // puts -1 in each field
		follow_links2(element_storage,storage,group,elements,0,group->nr_elements-1,steps,closed_nodes,rings);
		if(rings>0){
			unsigned int count=0;
			ring_storage=(unsigned int**)malloc(rings*sizeof(unsigned int*));
			for(unsigned int i=0; i<rings; i++){
				ring_storage[count]=(unsigned int*)malloc(3*sizeof(unsigned int));
				ring_storage[count][0]=2;
				bool found=false;
				for(unsigned int j=0; j<count; j++){
					if(closed_nodes[i]==ring_storage[j][1]){
						found=true;
						break;
					}
				}
				if(!found){
					ring_storage[count][1]=closed_nodes[i];
					count++;
				}
			}
			rings=count;
			unsigned int i=0;
			while(i<rings){
#if DEBUG_LEVEL>2
				cout << "\t-> Potential ring(s) around element " << ring_storage[i][1]+1 << ".\n";
#endif
				memset(back_storage,0xFF,group->nr_elements*sizeof(int)); // puts -1 in each field
				unsigned int temp=0;
//				follow_links2(element_storage,back_storage,group,elements,ring_storage[i][1],storage[ring_storage[i][1]],steps,closed_nodes,temp);
				follow_links2(element_storage,back_storage,group,elements,ring_storage[i][1],group->nr_elements-1,steps,closed_nodes,temp);
				unsigned int members=0;
				unsigned int current_count=0;
				if(temp>0){
#if DEBUG_LEVEL>3
					cout << "\t\t" << temp << "\n";
#endif
					unsigned int jj=0;
					unsigned int diff_count=0;
					while(jj<temp){
#if DEBUG_LEVEL>3
						cout << "\t\t\t" << jj << "\n";
#endif
						bool found=false;
						for(unsigned int j=jj; j<temp; j++){
							jj=j+1;
#if DEBUG_LEVEL>3
						cout << "\t\t\t\t" << storage[ring_storage[i][1]] << " > " << storage[closed_nodes[j]] << "? (" << closed_nodes[j] << ")\n";
#endif
							if(storage[ring_storage[i][1]]>storage[closed_nodes[j]]){ // choose closed loop which is closer to where we started
								bool different_end=true;
								for(unsigned int k=1; k<diff_count; k++){
									if(closed_nodes[j]==ring_storage[i-k][2]){
										different_end=false;
										break;
									}
								}
								if(different_end){
									diff_count++;
									ring_storage[i][2]=closed_nodes[j];
									found=true;
									break;
								}
							} else{
								if((group->nr_elements==3) && (storage[ring_storage[i][1]]==storage[closed_nodes[j]])){
									diff_count++;
									ring_storage[i][2]=closed_nodes[j];
									found=true;
								}
							}
						}
						if(found){
							rings++;
							ring_storage=(unsigned int**)realloc(ring_storage,rings*sizeof(unsigned int*));
							if(!ring_storage){
								cout << "ERROR: Not enough memory for ring storage.\n";
								exit(2);
							}
							for(unsigned int z=rings-1; z>i; z--) ring_storage[z]=ring_storage[z-1]; // move original entries
							// new rings need initialization
							ring_storage[i+1]=(unsigned int*)malloc(3*sizeof(unsigned int));
							ring_storage[i+1][1]=ring_storage[i][1];
							ring_storage[i+1][2]=ring_storage[i][2];
							unsigned int current=ring_storage[i][2];
							int distance=back_storage[current];
#if DEBUG_LEVEL>2
							cout << "\t\t-> other end is " << current+1 << " (" << distance << " steps away), now determining paths leading back ...\n";
#endif
							int* path=(int*)malloc(2*sizeof(int));
							unsigned int pathcount=0;
							for(unsigned int j=0; j<elements[group->elements[current]].nr_interactions; j++){
								if((back_storage[elements[group->elements[current]].interactions[j].partner]==distance-1) || (back_storage[elements[group->elements[current]].interactions[j].partner]==distance)){
									if(pathcount>=2){
										path=(int*)realloc(path,(pathcount+1)*sizeof(int));
										if(!path){
											cout << "ERROR: Not enough memory to store path starts.\n";
											exit(2);
										}
									}
									path[pathcount]=elements[group->elements[current]].interactions[j].partner;
									pathcount++;
								}
							}
							if(pathcount>2){
								unsigned int additional_rings=((pathcount*pathcount-pathcount)>>1)-1; // n^2-n is always even
								rings+=additional_rings;
								ring_storage=(unsigned int**)realloc(ring_storage,rings*sizeof(unsigned int*));
								if(!ring_storage){
									cout << "ERROR: Not enough memory for ring storage.\n";
									exit(2);
								}
								// move original entries
								for(unsigned int j=rings-1; j>i+additional_rings; j--) ring_storage[j]=ring_storage[j-additional_rings];
								for(unsigned int j=i+1; j<i+additional_rings+1; j++){
										// new rings need initialization
										ring_storage[j]=(unsigned int*)malloc(3*sizeof(unsigned int));
										ring_storage[j][1]=ring_storage[i][1];
										ring_storage[j][2]=ring_storage[i][2];
								}
							}
							unsigned int ii=i;
							for(unsigned int n=0; n<pathcount-1; n++){
								for(unsigned int m=n+1; m<pathcount; m++){
									members=2+back_storage[path[n]]+back_storage[path[m]];
									ring_storage[i]=(unsigned int*)realloc(ring_storage[i],(1+members)*sizeof(unsigned int));
									if(!ring_storage[i]){
										cout << "ERROR: Not enough memory to store ring members.\n";
										exit(2);
									}
									ring_storage[i][0]=members;
									ring_storage[i][3]=path[n];
									ring_storage[i][4]=path[m];
									count=5;
									current=path[n];
									distance=back_storage[path[n]]-1;
									for(int l=distance; l>0; l--){
										for(unsigned int k=0; k<elements[group->elements[current]].nr_interactions; k++){
											if(back_storage[elements[group->elements[current]].interactions[k].partner]==l){
												ring_storage[i][count]=elements[group->elements[current]].interactions[k].partner;
												current=ring_storage[i][count];
												count++;
												break;
											}
										}
									}
									current=path[m];
									distance=back_storage[path[m]]-1;
									bool invalid=false;
									for(int l=distance; l>0; l--){
										for(unsigned int k=0; k<elements[group->elements[current]].nr_interactions; k++){
											if(back_storage[elements[group->elements[current]].interactions[k].partner]==l){
												ring_storage[i][count]=elements[group->elements[current]].interactions[k].partner;
												current=ring_storage[i][count];
												for(unsigned int z=5; z<count; z++){
													if(current==ring_storage[i][z]){
														invalid=true;
														break;
													}
												}
												count++;
												break;
											}
											if(invalid) break;
										}
									}
									if(!invalid && (ii>0)){
										bool already_exists=false;
										for(unsigned int k=0; k<ii; k++){
											bool exists=true;
											if(ring_storage[k][0]==ring_storage[i][0]){
												for(unsigned int l=1; l<=ring_storage[k][0]; l++){
													bool foundit=false;
													for(unsigned int p=1; p<=ring_storage[i][0]; p++){
														if(ring_storage[k][l]==ring_storage[i][p]){
															foundit=true;
															break;
														}
													}
													if(!foundit){
														exists=false;
														break;
													}
												}
											} else exists=false;
											already_exists|=exists;
											if(already_exists) break;
										}
										invalid=already_exists;
									}
									if(invalid){ // get rid of invalid loops
#if DEBUG_LEVEL>2
										cout << "\t\t\t<- This is not the loop we're looking for (closes before startpoint).\n";
#endif
										free(ring_storage[i]);
										rings--;
										for(unsigned int z=i; z<rings; z++) ring_storage[z]=ring_storage[z+1];
										ring_storage=(unsigned int**)realloc(ring_storage,rings*sizeof(unsigned int*));
									} else{
										current_count++;
#if DEBUG_LEVEL>0
										cout << "\t\t-> Ring has " << ring_storage[i][0] << " members";
#if DEBUG_LEVEL>1
										cout << ":\n\t\t\t";
										for(unsigned int j=1; j<=ring_storage[i][0]; j++){
											cout << ring_storage[i][j]+1;
											if(j<ring_storage[i][0]){
												cout << ", ";
												if(j%10==0) cout << "\n\t\t\t";
											}
										}
										cout << "\n";
#else
										cout << ".\n";
#endif // DEBUG_LEVEL>1
#endif // DEBUG_LEVEL>0
										i++;
									}
								}
							}
							free(path);
						}
					}
					if(ring_storage[i]) free(ring_storage[i]);
					rings--;
					for(unsigned int z=i; z<rings; z++) ring_storage[z]=ring_storage[z+1];
					ring_storage=(unsigned int**)realloc(ring_storage,rings*sizeof(unsigned int*));
				} else{ // should not happen
					cout << "Something is badly wrong: No closed loop on way back ...\n";
					exit(42);
				}
				if(current_count==0){
					cout << "\t\t<- No new rings found.\n";
					if(rings>0){
						free(ring_storage[i]);
						rings--;
						for(unsigned int z=i; z<rings; z++) ring_storage[z]=ring_storage[z+1];
						ring_storage=(unsigned int**)realloc(ring_storage,rings*sizeof(unsigned int*));
					}
				}
			}
		}
		delete[] element_storage;
		delete[] storage;
		delete[] back_storage;
		delete[] closed_nodes;
	}
#if DEBUG_LEVEL>0
	cout << "<- Found " << rings << " closed loops (took " << steps << " steps).\n";
#endif
	group->Type->elements_in_ring=new bool[group->nr_elements];
	memset(group->Type->elements_in_ring,0,group->nr_elements*sizeof(bool));
	for(unsigned int i=0; i<rings; i++){
		for(unsigned int j=1; j<=ring_storage[i][0]; j++) group->Type->elements_in_ring[ring_storage[i][j]]=true;
	}
	// Take care of trivial rings (two links between two elements)
	unsigned int trivialcount=0;
#if DEBUG_LEVEL>0
	cout << "Determining trivial loops for group <" << group->Type->name << ">\n";
#endif
	unsigned int* multilinks=new unsigned int[group->nr_elements];
	for(unsigned int i=0; i<group->nr_elements; i++){
		Element* element=&elements[group->elements[i]];
		memset(multilinks,0,group->nr_elements*sizeof(unsigned int));
		for(unsigned int j=0; j<element->nr_interactions; j++){
			multilinks[element->interactions[j].partner]++;
			if((multilinks[element->interactions[j].partner]>1) && (i>element->interactions[j].partner)){ // at least two links to same partner -> ring
#if DEBUG_LEVEL>2
				cout << "\t-> Trivial ring between elements " << i+1 << " and " << element->interactions[j].partner+1 << ".\n";
#endif
				trivialcount++;
				group->Type->elements_in_ring[i]=true;
				group->Type->elements_in_ring[element->interactions[j].partner]=true;
			}
		}
	}
#if DEBUG_LEVEL>0
	cout << "<- Found " << trivialcount << " trivial loops.\n";
#endif
	delete[] multilinks;
	group->Type->nr_rings=rings;
	group->Type->rings=ring_storage;
}

/// Function that fills the bitfield of a group with the answer to the question: Is element i within bond_range # of bonds of k?
inline void SetGroupRangeField(Element_Group* group, Element* elements)
{
#if DEBUG_LEVEL>0
	cout << "Creating bond range bitfield for group <" << group->Type->name << "> (bond range: " << group->Type->bond_range << ")\n";
#endif
	unsigned int steps=0;
	// only create bitfield if bond_range is larger than one since nearest neighbors will always be there ...
	if((group->Type->bond_range>0) && (group->nr_elements>1)){
		// bitshift right by 5 is division by 32 (divide by two once more since we need 1/2(n^2-n), multiply by 4 for 4 bits; plus one for safety: 9>>(5-2)=1 ...)
		// - use 4 bit to encode special distances, key:
		//   0 ... less than bond_range
		//   1 ... more than bond_range, but not at a specific distance
		//   2-15 ... more than bond_range, at special distance according to distance table (entries with special_distance[2..15-2] != 0.0)
		unsigned int range_field_size=((group->nr_elements*(group->nr_elements-1))>>(6-SPECIAL_DISTANCES_BIT_EXPONENT))+1;
		group->Type->range_field=new __uint32_t[range_field_size];
		memset(group->Type->range_field,0,range_field_size<<2); // memset works on bytes: int32 is 4 bytes (hence the left shift by 2)
		
		int nb;
		unsigned int ij;
		int* element_storage = new int[group->nr_elements<<2];
		int* storage = new int[group->nr_elements+1];
		unsigned int cn_count;
		
		for(unsigned int i=1; i<group->nr_elements; i++){ // get distance from element i to j
			memset(storage,0xFF,group->nr_elements*sizeof(int)); // puts -1 in each field
			follow_links2(element_storage,storage,group,elements,i,group->nr_elements-1,steps,NULL,cn_count);
			for(unsigned int j=0; j<i; j++){
				nb=storage[j];
				// bitmatrix is [0..i-1][0..# group elements]
				// addressing is i*(i-1)/2+j for j<i
				ij=((i*(i-1))>>1)+j;
				// nb=0 means we did not get to the element
				// set bit (i modulo 32) in 32 bit integer i/32
				__uint32_t special=((nb>=(int)group->Type->bond_range));
				for(unsigned k=0; k<SPECIAL_DISTANCES; k++){
					if(nb==(int)group->Type->range_factors[k].distance){
						special=k+2;
						break;
					}
				}
				group->Type->range_field[ij>>(5-SPECIAL_DISTANCES_BIT_EXPONENT)] |= (special<<((ij&31)<<SPECIAL_DISTANCES_BIT_EXPONENT)); // for b being a power of 2: (a%b = a&(b-1))
#if DEBUG_LEVEL>2
				cout << j+1 << "<->" << i+1 << " closer than " << group->Type->bond_range << " bonds: ";
				if((nb<(int)group->Type->bond_range) && (nb>0)) cout << "yes"; else cout << "no";
				cout << " (" << storage[j] << " bonds away";
				if(special>1) cout << ", special scaling: " << group->Type->range_factors[special-2].factor;
				cout << ")\n";
#endif
			}
		}
		delete[] storage;
		delete[] element_storage;
	} else group->Type->range_field=NULL;
#if DEBUG_LEVEL>0
	cout << "<- Created bond range bitfield for group <" << group->Type->name << "> (" << group->nr_elements << " elements, took " << steps << " steps)\n";
#endif
}

/// Determine effectiv radius of sphere able to contain group (center is center of sphere)
inline double GroupRadius(Element_Group* group, Element* elements, Vec3 &center)
{
#if DEBUG_LEVEL>1
	cout << "Calculating effective sphere radius from current center.\n";
#endif
	double radius=0.0;
	double distance;
	Vec3 dist_vec;
	unsigned int start=0;
	for(unsigned int i=0; i<group->nr_elements; i++){
		dist_vec=center-elements[group->elements[i]].center;
		distance=dist_vec.V3Norm();
		// correct distance by maximum sigma*2^1/6 / 2 (average minimum energy distance/2)
		distance+=maxd(elements[group->elements[i]].MyType->saxes.vec,3)*two_to_onesixth; // use semi-axis and convert maximum semi-axis to rmin
		if(distance>radius){
			start=i;
			radius=distance;
		}
	}
#if DEBUG_LEVEL>1
	cout << "<- Done, distance from furthest element (" << start+1 << ") is: " << radius << " Angström\n";
	cout << "Searching for center with smaller bounding sphere radius ...\n";
#endif
	// optimize radius, use just found element as periphery starting point
	Vec3 A=elements[group->elements[start]].center;
	double newdist=0.0;
	unsigned int furthest=0;
	// get furthest element from starting point
	for(unsigned int i=0; i<group->nr_elements; i++){
		if(i!=start){
			dist_vec=A-elements[group->elements[i]].center;
			distance=dist_vec.V3Norm();
			// correct distances by maximum sigma*2^1/6 / 2 (average minimum energy distance/2)
			distance+=maxd(elements[group->elements[i]].MyType->saxes.vec,3)*two_to_onesixth; // use semi-axis and convert maximum semi-axis to rmin
			if(distance>newdist){
				furthest=i;
				newdist=distance;
			}
		}
	}
	Vec3 B=elements[group->elements[furthest]].center;
	// correct distances by maximum sigma*2^1/6 / 2 (average minimum energy distance/2)
	// calculate r_m_A first
	dist_vec=A-B; // points towards A
	// convert to unit vector
	distance=dist_vec.V3Norm();
	if(distance>1E-10){
		dist_vec/=distance; // unit vector
		// now correct with L-J sphere radius
		A+=dist_vec*(maxd(elements[group->elements[start]].MyType->saxes.vec,3))*two_to_onesixth; // use semi-axis and convert maximum semi-axis to rmin
		// same thing for r_m_B (with unit vector in other direction
		B-=dist_vec*(maxd(elements[group->elements[furthest]].MyType->saxes.vec,3))*two_to_onesixth; // use semi-axis and convert maximum semi-axis to rmin
	}
	Vec3 new_center=(A+B)/2.0;
	newdist=0.0;
	unsigned int furthest2=0;
	// get element which is furthest from center of just found ones
	for(unsigned int i=0; i<group->nr_elements; i++){
		if(i!=furthest){
			dist_vec=new_center-elements[group->elements[i]].center;
			distance=dist_vec.V3Norm();
			// correct distances by maximum sigma*2^1/6 / 2 (average minimum energy distance/2)
			distance+=maxd(elements[group->elements[i]].MyType->saxes.vec,3)*two_to_onesixth; // use semi-axis and convert maximum semi-axis to rmin
			if(distance>newdist){
				furthest2=i;
				newdist=distance;
			}
		}
	}
	Vec3 C=elements[group->elements[furthest2]].center;
	// calculate r_m_C
	dist_vec=C-new_center; // points towards C
	// convert to unit vector
	distance=dist_vec.V3Norm();
	if(distance>1E-10){
		dist_vec/=distance;
		// now correct with L-J sphere radius
		C+=dist_vec*(maxd(elements[group->elements[furthest2]].MyType->saxes.vec,3))*two_to_onesixth; // use semi-axis and convert maximum semi-axis to rmin
	}
	
	// now have three points A,B,C to base new minimum center on
	if((start==furthest) && (start==furthest2)){ // all elements are in one point
#if DEBUG_LEVEL>1
		cout << "-> all group element in same point.\n";
#endif
		new_center=elements[group->elements[start]].center;
	} else{
		if(furthest2==start){ // start and furthest are the only two points on periphery => new center is their middle
#if DEBUG_LEVEL>1
			cout << "-> new center determined from furthest apart elements (" << start+1 << " and " << furthest+1 << ").\n";
#endif
			new_center=(A+B)/2.0;
			// correct distances by maximum sigma*2^1/6 / 2 (average minimum energy distance/2)
			// calculate r_m_A first
			dist_vec=A-B; // points towards A
			// convert to unit vector
			distance=dist_vec.V3Norm();
			if(distance>1E-10){
				dist_vec/=distance;
				// now correct with L-J sphere radius
				A+=dist_vec*(maxd(elements[group->elements[start]].MyType->saxes.vec,3))*two_to_onesixth; // use semi-axis and convert maximum semi-axis to rmin
				// same thing for r_m_B (with unit vector in other direction
				B-=dist_vec*(maxd(elements[group->elements[furthest]].MyType->saxes.vec,3))*two_to_onesixth; // use semi-axis and convert maximum semi-axis to rmin
			}
			new_center=(A+B)/2.0;
		} else{ // triangle defines sphere, circumcenter is new center
			// calculation according to wikipedia (circumscribed circle)
			// a = A-C
			// b = B-C
			// circumcenter = (|a|²b-|b|²a)x(axb) / (2*|axb|²) + C
#if DEBUG_LEVEL>1
			cout << "-> triangulating new center (elements " << start+1 << ", " << furthest+1 << ", and " << furthest2+1 << ").\n";
#endif
			new_center=C;
			Vec3 a=A-C;
			Vec3 b=B-C;
			Vec3 c=a;
			c.V3Cross(b); // c=axb
			Vec3 d=b*(a*a)-(a*(b*b));
			new_center+=(d.V3Cross(c))/(2.0*(c*c));
		}
	}
	// calculate new radius
	start=0;
	double new_radius=0.0;
	for(unsigned int i=0; i<group->nr_elements; i++){
		dist_vec=new_center-elements[group->elements[i]].center;
		distance=dist_vec.V3Norm();
		// correct distance by maximum sigma*2^1/6 / 2 (average minimum energy distance/2)
		distance+=maxd(elements[group->elements[i]].MyType->saxes.vec,3)*two_to_onesixth; // use semi-axis and convert maximum semi-axis to rmin
		if(distance>new_radius){
			start=i;
			new_radius=distance;
		}
	}
#if DEBUG_LEVEL>1
	cout << "-> element furthest away from center: " << start+1 << " (radius: " << new_radius << ")\n";
#endif
	// use new radius and change center in case new one is smaller ...
	if(new_radius<radius){
		radius=new_radius;
		center=new_center;
#if DEBUG_LEVEL>1
		cout << "<- Done, new center is: (" << center.vec[0] << ", " << center.vec[1] << ", " << center.vec[2] << ").\n";
#endif
	} else{
#if DEBUG_LEVEL>1
		cout << "<- Done, could not improve center.\n";
#endif
	}
	
	return radius;
}

/// Determine (geometric) center of a group
inline Vec3 GroupCenter(Element_Group* group, Element* elements)
{
	Vec3 result(0.0);
	
	for(unsigned int i=0; i<group->nr_elements; i++){
		result+=elements[group->elements[i]].center;
	}
	result/=group->nr_elements;
	
	return result;
}

inline double MC_Elements::calc_lambda()
{
	double l=1.0;
	if(configuration->V>configuration->targetV){
		if(configuration->V<10.0*configuration->targetV){
			double sigma=9.0*configuration->targetV/sqrt(-2.0*log(EPS));
			l=exp(-0.5*qqpwr((configuration->V-configuration->targetV)/sigma,2));
		} else l=0.0;
	}
	return l;
}

inline double MC_Elements::calc_lambda(unsigned int step)
{
	double l=0.0;
	if(step>=configuration->transition_start){
		if(step<configuration->transition_start+configuration->transition_cycles){
			double sigma=(double)configuration->transition_cycles/sqrt(-2.0*log(EPS));
			l=1.0-exp(-0.5*qqpwr((step-configuration->transition_start)/sigma,2));
		} else l=1.0;
	}
	return l;
}

inline void CSPBC(const double bl, double &r)
{
	double bh=bl/2;
	if(r>bh){
		r-=bl;
		return;
	}
	if(r<-bh) r+=bl;
}

void MC_Elements::undo_PBCs(Vec3 &ref, Vec3 &rvec)
{
	Vec3 dv=rvec-ref;
	switch(configuration->latticetype){
		case 1:
		case 2: // Rectangular box (PBCs on arbitrary axes)
			// Apply PBC correction if PBCs are used along that axis
			if(configuration->PBCs[0]) CSPBC(configuration->boxlength[0],dv.vec[0]);
			if(configuration->PBCs[1]) CSPBC(configuration->boxlength[1],dv.vec[1]);
			if(configuration->PBCs[2]) CSPBC(configuration->boxlength[2],dv.vec[2]);
			break;
		case 3: break; // Spherical box (PBCs not yet supported)
		case 4: // Cylindrical box (optional PBCs on z-axis only)
			if(configuration->PBCs[2]) CSPBC(configuration->boxlength[2],dv.vec[2]);
			break;
		default: // Default condition
			cout << "Lattice type " << configuration->latticetype << " not recognized. Exiting.\n";
			exit(1);
	}
	rvec=ref+dv;
}

void MC_Elements::apply_PBCs(Vec3 &rvec)
{
	switch(configuration->latticetype){
		case 1:
		case 2: // Rectangular box (PBCs on arbitrary axes)
			// Apply PBC correction if PBCs are used along that axis
			if(configuration->PBCs[0]) CSPBC(configuration->boxlength[0],rvec.vec[0]);
			if(configuration->PBCs[1]) CSPBC(configuration->boxlength[1],rvec.vec[1]);
			if(configuration->PBCs[2]) CSPBC(configuration->boxlength[2],rvec.vec[2]);
			break;
		case 3: break; // Spherical box (PBCs not yet supported)
		case 4: // Cylindrical box (optional PBCs on z-axis only)
			if(configuration->PBCs[2]) CSPBC(configuration->boxlength[2],rvec.vec[2]);
			break;
		default: // Default condition
			cout << "Lattice type " << configuration->latticetype << " not recognized. Exiting.\n";
			exit(1);
	}
}

inline void MC_Elements::initialize_neighbors()
{
	if(configuration->cut){
		for(unsigned int i=0; i<nr_groups; i++){
			Groups[i]->nr_neighbors=0;
			Groups[i]->nr_neighbors_above=0;
		}
		for(unsigned int i=0; i<nr_individual_elements; i++){
			Elements[i+nr_group_elements].nr_neighbors=0;
			Elements[i+nr_group_elements].nr_neighbors_above=0;
		}
		find_neighbors();
	} else{
		unsigned int element_offset=nr_N_to_kk;
		unsigned int* i_nr_neighbors;
		unsigned int* j_nr_neighbors;
		unsigned int i_e=0;
		unsigned i_idx, j_idx;
		for(unsigned int i=0; i<nr_groups; i++){
			Groups[i]->nr_neighbors=0;
			Groups[i]->nr_neighbors_above=0;
		}
		for(unsigned int i=0; i<nr_individual_elements; i++){
			Elements[i+nr_group_elements].nr_neighbors=0;
			Elements[i+nr_group_elements].nr_neighbors_above=0;
		}
		for(unsigned int i=0; i<N-1; i++){ // go over all entities in box (groups and individual elements) - groups are first
			if(i<nr_groups){
				i_nr_neighbors=&(Groups[i]->nr_neighbors);
				Groups[i]->nr_neighbors_above=*i_nr_neighbors;
			} else{
				i_nr_neighbors=&(Elements[i+element_offset].nr_neighbors);
				Elements[i+element_offset].nr_neighbors_above=*i_nr_neighbors;
			}
			for(unsigned int j=i+1; j<N; j++){
				if(j<nr_groups) j_nr_neighbors=&(Groups[j]->nr_neighbors); else j_nr_neighbors=&(Elements[j+element_offset].nr_neighbors);
				i_idx=i*configuration->max_neighbors+*i_nr_neighbors;
				j_idx=j*configuration->max_neighbors+*j_nr_neighbors;
				neighbors[i_idx]=j;
				neighbors[j_idx]=i;
				energy_indices[i_idx]=i_e;
				energy_indices[j_idx]=i_e;
				pair_energies[i_e]=0.0;
				pair_energies[i_e+1]=0.0;
				(*i_nr_neighbors)++;
				(*j_nr_neighbors)++;
				i_e+=2;
			}
		}
	}
}

inline void MC_Elements::extend_neighbors()
{
	unsigned int new_nr_neighbors=configuration->max_neighbors+4;
	if(new_nr_neighbors>configuration->N-1) new_nr_neighbors=configuration->N-1;
	if(new_nr_neighbors>configuration->max_neighbors){
		unsigned int new_max_nr_neighbors=configuration->N*new_nr_neighbors;
		unsigned int size=new_max_nr_neighbors<<1;
		unsigned int* new_neighbor_store=new unsigned int[size];
		double* new_energy_store=new double[size];
		unsigned int* new_energy_index_store=new unsigned int[size];
		
		unsigned int* new_neighbors=&new_neighbor_store[(unsigned int)neighbor_switch*new_max_nr_neighbors];
		double* new_energies=&new_energy_store[(unsigned int)neighbor_switch*new_max_nr_neighbors];
		unsigned int* new_energy_indices=&new_energy_index_store[(unsigned int)neighbor_switch*new_max_nr_neighbors];
		
		// copy old content (not scaled to new limits)
		unsigned int idx, old_idx;
		idx=(unsigned int)(!neighbor_switch)*new_max_nr_neighbors;
		old_idx=(unsigned int)(!neighbor_switch)*max_nr_neighbors;
		memcpy(&new_neighbor_store[idx],&neighbor_store[old_idx],sizeof(unsigned int)*max_nr_neighbors);
		memcpy(&new_energy_store[idx],&pair_energy_store[old_idx],sizeof(double)*max_nr_neighbors);
		memcpy(&new_energy_index_store[idx],&energy_index_store[old_idx],sizeof(unsigned int)*max_nr_neighbors);
		
		// energy is linear and not adjusted to scaling lines
		memcpy(&new_energy_store[(unsigned int)(neighbor_switch)*new_max_nr_neighbors],&pair_energy_store[(unsigned int)(neighbor_switch)*max_nr_neighbors],sizeof(double)*max_nr_neighbors);
		
		// copy current content and adjust for new spacing limits
		for(unsigned int i=0; i<N; i++){
			idx=i*new_nr_neighbors;
			old_idx=i*configuration->max_neighbors;
			memcpy(&new_neighbors[idx],&neighbors[old_idx],sizeof(unsigned int)*configuration->max_neighbors);
			memcpy(&new_energy_indices[idx],&energy_indices[old_idx],sizeof(unsigned int)*configuration->max_neighbors);
		}
		delete[] neighbor_store;
		delete[] pair_energy_store;
		delete[] energy_index_store;
		
		max_nr_neighbors=new_max_nr_neighbors;
		configuration->max_neighbors=new_nr_neighbors;
		
		neighbor_store=new_neighbor_store;
		pair_energy_store=new_energy_store;
		energy_index_store=new_energy_index_store;
		
		neighbors=new_neighbors;
		pair_energies=new_energies;
		energy_indices=new_energy_indices;
#if DEBUG_LEVEL>2
		cout << "\nExpanded the neighborhood to include " << new_nr_neighbors << " neighbors.\n";
#endif
	}
}

inline void MC_Elements::find_neighbors()
{
	if(configuration->cut){
		unsigned int element_offset=nr_N_to_kk;
		Vec3 i_center, dist;
		unsigned int i_e=0;
		unsigned int* i_nr_neighbors;
		unsigned int* j_nr_neighbors;
		unsigned int old_max_neighbors=configuration->max_neighbors;
		unsigned int* old_neighbors=neighbors;
		double* old_energies=pair_energies;
		unsigned int* old_energy_indices=energy_indices;
		
		neighbor_switch=!neighbor_switch;
		neighbors=&neighbor_store[(unsigned int)(neighbor_switch)*max_nr_neighbors];
		pair_energies=&pair_energy_store[(unsigned int)(neighbor_switch)*max_nr_neighbors];
		energy_indices=&energy_index_store[(unsigned int)(neighbor_switch)*max_nr_neighbors];
		
		unsigned int old_neighbor_start, old_neighbor_idx;
		bool was_neighbor;
		unsigned i_idx, j_idx, o_idx;
		for(unsigned int i=0; i<nr_groups; i++){
			Groups[i]->old_nr_neighbors=Groups[i]->nr_neighbors;
			Groups[i]->nr_neighbors=0;
		}
		for(unsigned int i=0; i<nr_individual_elements; i++){
			Elements[i+nr_group_elements].old_nr_neighbors=Elements[i+nr_group_elements].nr_neighbors;
			Elements[i+nr_group_elements].nr_neighbors=0;
		}
		for(unsigned int i=0; i<N-1; i++){ // go over all entities in box (groups and individual elements) - groups are first
			if(i<nr_groups){
				i_nr_neighbors=&(Groups[i]->nr_neighbors);
				old_neighbor_idx=Groups[i]->nr_neighbors_above;
				Groups[i]->nr_neighbors_above=*i_nr_neighbors;
				i_center=Groups[i]->center;
			} else{
				i_nr_neighbors=&(Elements[i+element_offset].nr_neighbors);
				old_neighbor_idx=Elements[i+element_offset].nr_neighbors_above;
				Elements[i+element_offset].nr_neighbors_above=*i_nr_neighbors;
				i_center=Elements[i+element_offset].center;
			}
			old_neighbor_start=i*old_max_neighbors;
			for(unsigned int j=i+1; j<N; j++){
				if(j<nr_groups){
					was_neighbor=((old_neighbor_idx<Groups[i]->old_nr_neighbors) && (old_neighbors[old_neighbor_start+old_neighbor_idx]==j));
					j_nr_neighbors=&(Groups[j]->nr_neighbors);
					dist=i_center-Groups[j]->center;
				} else{
					was_neighbor=((old_neighbor_idx<Elements[i+element_offset].old_nr_neighbors) && (old_neighbors[old_neighbor_start+old_neighbor_idx]==j));
					j_nr_neighbors=&(Elements[j+element_offset].nr_neighbors);
					dist=i_center-Elements[j+element_offset].center;
				}
				apply_PBCs(dist);
				if((dist*dist)<configuration->cut2plus){
					if((*i_nr_neighbors>=configuration->max_neighbors) || (*j_nr_neighbors>=configuration->max_neighbors)){
						extend_neighbors();
						old_neighbors=&neighbor_store[(unsigned int)(!neighbor_switch)*max_nr_neighbors];
						old_energies=&pair_energy_store[(unsigned int)(!neighbor_switch)*max_nr_neighbors];
						old_energy_indices=&energy_index_store[(unsigned int)(!neighbor_switch)*max_nr_neighbors];
					}
					i_idx=i*configuration->max_neighbors+*i_nr_neighbors;
					j_idx=j*configuration->max_neighbors+*j_nr_neighbors;
					neighbors[i_idx]=j;
					neighbors[j_idx]=i;
					energy_indices[i_idx]=i_e;
					energy_indices[j_idx]=i_e;
					if(was_neighbor){ // still a neighbor, keep old energies
						o_idx=old_energy_indices[old_neighbor_start+old_neighbor_idx];
						pair_energies[i_e]=old_energies[o_idx];
						pair_energies[i_e+1]=old_energies[o_idx+1];
						old_neighbor_idx++;
					} else{ // new neighbor needs a new (zero) block of energies
						pair_energies[i_e]=0.0;
						pair_energies[i_e+1]=0.0;
					}
					(*i_nr_neighbors)++;
					(*j_nr_neighbors)++;
					i_e+=2;
				} else{
					if(was_neighbor) old_neighbor_idx++; // was a neighbor, but is not anymore
				}
			}
		}
	}
}

inline bool MC_Elements::get_distances(const unsigned int i)
{
	Rmus[i].vec[0]=0.0; Rmus[i].vec[1]=0.0; Rmus[i].vec[2]=0.0;
	Vec3 dist2center;
	Vec3 temp, gt;
	Vec3 element2group(0.0);
	unsigned int neighbor, j_curr, nr_neighbors, start;
	double t2;
	if(Elements[i].group){
		nr_neighbors=Elements[i].group->nr_neighbors;
		start=Elements[i].group->number*configuration->max_neighbors;
		dist2center=Elements[i].group->center;
		element2group=Elements[i].group->center-Elements[i].center;
	} else{
		nr_neighbors=Elements[i].nr_neighbors;
		start=(nr_groups+i-nr_group_elements)*configuration->max_neighbors;
		dist2center=Elements[i].center;
	}
	bool LJwallES=(configuration->LJwall_calc) && (fabs(dist2center.vec[0]-configuration->LJwall_xm)<configuration->rfcut);
	for(unsigned int j=start; j<start+nr_neighbors; ++j){
		neighbor=neighbors[j];
		if(neighbor<nr_groups){ // the neighbor is a group
			Element_Group* group=Groups[neighbor];
			temp=group->center-dist2center;
			apply_PBCs(temp);
			t2=temp.vec[0]*temp.vec[0]+temp.vec[1]*temp.vec[1]+temp.vec[2]*temp.vec[2];
			if(t2 < EPS) return false;
			gt=temp+element2group-group->center;
			for(unsigned int k=0; k<group->nr_elements; k++){
				j_curr=group->elements[k];
				Rmus[j_curr]=Elements[j_curr].center+gt;
				Rdist2[j_curr]=t2;
				Rdist2LJ[j_curr]=Rmus[j_curr].vec[0]*Rmus[j_curr].vec[0]+Rmus[j_curr].vec[1]*Rmus[j_curr].vec[1]+Rmus[j_curr].vec[2]*Rmus[j_curr].vec[2];
				LJwallMirror[j_curr]=false;
				if(LJwallES){
					if(Elements[j_curr].center.vec[0]<0.0){ // Rmirror is inverted
						Rmirror[j_curr]=2.0*configuration->LJwall_xm+Elements[j_curr].center.vec[0]-gt.vec[0];
					} else{
						Rmirror[j_curr]=-2.0*configuration->LJwall_xm+Elements[j_curr].center.vec[0]-gt.vec[0];
					}
					// note that this only works because both r_i_x and r_j_x are negative towards the negative wall:
					// +xm wall: 4*(r_i_x-xm)(r_j_x-xm),
					// -xm wall: 4*(r_i_x+xm)(r_j_m+xm) = 4*(-r_i_x-xm)(-r_j_x-xm) = 4*(|r_i_x|-xm)(|r_j_x|-xm
					LJwallMirror[j_curr]=(4.0*(fabs(dist2center.vec[0])-configuration->LJwall_xm)*(fabs(group->center.vec[0])-configuration->LJwall_xm)+t2<configuration->escut2);
				}
			}
		} else{ // individual element
			j_curr=neighbor+nr_N_to_kk;
			Element* element=&Elements[j_curr];
			temp=element->center-dist2center;
			apply_PBCs(temp);
			Rmus[j_curr]=temp+element2group;
			t2=temp.vec[0]*temp.vec[0]+temp.vec[1]*temp.vec[1]+temp.vec[2]*temp.vec[2];
			Rdist2[j_curr]=t2;
			Rdist2LJ[j_curr]=Rmus[j_curr].vec[0]*Rmus[j_curr].vec[0]+Rmus[j_curr].vec[1]*Rmus[j_curr].vec[1]+Rmus[j_curr].vec[2]*Rmus[j_curr].vec[2];
			// see comment above
			LJwallMirror[j_curr]=false;
			if(LJwallES){
				if(element->center.vec[0]<0.0){ // Rmirror is inverted
					Rmirror[j_curr]=2.0*configuration->LJwall_xm+element->center.vec[0]-temp.vec[0]-element2group.vec[0];
				} else{
					Rmirror[j_curr]=-2.0*configuration->LJwall_xm+element->center.vec[0]-temp.vec[0]-element2group.vec[0];
				}
				LJwallMirror[j_curr]=(4.0*(fabs(dist2center.vec[0])-configuration->LJwall_xm)*(fabs(element->center.vec[0])-configuration->LJwall_xm)+t2<configuration->escut2);
			}
			//Bomb out if a zero distance is detected to avoid dividing by zero
			if(Rdist2[j_curr] < EPS) return false;
		}
	}
	return true;
}

inline bool MC_Elements::get_distances_above_index(const unsigned int i)
{
	Rmus[i].vec[0]=0.0; Rmus[i].vec[1]=0.0; Rmus[i].vec[2]=0.0;
	Vec3 dist2center;
	Vec3 temp, gt;
	Vec3 element2group(0.0);
	unsigned int neighbor, j_curr, nr_neighbors, start, start_above;
	double t2;
	if(Elements[i].group){
		nr_neighbors=Elements[i].group->nr_neighbors;
		start_above=Elements[i].group->nr_neighbors_above;
		start=Elements[i].group->number*configuration->max_neighbors;
		dist2center=Elements[i].group->center;
		element2group=Elements[i].group->center-Elements[i].center;
	} else{
		nr_neighbors=Elements[i].nr_neighbors;
		start_above=Elements[i].nr_neighbors_above;
		start=(nr_groups+i-nr_group_elements)*configuration->max_neighbors;
		dist2center=Elements[i].center;
	}
	bool LJwallES=(configuration->LJwall_calc) && (fabs(dist2center.vec[0]-configuration->LJwall_xm)<configuration->rfcut);
	for(unsigned int j=start_above+start; j<start+nr_neighbors; ++j){
		neighbor=neighbors[j];
		if(neighbor<nr_groups){ // the neighbor is a group
			Element_Group* group=Groups[neighbor];
			temp=group->center-dist2center;
			apply_PBCs(temp);
			t2=temp.vec[0]*temp.vec[0]+temp.vec[1]*temp.vec[1]+temp.vec[2]*temp.vec[2];
			if(t2 < EPS) return false;
			gt=temp+element2group-group->center;
			for(unsigned int k=0; k<group->nr_elements; k++){
				j_curr=group->elements[k];
				Rmus[j_curr]=Elements[j_curr].center+gt;
				Rdist2[j_curr]=t2;
				Rdist2LJ[j_curr]=Rmus[j_curr].vec[0]*Rmus[j_curr].vec[0]+Rmus[j_curr].vec[1]*Rmus[j_curr].vec[1]+Rmus[j_curr].vec[2]*Rmus[j_curr].vec[2];
				LJwallMirror[j_curr]=false;
				if(LJwallES){
					if(Elements[j_curr].center.vec[0]<0.0){ // Rmirror is inverted
						Rmirror[j_curr]=2.0*configuration->LJwall_xm+Elements[j_curr].center.vec[0]-gt.vec[0];
					} else{
						Rmirror[j_curr]=-2.0*configuration->LJwall_xm+Elements[j_curr].center.vec[0]-gt.vec[0];
					}
					// note that this only works because both r_i_x and r_j_x are negative towards the negative wall:
					// +xm wall: 4*(r_i_x-xm)(r_j_x-xm),
					// -xm wall: 4*(r_i_x+xm)(r_j_m+xm) = 4*(-r_i_x-xm)(-r_j_x-xm) = 4*(|r_i_x|-xm)(|r_j_x|-xm
					LJwallMirror[j_curr]=(4.0*(fabs(dist2center.vec[0])-configuration->LJwall_xm)*(fabs(group->center.vec[0])-configuration->LJwall_xm)+t2<configuration->escut2);
				}
			}
		} else{ // individual element
			j_curr=neighbor+nr_N_to_kk;
			Element* element=&Elements[j_curr];
			temp=element->center-dist2center;
			apply_PBCs(temp);
			Rmus[j_curr]=temp+element2group;
			t2=temp.vec[0]*temp.vec[0]+temp.vec[1]*temp.vec[1]+temp.vec[2]*temp.vec[2];
			Rdist2[j_curr]=t2;
			Rdist2LJ[j_curr]=Rmus[j_curr].vec[0]*Rmus[j_curr].vec[0]+Rmus[j_curr].vec[1]*Rmus[j_curr].vec[1]+Rmus[j_curr].vec[2]*Rmus[j_curr].vec[2];
			LJwallMirror[j_curr]=false;
			if(LJwallES){
				if(element->center.vec[0]<0.0){ // Rmirror is inverted
					Rmirror[j_curr]=2.0*configuration->LJwall_xm+element->center.vec[0]-temp.vec[0]-element2group.vec[0];
				} else{
					Rmirror[j_curr]=-2.0*configuration->LJwall_xm+element->center.vec[0]-temp.vec[0]-element2group.vec[0];
				}
				// see comment above
				LJwallMirror[j_curr]=(4.0*(fabs(dist2center.vec[0])-configuration->LJwall_xm)*(fabs(element->center.vec[0])-configuration->LJwall_xm)+t2<configuration->escut2);
			}
			//Bomb out if a zero distance is detected to avoid dividing by zero
			if(Rdist2[j_curr] < EPS) return false;
		}
	}
	return true;
}

inline void MC_Elements::update_rf_correction()
{
	rf_correction = 2.0*(epsRF-configuration->n2)/(pwr3(configuration->rfcut)*(2.0*epsRF+configuration->n2));
	rf_neutralization = 3.0*epsRF/(configuration->rfcut*(2.0*epsRF+configuration->n2));
}

/*!
 * SCRF Update; called if configuration->dyneps = true
 * Updates reaction field dielectric constant. Convergence may be difficult for poorly-behaved systems like strongly polar Stockmayer fluids.
 */
inline void MC_Elements::update_eps(Vec3 &M, unsigned int step)
{
	double q;
	double n2 = configuration->n2;
	double n4 = configuration->n2*configuration->n2;
	
	// Variance calculation if no field
	if(configuration->noEfield){
		double M2 = M*M; // fallback in case no statistics is accumulated yet
		if((step>=configuration->dyneps_average_cycles) && (configuration->dyneps_average_cycles>1)){ // calculate averages or variance
			M2=0.0;
			Vec3 avgM(0.0);
			for(unsigned int i=step-configuration->dyneps_average_cycles; i<step; i++){
				M2+=Ms[i]*Ms[i];
				avgM+=Ms[i];
			}
			M2/=configuration->dyneps_average_cycles;
			if(configuration->dyneps_varm){
				avgM/=configuration->dyneps_average_cycles;
				M2-=avgM*avgM;
			}
		}
		// Calculate q = yg
		q = 4*pi*M2/(9*configuration->V*configuration->kT);
		// eps_RF = 1/4*(n^2+9*q +/- 3*sqrt(n^4+2*n^2*q+9*q^2))
		// - out of +/- sqrt solution only + is solution b/c for q->0 : eps_RF=1/4*(n^2+3*sqrt(n^4)) = 1/4*(4*n^2)
		epsRF = 0.25*n2+2.25*q+0.75*sqrt(n4+2.0*n2*q+9.0*q*q);
	} else{ // Calculation if E0z is present
		Vec3 avgM=M;
		if((step>=configuration->dyneps_average_cycles) && (configuration->dyneps_average_cycles>1)){ // calculate average M
			avgM.V3Zeros();
			for(unsigned int i=step-configuration->dyneps_average_cycles; i<step; i++) avgM+=Ms[i];
			avgM/=configuration->dyneps_average_cycles;
		}
		// Calculate epsilon
		if(configuration->realfield){
			epsRF = 4.0*pi*(avgM*Evec)/((Evec*Evec)*configuration->V) + n2;
		} else{
			q = (4.0*pi*(avgM*Evec))/(3.0*configuration->V*(Evec*Evec));
			epsRF = 0.25*n2+2.25*q+0.75*sqrt(n4+2.0*n2*q+9.0*q*q);
		}
	}
}

inline double MC_Elements::phi_q_mu(const unsigned int k, Vec3 &thermu, double &RFp, bool updateRF)
{
	Vec3 r=thermu; // point from charge on i to center of element k
	if(configuration->offctrmu) r+=(Elements[k].rot*Elements[k].MyType->mu_pos); // correction in case of offcenter dipoles
	double rmu=Elements[k].dipole*r;
	double ir2=1.0/(r*r);
	double ir=sqrt(ir2);
	// Add interacting charge to reaction field
	if(updateRF) RFp += rmu;
	// if dipole points in direction of charge to dipole vector then negative dipole charge is closer to charge q ...
	return -rmu*(ir2*ir); // attractive potentials are negative, which this one should be with a positive charge q
}

inline double MC_Elements::phi_q(const unsigned int k, Vec3 &thermu, double &RFp, bool updateRF)
{
	double result=0.0;
	RFp=0.0;
	if(Elements[k].MyType->hasmu) result+=phi_q_mu(k,thermu,RFp,updateRF);
	if(Elements[k].MyType->nr_charges>0){
		Vec3 r;
		double q, r2;
		for(unsigned int j=0; j<Elements[k].MyType->nr_charges; j++){
			// thermu points from m-th charge on i to center of k
			r=thermu+Elements[k].q_pos[j];
			r2=r*r;
			q=Elements[k].MyType->q[j];
			// Add interacting charge to reaction field
#ifndef FRIEDMAN_IMAGE
			if(updateRF) RFp += 0.5*q*r2;
#else
			if(updateRF) RFp += q*(1.0/epsRF+1.0/configuration->n2); // modified Friedman image method
#endif
			result+=q/sqrt(r2); // negative (attractive) for opposite charges
		}
	}
	return result;
}

inline Vec3 MC_Elements::E_mu_q(const unsigned int k, Vec3 &thermu, Vec3 &RF, bool updateRF)
{
	Vec3 result(0.0);
	Vec3 r;
	// element k is the one with charges
	// rmu points from dipole on i to center of k
	double q, ir2, ir3;
	for(unsigned int j=0; j<Elements[k].MyType->nr_charges; j++){
		// thermu points from dipole on i to center of k
		r=thermu+Elements[k].q_pos[j]; // vector from dipole on i to charge on k ...
		ir2=1.0/(r*r);
		ir3=sqrt(ir2)*ir2;
		q=Elements[k].MyType->q[j];
		// Add interacting charge to reaction field
		if(updateRF) RF += r*q;
		// if dipole points in direction of dipole to charge vector then positive dipole charge is closer to charge q ...
		result+=r*(q*ir3); // potential should be negative (attractive) when charge q is negative
	}
	return result;
}

inline Vec3 MC_Elements::E_mu(const unsigned int k, Vec3 &thermu, Vec3 &RF, bool updateRF)
{
	Vec3 result(0.0);
	RF=Vec3(0.0);
	if(Elements[k].MyType->hasmu){
		Vec3 r = thermu;
		if(configuration->offctrmu) r+=(Elements[k].rot*Elements[k].MyType->mu_pos); // correction in case of offcenter dipoles
		double ir2=1.0/(r*r);
		double ir3=sqrt(ir2)*ir2;
		
		Vec3* dipole=&Elements[k].dipole;
		
		// Add interacting dipole to reaction field
		if(updateRF) RF = *dipole;
		
		// Calculate the dipole-dipole energy (in picoerg)
		result = (*dipole-r*(3.0*ir2*((*dipole)*r)))*ir3; // will be negative (attractive) if both dipoles point in same direction congruent to distance vector
	}
	if(Elements[k].MyType->nr_charges>0) result+=E_mu_q(k,thermu,RF,updateRF);
	return result;
}

inline double MC_Elements::self_phi_q_mu(const unsigned int k, Vec3 &thermu, double &RFp, bool updateRF)
{
	Vec3 r=thermu; // point from charge on i to center of element k
	if(configuration->offctrmu) r+=(Elements[k].rot*Elements[k].MyType->mu_pos); // correction in case of offcenter dipoles
	double rmu=Elements[k].dipole*r;
	double ir2=1.0/(r*r);
	double ir=sqrt(ir2);
	// Add interacting charge to reaction field
	if(updateRF) RFp += rmu;
	// if dipole points in direction of charge to dipole vector then negative dipole charge is closer to charge q ...
	return -rmu*(ir2*ir); // attractive potentials are negative, which this one should be with a positive charge q
}

inline double MC_Elements::self_phi_q(const unsigned int k, Vec3 &thermu, double &RFp, bool updateRF)
{
	double result=0.0;
	if(Elements[k].MyType->hasmu) result+=self_phi_q_mu(k,thermu,RFp,updateRF);
	if(Elements[k].MyType->nr_charges>0){
		Vec3 r;
		double q;
		double r2;
		for(unsigned int j=0; j<Elements[k].MyType->nr_charges; j++){
			// thermu points from m-th charge on i to center of k
			r=thermu+Elements[k].q_pos[j];
			r2=r*r;
			q=Elements[k].MyType->q[j];
			// Add interacting charge to reaction field
			if(updateRF) RFp += 0.5*q*r2;
			result+=q/sqrt(r2); // negative (attractive) for opposite charges
		}
	}
	return result;
}

inline Vec3 MC_Elements::self_E_mu_q(const unsigned int k, Vec3 &thermu, Vec3 &RF, bool updateRF)
{
	Vec3 result(0.0);
	Vec3 r;
	// element k is the one with charges
	// rmu points from dipole on i to center of k
	double q, ir2, ir;
	for(unsigned int j=0; j<Elements[k].MyType->nr_charges; j++){
		// thermu points from dipole on i to center of k
		r=thermu+Elements[k].q_pos[j]; // vector from dipole on i to charge on k ...
		ir2=1.0/(r*r);
		ir=sqrt(ir2);
		q=Elements[k].MyType->q[j];
		// Add interacting charge to reaction field
		if(updateRF) RF += r*q;
		// if dipole points in direction of dipole to charge vector then positive dipole charge is closer to charge q ...
		result+=r*(q*ir2*ir); // potential should be negative (attractive) when charge q is negative
	}
	return result;
}

inline Vec3 MC_Elements::self_E_mu(const unsigned int k, Vec3 &thermu, Vec3 &RF, bool updateRF)
{
	Vec3 result(0.0);
	if(Elements[k].MyType->hasmu){
		Vec3 r = thermu;
		if(configuration->offctrmu) r+=(Elements[k].rot*Elements[k].MyType->mu_pos); // correction in case of offcenter dipoles
		double ir2=1.0/(r*r);
		double ir3=sqrt(ir2)*ir2;
		Vec3* dipole=&Elements[k].dipole;
		
		// Add interacting dipole to reaction field
		if(updateRF) RF += *dipole;
		
		// Calculate the dipole-dipole energy (in picoerg)
		result = (*dipole-r*(3.0*ir2*((*dipole)*r)))*ir3; // will be negative (attractive) if both dipoles point in same direction congruent to distance vector
	}
	if(Elements[k].MyType->nr_charges>0) result+=self_E_mu_q(k,thermu,RF,updateRF);
	return result;
}

/*!
 * Calculation of total moments and first three moments of cos(theta) wrt the z axis. LEJ 12/30/2008
 * Inputs:
 * 	BeOids - the array of all molecules in the system.
 * 	means - storage array for order parameters
 * 	M - storage array for total dipole moment
 */
inline void MC_Elements::get_means(Vec3 &means, Vec3 &M)
{
	double zproj, zp2;
	means.vec[0]=0.0; means.vec[1]=0.0; means.vec[2]=0.0;
	M.vec[0]=0.0; M.vec[1]=0.0; M.vec[2]=0.0;
	Vec3 dipole;
	
	unsigned int n_avg=0;
	Element* element;
	
	// loop over elements which are not in groups
	for(unsigned int i = nr_group_elements; i < nr_elements; i++){
		element=&Elements[i];
		dipole=element->dipole;
		for(unsigned int j=0; j<element->MyType->nr_charges; j++){
			Vec3 r=element->center+element->q_pos[j];
			dipole += r*element->MyType->q[j];
		}
		if(element->MyType->calculate_order){
			n_avg++;
			// get the moments of cos relative to the Z axis [0 0 1]. Use of any
			// other axis requires the extended means function.
			zproj = dipole.vec[2]/dipole.V3Norm(); // dipole z-component/norm(dipole)
			zp2 = zproj*zproj;
			means.vec[0] +=zproj;
			means.vec[1] +=zp2;
			means.vec[2] +=zproj*zp2;
		}
		//Add up dipoles and charges to obtain M
		M += dipole;
	}
	
	// loop over groups
	Element_Group* group;
	for(unsigned int i = 0; i < nr_groups; i++){
		group=Groups[i];
		dipole=Vec3(0.0);
		bool component_dipole=(group->Type->LOD!=NULL);
		unsigned int component=0;
		if(component_dipole) component_dipole&=(group->Type->LOD->order_dipole_component>0);
		for(unsigned int j=0; j<group->nr_elements; j++){
			element=&Elements[group->elements[j]];
			dipole+=element->dipole;
			for(unsigned int k=0; k<element->MyType->nr_charges; k++){
				Vec3 r=element->center+element->q_pos[k];
				if(group->Type->group_dipole){
					r-=group->center;
				}
				dipole += r*element->MyType->q[k];
			}
		}
		Vec3 order_dipole;
		if(component_dipole && group->Type->calculate_order){ // only calculate this if we'll use it ...
			component=group->Type->LOD->order_dipole_component-1;
			Vec3 component_center(0.0);
			double qsum=0.0;
			for(unsigned int j=0; j<group->Type->LOD->nr_elements[component]; j++){
				if(group->levelofdetail>0){
					unsigned int component_element=group->Type->LOD->component_elements[component][j];
					unsigned int ellipsoid_idx=group->Type->LOD->element_in_ellipsoid[group->levelofdetail-1][component_element];
					element=&Elements[group->elements[ellipsoid_idx]];
					Element* thiselement=&configuration->group_elements[configuration->groups[group->type]->elements[ellipsoid_idx]];
					Element* fullelement=&configuration->group_elements[configuration->groups[group->type]->Type->LOD->groups[0]->elements[component_element]];
					Mat33 delta_rot=element->rot*thiselement->rot.M3Transpose();
					Vec3 offset=delta_rot*(fullelement->center-thiselement->center);
					Mat33 element_rot=delta_rot*fullelement->rot;
					
					Vec3 position=element->center+offset;
					component_center+=position;
					order_dipole+=element_rot*fullelement->MyType->initial_dipole;
					for(unsigned int k=0; k<fullelement->MyType->nr_charges; k++){
						double q=fullelement->MyType->q[k];
						qsum+=q;
						order_dipole+=(position+element_rot*fullelement->MyType->q_pos[k])*q;
					}
				} else{
					element=&Elements[group->Type->LOD->component_elements[component][j]];
					component_center+=element->center;
					order_dipole+=element->dipole;
					for(unsigned int k=0; k<element->MyType->nr_charges; k++){
						Vec3 r=element->center+element->q_pos[k];
						order_dipole += r*element->MyType->q[k];
						qsum+=element->MyType->q[k];
					}
				}
			}
			component_center /= group->Type->LOD->nr_elements[component];
			order_dipole-=component_center*qsum;
		} else order_dipole=dipole;
		if(group->Type->calculate_order){
			n_avg++;
			// get the moments of cos relative to the Z axis [0 0 1]. Use of any
			// other axis requires the extended means function.
			zproj = order_dipole.vec[2]/order_dipole.V3Norm(); // dipole z-component/norm(dipole)
			zp2 = zproj*zproj;
			means.vec[0] +=zproj;
			means.vec[1] +=zp2;
			means.vec[2] +=zproj*zp2;
		}
		M += dipole;
	}
	
	//Calculate average
	if(n_avg>0) means /= n_avg;
}

inline double MC_Elements::Vwall(const unsigned int i, bool skip_charges, bool &success)
{
	// take care of Lennard-Jones contributions first
	double Rx;
	unsigned int ikk=configuration->num_element_types*configuration->num_element_types+2*Elements[i].mytype;
	double epsilon = configuration->pre_eps[ikk]*configuration->energy_scale[3];
	// wall normal is (1,0,0), multiplied by inverse element rotation matrix for min. distance normal, then squared elements to have (alpha^2,beta^2,gamma^2)
	Vec3 n2(Elements[i].rot.mat[0][0]*Elements[i].rot.mat[0][0],Elements[i].rot.mat[0][1]*Elements[i].rot.mat[0][1],Elements[i].rot.mat[0][2]*Elements[i].rot.mat[0][2]);
	double ns2=n2*Elements[i].MyType->saxes2; // alpha^2*a^2+beta^2*b^2+gamma^2*c^2
	double sigma=ns2/sqrt((n2.vec[0]+n2.vec[1]+n2.vec[2])*ns2); // minimum distance between ellipsoid and wall
	if(configuration->LJ_adjust_width){
		double wi=Elements[i].MyType->avg_width;
		if(wi<EPS){
			// Solve EllipsoidRmin explicitly (since direction is alway (1,0,0) in Lab frame)
			Vec3 dist(Elements[i].rot.mat[0][0]/Elements[i].MyType->saxes.vec[0],Elements[i].rot.mat[0][1]/Elements[i].MyType->saxes.vec[1],Elements[i].rot.mat[0][2]/Elements[i].MyType->saxes.vec[2]);
			// Now solve x^2/a^2+y^2/b^2+z^2/c^2=1/r^2 => r = sqrt(1/(x^2/a^2+y^2/b^2+z^2/c^2))
			wi=1.0/sqrt(dist*dist);
		}
		double s=wi+configuration->LJwall_width;
		double r=configuration->LJwall_xm-fabs(Elements[i].center.vec[0])+s-sigma; // denominator is Xm - (|x|-(width-sigma))
		if(r<0.0){
			r=0.0;
			success=false;
		}
		Rx=s/r;
	} else{
		Rx=sigma/(configuration->LJwall_xm-fabs(Elements[i].center.vec[0]));
		success&=((configuration->LJwall_xm-fabs(Elements[i].center.vec[0]))>=0.0); // treat wall as a hard wall
	}
	Rx*=Rx;
	if(configuration->LJ_interaction_area) epsilon*=IA(Elements[i].MyType,Vec3(Elements[i].rot.mat[0][0],Elements[i].rot.mat[0][1],Elements[i].rot.mat[0][2]));
	double disp = qqpwr(Rx,configuration->LJexp[0]>>1); // Since Rx is already squared, the small term must be halved (bitshift right by one is integer/2 - AT)
	double V=epsilon*(configuration->r*disp*(disp-LJlambda*configuration->pre_eps[ikk+1]));
	// take care of electrostatics (if within reach, interaction with mirror image and RF)
	if(!skip_charges){
		// make sure the mirror wall never overlaps with the ellipsoid
		double xm_eff = configuration->LJwall_xm;
		if(fabs(Elements[i].center.vec[0])+sigma>configuration->LJwall_xm) xm_eff = fabs(Elements[i].center.vec[0])+sigma;
		if(2.0*(xm_eff-fabs(Elements[i].center.vec[0]))<configuration->rfcut){
			double rmirror;
			if(Elements[i].center.vec[0]<0.0){ // rmu is inverted
				rmirror=2.0*(xm_eff+Elements[i].center.vec[0]);
			} else{
				rmirror=-2.0*(xm_eff-Elements[i].center.vec[0]);
			}
			double Ves=0.0;
			Vec3 rmu;
			Vec3 Emu, RF;
			double phiq, RFp;
			double i_charge=Elements[i].MyType->charge*rf_neutralization;
			bool updateRF=(fabs(rf_correction)>EPS);
			unsigned int nr_charges=Elements[i].MyType->nr_charges;
			if(Elements[i].MyType->hasmu){
				rmu=Elements[i].rot*Elements[i].MyType->mu_pos;
				rmu.vec[0]+=rmirror;
				rmu.vec[1]*=-1.0; rmu.vec[2]*=-1.0;
				Emu=E_mu(i,rmu,RF,updateRF); // rmu.x is inversed
				Emu.vec[0]*=-1.0;
				RF.vec[0]*=-1.0;
				Emu*=configuration->LJwall_mirror_factor;
				RF*=configuration->LJwall_mirror_factor;
				Ves += (Emu-RF*rf_correction)*Elements[i].dipole;
			}
			if(nr_charges){
				for(unsigned int m=0; m<nr_charges; m++){
					Ves-= i_charge*configuration->LJwall_mirror_factor*Elements[i].MyType->q[m];
					rmu.vec[0]=rmirror+Elements[i].q_pos[m].vec[0];
					rmu.vec[1]=-Elements[i].q_pos[m].vec[1]; rmu.vec[2]=-Elements[i].q_pos[m].vec[2];
					phiq=phi_q(i,rmu,RFp,updateRF)*configuration->LJwall_mirror_factor; // rmu.x is inversed
					RFp*=configuration->LJwall_mirror_factor;
					Ves += (phiq+RFp*rf_correction)*Elements[i].MyType->q[m];
				}
			}
			V+=Ves*configuration->in2;
		}
	}
	return V;
}

inline double MC_Elements::VmuE(const unsigned int i)
{
	//Calculate the muE energy
	Vec3 dipole=Elements[i].dipole;
	for(unsigned int j=0; j<Elements[i].MyType->nr_charges; j++) dipole += (Elements[i].center+Elements[i].q_pos[j])*Elements[i].MyType->q[j];
	double mE=dipole*Evec;
	if(configuration->VmuE_squared){
		double norm=dipole*dipole;
		norm*=Evec*Evec;
		mE*=mE*sqrt(1.0/norm);
	}
	return configuration->energy_scale[0]*mE; // energy_scale[0] contains factor -1
}

/*!
 * Touch, an algorithm designed for determining the closest contact distance between two ellipsoids at arbitrary angles,
 * reducing the three-dimensional problem to a one-dimensional problem in terms of a scalar variable lambda, then iteratively
 * solving a system of linear equations to optimize lambda. Initially written in Matlab by BHR, ported to C++ by RSB in 02/09, then
 * harmonized with current codebase by LEJ in 04/09, and optimized by LEJ in 01/10, finally fixed by AT in 02/11, optimized by AT in 2014
 * Inputs:
 * 	kk - the index of the currently active molecule
 * 	i - the index of the molecule it is interacting with
 * Return value in distance is the effective LJ sigma for the molecules at their current positions
 */
inline double MC_Elements::touch(const unsigned int kk, const unsigned int i)
{
	double matrices[15];
	// Minimization vectors
	Vec3 V, t, CI;
	
	// Create lab frame version of A and B matrices (both are symmetric matrices R * L_A/B * R^T)
	// A/B_ij=sum_k L_A/B_k*R_ik*R_jk
	// move into A ellipsoid frame of reference
	// doing so:
	// - replaces 36 multiplication and 12 additions with 36 multiplications and 24 additions once
	// - saves 4 multiplications and 6 additions in loop
	Vec3 z; // multply with kk's transposed (inverse) rotation matrix
	z.vec[0] = Rmus[i].vec[0]*Elements[kk].rot.mat[0][0]+Rmus[i].vec[1]*Elements[kk].rot.mat[1][0]+Rmus[i].vec[2]*Elements[kk].rot.mat[2][0];
	z.vec[1] = Rmus[i].vec[0]*Elements[kk].rot.mat[0][1]+Rmus[i].vec[1]*Elements[kk].rot.mat[1][1]+Rmus[i].vec[2]*Elements[kk].rot.mat[2][1];
	z.vec[2] = Rmus[i].vec[0]*Elements[kk].rot.mat[0][2]+Rmus[i].vec[1]*Elements[kk].rot.mat[1][2]+Rmus[i].vec[2]*Elements[kk].rot.mat[2][2];
	Mat33 irot=Elements[kk].rot.TransMulM3(Elements[i].rot); // new rotation matrix of i in kk's frame
	
	// A_LF is just the diagonal since we are in A's frame
	matrices[0]=Elements[kk].MyType->saxes2.vec[0];
	matrices[1]=Elements[kk].MyType->saxes2.vec[1];
	matrices[2]=Elements[kk].MyType->saxes2.vec[2];
	
	// B_LF upper diagonal elements
	matrices[3]=Elements[i].MyType->saxes2.vec[0]*irot.mat[0][0]*irot.mat[0][0]+Elements[i].MyType->saxes2.vec[1]*irot.mat[0][1]*irot.mat[0][1]+Elements[i].MyType->saxes2.vec[2]*irot.mat[0][2]*irot.mat[0][2];
	matrices[4]=Elements[i].MyType->saxes2.vec[0]*irot.mat[0][0]*irot.mat[1][0]+Elements[i].MyType->saxes2.vec[1]*irot.mat[0][1]*irot.mat[1][1]+Elements[i].MyType->saxes2.vec[2]*irot.mat[0][2]*irot.mat[1][2];
	matrices[5]=Elements[i].MyType->saxes2.vec[0]*irot.mat[0][0]*irot.mat[2][0]+Elements[i].MyType->saxes2.vec[1]*irot.mat[0][1]*irot.mat[2][1]+Elements[i].MyType->saxes2.vec[2]*irot.mat[0][2]*irot.mat[2][2];
	matrices[6]=Elements[i].MyType->saxes2.vec[0]*irot.mat[1][0]*irot.mat[1][0]+Elements[i].MyType->saxes2.vec[1]*irot.mat[1][1]*irot.mat[1][1]+Elements[i].MyType->saxes2.vec[2]*irot.mat[1][2]*irot.mat[1][2];
	matrices[7]=Elements[i].MyType->saxes2.vec[0]*irot.mat[1][0]*irot.mat[2][0]+Elements[i].MyType->saxes2.vec[1]*irot.mat[1][1]*irot.mat[2][1]+Elements[i].MyType->saxes2.vec[2]*irot.mat[1][2]*irot.mat[2][2];
	matrices[8]=Elements[i].MyType->saxes2.vec[0]*irot.mat[2][0]*irot.mat[2][0]+Elements[i].MyType->saxes2.vec[1]*irot.mat[2][1]*irot.mat[2][1]+Elements[i].MyType->saxes2.vec[2]*irot.mat[2][2]*irot.mat[2][2];
	
	if(configuration->vdwtype==6){ // modulate touch
		matrices[0]*=Elements[kk].delta;
		matrices[1]*=Elements[kk].delta;
		matrices[2]*=Elements[kk].delta;
		
		matrices[3]*=Elements[i].delta;
		matrices[4]*=Elements[i].delta;
		matrices[5]*=Elements[i].delta;
		matrices[6]*=Elements[i].delta;
		matrices[7]*=Elements[i].delta;
		matrices[8]*=Elements[i].delta;
	}
	// d = 1/det(B_LF) (uses the fact that rotations do not change the volume of another matrix)
	double d=Rdist2LJ[i]*Elements[i].MyType->invsaxesvolume;
	double e=d*d;
	if(e>1E10){
		d=1E5;
		e=1E10;
	}
	if(e>Rdist2LJ[i]){
		matrices[0]*=e;
		matrices[1]*=e;
		matrices[2]*=e;
		
		matrices[3]*=e;
		matrices[4]*=e;
		matrices[5]*=e;
		matrices[6]*=e;
		matrices[7]*=e;
		matrices[8]*=e;
		z*=d;
	}
	
	// save some calculation time in the loop (in order of calculation)
	matrices[9]=matrices[4]*matrices[5];
	matrices[10]=matrices[5]*matrices[7];
	matrices[11]=matrices[4]*matrices[7];
	matrices[12]=matrices[7]*matrices[7];
	matrices[13]=matrices[5]*matrices[5];
	matrices[14]=matrices[4]*matrices[4];
	
	double lambda=0.5;
	double xlx=1.0; // x=lambda/(1-lambda) -> do manual calculation with value above
	double xlx2=1.0;
	double Var = 1.0; // trying different forms found 6 loops gives 7 figs for Fx
	double VAV, VBV, det, Vz;
	
	while(Var > 1E-6){ // loop until variance in distance is sufficiently small
		Var = lambda; // keep lambda around but do CI matrix in terms of x (scale CI by 1/(1-lambda))
		//Populate CI matrix elements we need
		CI.vec[0] = xlx*matrices[3] + matrices[0];
		CI.vec[1] = xlx*matrices[6] + matrices[1];
		CI.vec[2] = xlx*matrices[8] + matrices[2];
		
		// Solve z = CI*V for V using inverse => V=CI^-1*z
		t.vec[0] = xlx2*matrices[9] - xlx*matrices[7]*CI.vec[0]; // a_12
		t.vec[1] = xlx2*matrices[10] - xlx*matrices[4]*CI.vec[2];
		t.vec[2] = xlx2*matrices[11] - xlx*matrices[5]*CI.vec[1];
		
		// do not need to calculate determinant because VAV/VBV will cancel it out anyway
//		det = CI.vec[0]*(CI.vec[1]*CI.vec[2]-xlx2*matrices[12]) + xlx*(matrices[4]*t.vec[1] + matrices[5]*t.vec[2]);
		
		// V = CI^-1*z (/det(CI) omitted)
		V.vec[0] = (CI.vec[1]*CI.vec[2]-xlx2*matrices[12])*z.vec[0]+t.vec[1]*z.vec[1]+t.vec[2]*z.vec[2];
		V.vec[1] = t.vec[1]*z.vec[0]+(CI.vec[0]*CI.vec[2]-xlx2*matrices[13])*z.vec[1]+t.vec[0]*z.vec[2];
		V.vec[2] = t.vec[2]*z.vec[0]+t.vec[0]*z.vec[1]+(CI.vec[0]*CI.vec[1]-xlx2*matrices[14])*z.vec[2];
		
		// numerator VAV=V*A_LF*V (uses fact that A_LF is symmetric)
		VAV = V.vec[0]*V.vec[0]*matrices[0]+V.vec[1]*V.vec[1]*matrices[1]+V.vec[2]*V.vec[2]*matrices[2];
		// denominator VBV=V*B_LF*V (uses fact that B_LF is symmetric)
		VBV = V.vec[0]*(V.vec[0]*matrices[3]+2.0*(V.vec[1]*matrices[4]+V.vec[2]*matrices[5]))
		     +V.vec[1]*(V.vec[1]*matrices[6]+2.0*V.vec[2]*matrices[7])
		     +V.vec[2]*V.vec[2]*matrices[8];
		//Calculate minimization parameter lambda
		if(VBV < EPS){
			cout << "ERROR: Denominator between oids " << kk << " and " << i << " in touch is too close to zero (" << VBV << ").\n";
			exit(3);
		}
		xlx2 = VAV/VBV;
		xlx = sqrt(xlx2); // independent of z-scaling (and determinant) -> also, interesting note: the sqrt is better than anything else in terms of speed and convergence
		lambda = xlx/(1.0+xlx);
		Var -= lambda;
		Var *= Var;
	}
	
	//Reconstruct CI and run a final iteration once converged
	CI.vec[0] = xlx*matrices[3] + matrices[0];
	CI.vec[1] = xlx*matrices[6] + matrices[1];
	CI.vec[2] = xlx*matrices[8] + matrices[2];
	
	t.vec[0] = CI.vec[1]*CI.vec[2] - xlx2*matrices[12];
	t.vec[1] = xlx*matrices[10] - matrices[4]*CI.vec[2];
	t.vec[2] = xlx*matrices[11] - matrices[5]*CI.vec[1];
	
	det = CI.vec[0]*t.vec[0] + xlx2*(matrices[4]*t.vec[1] + matrices[5]*t.vec[2]);
	if(fabs(det) < EPS){
		cout << "WARNING: Matrix is close to singular. det = " << det << "\n";
		cout << "Distance calculation failed. i = " << i << ", kk = " << kk << ", x = " << xlx << "\n";
		exit(2);
	}
	Vz = (z.vec[0]*t.vec[0]+2.0*xlx*(z.vec[1]*t.vec[1]+z.vec[2]*t.vec[2]))*z.vec[0]+
	     ((CI.vec[0]*CI.vec[2]-xlx2*matrices[13])*z.vec[1]+2.0*xlx*z.vec[2]*(xlx*matrices[9]-CI.vec[0]*matrices[7]))*z.vec[1]+
	     (CI.vec[0]*CI.vec[1]-xlx2*matrices[14])*z.vec[2]*z.vec[2];
	// return (sigma/r)^2
	return det/(lambda*Vz);
}

inline Mat33 MC_Elements::Rotate(Element* element, Vec3 &center, double &theta)
{
#ifdef USE_CMWC4096
	//Random number for picking an axis
	unsigned int choice=CMWC4096()%3;
	//Random number for the move itself
	theta = configuration->maxrot*(2.0*ranQ()-1.0);
#else
	//Random number for picking an axis
	unsigned int choice = (unsigned int)ran2_int(*configuration->idum)%3;
	//Random number for the move itself
	theta = configuration->maxrot*(2.0*ran2(*configuration->idum)-1.0);
#endif
	if(element->group) theta/=element->group->nr_elements;
	double co = cos(theta);
	double si = sin(theta);
	
	Mat33 rot;
	switch(choice){
		case 0:	/* yaw:
			 * (1   0   0)
			 * (0  co  si)
			 * (0 -si  co)
			 */
			rot.mat[0][0] = 1; rot.mat[0][1] = 0; rot.mat[0][2] = 0;
			rot.mat[1][0] = 0; rot.mat[1][1] = co; rot.mat[1][2] = si;
			rot.mat[2][0] = 0; rot.mat[2][1] = -si; rot.mat[2][2] = co;
			break;
		case 1:	/* pitch:
			 * ( co  0  si)
			 * ( 0   1   0)
			 * (-si  0  co)
			 */
			rot.mat[0][0] = co; rot.mat[0][1] = 0; rot.mat[0][2] = si;
			rot.mat[1][0] = 0; rot.mat[1][1] = 1; rot.mat[1][2] = 0;
			rot.mat[2][0] = -si; rot.mat[2][1] = 0; rot.mat[2][2] = co;
			break;
		case 2:	/* roll:
			 * ( co  si  0)
			 * (-si  co  0)
			 * ( 0   0   1)
			 */
			rot.mat[0][0] = co; rot.mat[0][1] = si; rot.mat[0][2] = 0;
			rot.mat[1][0] = -si; rot.mat[1][1] = co; rot.mat[1][2] = 0;
			rot.mat[2][0] = 0; rot.mat[2][1] = 0; rot.mat[2][2] = 1;
			break;
	}
	Rotate(element,center,rot);
	return rot;
}

/*!
 * Function provides randomization of individual group
 * Arguments:
 * 	group       ... group to randomize
 * 	kk          ... particular element to start rebuilding from (if independent is false)
 * 	independent ... if true:  rotation and translation of each group element is independent from the others
 * 			if false: rot/trans beginning from element kk (one way to think about it is "rebuilding" a molecule)
 */
inline double MC_Elements::GroupShuffle(Element_Group* group, const unsigned int kk)
{
	Element* element;
	double theta, co, si;
	double rmov=0.0;
	double inv_n=(double)1.0/group->Type->rand_elements;
	if(group->Type->rand_independent){
		// Now randomly rot/trans group_rand_frac*n individual elements from group
		// - currently does not use potential information to bias randomization
		Mat33 delta_rot;
		Vec3 delta_trans;
		unsigned int choice;
		for(unsigned int i=0; i<group->Type->rand_elements; i++){
#ifdef USE_CMWC4096
			unsigned int curr=(unsigned int)(ranQ()*group->nr_elements);
			element=&Elements[group->elements[curr]];
			// Random number for picking an axis
			choice=CMWC4096()%3;
			// Random number for the move itself
			theta = configuration->maxrot*(2.0*ranQ()-1.0)*inv_n;
#else
			unsigned int curr=(unsigned int)(ran2(*configuration->idum)*group->nr_elements);
			element=&Elements[group->elements[curr]];
			// Random number for picking an axis
			choice = (unsigned int)ran2_int(*configuration->idum)%3;
			// Random number for the move itself
			theta = configuration->maxrot*(2.0*ran2(*configuration->idum)-1.0)*inv_n;
#endif
			co = cos(theta); si = sin(theta);
			
			switch(choice){
				case 0:	/* yaw
					 * (1   0   0)
					 * (0  co  si)
					 * (0 -si  co)
					 */
					delta_rot.mat[0][0] = 1; delta_rot.mat[0][1] = 0; delta_rot.mat[0][2] = 0;
					delta_rot.mat[1][0] = 0; delta_rot.mat[1][1] = co; delta_rot.mat[1][2] = si;
					delta_rot.mat[2][0] = 0; delta_rot.mat[2][1] = -si; delta_rot.mat[2][2] = co;
					break;
				case 1:	/* pitch
					 * ( co  0  si)
					 * ( 0   1   0)
					 * (-si  0  co)
					 */
					delta_rot.mat[0][0] = co; delta_rot.mat[0][1] = 0; delta_rot.mat[0][2] = si;
					delta_rot.mat[1][0] = 0; delta_rot.mat[1][1] = 1; delta_rot.mat[1][2] = 0;
					delta_rot.mat[2][0] = -si; delta_rot.mat[2][1] = 0; delta_rot.mat[2][2] = co;
					break;
				case 2:	/* roll
					 * ( co  si  0)
					 * (-si  co  0)
					 * ( 0   0   1)
					 */
					delta_rot.mat[0][0] = co; delta_rot.mat[0][1] = si; delta_rot.mat[0][2] = 0;
					delta_rot.mat[1][0] = -si; delta_rot.mat[1][1] = co; delta_rot.mat[1][2] = 0;
					delta_rot.mat[2][0] = 0; delta_rot.mat[2][1] = 0; delta_rot.mat[2][2] = 1;
					break;
			}
#ifdef USE_CMWC4096
			delta_trans.vec[0]=(configuration->maxtrans)*(2.0*ranQ()-1.0)*inv_n;
			delta_trans.vec[1]=(configuration->maxtrans)*(2.0*ranQ()-1.0)*inv_n;
			delta_trans.vec[2]=(configuration->maxtrans)*(2.0*ranQ()-1.0)*inv_n;
#else
			delta_trans.vec[0]=(configuration->maxtrans)*(2.0*ran2(*configuration->idum)-1.0)*inv_n;
			delta_trans.vec[1]=(configuration->maxtrans)*(2.0*ran2(*configuration->idum)-1.0)*inv_n;
			delta_trans.vec[2]=(configuration->maxtrans)*(2.0*ran2(*configuration->idum)-1.0)*inv_n;
#endif
			
#ifdef LAB_FRAME_ROTATION
			element->rot=delta_rot*element->rot;
#else
			element->rot*=delta_rot;
#endif
			element->dipole=element->rot*element->MyType->initial_dipole;
			// charges
			for(unsigned int j=0; j<element->MyType->nr_charges; j++) element->q_pos[j]=element->rot*element->MyType->q_pos[j];
			// links
			for(unsigned int j=0; j<element->nr_interactions; j++){
				if(element->interactions[j].location) *element->interactions[j].location=element->rot*(*element->interactions[j].initial_location);
				if(element->interactions[j].normal) *element->interactions[j].normal=element->rot*(*element->interactions[j].initial_normal);
				if(element->interactions[j].tangent) *element->interactions[j].tangent=element->rot*(*element->interactions[j].initial_tangent);
			}
			rmov+=theta*theta;
			if(!element->MyType->rot_notrans) element->center+=delta_trans;
		}
	} else{ // "rebuild" group
		memset(rebuild_storage.visited,0,group->nr_elements*sizeof(unsigned int)); // initially no element has been visited
		unsigned int current, next, newcount;
		unsigned int group_start=0;
		for(unsigned int i=0; i<group->nr_elements; i++){
			// need to define where we start from in the group
			rebuild_storage.delta_trans[i].vec[0]=0.0; rebuild_storage.delta_trans[i].vec[1]=0.0; rebuild_storage.delta_trans[i].vec[2]=0.0;
			rebuild_storage.delta_rot[i].M3Eye();
			rebuild_storage.theta2[i]=0.0;
			if(group->elements[i]==kk){
				group_start=i;
#if DEBUG_LEVEL>3
				cout << "start from: " << group_start+1 << "\n";
#endif
			}
		}
		current=group_start; // the originating element counts as visited ...
		rebuild_storage.visited[current]=current+1; // the starting element has been visited by itself (to avoid going back to itself)
		unsigned int nr_elements=Elements[kk].nr_interactions;
		for(unsigned int i=0; i<nr_elements; i++){ // start from element kk
			next=Elements[kk].interactions[i].partner;
			rebuild_storage.visited[next]=group_start+1;
			rebuild_storage.to_link[next]=i;
			rebuild_storage.element_storage[i<<1]=next;
		}
		
		bool alternate=false;
		unsigned int nr_randomized=0;
#if DEBUG_LEVEL>3
		unsigned int statistics=0;
#endif
		unsigned int rand_per_round=0; // cumulative number of elements randomized per round
		while(nr_elements>0){
			newcount=0;
			for(unsigned int i=0; i<nr_elements; i++){
				current=rebuild_storage.element_storage[(i<<1)+alternate];
#if DEBUG_LEVEL>3
				cout << "current: " << current+1 << "\n";
#endif
				element=&Elements[group->elements[current]];
				// apply rot/trans of previous elements
				if(rand_per_round>0){
					// find out which of the randomized elements are en route to the current one
					memset(rebuild_storage.en_route,0,group->nr_elements*sizeof(bool)); // none en route by default
					unsigned int previous=current;
					int current_j=rand_per_round-1;
					do{
						previous=rebuild_storage.visited[previous]-1;
						unsigned int j=current_j;
						while((previous!=rebuild_storage.randomized[j]) && (j>0)) j--;
						if(previous==rebuild_storage.randomized[j]){
							rebuild_storage.en_route[j]=true;
							current_j=j-1;
						}
					} while((rebuild_storage.visited[previous]!=group_start+1) && (current_j>=0));
					// now apply rot/trans starting from first change in path to current element
					for(unsigned int j=0; j<rand_per_round; j++){
						if(rebuild_storage.en_route[j]){
#if DEBUG_LEVEL>3
							cout << "-> randomized en route: " << rebuild_storage.randomized[j]+1 << "\n";
#endif
							// rotate from current position of element where rotation originated
							element->center = Elements[group->elements[rebuild_storage.randomized[j]]].center+rebuild_storage.delta_rot[j]*(element->center+rebuild_storage.delta_trans[j]-Elements[group->elements[rebuild_storage.randomized[j]]].center);
							RotateElement(element,rebuild_storage.delta_rot[j]); // rotate element
							rmov+=rebuild_storage.theta2[j];
						}
					}
				}
#if DEBUG_LEVEL>3
				statistics++;
#endif
				Element* prev=&Elements[group->elements[rebuild_storage.visited[current]-1]];
				// get position of link on previous link
				Vec3 prevlink=prev->center;
				if(prev->interactions[rebuild_storage.to_link[current]].location) prevlink+=(*prev->interactions[rebuild_storage.to_link[current]].location);
				// enforce fixed bond length actively -- the math doesn't need it but numerical precision does
				if(!prev->interactions[rebuild_storage.to_link[current]].allow_bond_stretch){
					Vec3 v=element->center;
					if(element->interactions[prev->interactions[rebuild_storage.to_link[current]].back_link].location) v+=(*element->interactions[prev->interactions[rebuild_storage.to_link[current]].back_link].location);
					v-=prevlink;
					double bl=prev->interactions[rebuild_storage.to_link[current]].bond_length;
					double v2=v*v;
					if(fabs(v2-bl*bl)>EPS){ // correct only if needed (saves a sqrt)
#if DEBUG_LEVEL>3
						cout << "-> Adjusted bond distance\n";
#endif
						v*=bl/sqrt(v2);
						v+=prevlink; // now points to link on current element (element->center+location)
						if(element->interactions[prev->interactions[rebuild_storage.to_link[current]].back_link].location) v-=(*element->interactions[prev->interactions[rebuild_storage.to_link[current]].back_link].location);
						element->center=v;
					}
				}
				
#ifdef USE_CMWC4096
				if((nr_randomized<group->Type->rand_elements) && (ranQ()>0.5)){ // randomize only every other group element
#else
				if((nr_randomized<group->Type->rand_elements) && (ran2(*configuration->idum)>0.5)){ // randomize only every other group element
#endif
					// only randomize if the following conditions are met:
					// - current element *or* element we're coming from is not in a ring
					// - link between current and previous element is not fixed
					if(((!group->Type->elements_in_ring[current]) || (!group->Type->elements_in_ring[rebuild_storage.visited[current]-1])) && !element->interactions[prev->interactions[rebuild_storage.to_link[current]].back_link].fixed){
#if DEBUG_LEVEL>3
						cout << "randomize element: " << current+1;
#endif
						// choose what to randomize:
						// 0 ... bond rotation
						// 1 ... bond orientation (with respect to element connection we're coming from)
						// 2 ... bond orientation (with respect to current element connection)
						// 3 ... randomize bond length (if bond is not fixed length)
						unsigned int randmax=1;
						unsigned int addbend=0;
						if(element->interactions[prev->interactions[rebuild_storage.to_link[current]].back_link].allow_bond_bend){
							randmax+=2; // randmax=3
							addbend=2;
						}
						if(element->interactions[prev->interactions[rebuild_storage.to_link[current]].back_link].allow_bond_stretch) randmax++; // randmax = 2 or 4
						unsigned int randomize = 0;
#ifdef USE_CMWC4096
						if(randmax>1) randomize=CMWC4096()%randmax; // randomize now is random number between 0 and randmax-1
#else
						if(randmax>1) randomize=(unsigned int)ran2_int(*configuration->idum)%randmax;
#endif
						// store old position temporarily
						rebuild_storage.delta_trans[rand_per_round]=element->center;
						// do rotation around link on previous element and store rotational delta
						double thetatemp=0.0;
						
						if((addbend>0) && (randomize==1)){
							rebuild_storage.delta_rot[rand_per_round]=Rotate(element,prevlink,thetatemp);
						} else rebuild_storage.delta_rot[rand_per_round].M3Eye();
						double rmovtemp=thetatemp*thetatemp;
						// get position of link on current element (need to recalculate)
						Vec3 currlink=element->center;
						if(element->interactions[prev->interactions[rebuild_storage.to_link[current]].back_link].location) currlink+=(*element->interactions[prev->interactions[rebuild_storage.to_link[current]].back_link].location);
#if DEBUG_LEVEL>3
						// way too complicated (rebuild_storage.visited[current]) would do the trick -- if the complicated version works, however, then everything's good ...
						cout << " (came from: " << element->interactions[prev->interactions[rebuild_storage.to_link[current]].back_link].partner+1 << ", method " << (unsigned int)randomize << ")\n";
#endif
						// do rotation around link on current element and store rotational delta
						Mat33 rottemp;
						if((addbend>0) && (randomize==2)) rottemp=Rotate(element,currlink,thetatemp); // fully random
						if(randomize==0){ // only around link axis
#ifdef USE_CMWC4096
							thetatemp=configuration->maxrot*(2.0*ranQ()-1.0)*inv_n;
#else
							thetatemp=configuration->maxrot*(2.0*ran2(*configuration->idum)-1.0)*inv_n;
#endif
							Vec4 axistemp(currlink-prevlink,thetatemp);
							rottemp=AxisAngle2Rot(axistemp);
							RotateElement(element,currlink,rottemp);
						}
						rebuild_storage.delta_rot[rand_per_round]=rottemp*rebuild_storage.delta_rot[rand_per_round];
						rmovtemp+=thetatemp*thetatemp;
						rmov+=rmovtemp;
						// now determine overall (effective) rotation angle to use in successive rmov calculations
						rebuild_storage.theta2[rand_per_round]=rmovtemp;
						// change link stretch length by +/- 2% if bond is not fixed in length (in other words if bond length is variable)
						if(randomize==(1+addbend)){
#ifdef USE_CMWC4096
							element->center=prevlink+(currlink-prevlink)*(0.98+0.04*ranQ());
#else
							element->center=prevlink+(currlink-prevlink)*(0.98+0.04*ran2(*configuration->idum));
#endif
							// technically, element->center points to link on current element, correct by subtracting (link location = element center + link pos)
							if(element->interactions[prev->interactions[rebuild_storage.to_link[current]].back_link].location) element->center-=(*element->interactions[prev->interactions[rebuild_storage.to_link[current]].back_link].location);
						}
						// store translational change
						rebuild_storage.delta_trans[rand_per_round]=element->center-rebuild_storage.delta_trans[rand_per_round];
						
						nr_randomized++;
						rebuild_storage.randomized[rand_per_round]=current;
						rand_per_round++;
#if DEBUG_LEVEL>3
						if(nr_randomized>=group->Type->rand_elements) cout << "-> randomized: " << nr_randomized << ", randomized per round: " << rand_per_round << "\n";
#endif
					}
				}
#if DEBUG_LEVEL>3
				cout << "-> randomized: " << nr_randomized << ", randomized per round: " << rand_per_round << "\n";
#endif
				for(unsigned int j=0; j<element->nr_interactions; j++){ // prepare next round away
					next=element->interactions[j].partner;
					if(rebuild_storage.visited[next]==0){ // only follow link if other end is not were we've already been
						rebuild_storage.visited[next]=current+1;
						rebuild_storage.to_link[next]=j;
						rebuild_storage.element_storage[(newcount<<1)+(!alternate)]=next;
						newcount++;
					}
				}
			}
			if(((newcount==0) && (nr_randomized<group->Type->rand_elements)) && (nr_randomized>0)){ // we visited all elements of the group, randomized some, and want to randomize some more ...
#if DEBUG_LEVEL>3
				cout << "visited whole group, restart from " << group_start+1 << "\n";
#endif
				alternate=false;
				memset(rebuild_storage.visited,0,group->nr_elements*sizeof(unsigned int)); // initially no element has been visited
				current=group_start;
				rebuild_storage.visited[current]=current+1; // the originating element has been visited from itself
				nr_elements=Elements[kk].nr_interactions;
				for(unsigned int i=0; i<nr_elements; i++){
					next=Elements[kk].interactions[i].partner;
					rebuild_storage.visited[next]=group_start+1;
					rebuild_storage.to_link[next]=i;
					rebuild_storage.element_storage[i<<1]=next;
				}
				rand_per_round=0;
			} else{
				alternate=(!alternate); // flip-flop
				nr_elements=newcount;
			}
		}
#if DEBUG_LEVEL>3
		cout << "took " << statistics << " iterations\n";
#endif
	}
	return rmov;
}

inline double MC_Elements::ShrinkVolume(unsigned int step)
{
	double tmov=0.0;
	
	if(configuration->targetV<configuration->V && (!configuration->randsteps_novolume || (configuration->randsteps_novolume && (step>configuration->randsteps)))){
		unsigned int n0=(configuration->randsteps>>1);
		double newvolume;
		if(step-configuration->randsteps_novolume*configuration->randsteps>=n0){
			newvolume=configuration->targetV;
		} else newvolume=configuration->V-V_diff/sqrt((double)(step-configuration->randsteps_novolume*configuration->randsteps));
#if DEBUG_LEVEL>3
		cout << "step: " << step << ", " << configuration->V-newvolume << ", " << newvolume << "\n";
#endif
		double factor=cbrt(newvolume/configuration->V);
		// loop over elements which are not in groups
		for(unsigned int k = nr_group_elements; k < nr_elements; k++) Elements[k].center*=factor;
		// loop over groups
		Element_Group* group;
		for(unsigned int k = 0; k < nr_groups; k++){
			group=Groups[k];
			Vec3 delta_trans=group->center*factor-group->center; // new position - old position points to new position
			group->center+=delta_trans; // move group center but leave group elements undisturbed
			for(unsigned int i=0; i<group->nr_elements; i++){
				Elements[group->elements[i]].center+=delta_trans;
			}
		}
		// make changes "official"
		switch(configuration->latticetype){
			case 1:
			case 2: // Rectangular box
				configuration->boxlength[0]*=factor;
				configuration->inv_boxlength[0]=1.0/configuration->boxlength[0];
				configuration->boxlength[1]*=factor;
				configuration->inv_boxlength[1]=1.0/configuration->boxlength[1];
				configuration->boxlength[2]*=factor;
				configuration->inv_boxlength[2]=1.0/configuration->boxlength[2];
				break;
			case 3: // Spherical volume
				configuration->spherecylr*=factor;
				break;
			case 4: // Cylindrical volume
				configuration->spherecylr*=factor;
				configuration->boxlength[2]*=factor;
				configuration->inv_boxlength[2]=1.0/configuration->boxlength[2];
				break;
			default: // Default condition
				cout << "Cannot change volume. Lattice type " << configuration->latticetype << " not recognized. Exiting.\n";
				exit(1);
		}
		configuration->nndist*=factor;
		configuration->V=newvolume;
		if(fabs(newvolume-configuration->targetV)<EPS){
			cout << "\n\tTarget volume " << configuration->targetV << " Angström³ reached.\n";
			changevolume=false;
		}
	}
	return tmov;
}

#endif

