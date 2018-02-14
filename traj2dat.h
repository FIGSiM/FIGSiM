/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * Robinson group Monte Carlo simulation code trajectory file analysis tool (e.g. pair correlation functions)
 * Currently requires the entire simulation code for configuration trajectory
 * parsing, as well as the simulation configuration (.conf) file. Replaces and
 * expands upon functionality of gofr.cpp
 * Command line utility with two arguments, (1) the trajectory file, and, optionally, (2) the starting step
 * Usage example: traj2dat test0.traj 1000
 * Created by Andreas Tillack and Lewis Johnson, August 2011
 */

#include "ConfigReader.h"
#include "MC_Config.h"
#include "ScalarMat.h"
#include "VecMat.h"
#include "setup.h"

const unsigned int nr_additional_types=5;
string additional_types[nr_additional_types] = {"x-Axis", "y-Axis", "z-Axis", "Efield", "totalM"};
const unsigned int nr_additional_subtypes=2;
string additional_subtypes[nr_additional_subtypes] = {"Center", "Dipole"};
const unsigned nr_corr_types=12;
const unsigned nr_calc_types=1;
const unsigned nr_corr_base_types=5; // the first 5 correlation types are base types (no number after type)
string corr_types[nr_corr_types]={"gr","sdf","angle","alphabetagamma","gammaalpha","Q","distances","gmu","galpha","cosine","VLJ","dihedral"};
string output_corr_types[nr_corr_types]={"g","sdf","angle","Euler angles","gammaalpha","Q-tensor","distances","g_{mu}","g_{alpha}","cos","V_{LJ}","dihedral"};
int corr_types_items_involved[nr_corr_types]={2,2,2,1,2,-1,2,2,2,-1,2,4}; // items to be correlated (Q can have any number of items)

// Enumerated list of correlation function types
enum cftypes_t{
	gr=0,
	sdf,
	angle,
	alphabetagamma,
	gammaalpha,
	Q,
	distances,
	gmu,
	galpha,
	cosine,
	VLJ,
	dihedral
};

typedef struct _Corr_Item{
	bool is_element, is_group, is_LOD, is_component;
	int element_idx, group_idx, type_idx, traj_group;
	unsigned int LOD_nr;
} Corr_Item;

///Struct for passing correlation function parameters
typedef struct _cfparams_t {
	Config_Data* configuration;
	Traj_EP** trajectory;
	double rmax;
	double* V;
	double* Xm;
	double* time;
	
	unsigned nr_correlations;
	cftypes_t* types;
	unsigned int* typenr;
	double** additional_parameters;
	
	double*** results;
	
	bool group_internal;
	int groupnr;
	
	unsigned int items_involved;
	bool any_number_of_items;
	Corr_Item* items;
	double* weights;
	string* item_names;
	
	Vec3** positions;
	Vec3** vectors;
	Mat33** rots;
	Vec3** max_vectors;
	Vec3 rotatedE_direction;
	unsigned int* total_nr;
	
	unsigned int start_frame_in_data;
	unsigned int startstep, endstep, startframe, endframe, n_frames, nr_bins, average_frames;
	
	bool individual_output; //if true, do old line-by-line output
	bool find_max_item;
	bool find_max_abs_item;
	bool report;
} cfparams_t;

inline void apply_PBCs(Vec3 &rvec, Config_Data* configuration) // copied from MC_Elements.h (not included for trajectory tool)
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

unsigned int GetPosDipole(cfparams_t &params, unsigned int item_nr, unsigned int frame, Vec3* &positions, Vec3* &vectors, Vec3* &max_vectors, Mat33* &rots, unsigned int &total_nr, unsigned int start_idx, bool just_one) // <- note the "Vec3* &" -- it is needed in order to change the actual pointer
{
	unsigned int nr=0;
	unsigned int i;
	double q;
	
	int group_t, element_t;
	unsigned int idx;
	
	Element_Type* e_t;
	Mat33 e_r;
	
	if(item_nr<params.items_involved){
		group_t=params.items[item_nr].group_idx;
		element_t=params.items[item_nr].element_idx;
		if(params.items[item_nr].is_element || params.items[item_nr].is_group || params.items[item_nr].is_component || params.items[item_nr].is_LOD || (params.group_internal && (params.items[item_nr].type_idx<0))){
			idx=start_idx;
			while((idx<params.configuration->n_oids) && !(just_one && idx>start_idx)){
				if(params.items[item_nr].is_element){ // single element
					if((int)params.trajectory[idx][frame].element_type==element_t){
						nr++;
						if(nr>total_nr){
							total_nr=nr;
							positions=(Vec3*)realloc(positions,sizeof(Vec3)*total_nr);
							vectors=(Vec3*)realloc(vectors,sizeof(Vec3)*total_nr);
							max_vectors=(Vec3*)realloc(max_vectors,sizeof(Vec3)*total_nr);
							rots=(Mat33*)realloc(rots,sizeof(Mat33)*total_nr);
							if((!positions || !vectors) || !max_vectors || !rots){
								cout << "Not enough memory.\n";
								exit(3);
							}
						}
						e_t=params.configuration->element_types[element_t];
						e_r=AxisAngle2Rot(params.trajectory[idx][frame].rotation_vector);
						positions[nr-1]=params.trajectory[idx][frame].position;
						if(params.configuration->offctrmu) positions[nr-1]+=e_r*e_t->mu_pos;
						vectors[nr-1]=e_r*e_t->initial_dipole;
						max_vectors[nr-1]=e_t->initial_dipole;
						for(i=0; i<e_t->nr_charges; i++){
							q=e_t->q[i];
							vectors[nr-1]+=(e_r*e_t->q_pos[i])*q;
							max_vectors[nr-1]+=(e_t->q_pos[i])*q;
						}
						rots[nr-1]=e_r; // orientation is in the element's frame
						// return element's z-direction in case no dipoles exist
						if(vectors[nr-1]*vectors[nr-1]<EPS_CL) vectors[nr-1]=e_r*Vec3(0,0,1);
						if(max_vectors[nr-1]*max_vectors[nr-1]<EPS_CL) max_vectors[nr-1]=Vec3(0,0,1);
					}
					idx++;
				} else{ // not a single element
					int search_group=-2; // -2 does not exist in trajectory file and hence can safely be used as no-search fallback
					// If user specified fully atomistic model but we only got LOD version in trajectory file make sure to find that as well
					if(params.configuration->groups[params.trajectory[idx][frame].group_type]->levelofdetail>0) search_group=params.configuration->groups[params.trajectory[idx][frame].group_type]->Type->LOD->groups[0]->type;
					if((group_t==search_group) || (params.trajectory[idx][frame].group_type==group_t)){
						params.items[item_nr].traj_group=params.trajectory[idx][frame].group_type;
						nr++;
						if(nr>total_nr){
							total_nr=nr;
							positions=(Vec3*)realloc(positions,sizeof(Vec3)*total_nr);
							vectors=(Vec3*)realloc(vectors,sizeof(Vec3)*total_nr);
							max_vectors=(Vec3*)realloc(max_vectors,sizeof(Vec3)*total_nr);
							rots=(Mat33*)realloc(rots,sizeof(Mat33)*total_nr);
							if((!positions || !vectors) || !max_vectors || !rots){
								cout << "Not enough memory.\n";
								exit(3);
							}
						}
						// calculate charge center for whole (fully-atomistic) group
						Vec3 charge_center(0.0);
						double charge_wsum=0.0;
						if(params.configuration->groups[params.trajectory[idx][frame].group_type]->levelofdetail>0){
							for(i=0; i<params.configuration->groups[params.trajectory[idx][frame].group_type]->Type->LOD->groups[0]->nr_elements; i++){
								unsigned int ellipsoid_idx=params.configuration->groups[group_t]->Type->LOD->element_in_ellipsoid[params.configuration->groups[params.trajectory[idx][frame].group_type]->levelofdetail-1][i];
								
								Element* thiselement=&params.configuration->group_elements[params.configuration->groups[params.trajectory[idx][frame].group_type]->elements[ellipsoid_idx]];
								Element* fullelement=&params.configuration->group_elements[params.configuration->groups[params.trajectory[idx][frame].group_type]->Type->LOD->groups[0]->elements[i]];
								// rotation difference for fully atomistic element
								Mat33 delta_rot=AxisAngle2Rot(params.trajectory[idx+ellipsoid_idx][frame].rotation_vector)*thiselement->rot.M3Transpose();
								Vec3 offset=delta_rot*(fullelement->center-thiselement->center);
								Mat33 element_rot=delta_rot*fullelement->rot;
								
								Vec3 center=params.trajectory[idx+ellipsoid_idx][frame].position+offset;
								for(unsigned int j=0; j<fullelement->MyType->nr_charges; j++){
									charge_wsum+=fabs(fullelement->MyType->q[j]);
									charge_center+=(center+element_rot*fullelement->MyType->q_pos[j])*fabs(fullelement->MyType->q[j]);
								}
							}
						} else{
							for(i=0; i<params.configuration->groups[params.trajectory[idx][frame].group_type]->nr_elements; i++){
								Vec3 center=params.trajectory[idx+i][frame].position;
								Mat33 rot=AxisAngle2Rot(params.trajectory[idx+i][frame].rotation_vector);
								Element_Type* type=params.configuration->element_types[params.trajectory[idx+i][frame].element_type];
								for(unsigned int j=0; j<type->nr_charges; j++){
									charge_wsum+=fabs(type->q[j]);
									charge_center+=(center+rot*type->q_pos[j])*fabs(type->q[j]);
								}
							}
						}
						if(charge_wsum>EPS) charge_center/=charge_wsum; // sanity check
						if(params.items[item_nr].is_component){ // component inside a group (includes LOD groups)
							positions[nr-1]=Vec3(0.0);
							Vec3 max_position=Vec3(0.0);
							vectors[nr-1]=Vec3(0.0);
							max_vectors[nr-1]=Vec3(0.0);
							charge_center=Vec3(0.0);
							charge_wsum=0.0;
							rots[nr-1].M3Eye(); // for the moment, components have no particular rotation
							unsigned int nr_component_elements=params.configuration->groups[params.trajectory[idx][frame].group_type]->Type->LOD->nr_elements[element_t];
							double qsum=0.0;
							for(i=0; i<nr_component_elements; i++){
								Vec3 position;
								if(params.configuration->groups[params.trajectory[idx][frame].group_type]->levelofdetail>0){
									// fully atomistic component inside LOD group
									unsigned int component_element=params.configuration->groups[params.trajectory[idx][frame].group_type]->Type->LOD->component_elements[element_t][i];
									unsigned int ellipsoid_idx=params.configuration->groups[group_t]->Type->LOD->element_in_ellipsoid[params.configuration->groups[params.trajectory[idx][frame].group_type]->levelofdetail-1][component_element];
									Element* thiselement=&params.configuration->group_elements[params.configuration->groups[params.trajectory[idx][frame].group_type]->elements[ellipsoid_idx]];
#if DEBUG_LEVEL>3
									cout << "\n" << thiselement->MyType->name;
#endif
									Element* fullelement=&params.configuration->group_elements[params.configuration->groups[group_t]->Type->LOD->groups[0]->elements[component_element]];
#if DEBUG_LEVEL>3
									cout << "\n" << fullelement->MyType->name << "\n";
#endif
									// rotation difference for fully atomistic element
									Mat33 delta_rot=AxisAngle2Rot(params.trajectory[idx+ellipsoid_idx][frame].rotation_vector)*thiselement->rot.M3Transpose();
									Vec3 offset=delta_rot*(fullelement->center-thiselement->center);
									Mat33 element_rot=delta_rot*fullelement->rot;
									position=params.trajectory[idx+ellipsoid_idx][frame].position+offset;
									vectors[nr-1]+=element_rot*fullelement->MyType->initial_dipole;
									max_position+=fullelement->center-thiselement->center;
									max_vectors[nr-1]+=fullelement->rot*fullelement->MyType->initial_dipole;
									for(unsigned int j=0; j<fullelement->MyType->nr_charges; j++){
										q=fullelement->MyType->q[j];
										qsum+=q;
										vectors[nr-1]+=(position+element_rot*fullelement->MyType->q_pos[j])*q;
										max_vectors[nr-1]+=(fullelement->center-thiselement->center+fullelement->rot*fullelement->MyType->q_pos[j])*q;
										charge_wsum+=fabs(q);
										charge_center+=(position+element_rot*fullelement->MyType->q_pos[j])*fabs(q);
									}
								} else{ // component inside fully atomistic group
									unsigned int component_element=params.configuration->groups[params.trajectory[idx][frame].group_type]->Type->LOD->component_elements[element_t][i];
									Element* fullelement=&params.configuration->group_elements[params.configuration->groups[group_t]->elements[component_element]];
									Element_Type* type=params.configuration->element_types[params.trajectory[idx+i][frame].element_type];
									Mat33 element_rot=AxisAngle2Rot(params.trajectory[idx+i][frame].rotation_vector);
									
									position=params.trajectory[idx+i][frame].position;
									vectors[nr-1]+=element_rot*type->initial_dipole;
									max_position+=fullelement->center;
									max_vectors[nr-1]+=fullelement->rot*type->initial_dipole;
									for(unsigned int j=0; j<type->nr_charges; j++){
										q=type->q[j];
										qsum+=q;
										vectors[nr-1]+=(position+element_rot*type->q_pos[j])*q;
										max_vectors[nr-1]+=(fullelement->rot*type->q_pos[j])*q;
										charge_wsum+=fabs(q);
										charge_center+=(position+element_rot*type->q_pos[j])*fabs(q);
									}
								}
								positions[nr-1]+=position;
							}
							if(charge_wsum>EPS) charge_center/=charge_wsum; // sanity check
							if(nr_component_elements>0){
								positions[nr-1]/=nr_component_elements; // component center
//								vectors[nr-1]-=positions[nr-1]*qsum; // overall dipole moment with respect to component center
								vectors[nr-1]-=charge_center*qsum; // overall dipole moment with respect to component center
								max_position/=nr_component_elements;
								max_vectors[nr-1]-=max_position*qsum;
							}
						} else{ // group
							if(element_t>=0){ // element in group
								if(params.trajectory[idx][frame].group_type==group_t){ // element in normal group
									e_t=params.configuration->element_types[params.trajectory[idx+element_t][frame].element_type];
									e_r=AxisAngle2Rot(params.trajectory[idx+element_t][frame].rotation_vector);
									positions[nr-1]=params.trajectory[idx+element_t][frame].position;
									vectors[nr-1]=e_r*e_t->initial_dipole;
									max_vectors[nr-1]=e_t->initial_dipole;
									for(i=0; i<e_t->nr_charges; i++){
										q=e_t->q[i];
										vectors[nr-1]+=(e_r*e_t->q_pos[i])*q;
										max_vectors[nr-1]+=(e_t->q_pos[i])*q;
									}
									rots[nr-1]=e_r*params.configuration->group_elements[params.configuration->groups[params.trajectory[idx][frame].group_type]->elements[element_t]].rot.M3Transpose();
									if(max_vectors[nr-1]*max_vectors[nr-1]<EPS_CL) max_vectors[nr-1]=Vec3(0,0,1);
								} else{ // fully atomistic element in LOD group
									if(params.configuration->groups[params.trajectory[idx][frame].group_type]->levelofdetail>0){
										unsigned int ellipsoid_idx=params.configuration->groups[group_t]->Type->LOD->element_in_ellipsoid[params.configuration->groups[params.trajectory[idx][frame].group_type]->levelofdetail-1][element_t];
										
										Element* thiselement=&params.configuration->group_elements[params.configuration->groups[params.trajectory[idx][frame].group_type]->elements[ellipsoid_idx]];
										Element* fullelement=&params.configuration->group_elements[params.configuration->groups[params.trajectory[idx][frame].group_type]->Type->LOD->groups[0]->elements[element_t]];
#if DEBUG_LEVEL>3
										cout << "\n" << thiselement->MyType->name << "\n" << fullelement->MyType->name << "\n";
#endif
										// rotation difference for fully atomistic element
										Mat33 delta_rot=AxisAngle2Rot(params.trajectory[idx+ellipsoid_idx][frame].rotation_vector)*thiselement->rot.M3Transpose();
										Vec3 offset=delta_rot*(fullelement->center-thiselement->center);
										Mat33 element_rot=delta_rot*fullelement->rot;
										
										positions[nr-1]=params.trajectory[idx+ellipsoid_idx][frame].position+offset;
										vectors[nr-1]=element_rot*fullelement->MyType->initial_dipole;
										max_vectors[nr-1]=fullelement->rot*fullelement->MyType->initial_dipole;
										rots[nr-1]=element_rot;
										for(i=0; i<fullelement->MyType->nr_charges; i++){
											q=fullelement->MyType->q[i];
											vectors[nr-1]+=(element_rot*fullelement->MyType->q_pos[i])*q;
											max_vectors[nr-1]+=(fullelement->rot*fullelement->MyType->q_pos[i])*q;
										}
										if(max_vectors[nr-1]*max_vectors[nr-1]<EPS_CL) max_vectors[nr-1]=Vec3(0,0,1);
									} else{ // This should not happen
										cout << "Strange days: Oh Ellipsoid, Where art thou?\n";
										exit(42);
									}
								}
								if(vectors[nr-1]*vectors[nr-1]<EPS_CL){ // in case no dipole inside the chosen element, return whole group's dipole (but keep element location)
									vectors[nr-1]=Vec3(0.0);
									max_vectors[nr-1]=Vec3(0.0);
									unsigned int nr_elements=params.configuration->groups[params.trajectory[idx][frame].group_type]->nr_elements;
									double sigma=0.0;
									for(i=0; i<nr_elements; i++){
										Element* fullelement=&params.configuration->group_elements[params.configuration->groups[params.trajectory[idx][frame].group_type]->elements[i]];
										Vec3 center=params.trajectory[idx+i][frame].position;
										Mat33 rot=AxisAngle2Rot(params.trajectory[idx+i][frame].rotation_vector);
										Element_Type* type=params.configuration->element_types[params.trajectory[idx+i][frame].element_type];
										vectors[nr-1]+=rot*type->initial_dipole;
										max_vectors[nr-1]+=fullelement->rot*type->initial_dipole;
										// since the sum of charges per group is zero the following works
										for(unsigned int j=0; j<type->nr_charges; j++){
											sigma+=type->q[j];
											vectors[nr-1]+=(center+rot*type->q_pos[j])*type->q[j];
											max_vectors[nr-1]+=(fullelement->center+fullelement->rot*type->q_pos[j])*type->q[j];
										}
									}
									vectors[nr-1]-=charge_center*sigma;
									max_vectors[nr-1]-=positions[nr-1]*sigma;
								}
							} else{ // get center and dipole for the whole group -- the good news is that's independent of the level of detail, the center position and dipole moment are conserved
								positions[nr-1]=Vec3(0.0);
								vectors[nr-1]=Vec3(0.0);
								max_vectors[nr-1]=Vec3(0.0);
								rots[nr-1].M3Eye(); // for the moment, no particular rotation
								Element_Group* group=params.configuration->groups[params.trajectory[idx][frame].group_type];
								unsigned int nr_elements=group->nr_elements;
								double sigma=0.0;
								double w;
								double wsum=0.0;
								if(group->levelofdetail<=0){ // all-atom group
									for(i=0; i<nr_elements; i++){
										Element* fullelement=&params.configuration->group_elements[group->elements[i]];
										Vec3 center=params.trajectory[idx+i][frame].position;
										Mat33 rot=AxisAngle2Rot(params.trajectory[idx+i][frame].rotation_vector);
										Element_Type* type=params.configuration->element_types[params.trajectory[idx+i][frame].element_type];
										vectors[nr-1]+=rot*type->initial_dipole;
										w=type->initial_dipole.V3Norm();
										max_vectors[nr-1]+=fullelement->rot*type->initial_dipole;
										// since the sum of charges per group is zero the following works
										for(unsigned int j=0; j<type->nr_charges; j++){
											sigma+=type->q[j];
											w+=fabs(type->q[j]);
											vectors[nr-1]+=(center+rot*type->q_pos[j])*type->q[j];
											max_vectors[nr-1]+=(fullelement->center+fullelement->rot*type->q_pos[j])*type->q[j];
										}
										
										if(w<EPS) w=1.0;
										positions[nr-1]+=center*w;
										wsum+=w;
									}
								} else{
									nr_elements=group->Type->LOD->groups[0]->nr_elements; // go over all elements in underlying fully-atomistic group
									for(i=0; i<nr_elements; i++){
										unsigned int ellipsoid_idx=group->Type->LOD->element_in_ellipsoid[group->levelofdetail-1][i];
										Element* lodelement=&params.configuration->group_elements[group->elements[ellipsoid_idx]];
										Element* fullelement=&params.configuration->group_elements[group->Type->LOD->groups[0]->elements[i]];
										
										// rotation difference for fully atomistic element
										Mat33 delta_rot=AxisAngle2Rot(params.trajectory[idx+ellipsoid_idx][frame].rotation_vector)*lodelement->rot.M3Transpose();
										Vec3 offset=delta_rot*(fullelement->center-lodelement->center);
										Mat33 element_rot=delta_rot*fullelement->rot;
										Vec3 center=params.trajectory[idx+ellipsoid_idx][frame].position+offset;
										
										vectors[nr-1]+=element_rot*fullelement->MyType->initial_dipole;
										w=fullelement->MyType->initial_dipole.V3Norm();
										max_vectors[nr-1]+=fullelement->rot*fullelement->MyType->initial_dipole;
										for(unsigned int j=0; j<fullelement->MyType->nr_charges; j++){
											q=fullelement->MyType->q[j];
											w+=fabs(q);
											vectors[nr-1]+=(center+element_rot*fullelement->MyType->q_pos[j])*q;
											max_vectors[nr-1]+=(fullelement->center-lodelement->center+fullelement->rot*fullelement->MyType->q_pos[j])*q;
										}
										
										if(w<EPS) w=1.0;
										positions[nr-1]+=center*w;
										wsum+=w;
									}
								}
								positions[nr-1]/=wsum; // group center
//								vectors[nr-1]-=positions[nr-1]*sigma;
								vectors[nr-1]-=charge_center*sigma;
								if(vectors[nr-1]*vectors[nr-1]<EPS_CL){ // things get tricky here, no dipole, no unique direction for the whole group...
									// idea now is to use the z-semi-axis vector of the gyration tensor
									// Determine gyration tensor
									Vec3 points[6];
									Mat33 S;
									S.M3Zeros();
									Vec3 zdir(0.0);
									double gcount=0.0;
									for(i=0; i<nr_elements; i++){
										Element* curr;
										Vec3 center;
										Mat33 rot;
										if(group->levelofdetail<=0){
											curr=&params.configuration->group_elements[params.configuration->groups[params.trajectory[idx][frame].group_type]->elements[i]];
											center=params.trajectory[idx+i][frame].position-positions[nr-1];
											rot=AxisAngle2Rot(params.trajectory[idx+i][frame].rotation_vector);
										} else{
											unsigned int ellipsoid_idx=group->Type->LOD->element_in_ellipsoid[group->levelofdetail-1][i];
											Element* lodelement=&params.configuration->group_elements[group->elements[ellipsoid_idx]];
											curr=&params.configuration->group_elements[group->Type->LOD->groups[0]->elements[i]]; // underlying all-atom element
											
											// rotation difference for fully atomistic element
											Mat33 delta_rot=AxisAngle2Rot(params.trajectory[idx+ellipsoid_idx][frame].rotation_vector)*lodelement->rot.M3Transpose();
											Vec3 offset=delta_rot*(curr->center-lodelement->center);
											
											rot=delta_rot*curr->rot;
											center=params.trajectory[idx+ellipsoid_idx][frame].position+offset-positions[nr-1];
										}
										
										zdir+=rot.ColumnVec3(2);
										for(unsigned int l=0; l<3; l++){
											points[l<<1].vec[0]=(center.vec[0]+rot.mat[l][0]*curr->MyType->saxes.vec[0]);
											points[(l<<1)+1].vec[0]=(center.vec[0]-rot.mat[l][0]*curr->MyType->saxes.vec[0]);
											
											points[l<<1].vec[1]=(center.vec[1]+rot.mat[l][1]*curr->MyType->saxes.vec[1]);
											points[(l<<1)+1].vec[1]=(center.vec[1]-rot.mat[l][1]*curr->MyType->saxes.vec[1]);
											
											points[l<<1].vec[2]=(center.vec[2]+rot.mat[l][2]*curr->MyType->saxes.vec[2]);
											points[(l<<1)+1].vec[2]=(center.vec[2]-rot.mat[l][2]*curr->MyType->saxes.vec[2]);
										}
										for(unsigned int k=0; k<6; k++){
											for(unsigned int m=0; m<3; m++){
												for(unsigned int n=0; n<3; n++){
													// S_mn = sum_i r_m(i)*r_n(i) (corrected for center)
													S.mat[m][n]+=points[k].vec[m]*points[k].vec[n];
												}
											}
#if DEBUG_LEVEL>3
											cout << points[k].V3Str(',') << "\n";
#endif
											gcount+=1.0;
										}
#if DEBUG_LEVEL>3
										cout << "---\n";
#endif
									}
									S*=3.0/gcount; // since this gets square rooted later -> *3^(1/2) gives effective radius and principle axes (is needed to pass sanity check of single ellipsoid points)
									// get eigenvalues for ellipse semiaxis
									CVec3 eigenvalues=S.Eigenvalues();
									// in order for that to work, the eigenvalues need to be real
									if(eigenvalues.Im()*eigenvalues.Im()>EPS*EPS){
										cout << "Can only find imaginary ellipsoid semiaxes but we need something real here.\n";
										exit(2);
									}
									// semi axis are sqrts of the eigenvalues
									Vec3 semiaxis=eigenvalues.Re(); // lambda^2's
									Mat33 ev=S.Eigenvectors(semiaxis);
									Vec3 z=ev.ColumnVec3(2);
									vectors[nr-1]=z*(zdir*z); // make sure we always point in the same relative direction
								}
								max_vectors[nr-1]-=positions[nr-1]*sigma;
							}
						}
						idx+=params.configuration->groups[params.trajectory[idx][frame].group_type]->nr_elements; // works b/c group elements are grouped together (and we just hit the first one)
					} else idx++;
				}
			}
		} else{
			if(params.items[item_nr].type_idx>=0){ // special type
				if(params.items[item_nr].type_idx>=(int)nr_additional_subtypes){
					nr++;
					if(nr>total_nr){
						total_nr=nr;
						positions=(Vec3*)realloc(positions,sizeof(Vec3)*total_nr);
						vectors=(Vec3*)realloc(vectors,sizeof(Vec3)*total_nr);
						max_vectors=(Vec3*)realloc(max_vectors,sizeof(Vec3)*total_nr);
						rots=(Mat33*)realloc(rots,sizeof(Mat33)*total_nr);
						if((!positions || !vectors) || !max_vectors || !rots){
							cout << "Not enough memory.\n";
							exit(3);
						}
					}
					switch((unsigned int)(params.items[item_nr].type_idx-nr_additional_subtypes)){
						case 0: // x-Axis
							positions[nr-1]=Vec3(0.0);
							vectors[nr-1]=Vec3(1.0,0.0,0.0);
							break;
						case 1: // y-Axis
							positions[nr-1]=Vec3(0.0);
							vectors[nr-1]=Vec3(0.0,1.0,0.0);
							break;
						default:
						case 2: // z-Axis
							positions[nr-1]=Vec3(0.0);
							vectors[nr-1]=Vec3(0.0,0.0,1.0);
							break;
						case 3: // Efield
							positions[nr-1]=Vec3(0.0);
							vectors[nr-1]=Vec3(0.0,0.0,1.0); // in z-direction by default, except when rotated (or negated)
							if(params.configuration->Efield.V3Norm()>EPS){
								vectors[nr-1]=params.configuration->Efield;
								if((params.startstep>params.configuration->rotate_Efield_steps) && (params.configuration->rotate_Efield_steps>0)){
									vectors[nr-1]=params.rotatedE_direction;
								}
							}
							break;
						case 4: // totalM
							vectors[nr-1]=Vec3(0.0);
							// loop over all dipoles in system
							idx=0;
							while(idx<params.configuration->n_oids){
								e_t=params.configuration->element_types[params.trajectory[idx][frame].element_type];
								e_r=AxisAngle2Rot(params.trajectory[idx][frame].rotation_vector);
								positions[nr-1]=params.trajectory[idx][frame].position;
								if(params.configuration->offctrmu) positions[nr-1]+=e_r*e_t->mu_pos;
								vectors[nr-1]+=e_r*e_t->initial_dipole;// /e_t->initial_dipole.V3Norm();
								for(i=0; i<e_t->nr_charges; i++) vectors[nr-1]+=(e_r*e_t->q_pos[i]+positions[nr-1])*e_t->q[i];
								idx++;
							}
							positions[nr-1]=Vec3(0.0);
							break;
					}
					Vec3 zAxis(0.0,0.0,1.0);
					rots[nr-1]=RotAtoB(zAxis,vectors[nr-1]);
				} else{ // should not happen (subtypes are only allowed with associated group and we already took care of that)
					cout << "Somebody forgot to raise an error earlier: ERROR!\n";
					exit(42);
				}
			} else{ // single element out of trajectory file
				nr++;
				if(nr>total_nr){
					total_nr=nr;
					positions=(Vec3*)realloc(positions,sizeof(Vec3)*total_nr);
					vectors=(Vec3*)realloc(vectors,sizeof(Vec3)*total_nr);
					max_vectors=(Vec3*)realloc(max_vectors,sizeof(Vec3)*total_nr);
					rots=(Mat33*)realloc(rots,sizeof(Mat33)*total_nr);
					if((!positions || !vectors) || !max_vectors | !rots){
						cout << "Not enough memory.\n";
						exit(3);
					}
				}
				e_t=params.configuration->element_types[params.trajectory[params.items[item_nr].element_idx][frame].element_type];
				e_r=AxisAngle2Rot(params.trajectory[params.items[item_nr].element_idx][frame].rotation_vector);
				positions[nr-1]=params.trajectory[params.items[item_nr].element_idx][frame].position;
				if(params.configuration->offctrmu) positions[nr-1]+=e_r*e_t->mu_pos;
				vectors[nr-1]=e_r*e_t->initial_dipole;
				max_vectors[nr-1]=e_t->initial_dipole;
				rots[nr-1]=e_r;
				for(i=0; i<e_t->nr_charges; i++){
					q=e_t->q[i];
					vectors[nr-1]+=(e_r*e_t->q_pos[i])*q;
					max_vectors[nr-1]+=(e_t->q_pos[i])*q;
				}
				// return element's z-direction in case no dipoles exist
				if(vectors[nr-1]*vectors[nr-1]<EPS_CL) vectors[nr-1]=e_r*Vec3(0,0,1);
				if(max_vectors[nr-1]*max_vectors[nr-1]<EPS_CL) max_vectors[nr-1]=Vec3(0,0,1);
			}
		}
	}
	total_nr=nr;
	positions=(Vec3*)realloc(positions,sizeof(Vec3)*total_nr); // we're shrinking (in the worst case) here -- there's enough memory for that
	vectors=(Vec3*)realloc(vectors,sizeof(Vec3)*total_nr);
	max_vectors=(Vec3*)realloc(max_vectors,sizeof(Vec3)*total_nr);
	return nr;
}

unsigned int GetPosDipole(cfparams_t &params, unsigned int item_nr, unsigned int frame, Vec3* &positions, Vec3* &vectors, Vec3* &max_vectors, Mat33* &rots, unsigned int &total_nr)
{
	return GetPosDipole(params,item_nr,frame,positions,vectors,max_vectors,rots,total_nr,0,false);
}

inline double Vshell(double r, double dr)
{
	// 4/3*pi*((r+dr)^3-r^3) = 4/3*pi*(r*(r+dr)^2+dr*(r+dr)^2-r^3)
	// = 4/3*pi*(r*(r^2+dr^2+2*r*dr)+dr*(r^2+dr^2+2*r*dr)-r^3) = 4/3*pi*(r^3+r*dr^2+2*r^2*dr+dr*r^2+dr^3+2*r*dr^2-r^3)
	// = 4/3*pi*(3*r*dr^2+3*r^2*dr+dr^3) = 4*pi*r^2*dr+4*pi*r*dr^2+4/3*pi*dr^3 = 4*pi*(r^2*dr+r*dr^2+dr^3/3)
	return 4.0*pi*(r*r*dr+r*dr*dr+dr*dr*dr/3.0);
}

inline double Vvoxel(double r, double dr, double dcos, double dphi)
{
	// Vvoxel = int r^2*dr*sin(theta)dtheta*dphi = [1/3*r^3*cos(theta)*phi]_r^r+dr = (r^2*dr+r*dr^2+dr^3/3)*dcos*dphi
	return (r*r*dr+r*dr*dr+dr*dr*dr/3.0)*dcos*dphi;
}

