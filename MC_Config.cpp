/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

#include "MC_Config.h"
#define IN_MC_CONFIG
#include "MC_Elements.h"

/*!
Nothing to see here but a default constructor and destructor. Move along.
Not anymore, look, there's something going on :-) -- AT, Feb 8, 2011
*/
MC_Config::MC_Config()
{
	parameters.configfile="";
	parameters.Rand_Time=0.0;
	parameters.LJlambda=1.0;
	
	nr_group_elements=0;
	
	parameters.group_elements=NULL; // gets altered dynamically by realloc and hence needs to be NULL before
	max_nr_components=0;
	max_levelofdetail=0;
	max_LOD_elements=0;
	
	parameters.fit2lod=false;
}

MC_Config::~MC_Config()
{
}

/*!
 * Calculate group-internal "LJ dispersive energy" if elements after element nr <divider> are rotated around <axis> by <theta> at <location>
 * - ignores elements marked for deletion
 * - only calculates energy between element divider (one subgroup vs other subgroup)
 * - value is proportional and shifted compared to true energy, so only use for energy minimization ...
 */
double MC_Config::GroupLJdisp(Element_Group* group, unsigned int* deleted, unsigned int nr_deleted, bool* connection_site, unsigned int start, unsigned int divider, double theta, Vec3 axis, Vec3 location)
{
	double E=0.0;
	Mat33 rot=AxisAngle2Rot(axis,theta);
	for(unsigned int i=start; i<divider; i++){
		bool present=true;
		for(unsigned int k=0; k<nr_deleted; k++){
			if(i==deleted[k]){
				present=false;
				break;
			}
		}
		if(present){
			Element* A=&parameters.group_elements[group->elements[i]];
			for(unsigned int j=divider; j<group->nr_elements; j++){
				bool Bpresent=true;
				for(unsigned int k=0; k<nr_deleted; k++){
					if(j==deleted[k]){
						Bpresent=false;
						break;
					}
				}
				if(Bpresent){
					Element* B=&parameters.group_elements[group->elements[j]];
					double eps=sqrt(A->MyType->Vvdw*B->MyType->Vvdw);
					if(connection_site[j-divider]) eps*=8;
					double sigma=average(A->MyType->saxes.vec,3)+average(B->MyType->saxes.vec,3);
					double distance=(A->center-(location+rot*(B->center-location))).V3Norm();
					E+=eps*qqpwr(sigma/distance,parameters.LJexp[0]<<1); // bitshift left by one is multiplication by 2;
				}
			}
		}
	}
	return E;
}

Element_Type* MC_Config::NewType(Element_Type* source)
{
	// Need new element type for the new element
	parameters.num_element_types++;
	parameters.element_types = (Element_Type**)realloc(parameters.element_types,parameters.num_element_types*sizeof(Element_Type*));
	if(parameters.element_types==NULL){
		cout << "Could not find a memory block large enough to hold more element types.\n";
		exit(1);
	}
	parameters.element_types[parameters.num_element_types-1]=new Element_Type;
	Element_Type* element_type=parameters.element_types[parameters.num_element_types-1];
	
	element_type->thistype=parameters.num_element_types-1;
	
	element_type->number=source->number;
	element_type->dof=source->dof;
	element_type->archetype=source->archetype;
	element_type->saxes=source->saxes;
	element_type->saxes2=source->saxes2;
	element_type->invsaxes2=source->invsaxes2;
	element_type->invsaxesvolume=source->invsaxesvolume;
	element_type->mu_pos=source->mu_pos;
	element_type->dipole=source->dipole;
	element_type->initial_dipole=source->initial_dipole;
	element_type->Vvdw=source->Vvdw;
	element_type->sqrtVvdw=source->sqrtVvdw;
	element_type->inv_avg_area=source->inv_avg_area;
	element_type->sphericity=source->sphericity;
	element_type->IA_average=source->IA_average;
	element_type->avg_width=source->avg_width;
	element_type->mass=source->mass;
	element_type->MOI=source->MOI;
	element_type->hasmu=source->hasmu;
	element_type->still=source->still;
	element_type->rot_notrans=source->rot_notrans;
	element_type->calculate_order=source->calculate_order;
	element_type->color=source->color;
	element_type->transparency=source->transparency;
	element_type->label_elements=source->label_elements;
	element_type->texture_is_dipole=source->texture_is_dipole;
	element_type->texture=source->texture;
	element_type->name=source->name;
	element_type->nr_Vvdw_coefficients=source->nr_Vvdw_coefficients;
	if(element_type->nr_Vvdw_coefficients>0){
		element_type->Vvdw_coefficients=new double[element_type->nr_Vvdw_coefficients];
		for(unsigned int i=0; i<element_type->nr_Vvdw_coefficients; i++) element_type->Vvdw_coefficients[i]=source->Vvdw_coefficients[i];
	} else element_type->Vvdw_coefficients=NULL;
	if(source->IA_coefficients){
		element_type->IA_coefficients=new double[7];
		for(unsigned int i=0; i<7; i++) element_type->IA_coefficients[i]=source->IA_coefficients[i];
	} else element_type->IA_coefficients=NULL;
	element_type->nr_charges=source->nr_charges;
	if(element_type->nr_charges>0){
		element_type->q_pos=(Vec3*)malloc(element_type->nr_charges*sizeof(Vec3));
		element_type->q=(double*)malloc(element_type->nr_charges*sizeof(double));
		if(!element_type->q_pos || !element_type->q){
			cout << "ERROR: Not enough memory to store charges for type copy.\n";
			exit(11);
		}
		for(unsigned int i=0; i<element_type->nr_charges; i++){
			element_type->q_pos[i]=source->q_pos[i];
			element_type->q[i]=source->q[i];
		}
	} else{
		element_type->q_pos=NULL;
		element_type->q=NULL;
	}
	element_type->eps_texture=NULL;
	if(source->eps_texture){
		element_type->eps_texture=new double[phi_res*theta_res+4];
		for(unsigned int i=0; i<phi_res*theta_res+4; i++) element_type->eps_texture[i]=source->eps_texture[i];
	}
	element_type->r_diff=NULL;
	if(source->r_diff){
		element_type->r_diff=new double[phi_res*theta_res+2];
		for(unsigned int i=0; i<phi_res*theta_res+2; i++) element_type->r_diff[i]=source->r_diff[i];
	}
	element_type->width_diff=NULL;
	if(source->width_diff){
		element_type->width_diff=new double[phi_res*theta_res+2];
		for(unsigned int i=0; i<phi_res*theta_res+2; i++) element_type->width_diff[i]=source->width_diff[i];
	}
	
	return element_type;
}

void MC_Config::CopyElement(Element* source, Element* dest, int offset, bool copyproperties)
{
	if(copyproperties){
		dest->dipole=source->dipole;
		dest->rot=source->rot;
		dest->center=source->center;
		dest->group=NULL; // needs to be set later
		dest->q_pos=NULL;
		if(source->MyType){
			dest->MyType=NewType(source->MyType);
			dest->mytype=dest->MyType->thistype;
			if(dest->MyType->nr_charges>0){
				dest->q_pos=(Vec3*)malloc(dest->MyType->nr_charges*sizeof(Vec3));
				if(!dest->q_pos){
					cout << "ERROR: Can not allocate memory for charge position copy.\n";
					exit(11);
				}
			}
		} else{
			dest->mytype=0;
			dest->MyType=NULL;
		}
	}
	dest->nr_interactions=source->nr_interactions;
	dest->interactions=(Interaction*)malloc(dest->nr_interactions*sizeof(Interaction));
	if(!dest->interactions){
		cout << "ERROR: Can not allocate memory for interaction copy.\n";
		exit(11);
	}
	for(unsigned int i=0; i<source->nr_interactions; i++){
		dest->interactions[i].fixed=source->interactions[i].fixed;
		dest->interactions[i].allow_bond_stretch=source->interactions[i].allow_bond_stretch;
		dest->interactions[i].allow_bond_bend=source->interactions[i].allow_bond_bend;
		dest->interactions[i].bond_length=source->interactions[i].bond_length;
		dest->interactions[i].bond_order=source->interactions[i].bond_order;
		dest->interactions[i].nr_potentials=source->interactions[i].nr_potentials;
		if(dest->interactions[i].nr_potentials>0){
			unsigned int start_potentials=parameters.num_potentials;
			parameters.num_potentials+=dest->interactions[i].nr_potentials;
			parameters.potentials=(Interaction_Potential**)realloc(parameters.potentials,parameters.num_potentials*sizeof(Interaction_Potential*));
			if(!parameters.potentials){
				cout << "ERROR: Can not allocate memory for interaction potential copy.\n";
				exit(11);
			}
			dest->interactions[i].potentials=(Interaction_Potential**)malloc(dest->interactions[i].nr_potentials*sizeof(Interaction_Potential*));
			if(!dest->interactions[i].potentials){
				cout << "ERROR: Can not allocate memory for interaction potential copy.\n";
				exit(11);
			}
			for(unsigned int j=0; j<dest->interactions[i].nr_potentials; j++){
				parameters.potentials[start_potentials+j]=new Interaction_Potential;
				dest->interactions[i].potentials[j]=parameters.potentials[start_potentials+j];
				dest->interactions[i].potentials[j]->name=source->interactions[i].potentials[j]->name; // pointer
				dest->interactions[i].potentials[j]->type=source->interactions[i].potentials[j]->type;
				dest->interactions[i].potentials[j]->n=source->interactions[i].potentials[j]->n;
				dest->interactions[i].potentials[j]->nr_parameters=source->interactions[i].potentials[j]->nr_parameters;
				dest->interactions[i].potentials[j]->parameters=source->interactions[i].potentials[j]->parameters; // pointer
				dest->interactions[i].potentials[j]->links=source->interactions[i].potentials[j]->links; // pointer (stays the same, only relative partner numbers change)
				dest->interactions[i].potentials[j]->partners=new unsigned int[dest->interactions[i].potentials[j]->n];
				for(unsigned int k=0; k<dest->interactions[i].potentials[j]->n; k++){
					if((int)source->interactions[i].potentials[j]->partners[k]+offset>=0){
						dest->interactions[i].potentials[j]->partners[k]=(unsigned int)(source->interactions[i].potentials[j]->partners[k]+offset);
					} else{ // should not happen, but one never knows ...
						cout << "ERROR: Cannot copy element due to negative offset being too large.\n";
						exit(12);
					}
				}
			}
		}
		if((int)source->interactions[i].partner+offset>=0){
			dest->interactions[i].partner=(unsigned int)(source->interactions[i].partner+offset);
		} else{ // also should not happen, but better be safe ...
			cout << "ERROR: Cannot copy element due to negative offset being too large.\n";
			exit(12);
		}
		dest->interactions[i].back_link=source->interactions[i].back_link;
		
		dest->interactions[i].location=NULL;
		dest->interactions[i].initial_location=new Vec3;
		*dest->interactions[i].initial_location=*source->interactions[i].initial_location;
		if(source->interactions[i].location){
			dest->interactions[i].location=new Vec3;
			*dest->interactions[i].location=*source->interactions[i].location;
		}
		dest->interactions[i].normal=NULL;
		dest->interactions[i].initial_normal=new Vec3;
		*dest->interactions[i].initial_normal=*source->interactions[i].initial_normal;
		if(source->interactions[i].normal){
			dest->interactions[i].normal=new Vec3;
			*dest->interactions[i].normal=*source->interactions[i].normal;
		}
		dest->interactions[i].tangent=NULL;
		dest->interactions[i].initial_tangent=new Vec3;
		*dest->interactions[i].initial_tangent=*source->interactions[i].initial_tangent;
		if(source->interactions[i].tangent){
			dest->interactions[i].tangent=new Vec3;
			*dest->interactions[i].tangent=*source->interactions[i].location;
		}
	}
}

/*!
 * Adds group elements to another group and adjusts connections (also takes care of LOD meta-information)
 */
void MC_Config::AddGroup2Group(Element_Group* group, Element_Group* add_group, unsigned int basegroup_start)
{
	unsigned int nr_elements=group->nr_elements;
	unsigned int start_group_nr=nr_group_elements;
	
	group->nr_elements+=add_group->nr_elements;
	group->elements = (unsigned int*)realloc(group->elements,group->nr_elements*sizeof(unsigned int));
	if(group->elements==NULL){
		cout << "ERROR: Not enough memory to hold elements numbers of group " << group->Type->name << "\n";
		exit(1);
	}
	nr_group_elements+=add_group->nr_elements; // Adding one more element, increasing running total
	parameters.group_elements = (Element*)realloc(parameters.group_elements,nr_group_elements*sizeof(Element));
	if(parameters.group_elements==NULL){ // Danger, danger
		cout << "ERROR: Could not find a memory block large enough to hold a total of " << nr_group_elements << " group associated elements.\n";
		exit(1);
	}
	// for LOD level groups, obtain size of "element_groups" array (each ellipsoid has 1 entry for nr base elements plus the nr of base elements)
	unsigned int lod_eg_nr=0;
	if(group->levelofdetail>0){
		for(unsigned int i=0; i<group->Type->LOD->nr_ellipsoids[group->levelofdetail-1]; i++){
			lod_eg_nr+=group->Type->LOD->element_groups[group->levelofdetail-1][group->Type->LOD->ellipsoid_counts[group->levelofdetail-1][i]]+1;
		}
		group->Type->LOD->element_in_ellipsoid[group->levelofdetail-1]=(int*)realloc(group->Type->LOD->element_in_ellipsoid[group->levelofdetail-1],group->Type->LOD->groups[0]->nr_elements*sizeof(int));
		if(!group->Type->LOD->element_in_ellipsoid[group->levelofdetail-1]){
			cout << "ERROR: Not enough memory to store base model ellipsoid affiliations.\n";
			exit(6);
		}
	}
	for(unsigned int i=0; i<add_group->nr_elements; i++){
		group->elements[i+nr_elements]=start_group_nr+i;
		CopyElement(&parameters.group_elements[add_group->elements[i]],&parameters.group_elements[group->elements[i+nr_elements]],nr_elements);
		parameters.group_elements[group->elements[i+nr_elements]].group=group;
		if(group->levelofdetail>0){ // take care of LOD meta-information
			group->Type->LOD->nr_ellipsoids[group->levelofdetail-1]++; // we're adding one ellipsoid (=element) to it each iteration of this loop
			group->Type->LOD->ellipsoid_counts[group->levelofdetail-1]=(unsigned int*)realloc(group->Type->LOD->ellipsoid_counts[group->levelofdetail-1],group->Type->LOD->nr_ellipsoids[group->levelofdetail-1]*sizeof(unsigned int));
			if(!group->Type->LOD->ellipsoid_counts[group->levelofdetail-1]){
				cout << "ERROR: Out of memory adding an ellipsoid.\n";
				exit(7);
			}
			group->Type->LOD->ellipsoid_counts[group->levelofdetail-1][group->Type->LOD->nr_ellipsoids[group->levelofdetail-1]-1]=lod_eg_nr;
			group->Type->LOD->reduce_electrostatics[group->levelofdetail-1]=(double*)realloc(group->Type->LOD->reduce_electrostatics[group->levelofdetail-1],group->Type->LOD->nr_ellipsoids[group->levelofdetail-1]*sizeof(double));
			if(!group->Type->LOD->reduce_electrostatics[group->levelofdetail-1]){
				cout << "ERROR: Out of memory adding ellipsoid charge reduction parameters.\n";
				exit(8);
			}
			group->Type->LOD->reduce_electrostatics[group->levelofdetail-1][group->Type->LOD->nr_ellipsoids[group->levelofdetail-1]-1]=-1.0; // default is to keep original charges
			if(add_group->levelofdetail>0){ // group added has its own LOD meta-information
				group->Type->LOD->reduce_electrostatics[group->levelofdetail-1][group->Type->LOD->nr_ellipsoids[group->levelofdetail-1]-1]=add_group->Type->LOD->reduce_electrostatics[add_group->levelofdetail-1][i];
				unsigned int agec=add_group->Type->LOD->ellipsoid_counts[add_group->levelofdetail-1][i];
				unsigned int nr_base_elements=add_group->Type->LOD->element_groups[add_group->levelofdetail-1][agec];
				group->Type->LOD->element_groups[group->levelofdetail-1]=(int*)realloc(group->Type->LOD->element_groups[group->levelofdetail-1],(nr_base_elements+1+lod_eg_nr)*sizeof(int));
				if(!group->Type->LOD->element_groups[group->levelofdetail-1]){
					cout << "ERROR: Not enough memory to store base model dependencies.\n";
					exit(9);
				}
				group->Type->LOD->element_groups[group->levelofdetail-1][lod_eg_nr]=nr_base_elements;
				for(unsigned j=0; j<nr_base_elements; j++){
					unsigned int benr=add_group->Type->LOD->element_groups[add_group->levelofdetail-1][agec+j+1]+basegroup_start; // elements copied to basegroup are shifted
					group->Type->LOD->element_groups[group->levelofdetail-1][lod_eg_nr+j+1]=benr;
					group->Type->LOD->element_in_ellipsoid[group->levelofdetail-1][benr-1]=group->Type->LOD->nr_ellipsoids[group->levelofdetail-1]-1;
				}
				lod_eg_nr+=nr_base_elements+1;
			} else{ // just add single ellipsoids
				group->Type->LOD->element_groups[group->levelofdetail-1]=(int*)realloc(group->Type->LOD->element_groups[group->levelofdetail-1],(lod_eg_nr+2)*sizeof(int));
				if(!group->Type->LOD->element_groups[group->levelofdetail-1]){
					cout << "ERROR: Not enough memory to store base model dependencies.\n";
					exit(9);
				}
				group->Type->LOD->element_groups[group->levelofdetail-1][lod_eg_nr]=1;
				unsigned int benr=i+basegroup_start+1; // elements copied to basegroup are shifted
				group->Type->LOD->element_groups[group->levelofdetail-1][lod_eg_nr+1]=benr;
				group->Type->LOD->element_in_ellipsoid[group->levelofdetail-1][benr-1]=group->Type->LOD->nr_ellipsoids[group->levelofdetail-1]-1;
				lod_eg_nr+=2;
			}
		}
	}
}

char* MC_Config::IncludeFile(readfile* file, char* include_type, unsigned int &pos, unsigned int &include_position, string name, string &subname, bool combine_sections, string* inc_files)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	char* include=NULL;
	if(compare_strings(include_type,"CONF")) include=ConfigReader::IncludeFile(file,include_type,pos,include_position,name,subname,combine_sections,inc_files); // Base type handles config file including
		else{
			cout << "WARNING: Can not include configuration file \"" << file->filename << "\" as type <" << include_type << ">.\n";
		}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return include;
}

void MC_Config::SetSystemProperties()
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	cout << "Calculating system properties from configuration.\n";
	parameters.max_dipole=0.0;
	parameters.Mass=0.0;
	bool mass_defined=true;
	bool adjust_dipoleQ=false;
	if(fabs(parameters.dipoleQ)<EPS){
		parameters.dipoleQ=0.0;
		adjust_dipoleQ=true;
	}
	Mat33 ident;
	for(unsigned int i=0; i<parameters.num_element_types; i++){
		// set saxesbox while here
		parameters.element_types[i]->invsaxesvolume=1.0/(parameters.element_types[i]->saxes.vec[0]*parameters.element_types[i]->saxes.vec[1]*parameters.element_types[i]->saxes.vec[2]);
		parameters.element_types[i]->saxes2.vec[0]=parameters.element_types[i]->saxes.vec[0]*parameters.element_types[i]->saxes.vec[0];
		parameters.element_types[i]->saxes2.vec[1]=parameters.element_types[i]->saxes.vec[1]*parameters.element_types[i]->saxes.vec[1];
		parameters.element_types[i]->saxes2.vec[2]=parameters.element_types[i]->saxes.vec[2]*parameters.element_types[i]->saxes.vec[2];
		parameters.element_types[i]->invsaxes2.vec[0]=1.0/(parameters.element_types[i]->saxes.vec[0]*parameters.element_types[i]->saxes.vec[0]);
		parameters.element_types[i]->invsaxes2.vec[1]=1.0/(parameters.element_types[i]->saxes.vec[1]*parameters.element_types[i]->saxes.vec[1]);
		parameters.element_types[i]->invsaxes2.vec[2]=1.0/(parameters.element_types[i]->saxes.vec[2]*parameters.element_types[i]->saxes.vec[2]);
		parameters.element_types[i]->sqrtVvdw=sqrt(parameters.element_types[i]->Vvdw);
		// test sphere size is set to radius of sphere with same surface area
		parameters.element_types[i]->inv_avg_area=1.0;
		parameters.element_types[i]->IA_average=0.0;
		for(unsigned int ct=0; ct<theta_res*10; ct++){
			double cost=1.0-2.0*(ct+0.5)/(theta_res*10);
			double sint=sqrt(1.0-cost*cost);
			for(unsigned int j=0; j<phi_res*10; j++){
				double phi=2.0*pi*(j+0.5)/(phi_res*10);
				parameters.element_types[i]->IA_average+=IA0(parameters.element_types[i],Vec3(sint*cos(phi),sint*sin(phi),cost));
			}
		}
		parameters.element_types[i]->IA_average+=2.0*IA0(parameters.element_types[i],Vec3(0.0,0.0,1.0));
		parameters.element_types[i]->IA_average/=theta_res*phi_res*100+2;
//		parameters.element_types[i]->rT=sqrt(EllipsoidSurface(parameters.element_types[i]->saxes.vec[0],parameters.element_types[i]->saxes.vec[1],parameters.element_types[i]->saxes.vec[2])/(4.0*pi));
		parameters.element_types[i]->radius=sqrt(parameters.element_types[i]->IA_average/pi);
		if(isnan(parameters.element_types[i]->radius)) parameters.element_types[i]->radius=0.0; // Interaction Area uses inverse semi-axis which of course NaN's when zero axis are used
#if DEBUG_LEVEL>3
		cout << i << ": " << "sqrt(" << parameters.element_types[i]->IA_average << ") = " << parameters.element_types[i]->radius << "\n";
#endif
		parameters.element_types[i]->rT=parameters.element_types[i]->radius;
		if(parameters.element_types[i]->avg_width>0.0) parameters.element_types[i]->rT=parameters.element_types[i]->avg_width;
//		parameters.element_types[i]->inv_avg_area=1.0/(pi*qqpwr(parameters.element_types[i]->rT,2));
		parameters.element_types[i]->inv_avg_area/=parameters.element_types[i]->IA_average;
		if(isnan(parameters.element_types[i]->inv_avg_area)) parameters.element_types[i]->inv_avg_area=0.0;
		// total dipole strength of individual element types created
		Vec3 dipole=parameters.element_types[i]->initial_dipole;
		for(unsigned int j=0; j<parameters.element_types[i]->nr_charges; j++) dipole += parameters.element_types[i]->q_pos[j]*parameters.element_types[i]->q[j];
		parameters.max_dipole += (double)parameters.element_types[i]->number*dipole.V3Norm();
		if(adjust_dipoleQ){
			double R=EllipsoidRmin(dipole,parameters.element_types[i]->saxes,ident);
			double Q=dipole.V3Norm()/R;
			if(Q>parameters.dipoleQ) parameters.dipoleQ=Q;
		}
		if((parameters.element_types[i]->mass<=EPS) && (parameters.element_types[i]->number>0)) mass_defined=false;
		parameters.Mass += (double)parameters.element_types[i]->number*parameters.element_types[i]->mass;
	}
	// total dipole strength of individual groups
	Vec3 dipole;
	parameters.max_group_dipoles=new double[parameters.num_groups];
	parameters.group_max_centerdist=0.0;
	double maxradius=0.0;
	double d=0.0;
	for(unsigned int i=0; i<parameters.num_groups; i++){
		if((parameters.groups[i]->number>0) && (parameters.groups[i]->Type->Volume<EPS)){
#ifdef USE_CMWC4096
			init_CMWC4096(8);
			zigset();
#else
			__int32_t idum=-8;
#endif
			DetermineGroupVolume(parameters.groups[i],parameters.group_elements);
		}
		dipole.V3Zeros();
		parameters.max_group_dipoles[i]=0.0;
		Vec3 dist;
		if((parameters.groups[i]->number>0) && (parameters.group_radii[i]>maxradius)) maxradius=parameters.group_radii[i];
		for(unsigned int j=0; j<parameters.groups[i]->nr_elements; j++){
			Element* element = &parameters.group_elements[parameters.groups[i]->elements[j]];
			if(parameters.groups[i]->number>0){
				dist=parameters.group_centers[i]-element->center;
				d=dist*dist;
				if(d>parameters.group_max_centerdist) parameters.group_max_centerdist=d;
			}
			if(parameters.groups[i]->Type->group_dipole){
				dipole+=element->dipole;
			} else dipole=element->dipole;
			for(unsigned int k=0; k<element->MyType->nr_charges; k++){
				dipole+=(element->center+element->rot*element->MyType->q_pos[k])*element->MyType->q[k];
			}
			// if group_dipole is false, add individual dipole strength of each element to overall maximum dipole
			if(!parameters.groups[i]->Type->group_dipole) parameters.max_group_dipoles[i] += dipole.V3Norm();
			if((element->MyType->mass<=EPS) && (parameters.groups[i]->number>0)) mass_defined=false;
			parameters.Mass += (double)parameters.groups[i]->number*element->MyType->mass;
		}
		// if group_dipole is true, add dipole strength of whole group to overall maximum dipole
		if(parameters.groups[i]->Type->group_dipole) parameters.max_group_dipoles[i]=dipole.V3Norm();
		parameters.max_dipole += (double)parameters.groups[i]->number*parameters.max_group_dipoles[i];
		if(adjust_dipoleQ){
			double R=EllipsoidRmin(dipole,parameters.element_types[i]->saxes,ident);
			double Q=dipole.V3Norm()/R;
			if(Q>parameters.dipoleQ) parameters.dipoleQ=Q;
		}
	}
	if(parameters.dipoleQ>EPS) cout << "\t-> Charge for dipole expansion: " << parameters.dipoleQ/e_in_esu << " e\n";
	parameters.group_max_centerdist=sqrt(parameters.group_max_centerdist);
	if(maxradius>parameters.group_max_centerdist) parameters.group_max_centerdist=maxradius;
	parameters.max_neighbor_trans=parameters.maxtrans*2*NEIGHBOR_UPDATE_FREQ+parameters.group_max_centerdist;
	cout << "\t-> Maximum translation of all elements and groups: " << parameters.max_neighbor_trans << " Angström (contains group contribution of " << parameters.group_max_centerdist << " Angström)\n";
	if(!mass_defined && (parameters.Mass>EPS)){
		cout << "WARNING: Overall system mass undefined because some elements contributing to the system do not have an associated mass.\n";
		parameters.Mass=0.0;
	}
	cout << "\t-> Maximum total dipole moment is " << parameters.max_dipole << " Debye\n";
	if(parameters.Mass>EPS){
		cout << "\t-> Overall system mass is " << parameters.Mass << " g/mol";
		if(!parameters.NpT) cout << " (target density: " << parameters.Mass/(NA*1E-24*parameters.targetV) << " g/cm^3)";
		cout << "\n";
	}
	cout << "<- Done.\n";
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void MC_Config::SetSystemVolume()
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	cout << "Adjusting simulation volume based on configuration:\n";
	// Calculate volume for sphere, cylinder, and rectangular boxes, respectively;
	if (parameters.latticetype == 3) parameters.V = (4.0/3.0)*pi*qqpwr(parameters.spherecylr,3);
		else if (parameters.latticetype == 4) parameters.V = 2.0*pi*qqpwr(parameters.spherecylr,2)*parameters.boxlength[2];
			else parameters.V=parameters.boxlength[0]*parameters.boxlength[1]*parameters.boxlength[2];
	
	// Check to make sure box has finite volume
	if(parameters.V<0.0){
		cout << "Box volume must be positive.\n";
		exit(1);
	}
	
	parameters.targetV=parameters.V; // user provided value is target volume
	
	// Calculate overall VdW volume of things in system
	double VdW_volume=0.0;
	// elements first
	for(unsigned int i=0; i<parameters.num_element_types; i++){
		double V=parameters.element_types[i]->number*parameters.element_types[i]->saxes.vec[0]*parameters.element_types[i]->saxes.vec[1]*parameters.element_types[i]->saxes.vec[2];
		VdW_volume+=V;
	}
	// groups next
	for(unsigned int i=0; i<parameters.num_groups; i++){
		for(unsigned int j=0; j<parameters.groups[i]->nr_elements; j++){
			double V=parameters.groups[i]->number*parameters.group_elements[parameters.groups[i]->elements[j]].MyType->saxes.vec[0]*parameters.group_elements[parameters.groups[i]->elements[j]].MyType->saxes.vec[1]*parameters.group_elements[parameters.groups[i]->elements[j]].MyType->saxes.vec[2];
			VdW_volume+=V;
		}
	}
	
	VdW_volume*=4.0/3.0*pi;
#if DEBUG_LEVEL>1
		cout << "-> VdW volume of elements in simulation volume: " << VdW_volume << " Angström³ (target simulation volume: " << parameters.targetV << " Angström³)\n";
#endif
	if(VdW_volume>1000000*parameters.targetV){
		cout << "Creation-of-universe-simulation not supported on current hardware.\n";
		exit(42);
	} else{
		if(VdW_volume>1000*parameters.targetV){
			cout << "WARNING: Attempting to simulate black hole.\n";
		} else{
			if(VdW_volume>2*parameters.targetV){
				cout << "WARNING: Attempting to simulate white dwarf.\n";
			} else if(VdW_volume>parameters.targetV){
					cout << "\nWARNING: VdW volume of simulated elements is larger than simulation volume.\n\n";
				}
		}
	}
	
	unsigned int nr_elements=parameters.n_oids-parameters.n_group_oids;
	
	unsigned int largest;
	double box[3];
	
	if(parameters.n_group_oids>0){ // Houston control, we have groups
		// First, find out how big the problem is going to be (literally)
#if DEBUG_LEVEL>0
		cout << "Arranging " << parameters.n_groups << " groups and " << nr_elements << " elements in " << parameters.V << " Angström³ simulation volume ...\n";
#endif
		double effective_group_volume=0.0;
		double V;
		parameters.rmax_g=0.0;
#if DEBUG_LEVEL>1
		cout << "-> Groups:\n";
#endif
		for(unsigned int i=0; i<parameters.num_groups; i++){
			if(parameters.groups[i]->number>0){
				V=parameters.group_radii[i];
				if(V>parameters.rmax_g) parameters.rmax_g=V;
				V=parameters.groups[i]->number*qqpwr(V,3);
				V*=4.0*pi/3.0;
				effective_group_volume+=V;
#if DEBUG_LEVEL>1
				cout << "\tOverall bounding sphere volume of " << parameters.groups[i]->number << " <" << parameters.groups[i]->Type->name << ">: " << V << " Angström³\n";
#endif
			}
		}
		parameters.a_group=sqrt(8.0)*parameters.rmax_g;
#if DEBUG_LEVEL>1
		cout << "-> Total volume of group bounding spheres: " << effective_group_volume << " Angström³\n";
		cout << "-> Maximum group bounding sphere radius: " << parameters.rmax_g << " Angström\n";
		cout << "=> fcc lattice constant a: " << parameters.a_group << " Angström\n";
#endif
		// Now adjust fcc unit cell number to fit in provided space
		unsigned int need_fcc_cells = parameters.n_groups/4;
		if(parameters.n_groups%4>0) need_fcc_cells+=1; // cannot distribute fractional fcc cells
		// number of fcc cells fitting in user defined box
		parameters.fcc[0] = (unsigned int)(parameters.boxlength[0]/parameters.a_group);
		parameters.fcc[1] = (unsigned int)(parameters.boxlength[1]/parameters.a_group);
		parameters.fcc[2] = (unsigned int)(parameters.boxlength[2]/parameters.a_group);
		// simulation volume needs to be symetrically spaced with unit cells (lin_dim takes care of that)
		largest=0;
		for(unsigned int i=0; i<3; i++){
			box[i]=parameters.boxlength[i]; // initialize with current size
			if(box[i]>box[largest]) largest=i;
		}
		parameters.nr_fcc_cells = parameters.fcc[0]*parameters.fcc[1]*parameters.fcc[2];
		if(need_fcc_cells>parameters.nr_fcc_cells){ // more cells needed than space available
#if DEBUG_LEVEL>0
			cout << "Need " << need_fcc_cells << " fcc unit cells but specified system volume only holds " << parameters.nr_fcc_cells << "\n";
#endif
			if(!parameters.auto_volume) exit(3);
#if DEBUG_LEVEL>0
			cout << "-> Adjusting system volume, specified volume will be reached during randomization.\n";
#endif
			// initialize with current volume shrinked to incorporate currently fitting fcc elements
			for(unsigned int i=0; i<3; i++) box[i]=(double)parameters.fcc[i]*parameters.a_group;
			while(parameters.nr_fcc_cells<need_fcc_cells){ // keep volume shape identical while growing
				parameters.nr_fcc_cells=1;
				for(unsigned int i=0; i<3; i++){ // scale up by adding one fcc unit to the largest dimension each round
					box[i]+=parameters.a_group*parameters.boxlength[i]/parameters.boxlength[largest];
					parameters.fcc[i]=(unsigned int)(box[i]/parameters.a_group);
					parameters.nr_fcc_cells*=parameters.fcc[i];
				}
			}
			parameters.boxlength[0]=box[0]; parameters.boxlength[1]=box[1]; parameters.boxlength[2]=box[2];
			parameters.inv_boxlength[0]=1.0/box[0]; parameters.inv_boxlength[1]=1.0/box[1]; parameters.inv_boxlength[2]=1.0/box[2];
			parameters.V=box[0]*box[1]*box[2]; // still need to check if regular elements will fit ...
#if DEBUG_LEVEL>1
			cout << "-> System volume adjusted to " << parameters.V << " Angström³.\n";;
#endif
		} else{ // also adjust if too many unit cells ...
			unsigned int cells[3];
			unsigned int solution[3];
			largest=0;
			for(unsigned int i=0; i<3; i++){
				cells[i]=parameters.fcc[i]; // initialize with current size
				solution[i]=parameters.fcc[i];
				if(parameters.fcc[i]>parameters.fcc[largest]) largest=i;
			}
			unsigned int minus=0;
			while(cells[0]*cells[1]*cells[2]>need_fcc_cells){
				minus++; // slightly paradox ...
				for(unsigned int i=0; i<3; i++) cells[i]=parameters.fcc[i]-(unsigned int)(minus*double(parameters.fcc[i]/parameters.fcc[largest])); // take of from largest first ...
				if(cells[0]*cells[1]*cells[2]>=need_fcc_cells) for(unsigned int i=0; i<3; i++) solution[i]=cells[i];
			}
			parameters.nr_fcc_cells=1;
			for(unsigned int i=0; i<3; i++){
				parameters.fcc[i]=solution[i];
				parameters.nr_fcc_cells*=parameters.fcc[i];
			}
		}
		double fcc_volume=parameters.nr_fcc_cells*qqpwr(parameters.a_group,3);
		// initial fcc unit cell scaling is 1.0 in each direction
		parameters.scale[0]=1.0; parameters.scale[1]=1.0; parameters.scale[2]=1.0;
		unsigned int min_scale=0;
		if(parameters.V>fcc_volume){ // available volume bigger than minimally needed for groups - adjust and use :-)
#if DEBUG_LEVEL>1
			cout << "Using available system volume to equally distribute groups:\n";
#endif
			fcc_volume=1.0; // multiplication's neutral element
			for(unsigned int i=0; i<3; i++){
				parameters.scale[i]=box[i]/(double(parameters.fcc[i])*parameters.a_group);
				fcc_volume*=parameters.scale[i];
				if(parameters.scale[i]<parameters.scale[min_scale]) min_scale=i;
			}
#if DEBUG_LEVEL>1
			cout << "=> scale fcc lattice constant by (" << parameters.scale[0] << ", " << parameters.scale[1] << ", " << parameters.scale[2] << ")\n";
#endif
			fcc_volume*=parameters.nr_fcc_cells*qqpwr(parameters.a_group,3); // scaled volume: "*=" is not "="
		}
#if DEBUG_LEVEL>1
		cout << "<- Total volume of " << parameters.nr_fcc_cells << " scaled fcc unit cells: " << fcc_volume << " Angström³ (available volume)\n";
#endif
		// calculate empty radius in group fcc (12*1/4 = 3 such spaces are available ...)
		for(unsigned int i=0; i<3; i++) parameters.c_group[i]=(parameters.scale[i]*parameters.a_group/2.0)-parameters.rmax_g; // scaled versions
#if DEBUG_LEVEL>1
		cout << "Available empty space per scaled fcc unit: [" << parameters.c_group[0] << ", " << parameters.c_group[1] << ", " << parameters.c_group[2] << "] Angström (corresponds to volume: " << 3*parameters.n_groups*parameters.c_group[0]*parameters.c_group[1]*parameters.c_group[2] << " Angström³)\n";
#endif
		// Find out how much space all elements rotation unperturbed occupy (largest ellipsoid axis determines bounding sphere)
		parameters.rmax_e=0.0;
		double effective_element_volume=0.0;
#if DEBUG_LEVEL>1
		cout << "-> Elements:\n";
#endif
		for(unsigned int i=0; i<parameters.num_element_types; i++){
			if(parameters.element_types[i]->number>0){
				V=maxd(parameters.element_types[i]->saxes.vec,3)*two_to_onesixth; // use semi-axis and convert maximum semi-axis to rmin
				if(V>parameters.rmax_e) parameters.rmax_e=V;
				V=parameters.element_types[i]->number*qqpwr(V,3);
				V*=4.0*pi/3.0;
				effective_element_volume+=V;
#if DEBUG_LEVEL>1
				cout << "\tOverall bounding sphere volume of " << parameters.element_types[i]->number << " <" << parameters.element_types[i]->name << ">: " << V << " Angström³\n";
#endif
			}
		}
#if DEBUG_LEVEL>1
		cout << "-> Total volume of element bounding spheres: " << effective_element_volume << " Angström³\n";
		cout << "-> Maximum effective element type radius: " << parameters.rmax_e << " Angström\n";
#endif
		parameters.elements_per_fcc=3;
		for(unsigned int i=0; i<3; i++){
			parameters.sc[i]=(unsigned int)(parameters.c_group[i]/parameters.rmax_e);
			parameters.elements_per_fcc*=parameters.sc[i];
		}
		if(parameters.elements_per_fcc*parameters.nr_fcc_cells<nr_elements){ // not enough space to squeeze in elements or not enough spots
#if DEBUG_LEVEL>0
			cout << "Not enough space left in fcc lattice for regular elements";
#endif
#if DEBUG_LEVEL>1
			cout << " (space for " << parameters.elements_per_fcc << " elements per unit cell).";
#endif
#if DEBUG_LEVEL>0
			cout << "\n";
#endif
			if(!parameters.auto_volume) exit(3);
#if DEBUG_LEVEL>0
			if(!parameters.NpT) cout << "-> Adjusting system volume, specified volume will be reached during randomization.\n";
#endif
			while(parameters.elements_per_fcc*parameters.nr_fcc_cells<nr_elements){
				parameters.elements_per_fcc=3;
				fcc_volume=1.0; // multiplication's neutral element
				for(unsigned int i=0; i<3; i++){
					box[i]*=1.001; // each round increase system volume by 0.1% in each dimension (to keep volume shape constant)
					parameters.scale[i]=box[i]/(double(parameters.fcc[i])*parameters.a_group);
					fcc_volume*=parameters.scale[i];
					parameters.c_group[i]=(parameters.scale[i]*parameters.a_group/2.0)-parameters.rmax_g;
					parameters.sc[i]=(unsigned int)(parameters.c_group[i]/parameters.rmax_e);
					parameters.elements_per_fcc*=parameters.sc[i];
				}
				fcc_volume*=parameters.nr_fcc_cells*qqpwr(parameters.a_group,3); // scaled volume: "*=" is not "="
			}
			update_volume(&parameters,fcc_volume);
			parameters.nndist=cbrt(parameters.V/parameters.n_oids);
#if DEBUG_LEVEL>1
			cout << "-> System volume adjusted to " << parameters.V << " Angström³ (space for " << parameters.elements_per_fcc << " elements per unit cell).\n";
			cout << "=> scale fcc lattice constant by: (" << parameters.scale[0] << ", " << parameters.scale[1] << ", " << parameters.scale[2] << ")\n";
			cout << "Available empty space per scaled fcc unit: [" << parameters.c_group[0] << ", " << parameters.c_group[1] << ", " << parameters.c_group[2] << "] Angström (corresponds to volume: " << 3*parameters.n_groups*parameters.c_group[0]*parameters.c_group[1]*parameters.c_group[2] << " Angström³)\n";
#endif
		}
		if(parameters.LJwall_calc){
			double additional=((parameters.fcc[0]+1)*(parameters.fcc[1]+1)*(parameters.fcc[2]+1)-parameters.fcc[0]*parameters.fcc[1]*parameters.fcc[2])*qqpwr(parameters.a_group,3);
			cout << "-> Making additional room (+" << additional << " Angström³) so freely rotated groups and elements cannot enter the LJ wall.\n";
			update_volume(&parameters,parameters.V+additional);
			parameters.nndist=cbrt(parameters.V/parameters.n_oids);
		}
	} else{ // placing only elements
#if DEBUG_LEVEL>0
		cout << "Arranging " << nr_elements << " elements in " << parameters.V << " Angström³ simulation volume ...\n";
#endif
		parameters.rmax_e=0.0;
		double effective_element_volume=0.0;
#if DEBUG_LEVEL>1
		cout << "-> Elements:\n";
#endif
		double V;
		for(unsigned int i=0; i<parameters.num_element_types; i++){
			if(parameters.element_types[i]->number>0){
				V=maxd(parameters.element_types[i]->saxes.vec,3)*two_to_onesixth; // use semi-axis and convert maximum semi-axis to rmin
				if(V>parameters.rmax_e) parameters.rmax_e=V;
				V=parameters.element_types[i]->number*qqpwr(V,3);
				V*=4.0*pi/3.0;
				effective_element_volume+=V;
#if DEBUG_LEVEL>1
				cout << "\tOverall bounding sphere volume of " << parameters.element_types[i]->number << " <" << parameters.element_types[i]->name << ">: " << V << " Angström³\n";
#endif
			}
		}
#if DEBUG_LEVEL>1
		cout << "-> Total volume of element bounding spheres: " << effective_element_volume << " Angström³\n";
		cout << "-> Maximum effective element type radius: " << parameters.rmax_e << " Angström\n";
#endif
		largest=0;
		for(unsigned int i=0; i<3; i++){
			box[i]=parameters.boxlength[i]; // initialize with current size
			if(box[i]>box[largest]) largest=i;
		}
		parameters.a_element=2.0*parameters.rmax_e;
		for(unsigned int i=0; i<3; i++) parameters.sc[i]=(unsigned int)(parameters.boxlength[i]/parameters.a_element);
		parameters.nr_sc_cells=parameters.sc[0]*parameters.sc[1]*parameters.sc[2];
		if(parameters.nr_sc_cells<nr_elements){
#if DEBUG_LEVEL>0
			cout << "Not enough space in specified system volume to place all elements.\n";
#endif
			if(parameters.auto_volume){
#if DEBUG_LEVEL>0
				cout << "-> Adjusting system volume, specified volume will be reached during randomization.\n";
#endif
				for(unsigned int i=0; i<3; i++) box[i]=(double)parameters.sc[i]*parameters.a_element; // initialize with current size
				while(parameters.nr_sc_cells<nr_elements){ // keep volume shape identical while growing
					parameters.nr_sc_cells=1;
					for(unsigned int i=0; i<3; i++){ // scale up by adding one sc unit to the largest dimension each round
						box[i]+=parameters.a_element*parameters.boxlength[i]/parameters.boxlength[largest];
						parameters.sc[i] = (unsigned int)(box[i]/parameters.a_element);
						parameters.nr_sc_cells*=parameters.sc[i];
					}
				}
				update_volume(&parameters,box[0]*box[1]*box[2]);
				parameters.nndist=cbrt(parameters.V/parameters.n_oids);
#if DEBUG_LEVEL>1
				cout << "-> System volume adjusted to " << parameters.V << " Angström³.\n";
#endif
			} else{
#if DEBUG_LEVEL>0
				cout << "-> WARNING: Adjusting sc lattice constant, element overlap may occur.\n";
#endif
				// Adjust volume first to find out how many sc unit cells in each direction are needed to maintain aspect ratios
				for(unsigned int i=0; i<3; i++) box[i]=(double)parameters.sc[i]*parameters.a_element; // initialize with current size
				while(parameters.nr_sc_cells<nr_elements){ // keep volume shape identical while growing
					parameters.nr_sc_cells=1;
					for(unsigned int i=0; i<3; i++){ // scale up by adding one sc unit to the largest dimension each round
						box[i]+=parameters.a_element*parameters.boxlength[i]/parameters.boxlength[largest];
						parameters.sc[i] = (unsigned int)(box[i]/parameters.a_element);
						parameters.nr_sc_cells*=parameters.sc[i];
					}
				}
				// Calculate new lattice constant
				parameters.a_element=parameters.boxlength[largest]/parameters.sc[largest];
				parameters.rmax_e=parameters.a_element/2.0;
				parameters.nndist=parameters.a_element;
#if DEBUG_LEVEL>1
				cout << "-> sc lattice constant adjusted to " << parameters.a_element << " Angström.\n";
#endif
			}
		} else{ // also adjust if too many unit cells ...
			unsigned int cells[3];
			unsigned int solution[3];
			largest=0;
			for(unsigned int i=0; i<3; i++){
				cells[i]=parameters.sc[i]; // initialize with current size
				solution[i]=parameters.sc[i];
				if(parameters.sc[i]>parameters.sc[largest]) largest=i;
			}
			unsigned int minus=0;
			while(cells[0]*cells[1]*cells[2]>nr_elements){
				minus++; // slightly paradox ...
				for(unsigned int i=0; i<3; i++) cells[i]=parameters.sc[i]-(unsigned int)(minus*double(parameters.sc[i]/parameters.sc[largest])); // take of from largest first ...
				if(cells[0]*cells[1]*cells[2]>=nr_elements) for(unsigned int i=0; i<3; i++) solution[i]=cells[i];
			}
			
			parameters.nr_sc_cells=1;
			for(unsigned int i=0; i<3; i++){
				parameters.sc[i]=solution[i];
				parameters.nr_sc_cells*=parameters.sc[i];
			}
		}
		// initial unit cell scaling is 1.0 in each direction
		parameters.scale[0]=1.0; parameters.scale[1]=1.0; parameters.scale[2]=1.0;
		if(parameters.V>parameters.nr_sc_cells*qqpwr(parameters.a_element,3)){ // available volume bigger than minimally needed for elements - adjust and use :-)
#if DEBUG_LEVEL>1
			cout << "Using available system volume to equally distribute elements:\n";
#endif
			for(unsigned int i=0; i<3; i++) parameters.scale[i]=parameters.boxlength[i]/(double(parameters.sc[i])*parameters.a_element);
#if DEBUG_LEVEL>1
			cout << "=> scale sc lattice constant by (" << parameters.scale[0] << ", " << parameters.scale[1] << ", " << parameters.scale[2] << ")\n";
#endif
		}
		if(parameters.LJwall_calc){
			double additional=((parameters.sc[0]+1)*(parameters.sc[1]+1)*(parameters.sc[2]+1)-parameters.sc[0]*parameters.sc[1]*parameters.sc[2])*qqpwr(parameters.a_element,3);
			cout << "-> Making additional room (+" << additional << " Angström³) so freely rotated elements cannot enter the LJ wall.\n";
			update_volume(&parameters,parameters.V+additional);
			parameters.nndist=cbrt(parameters.V/parameters.n_oids);
		}
	}
	cout << "<- Done, simulation volume is " << parameters.V << " Angström³";
	if((fabs(parameters.V-parameters.targetV)>EPS) && !parameters.NpT) cout << " (will be adjusted to: " << parameters.targetV << " Angström³)."; else cout << ".";
	cout << "\n";
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void MC_Config::AssignFixedElements(Element_Group* group)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	group->Type->nr_movable=group->nr_elements;
	group->Type->fixed_elements=new bool[group->nr_elements];
	for(unsigned int i=0; i<group->nr_elements; i++){
		Element* element=&parameters.group_elements[group->elements[i]];
		group->Type->fixed_elements[i]=false;
		bool fixed=true;
		// element in group is fixed if all interactions to and with it are fixed
		for(unsigned int j=0; j<element->nr_interactions; j++) fixed&=element->interactions[j].fixed;
		if(fixed) group->Type->nr_movable--;
		group->Type->fixed_elements[i]=fixed;
	}
#if DEBUG_LEVEL>1
	cout << "-> Found " << group->nr_elements-group->Type->nr_movable << " elements with fixed bonds (" << group->Type->nr_movable << " movable elements in group)\n";
#endif
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void MC_Config::AutoPotentials(Element_Group* group)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	cout << "Automatically assigning bond potentials for group <" << group->Type->name << "> ...\n";
	unsigned int count=0;
	unsigned int bondcount=0;
	unsigned int nr_statics=0;
	
	for(unsigned int i=0; i<group->nr_elements; i++){
		Element* element=&parameters.group_elements[group->elements[i]];
		for(unsigned int j=0; j<element->nr_interactions; j++){
			bool hasbondstretchpotential=false;
			bool hasbondbendpotential=false;
			bool hasspringpotential=false;
			Interaction* interaction=&element->interactions[j];
			// If this is an LOD group take care of fully-atomistic potentials which may exist for this link
			if(group->levelofdetail>0){
				if(group->Type->LOD->keep_original_potentials[group->levelofdetail-1]){
					Element_Group* basegroup=group->Type->LOD->groups[0];
					Vec3 link1_pos=element->center+element->rot*(*interaction->initial_location);
#if DEBUG_LEVEL>2
					cout << "link1_pos = (" << link1_pos.V3Str(',') << "\n";
#endif
					Element* partner=&parameters.group_elements[group->elements[interaction->partner]];
					Vec3 link2_pos=partner->center+partner->rot*(*(partner->interactions[interaction->back_link].initial_location));
#if DEBUG_LEVEL>2
					cout << "link2_pos = (" << link2_pos.V3Str(',') << "\n";
#endif
					// find associated links
					unsigned int partner_id=interaction->partner;
					for(unsigned int k=0; k<basegroup->nr_potentials; k++){
						Interaction_Potential* orig_potential=basegroup->potentials[k];
#if DEBUG_LEVEL>2
						cout << "potential: " << k << " (" << orig_potential->n << ")\n";
#endif
						unsigned int match_A=0;
						unsigned int match_B=0;
						unsigned int link1=0;
						unsigned int link2=0;
						int link_id=0;
						for(unsigned int l=0; l<orig_potential->n; l++){
							unsigned int e_idx=group->Type->LOD->element_in_ellipsoid[group->levelofdetail-1][orig_potential->partners[l]];
							Element* p=&parameters.group_elements[basegroup->elements[orig_potential->partners[l]]];
#if DEBUG_LEVEL>2
							cout << "partner: " << orig_potential->partners[l] << "\n";
#endif
							for(unsigned int m=0; m<orig_potential->n; m++){
								unsigned int link=orig_potential->links[l*orig_potential->n+m];
#if DEBUG_LEVEL>2
								cout << "\t-> link: " << link << " (" << p->nr_interactions << ")\n";
#endif
								if(l!=m){ // links to self not allowed ...
									Vec3 pos=p->center+p->rot*(*p->interactions[link].initial_location);
#if DEBUG_LEVEL>2
									cout << pos.V3Str(',') << "\n";
#endif
									Vec3 delta=link1_pos-pos;
									if(delta.V3Norm()<EPS){
#if DEBUG_LEVEL>2
										cout << "l1\n";
#endif
										link1=l+1;
										if(link_id==0) link_id=link+1;
										if(link2!=0) break;
									}
									delta=link2_pos-pos;
									if(delta.V3Norm()<EPS){
#if DEBUG_LEVEL>2
										cout << "l2\n";
#endif
										link2=l+1;
										if(link_id==0) link_id=-(link+1);
										if(link1!=0) break;
									}
								}
							}
							match_A+=(e_idx==i);
							match_B+=(e_idx==partner_id);
						}
						if((match_A && match_B) && (link_id!=0)){ // we have a potential match
							if(i<partner_id){ // only set up once
								Interaction_Potential* new_potential;
								unsigned int shift=0;
								switch(orig_potential->n){ // for the moment focus on two-partner potentials only
									case 2:	cout << "\t-> Setting up potential between " << i+1 << " and " << partner_id+1 << "\n";
										// create new potential
										parameters.num_potentials++;
										parameters.potentials=(Interaction_Potential**)realloc(parameters.potentials,parameters.num_potentials*sizeof(Interaction_Potential*));
										if(!parameters.potentials){
											cout << "ERROR: Not enough memory to create additional potential.\n";
											exit(11);
										}
										new_potential=new Interaction_Potential;
										new_potential->name=orig_potential->name;
										new_potential->type=orig_potential->type;
										new_potential->n=orig_potential->n;
										new_potential->nr_parameters=orig_potential->nr_parameters;
										new_potential->parameters=orig_potential->parameters; // copy just pointer here (parameters don't change)
										new_potential->partners=new unsigned int[2];
										if(link_id<0) shift=1;
										new_potential->partners[shift]=i;
										new_potential->partners[(1+shift)%2]=partner_id;
										new_potential->links=new unsigned int[4];
										new_potential->links[0]=0; new_potential->links[3]=0;
										new_potential->links[shift*2+(1+shift)%2]=j; // link nr from CA->CB
										new_potential->links[((1+shift)%2)*2+shift]=interaction->back_link; // link nr from CB->CA
										// Store new potential
										parameters.potentials[parameters.num_potentials-1]=new_potential;
										// Attach to group
										group->nr_potentials++;
										group->potentials=(Interaction_Potential**)realloc(group->potentials,group->nr_potentials*sizeof(Interaction_Potential*));
										if(!group->potentials){
											cout << "ERROR: Not enough memory to create additional potential.\n";
											exit(11);
										}
										group->potentials[group->nr_potentials-1]=new_potential;
										// Attach to both links (forward and backward)
										if(interaction->nr_potentials==0) interaction->potentials=NULL;
										interaction->nr_potentials++;
										interaction->potentials=(Interaction_Potential**)realloc(interaction->potentials,interaction->nr_potentials*sizeof(Interaction_Potential*));
										if(!interaction->potentials){
											cout << "ERROR: Not enough memory to create additional potential.\n";
											exit(11);
										}
										interaction->potentials[interaction->nr_potentials-1]=new_potential;
										if(partner->interactions[interaction->back_link].nr_potentials==0) partner->interactions[interaction->back_link].potentials=NULL;
										partner->interactions[interaction->back_link].nr_potentials++;
										partner->interactions[interaction->back_link].potentials=(Interaction_Potential**)realloc(partner->interactions[interaction->back_link].potentials,partner->interactions[interaction->back_link].nr_potentials*sizeof(Interaction_Potential*));
										if(!partner->interactions[interaction->back_link].potentials){
											cout << "ERROR: Not enough memory to create additional potential.\n";
											exit(11);
										}
										partner->interactions[interaction->back_link].potentials[partner->interactions[interaction->back_link].nr_potentials-1]=new_potential;
										break;
									default:
										cout << "\t\t-> " << orig_potential->n << " partner potential transcription not implemented currently.\n";
								}
							}
						}
					}
				}
			}
			for(unsigned int k=0; k<interaction->nr_potentials; k++){
				switch(interaction->potentials[k]->type){
					case stretch_potential:
								hasbondstretchpotential=true;
								break;
					case bend_potential:
								hasbondbendpotential=true;
								break;
					case chainspring_potential:
					case spring_potential:
								hasspringpotential=true;
								break;
				}
			}
			if(interaction->fixed){
#if DEBUG_LEVEL>2
				cout << "\t-> Static link between: " << i+1 << " and " << interaction->partner+1 << "\n";
#endif
				nr_statics++;
			}
			bondcount++;
			if(hasbondstretchpotential){
				interaction->allow_bond_stretch=true;
			} else{
				interaction->allow_bond_stretch=false;
				interaction->bond_length=(element->center+(element->rot*(*interaction->initial_location))-parameters.group_elements[group->elements[interaction->partner]].center-(parameters.group_elements[group->elements[interaction->partner]].rot*(*parameters.group_elements[group->elements[interaction->partner]].interactions[interaction->back_link].initial_location))).V3Norm();
				group->Type->rand_independent=false;
#if DEBUG_LEVEL>2
				cout << "\t-> Bond distance between " << i+1 << " and " << interaction->partner+1 << ": " << interaction->bond_length << " Angström\n";
#endif
				count++;
			}
			if(hasbondbendpotential){
				interaction->allow_bond_bend=true;
			} else{
				interaction->allow_bond_bend=group->Type->allow_bond_bend;
			}
			if(hasspringpotential){
				interaction->allow_bond_stretch=true;
				interaction->allow_bond_bend=true;
				// calculate bond vector (points from current element to partner)
				Vec3 bond=parameters.group_elements[group->elements[interaction->partner]].center+(parameters.group_elements[group->elements[interaction->partner]].rot*(*parameters.group_elements[group->elements[interaction->partner]].interactions[interaction->back_link].initial_location))-element->center-(element->rot*(*interaction->initial_location));
				interaction->bond_length=bond.V3Norm();
				// use current bond direction as bond normal for bend potential
				*interaction->initial_normal=element->rot.M3Transpose()*bond;
				*interaction->initial_normal/=interaction->initial_normal->V3Norm();
				if(!interaction->normal) interaction->normal=new Vec3;
				*interaction->normal=element->rot*(*interaction->initial_normal);
				group->Type->rand_independent&=true;
#if DEBUG_LEVEL>1
				cout << "\t-> Bond distance between " << i+1 << " and " << interaction->partner+1 << ": " << interaction->bond_length << " Angström\n";
#endif
			}
		}
	}
#if DEBUG_LEVEL>1
	cout << "\t-> " << bondcount << " links found: " << count << " fixed bond potentials assigned (" << nr_statics << " static links found).\n";
#endif
	if(bondcount-nr_statics<1){
		cout << "\t-> Group is rigid, setting number of elements to be randomized to zero.\n";
		group->Type->rand_elements=0;
	}
	cout << "<- Done\n";
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void MC_Config::LoadSimParams(char* conf)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	cout << "\nReading simulation parameters\n";
	SetParam(parameters.output_code,"output_code",conf,true);
	SetParam(parameters.runid,"runid",conf,0);
	SetParam(parameters.fileout,"FileOut",conf,"test");
	SetParam(parameters.test,"test",conf,false);
	SetParam(parameters.dogr,"dogr",conf,false);
	SetParam(parameters.doxm,"doxm",conf,false);
	SetParam(parameters.chkcoord,"chkcoord",conf,false);
	SetParam(parameters.latticetype,"latticetype",conf,1); // default is simple cubic
	SetParam(parameters.PBCs,"PBCs",conf,Vec3(1.0)); // default is PBCs in all directions
	SetParam(parameters.cut,"cut",conf,true);
	SetParam(parameters.cut_internal,"cut_internal",conf,parameters.cut);
	SetParam(parameters.dyneps,"dyneps",conf,true);
	SetParam(parameters.polarize,"polarize",conf,false);
	SetParam(parameters.dyneps_varm,"dyneps_varm",conf,false);
	SetParam(parameters.dyneps_average_cycles,"dyneps_average_cycles",conf,0);
	SetParam(parameters.NpT,"NpT",conf,false);
	SetParam(parameters.adjust_internal,"adjust_internal",conf,false);
	SetParam(parameters.collective_moves,"collective_moves",conf,false);
	SetParam(parameters.VmuE_squared,"VmuE_squared",conf,false);
	SetParam(parameters.pext,"pext",conf,1.0);
	parameters.pext*=atm_in_pErg_per_Ang3;
	SetParam(parameters.selfrf,"selfrf",conf,false);
	SetParam(parameters.vdwtype,"vdwtype",conf,4);
	SetParam(parameters.loader,"loader",conf,2);
	SetParam(parameters.realfield,"realfield",conf,false);
	SetParam(parameters.anneal,"anneal",conf,false);
	SetParam(parameters.use_gpu,"use_gpu",conf,false);
	SetParam(parameters.xmfreq,"xmfreq",conf,100);
	SetParam(parameters.dipoleQ,"dipole_charge",conf,0.0);
	parameters.dipoleQ*=e_in_esu;
#ifdef USE_CMWC4096
	SetParam(parameters.time,"time",conf,false);
#else
	parameters.time=false;
#endif
	SetParam(parameters.randsteps_nointernal,"randsteps_nointernal",conf,false);
	SetParam(parameters.randsteps_novolume,"randsteps_novolume",conf,false);
	SetParam(parameters.steps,"steps",conf);
	SetParam(parameters.group_rand_frac,"group_rand_frac",conf,1.0);
	if(parameters.group_rand_frac<0.0){
		cout << "WARNING: Fraction of elements to be randomized in a group needs to >=0. Adjusted to 0.\n";
		parameters.group_rand_frac=0.0;
	}
	SetParam(parameters.randsteps,"randsteps",conf,1000);
	SetParam(parameters.laststep,"laststep",conf,parameters.steps/2);
	SetParam(parameters.stepsize_average,"stepsize_average",conf,0);
	SetParam(parameters.transition,"transition",conf,false);
	SetParam(parameters.transition_start,"transition_start",conf,parameters.randsteps);
	SetParam(parameters.transition_cycles,"transition_cycles",conf,parameters.randsteps);
	SetParam(parameters.transition_delta_n2,"transition_delta_n2",conf,0.0);
	SetParam(parameters.transition_LJ_start,"transition_LJ_start",conf,0.0);
	SetParam(parameters.transition_LJ_end,"transition_LJ_end",conf,1.0);
	parameters.transition_LJ_delta=parameters.transition_LJ_end-parameters.transition_LJ_start;
	SetParam(parameters.add_slowly,"add_slowly",conf,0);
	if(parameters.add_slowly>0){ // if element/groups are added slowly a couple of things are not needed
		parameters.randsteps=0;
		parameters.transition=false;
	}
	SetParam(parameters.nndist,"nndist",conf);
	parameters.Efield=Vec3(0.0);
	double* E=get_flex_tupel("Efield",conf);
	if(E){
		if((int)E[0]==3){
			parameters.Efield.vec[0]=E[1];
			parameters.Efield.vec[1]=E[2];
			parameters.Efield.vec[2]=E[3];
			parameters.noEfield=false;
		} else{
			if((int)E[0]==1){
				parameters.Efield.vec[2]=E[1];
				parameters.noEfield=false;
			} else{
				cout << "ERROR: Please specify <Efield> parameter as either a scalar (field in z-direction) or a 3-vector.\n";
				exit(2);
			}
		}
		delete[] E;
	}
	if(parameters.Efield.V3Norm()>EPS) parameters.noEfield=false; else parameters.noEfield=true;
	SetParam(parameters.rand_Efield,"rand_Efield",conf,false);
	SetParam(parameters.rotate_Efield_steps,"rotate_Efield_steps",conf,0);
	SetParam(parameters.T,"T",conf,298.0);
	SetParam(parameters.epsilon,"epsilon",conf,1.0);
	SetParam(parameters.n2,"n2",conf,1.0);
	SetParam(parameters.boxlength,"boxlength",conf);
	parameters.inv_boxlength[0]=1.0/parameters.boxlength[0]; parameters.inv_boxlength[1]=1.0/parameters.boxlength[1]; parameters.inv_boxlength[2]=1.0/parameters.boxlength[2];
	SetParam(parameters.auto_volume,"auto_volume",conf,true);
	SetParam(parameters.maxtrans,"maxtrans",conf);
	SetParam(parameters.maxrot,"maxrot",conf);
	SetParam(parameters.rT,"rT",conf,1.65); // if not specified use default (needs to be specified)
	if(parameters.rT<0.0){
		cout << "ERROR: Test sphere radius rT needs to be positive value.\n";
		exit(3);
	}
	SetParam(parameters.LJwall_a,"LJwall_a",conf,0.0);
	SetParam(parameters.LJwall_b,"LJwall_b",conf,parameters.LJwall_a);
	SetParam(parameters.LJwall_width,"LJwall_width",conf,0.0);
	SetParam(parameters.LJwall_xm,"LJwall_xm",conf,parameters.boxlength[0]*0.5);
	SetParam(parameters.LJwall_epsilon,"LJwall_epsilon",conf,1.0);
	parameters.LJwall_mirror_factor=(parameters.n2-parameters.LJwall_epsilon)/(parameters.n2+parameters.LJwall_epsilon);
	parameters.LJwall_x=parameters.LJwall_xm;
	SetParam(parameters.LJwall_fixed,"LJwall_fixed",conf,false);
	if((fabs(parameters.LJwall_a)>EPS) || (fabs(parameters.LJwall_b)>EPS)){
		parameters.PBCs[0]=false;
		parameters.LJwall_calc=true;
	} else parameters.LJwall_calc=false;
	if(parameters.NpT){
		SetParam(parameters.NpT_lnV_change,"NpT_lnV_change",conf,true);
		SetParam(parameters.NpT_relative_coords,"NpT_relative_coords",conf,true);
		if(parameters.NpT && !parameters.NpT_relative_coords && parameters.LJwall_calc){
			cout << "WARNING: Need NpT_relative_coords = 1 for system with a LJ wall. Changed NpT_relative_coords to 1.\n";
			parameters.NpT_relative_coords=true;
		}
		if(parameters.NpT && !parameters.NpT_relative_coords && ((fabs(parameters.PBCs[0])<EPS) || (fabs(parameters.PBCs[1])<EPS) || (fabs(parameters.PBCs[2])<EPS))){
			cout << "WARNING: Need NpT_relative_coords = 1 for system without periodic boundary conditions in all directions. Changed NpT_relative_coords to 1.\n";
			parameters.NpT_relative_coords=true;
		}
		if(parameters.NpT_lnV_change){
			SetParam(parameters.maxdV,"maxdV",conf,0.01);
			if(parameters.maxdV>log(parameters.boxlength[0]*parameters.boxlength[1]*parameters.boxlength[2])){
				cout << "ERROR: Specified maximum logarithmic change in NpT volume (" << parameters.maxdV << ") is larger then specified logarithm of simulation volume (" << log(parameters.boxlength[0]*parameters.boxlength[1]*parameters.boxlength[2]) <<").\n";
				exit(4);
			}
		} else{
			SetParam(parameters.maxdV,"maxdV",conf,50.0);
			if(parameters.maxdV>parameters.boxlength[0]*parameters.boxlength[1]*parameters.boxlength[2]){
				cout << "ERROR: Specified maximum change in NpT volume (" << parameters.maxdV << ") is larger then specified of simulation volume (" << parameters.boxlength[0]*parameters.boxlength[1]*parameters.boxlength[2] <<").\n";
				exit(4);
			}
		}
		SetParam(parameters.NpT_adjust_EScut,"NpT_adjust_EScut",conf,parameters.cut);
	}
	if((parameters.latticetype==3) || (parameters.latticetype==4)) SetParam(parameters.spherecylr,"spherecylr",conf);
	SetParam(parameters.LJ_interaction_area,"LJ_interaction_area",conf,false);
	SetParam(parameters.LJ_interaction_area_fit,"LJ_interaction_area_fit",conf,false);
	SetParam(parameters.LJ_adjust_width,"LJ_adjust_width",conf,false);
	SetParam(parameters.LJ_epsilon_sixthpower,"LJ_epsilon_sixthpower",conf,false);
	SetParam(parameters.LJ_radius_mixing,"LJ_radius_mixing",conf,0);
	SetParam(parameters.ljcut,"ljcut",conf,4*parameters.nndist);
	parameters.ljcut2=parameters.ljcut*parameters.ljcut;
	SetParam(parameters.escut,"escut",conf,parameters.boxlength[0]/2.0);
	parameters.escut2=parameters.escut*parameters.escut;
	SetParam(parameters.ES_cutoff_precision,"ES_cutoff_precision",conf,EPS);
	SetParam(parameters.LJ_cutoff_precision,"LJ_cutoff_precision",conf,EPS);
	parameters.LJexp[0]=6; parameters.LJexp[1]=12; // default value
	SetParam(parameters.LJexp,2,"LJexp",conf,parameters.LJexp);
	SetParam(parameters.rmax,"rmax",conf,parameters.boxlength[0]/2.0);
	SetParam(parameters.Solvent,"Solvent",conf,Vec3(-1.01,0.1,12.43));
	SetParam(parameters.xmfreq,"xmfreq",conf,100);
	SetParam(parameters.grfreq,"grfreq",conf,100);
	SetParam(parameters.correlation_length,"correlation_length",conf,0);
	SetParam(parameters.running_average,"running_average",conf,4);
	SetParam(parameters.grinc,"grinc",conf,200);
	SetParam(parameters.rmax,"rmax",conf,parameters.boxlength[0]/2.0);
	SetParam(parameters.refdr,"refdr",conf,parameters.rmax/parameters.grinc);
	SetParam(parameters.stillrot,"stillrot",conf,Mat33());
	SetParam(parameters.stilleuler,"stilleuler",conf,Vec3(0.0));
	SetParam(parameters.smooth_animation,"smooth_animation",conf,false);
	SetParam(parameters.quadrupole_charge_factor,"quadrupole_charge_factor",conf,8.0);
	if(parameters.quadrupole_charge_factor<1.0){
		cout << "WARNING: Quadrupole charge factor < 1. Is this really what you want?.\n";
	}
	parameters.quadrupole_charge_factor*=e_in_esu;
	SetParam(parameters.energy_scale,4,"energy_scale",conf,Vec4(1.0,1.0,1.0,1.0).vec);
	parameters.energy_scale[0]*=-1.0; // take care of negative prefactor
/*	SetParam(parameters.Vvdwboundary,"Vvdwboundary",conf,0.5);
	SetParam(parameters.stilleuler,"uBoundary",conf,12.0);
	SetParam(parameters.stilleuler,"rBoundary",conf,5.55);
	SetParam(parameters.stilleuler,"uSigma",conf,14.1);*/
#if DEBUG_LEVEL>0
	cout << "<- Finished reading simulation parameters.\n";
#endif
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
} // general simulation parameters initialized - next is element type creation ...

bool MC_Config::VLJatomisticZero(Vec3 &center, unsigned int* combine_elements, unsigned int count, Vec3 direction, Vec3 &r, double rT)
{
	bool success=true;
	r.vec[2]=(r.vec[0]+r.vec[1])/2.0; // r[0] is in positive region, r[1] negative, and r[2] in middle between both
	Vec3 VLJatomistic;
	while(fabs(r.vec[0]-r.vec[1])*2.0>EPS){
		VLJatomistic.V3Zeros();
		for(unsigned int i=0; i<count; i++){
			Element* curr=&parameters.group_elements[combine_elements[i]];
			double curr_eps=sqrt(curr->MyType->Vvdw);
			for(unsigned int k=0; k<3; k++){
				Vec3 pv=center+direction*r.vec[k]-curr->center; // test position - current center (points from current center to test position)
				double rc=pv.V3Norm();
				double rminc=EllipsoidRmin(pv,curr->MyType->saxes,curr->rot)+rT;
				double disp_i;
				switch(parameters.vdwtype){
					default:
						disp_i = qqpwr(rminc/rc,parameters.LJexp[0]);
						VLJatomistic.vec[k] += curr_eps*disp_i*(disp_i-1.0);
						break;
					case 2:
						disp_i = rminc/rc;
						//LJ + Bruce correction around r=rmin. Tanh functions have been replaced with a Gaussian for speed
						VLJatomistic.vec[k] += curr_eps*(qqpwr(disp_i,parameters.LJexp[1])-qqpwr(disp_i,parameters.LJexp[0])) + parameters.Solvent[1]*exp(parameters.Solvent[0]*(rc - parameters.Solvent[2])*(rc - parameters.Solvent[2]));
						break;
				}
			}
		}
		if(VLJatomistic.vec[1]>0.0){ // safety check that r[1] in outside (if it isn't move towards outside)
			r.vec[1]+=fabs(r.vec[0]-r.vec[1]);
		} else if(VLJatomistic.vec[2]>0.0) r.vec[0]=r.vec[2]; else r.vec[1]=r.vec[2];
		r.vec[2]=(r.vec[0]+r.vec[1])/2.0;
	}
	return success;
}

bool MC_Config::dVLJatomisticZero(Vec3 &center, unsigned int* combine_elements, unsigned int count, Vec3 direction, Vec3 &r, double rT)
{
	bool success=true;
	r.vec[2]=(r.vec[0]+r.vec[1])/2.0; // r[0] is in positive region, r[1] negative, and r[2] in middle between both
	Vec3 VLJatomistic;
	while(fabs(r.vec[0]-r.vec[1])*2.0>EPS){
		VLJatomistic.V3Zeros();
		for(unsigned int i=0; i<count; i++){
			Element* curr=&parameters.group_elements[combine_elements[i]];
			double curr_eps=sqrt(curr->MyType->Vvdw);
			for(unsigned int k=0; k<3; k++){
				Vec3 pv=center+direction*r.vec[k]-curr->center; // test position - current center (points from current center to test position)
				double rc=pv.V3Norm();
				double rminc=EllipsoidRmin(pv,curr->MyType->saxes,curr->rot)+rT;
				double disp_i;
				switch(parameters.vdwtype){
					default:
						disp_i = qqpwr(rminc/rc,parameters.LJexp[0]);
						VLJatomistic.vec[k] += curr_eps*disp_i*(1.0-2.0*disp_i)*parameters.LJexp[0]/rc;
						break;
					case 2:
						disp_i = rminc/rc;
						//LJ + Bruce correction around r=rmin. Tanh functions have been replaced with a Gaussian for speed
						VLJatomistic.vec[k] += curr_eps*(parameters.LJexp[0]*qqpwr(disp_i,parameters.LJexp[0])-parameters.LJexp[1]*qqpwr(disp_i,parameters.LJexp[1])) + 2.0*parameters.Solvent[0]*(rc - parameters.Solvent[2])*parameters.Solvent[1]*exp(parameters.Solvent[0]*(rc - parameters.Solvent[2])*(rc - parameters.Solvent[2]));
						break;
				}
			}
		}
		if(VLJatomistic.vec[1]<0.0){ // safety check that r[1] in outside (if it isn't move towards outside)
			r.vec[1]+=fabs(r.vec[0]-r.vec[1]);
		} else if(VLJatomistic.vec[2]<0.0) r.vec[0]=r.vec[2]; else r.vec[1]=r.vec[2];
		r.vec[2]=(r.vec[0]+r.vec[1])/2.0;
	}
	return success;
}

void MC_Config::CombineElements(Element_Group* group, Element_Group* newgroup, unsigned int* combine_elements, unsigned int count, string specifier, double origV, double origEps, double rT)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	cout << "\t-> Creating ellipsoid <" << newgroup->Type->name+specifier << ">\n";
	nr_group_elements++;
	parameters.group_elements = (Element*)realloc(parameters.group_elements,nr_group_elements*sizeof(Element));
	if(parameters.group_elements==NULL){ // Danger, danger
		cout << "Could not find a memory block large enough to hold a total of " << nr_group_elements << " group associated elements.\n";
		exit(1);
	}
	
	newgroup->nr_elements++;
	newgroup->elements = (unsigned int*)realloc(newgroup->elements,newgroup->nr_elements*sizeof(unsigned int));
	if(newgroup->elements==NULL){ // More danger ...
		cout << "Not enough memory to hold element numbers of group " << newgroup->Type->name << "\n";
		exit(1);
	}
	newgroup->elements[newgroup->nr_elements-1]=nr_group_elements-1;
	Element* element=&parameters.group_elements[nr_group_elements-1];
	// Need new element type for the new element
	parameters.num_element_types++;
	parameters.element_types = (Element_Type**)realloc(parameters.element_types,parameters.num_element_types*sizeof(Element_Type*));
	if(parameters.element_types==NULL){
		cout << "Could not find a memory block large enough to hold more element types.\n";
		exit(1);
	}
	parameters.element_types[parameters.num_element_types-1]=new Element_Type;
	Element_Type* element_type=parameters.element_types[parameters.num_element_types-1];
	
	element->MyType=element_type;
	element->mytype=parameters.num_element_types-1;
	element_type->eps_texture=NULL;
	element_type->name=newgroup->Type->name+specifier;
	element_type->archetype=-1;
	element_type->thistype=element->mytype;
	element_type->number=0;
	element_type->dof=0;
	element_type->initial_dipole=Vec3(0.0);
	element_type->mu_pos=Vec3(0.0);
	element_type->dipole=0.0;
	element_type->nr_charges=0;
	element_type->calculate_order=group->Type->calculate_order; // TODO: check if necessary (I don't think so at the moment, but it can't hurt either)
	element_type->q=NULL;
	element_type->q_pos=NULL;
	element_type->saxes=Vec3(0.0);
	element_type->rot_notrans=false;
	element_type->still=false;
	element->q_pos=NULL;
	element->nr_interactions=0;
	element->interactions=NULL;
	element->group=newgroup;
	element->dipole=Vec3(0.0);
	
	// Get arithmetic center of elements as new center
	element->center=Vec3(0.0);
	element_type->color=Vec3(0.0);
	element_type->transparency=0.0;
	element_type->label_elements=true;
	element_type->mass=0.0;
	element_type->MOI.M3Zeros(); // moment of inertia tensor is calculated once rotation of ellipsoid is known
	double wsum=0;
	unsigned int group_start_element=group->elements[0];
	bool mass_defined=true;
	for(unsigned int i=0; i<count; i++){
		Element* curr=&parameters.group_elements[combine_elements[i]];
#ifdef LOD_CENTER_WEIGHTED
		double w=sqrt(curr->MyType->Vvdw);
#else
		double w=1.0;
#endif
		wsum+=w;
		element->center+=curr->center*w;
		if((group->Type->LOD->nr_components>0) && (group->Type->LOD->element_in_component[combine_elements[i]-group_start_element]>=0)){
			element_type->color+=group->Type->LOD->component_color[group->Type->LOD->element_in_component[combine_elements[i]-group_start_element]];
			if(group->Type->LOD->component_texture[group->Type->LOD->element_in_component[combine_elements[i]-group_start_element]]!="") element_type->texture=group->Type->LOD->component_texture[group->Type->LOD->element_in_component[combine_elements[i]-group_start_element]];
			if(group->Type->LOD->component_texture_is_dipole[group->Type->LOD->element_in_component[combine_elements[i]-group_start_element]]) element_type->texture_is_dipole=true; else element_type->texture_is_dipole=false;
			if(group->Type->LOD->component_transparency[group->Type->LOD->element_in_component[combine_elements[i]-group_start_element]]>=0.0)
				element_type->transparency+=group->Type->LOD->component_transparency[group->Type->LOD->element_in_component[combine_elements[i]-group_start_element]];
		} else{
			element_type->color+=curr->MyType->color;
			element_type->transparency+=curr->MyType->transparency;
		}
		if(curr->MyType->mass<=EPS) mass_defined=false;
		element_type->mass+=curr->MyType->mass;
	}
	if(!mass_defined && (element_type->mass>EPS)){
		cout << "WARNING: Ellipsoid mass is not defined.\n";
		element_type->mass=0.0;
	}
	element->center/=wsum;
	element_type->color/=count;
	element_type->transparency/=count;
	
	double cr = (double)(parameters.LJexp[1])/(double)(parameters.LJexp[0]);
	parameters.r=pow(cr,cr/(cr-1.0))/(cr-1.0);
	
	double Rtest=1E6;
	if(newgroup->Type->LOD->symmetry_center){
		bool success=true;
		Vec3 newcenter=element->center;
		Vec3 centerdiff;
		unsigned int N=0;
		cout << "\t\t-> Determining maximum symmetry center\n";
		do{
			centerdiff.V3Zeros();
			wsum=0.0;
			// calculate zero crossing of fully atomistic potential in +/- direction defined through angles theta and phi
			for(unsigned int ct=0; ct<=theta_res/2; ct++){ // theta shall only go from 0 to pi/2
				double costheta=1.0-2.0*(double)ct/(double)theta_res; // cos(theta) ranges from 1.0 to 0.0
				double sint=sqrt(1.0-costheta*costheta); // cos²(theta) + sin²(theta) = 1
				for(unsigned int p=0; p<phi_res; p++){ // phi ranges from 0 to 2*pi
					double phi=(double)p/(double)phi_res*2.0*pi;
					
					Vec3 direction(sint*cos(phi),sint*sin(phi),costheta);
					Vec3 r;
					Vec3 currcenter(0.0);
					for(unsigned int s=0; s<2; s++){
						r.vec[0]=(2.0*s-1.0)*pow(EPS,1.0/(double)parameters.LJexp[0]);
						r.vec[1]=(2.0*s-1.0)*(parameters.ljcut+Rtest);
						success&=VLJatomisticZero(newcenter,combine_elements,count,direction,r,Rtest);
						double w=r.vec[2]*r.vec[2]; // weight accelerates convergence
						currcenter+=direction*r.vec[2]*w;
						wsum+=w;
					}
					centerdiff+=currcenter;
				}
			}
			centerdiff/=wsum;
			N++;
			cout << "\t\t\t-> center shift " << N << ": (" << centerdiff.V3Str(',') << ")\n";
			newcenter+=centerdiff;
		} while((centerdiff*centerdiff>EPS*EPS) && success);
		
		if(success){
			cout << "\t\t-> New center: (" << newcenter.V3Str(',') << "), ";
		} else{
			cout << "\t\t-> Error calculating new center, using ";
		}
#ifdef LOD_CENTER_WEIGHTED
			cout << "LJ epsilon weighted ";
#endif
			cout << "geometric center: (" << element->center.V3Str(',') << ")\n";
			element->center=newcenter;
	}
	if(rT>=0.0) Rtest=rT;
	// Determine gyration tensor
	Vec3 points[6];
	Mat33 S;
	S.M3Zeros();
	double gcount=0.0;
	for(unsigned int i=0; i<count; i++){
		Element* curr=&parameters.group_elements[combine_elements[i]];
#if DEBUG_LEVEL>3
		cout << curr->MyType->name << "\n";
#endif
		Vec3 center=curr->center-element->center;
		for(unsigned int l=0; l<3; l++){
			points[l<<1].vec[0]=(center.vec[0]+curr->rot.mat[l][0]*curr->MyType->saxes.vec[0]);
			points[(l<<1)+1].vec[0]=(center.vec[0]-curr->rot.mat[l][0]*curr->MyType->saxes.vec[0]);
			
			points[l<<1].vec[1]=(center.vec[1]+curr->rot.mat[l][1]*curr->MyType->saxes.vec[1]);
			points[(l<<1)+1].vec[1]=(center.vec[1]-curr->rot.mat[l][1]*curr->MyType->saxes.vec[1]);
			
			points[l<<1].vec[2]=(center.vec[2]+curr->rot.mat[l][2]*curr->MyType->saxes.vec[2]);
			points[(l<<1)+1].vec[2]=(center.vec[2]-curr->rot.mat[l][2]*curr->MyType->saxes.vec[2]);
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
	// Gyration radius Rg^2 = lambda_x^2+lambda_y^2+lambda_z^2
	double Rg=sqrt(semiaxis.vec[0]+semiaxis.vec[1]+semiaxis.vec[2]);
#if DEBUG_LEVEL>1
	cout << "\t\t-> Radius of gyration is " << Rg << " Angström.\n";
#endif
	// use smallest dimension as test sphere radius in case no value is specified
	if(rT<0.0) Rtest=mind(semiaxis.vec,3);
	if(newgroup->Type->LOD->gyration_sphere){
#if DEBUG_LEVEL>1
		cout << "\t\t-> Constructing sphere from radius of gyration.\n";
#endif
		element_type->saxes.vec[0]=Rg; element_type->saxes.vec[1]=Rg; element_type->saxes.vec[2]=Rg;
	} else{
		element_type->saxes.vec[0]=sqrt(semiaxis.vec[0]); element_type->saxes.vec[1]=sqrt(semiaxis.vec[1]); element_type->saxes.vec[2]=sqrt(semiaxis.vec[2]);
	}
#if DEBUG_LEVEL>1
	cout << "\t\t-> Gyration tensor ellipsoid semiaxes: " << element_type->saxes.V3Str(',') << "\n";
#endif
	Mat33 ev=S.Eigenvectors(semiaxis);
	element->rot=RotFromEigenvectors(ev);
	// Done with rotation matrix, can now define moment of inertia tensor
	for(unsigned int i=0; i<count; i++){
		Element* curr=&parameters.group_elements[combine_elements[i]];
		Vec3 center=element->rot.M3Transpose()*(curr->center-element->center); // bring coordinates into ellipsoid coordinate system
		element_type->MOI.mat[0][0]+=curr->MyType->mass*(center.vec[1]*center.vec[1]+center.vec[2]*center.vec[2]);
		element_type->MOI.mat[1][1]+=curr->MyType->mass*(center.vec[0]*center.vec[0]+center.vec[2]*center.vec[2]);
		element_type->MOI.mat[2][2]+=curr->MyType->mass*(center.vec[0]*center.vec[0]+center.vec[1]*center.vec[1]);
		element_type->MOI.mat[0][1]-=curr->MyType->mass*center.vec[0]*center.vec[1];
		element_type->MOI.mat[0][2]-=curr->MyType->mass*center.vec[0]*center.vec[2];
		element_type->MOI.mat[1][2]-=curr->MyType->mass*center.vec[1]*center.vec[2];
	}
	element_type->MOI.mat[1][0]=element_type->MOI.mat[0][1];
	element_type->MOI.mat[2][0]=element_type->MOI.mat[0][2];
	element_type->MOI.mat[2][1]=element_type->MOI.mat[1][2];
	
#ifdef RT_TEST_LOOP
	ofstream datout;
	datout.open("rT_dependence.dat");
	if (datout.fail()){
		cout << "Unable to open output file.\n";
		exit(1);
	}
	datout << "# gyration tensor semi-axes: " << element_type->saxes.V3Str('\t') << "\n";
	datout << "# rT\ta\tb\tc\tvolume\tavg.width\tsphericity\n";
	Rtest=0;
	do{
#endif
	// Take care of semi-axes
	double newV=element_type->saxes.vec[0]*element_type->saxes.vec[1]*element_type->saxes.vec[2]*4.0/3.0*pi;
	if(newgroup->Type->LOD->match_volume){ // adjust volume to match originally excluded Lennard-Jones volume
		cout << "\t\t-> Adjusting ellipsoid volume (" << newV << " Angström³) to match";
		if(origV<0.0) cout << " originally excluded"; else cout << " user-specified";
		cout << " Lennard-Jones volume\n";
		if(origV<0.0){
			cout << "\t\t\t-> No volume specified, determining volume using Monte-Carlo routine (may take a while)\n";
#if DEBUG_LEVEL>1
			double tstart = (double)clock();
#endif
			// Find bounding box first using rotation matrix of new element
			Vec3 box(0.0);
			Vec3 box_corner=element->center;
			for(unsigned int i=0; i<count; i++){
				Element* curr=&parameters.group_elements[combine_elements[i]];
				Vec3 dist=element->rot.M3Transpose()*(curr->center-element->center);
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
#if DEBUG_LEVEL>3
			cout << box_corner.V3Str(',') << " | " << box.V3Str(',') << " | " << Rot2AxisAngle(element->rot).V4Str(',') << " | " << element->center.V3Str(',') << "\n";
#endif
			double inside;
#if DEBUG_LEVEL>1
			double var;
#endif
#ifdef USE_CMWC4096
			init_CMWC4096(8);
			zigset();
#else
			__int32_t idum=-8;
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
				point=element->rot*point+element->center;
				for(unsigned int i=0; i<count; i++){
					Element* curr=&parameters.group_elements[combine_elements[i]];
					if(Point_in_Ellipsoid(point,curr->MyType->saxes,curr->center,curr->rot)){
						Nin++;
						break; // only single-count
					}
				}
				N++;
				inside=(double)Nin/N;
			} while(((1-inside)/(N*inside)>1E-6) || (Nin<1000*count)); // first statement is testing for relative error of 1E-3 (0.1%)
			origV=inside*boundingV;
#if DEBUG_LEVEL>1
			var=boundingV*boundingV*inside*(1-inside)/N;
			double tend = (double)clock();
			cout << "\t\t\t\t-> Originally excluded volume is " << origV << " +/- " << sqrt(var) << " Angström³ (took " << (unsigned int)N << " iterations, took " << (tend-tstart)/CLOCKS_PER_SEC << " s)\n";
#endif
		}
		double scale=origV/newV;
		if((scale<0.25) || (scale>4.0)){
			cout << "WARNING: Adjusting the gyration tensor determined Lennard-Jones volume by more than a factor of four to " << origV << " Angström³ is not recommended, please choose a different set of elements.\n";
		}
		scale=cbrt(scale);
		cout << "\t\t\t-> Scaling factor: " << scale*100.0 << "% in each direction\n";
		element_type->saxes*=scale;
#if DEBUG_LEVEL>1
		cout << "\t\t\t-> Adjusted ellipsoid semiaxes: " << element_type->saxes.V3Str(',') << "\n";
#endif
	} else{
		if(newgroup->Type->LOD->symmetry_axes){
			cout << "\t\t-> Adjusting ellipsoid semiaxes- through average zero-crossing of fully-atomistic potential.\n";
			Vec3 curr_axes=element_type->saxes;
			Vec3 old_axes(0.0);
			bool ff=false;
			while((curr_axes-old_axes).V3Norm()>sqrt(EPS)){ // stop when difference between running average of last two
				if(ff) old_axes=curr_axes;
				ff=!ff;
				Vec3 axes(0.0);
				Vec3 variances(0.0);
				bool success=true;
				// x-axis semi-axis
				Vec3 N(0.0);
				Mat33 rotation;
				// calculate zero crossing of fully atomistic potential in +/- direction defined through angles theta and phi
				for(unsigned int ct=0; ct<=theta_res*pi_quarter_fraction; ct++){ // theta shall only go from 0 to pi/4
					double costheta=1.0-2.0*(double)ct/(double)theta_res;
					double sint=sqrt(1.0-costheta*costheta); // cos²(theta) + sin²(theta) = 1
					for(unsigned int p=0; p<phi_res; p++){ // 0 < phi < 2*pi (yes, it is correct to have phi_res points at theta=0 - otherwise it ain't uniform sampling)
						double phi=(double)p/(double)phi_res*2.0*pi;
						Vec3 zdir=Vec3(sint*cos(phi),sint*sin(phi),costheta); // in z-direction
						for(unsigned int xyz=0; xyz<3; xyz++){
							rotation.M3Eye();
							Vec3 ca(0,0,1);
							if(xyz==0){ // x-direction: rotate z 90 degree around y
								ca.vec[0]=1.0; ca.vec[2]=0.0;
								rotation.M3Zeros();
								rotation.mat[1][1] = 1.0; rotation.mat[0][2] = 1.0; rotation.mat[2][0] = -1.0;
							}
							if(xyz==1){ // y-direction: rotate z 90 degree around x
								ca.vec[1]=1.0; ca.vec[2]=0.0;
								rotation.M3Zeros();
								rotation.mat[0][0] = 1.0; rotation.mat[1][2] = 1.0; rotation.mat[2][1] = -1.0;
							}
							Vec3 dirvec=(rotation*zdir);
							dirvec.vec[0]*=curr_axes.vec[0]; dirvec.vec[1]*=curr_axes.vec[1]; dirvec.vec[2]*=curr_axes.vec[2];
							dirvec/=dirvec.V3Norm();
							Vec3 direction=element->rot*dirvec;
							double rmin=EllipsoidRmin(direction,curr_axes,element->rot);
							double rt=touch_sphere_sigma(direction,curr_axes,element->rot,Rtest)-rmin;
//							double rt=curr_axes.vec[xyz]/touch_sphere_sigma(direction,curr_axes,element->rot,Rtest);
							double w=qqpwr(dirvec*ca,2);
//							if(w<0.0) w=0.0; // else 
//							double w=1.0/sqrt(qqpwr(dirvec.vec[0]*curr_axes.vec[0],2)+qqpwr(dirvec.vec[1]*curr_axes.vec[1],2)+qqpwr(dirvec.vec[2]*curr_axes.vec[2],2));
							N.vec[xyz]+=2.0*w;
							Vec3 r;
							for(unsigned int s=0; s<2; s++){
								r.vec[0]=(2.0*s-1.0)*pow(EPS,1.0/(double)parameters.LJexp[0]);
								r.vec[1]=(2.0*s-1.0)*(parameters.ljcut+Rtest);
								success&=VLJatomisticZero(element->center,combine_elements,count,direction,r,Rtest);
								double axis=(fabs(r.vec[2])-rt)/rmin*curr_axes.vec[xyz];
//								double axis=fabs(r.vec[2])*rt;
								axes.vec[xyz]+=w*axis;
								variances.vec[xyz]+=w*axis*axis;
							}
						}
					}
				}
				for(unsigned int xyz=0; xyz<3; xyz++){
					axes.vec[xyz]/=N.vec[xyz];
					variances.vec[xyz]/=N.vec[xyz];
					variances.vec[xyz]-=axes.vec[xyz]*axes.vec[xyz];
					variances.vec[xyz]/=N.vec[xyz];
				}
				if(success){
					// t-distribution critical value taken from: www.itl.nist.gov/div898/handbook/eda/section3/eda3672.htm
					// t<|x-y|/sqrt(var(x)/n_x+var(y)/n_y) => |x-y|>t*sqrt(var(x)/n_x+var(y)/n_y)
					double tcritical=2.576; // 0.01 (two-sided) level, infinite degrees of freedom
//					double tcritical=1.96; // 0.05 (two-sided) level, infinite degrees of freedom
					bool sig12=((sqrt(variances.vec[0]+variances.vec[1])>EPS) && (fabs(axes.vec[0]-axes.vec[1])>tcritical*sqrt(variances.vec[0]+variances.vec[1])));
					bool sig13=((sqrt(variances.vec[0]+variances.vec[2])>EPS) && (fabs(axes.vec[0]-axes.vec[2])>tcritical*sqrt(variances.vec[0]+variances.vec[2])));
					bool sig23=((sqrt(variances.vec[1]+variances.vec[2])>EPS) && (fabs(axes.vec[1]-axes.vec[2])>tcritical*sqrt(variances.vec[1]+variances.vec[2])));
					if(!sig12 && !sig13 && !sig23){ // all axes differences not statistically significant
#if DEBUG_LEVEL>3
						cout << "\t\t\t-> Differences in axis lengths are not statistically significant (at 0.01 level).\n";
#endif
						double avg123=(axes.vec[0]+axes.vec[1]+axes.vec[2])/3.0;
						curr_axes.vec[0]=avg123;
						curr_axes.vec[1]=avg123;
						curr_axes.vec[2]=avg123;
					} else{ // check for two axes being equal
						if(!sig12){
#if DEBUG_LEVEL>3
							cout << "\t\t\t-> Length difference between axis #1 and #2 is not statistically significant (at 0.01 level).\n";
#endif
							double avg12=(axes.vec[0]+axes.vec[1])/2.0;
							curr_axes.vec[0]=avg12;
							curr_axes.vec[1]=avg12;
							curr_axes.vec[2]=axes.vec[2];
						} else{
							if(!sig13){
#if DEBUG_LEVEL>3
								cout << "\t\t\t-> Length difference between axis #1 and #3 is not statistically significant (at 0.01 level).\n";
#endif
								double avg13=(axes.vec[0]+axes.vec[2])/2.0;
								curr_axes.vec[0]=avg13;
								curr_axes.vec[1]=axes.vec[1];
								curr_axes.vec[2]=avg13;
							} else{
								if(!sig23){
#if DEBUG_LEVEL>3
									cout << "\t\t\t-> Length difference between axis #2 and #3 is not statistically significant (at 0.01 level).\n";
#endif
									double avg23=(axes.vec[1]+axes.vec[2])/2.0;
									curr_axes.vec[0]=axes.vec[0];
									curr_axes.vec[1]=avg23;
									curr_axes.vec[2]=avg23;
								} else{
									curr_axes.vec[0]=axes.vec[0];
									curr_axes.vec[1]=axes.vec[1];
									curr_axes.vec[2]=axes.vec[2];
								}
							}
						}
					}
				}
#if DEBUG_LEVEL>1
				cout << "\t\t\t-> Adjusted ellipsoid semiaxes, rT = " << Rtest << " Angström: " << curr_axes.V3Str(',') << "\n";
#endif
			}
			element_type->saxes=curr_axes;
		}
	}
	newV=element_type->saxes.vec[0]*element_type->saxes.vec[1]*element_type->saxes.vec[2]*4.0/3.0*pi;
#if DEBUG_LEVEL>1
	cout << "\t\t-> Volume: " << newV << "\n";
#endif
	// test sphere size is set to radius of sphere with same surface area
//	element_type->rT=sqrt(EllipsoidSurface(element_type->saxes.vec[0],element_type->saxes.vec[1],element_type->saxes.vec[2])/(4.0*pi));
//	element_type->rT=cbrt(element_type->saxes.vec[0]*element_type->saxes.vec[1]*element_type->saxes.vec[2]);
	// set up things for IA so we can use it right away
//	element_type->inv_avg_area=1.0/(pi*qqpwr(element_type->rT,2));
	element_type->inv_avg_area=1.0;
	element_type->IA_coefficients=NULL;
	element_type->invsaxes2.vec[0]=1.0/(element_type->saxes.vec[0]*element_type->saxes.vec[0]);
	element_type->invsaxes2.vec[1]=1.0/(element_type->saxes.vec[1]*element_type->saxes.vec[1]);
	element_type->invsaxes2.vec[2]=1.0/(element_type->saxes.vec[2]*element_type->saxes.vec[2]);
//	Rtest=element_type->rT;
	// sphericity is ratio of volume-corresponding sphere surface area to ellipsoid surface area
	element_type->sphericity=4*pi*qqpwr(cbrt(element_type->saxes.vec[0]*element_type->saxes.vec[1]*element_type->saxes.vec[2]),2)/EllipsoidSurface(element_type->saxes.vec[0],element_type->saxes.vec[1],element_type->saxes.vec[2]);
	element_type->IA_average=0.0;
	for(unsigned int ct=0; ct<theta_res*10; ct++){
		double cost=1.0-2.0*(ct+0.5)/(theta_res*10);
		double sint=sqrt(1.0-cost*cost);
		for(unsigned int j=0; j<phi_res*10; j++){
			double phi=2.0*pi*(j+0.5)/(phi_res*10);
			element_type->IA_average+=IA(element_type,Vec3(sint*cos(phi),sint*sin(phi),cost));
		}
	}
	element_type->IA_average+=2.0*IA(element->MyType,Vec3(0.0,0.0,1.0));
	element_type->IA_average/=theta_res*phi_res*100+2;
//	cout << "\t\t-> <IA> = " << element_type->IA_average << ", sphericity^(3/2) = " << element_type->sphericity*sqrt(element_type->sphericity) << " (" << 100*(element_type->IA_average-element_type->sphericity*sqrt(element_type->sphericity))/element_type->sphericity*sqrt(element_type->sphericity) << "% error)\n";
	cout << "\t\t-> Sphericity = " << element_type->sphericity << ", <IA> = " << element_type->IA_average << "\n";
	element_type->inv_avg_area/=element_type->IA_average;
	element_type->rT=sqrt(element_type->IA_average/pi);
//	Rtest=element_type->rT;
//	cout << "\t\t-> <IA> = " << element_type->IA_average << " => rT = " << element_type->rT << " Angström\n";
	element_type->r_diff=NULL;
	element_type->width_diff=NULL;
	element_type->avg_width=0.0;
	if(newgroup->Type->LOD->match_epsilon || parameters.LJ_adjust_width){
		double phi;
		Vec3 direction;
		bool success=true;
		Vec3 r;
		element_type->r_diff=new double[theta_res*phi_res+2]; // treat poles separately
		element_type->width_diff=new double[theta_res*phi_res+2];
		double avg_w=0.0;
		double avg_w2=0.0;
		double w;
		cout << "\t\t-> Calculating ellipsoid deviations from fully-atomistic shape and potential width: 0%";
		for(unsigned int t=0; t<theta_res; t++){
			if((t+1)%(theta_res/16)==0){
				if((t+1)%(theta_res/4)==0) cout << (double)(t+1)/theta_res*100.0 << "%"; else cout << ".";
				cout.flush();
			}
			for(unsigned int p=0; p<phi_res; p++){
				phi=2.0*pi*(double)(p+0.5)/(double)phi_res;
				double cost=1.0-2.0*(double)(t+0.5)/(double)theta_res;
				double sint=sqrt(1.0-cost*cost);
				direction.vec[0]=sint*cos(phi);
				direction.vec[1]=sint*sin(phi);
				direction.vec[2]=cost;
				direction=element->rot*direction; // get into ellipsoid coordinate system
				double r_ellipsoid=touch_sphere_sigma(direction,element_type->saxes,element->rot,Rtest);
				r.vec[0]=0.0;
				r.vec[1]=r_ellipsoid*2.0;
				success&=VLJatomisticZero(element->center,combine_elements,count,direction,r,Rtest);
				double r_atomistic=r.vec[2];
				element_type->r_diff[t*phi_res+p]=r_atomistic-r_ellipsoid;
				
				r.vec[0]=r_atomistic;
				r.vec[1]=r_atomistic*2.0;
				success&=dVLJatomisticZero(element->center,combine_elements,count,direction,r,Rtest);
				w=(r.vec[2]-r_atomistic)/(pow(2.0,1.0/(parameters.LJexp[1]-parameters.LJexp[0]))-1.0)-Rtest;
				avg_w += w;
				avg_w2 += w*w;
				element_type->width_diff[t*phi_res+p]=w;
#if DEBUG_LEVEL>4
				cout << phi << "\t" << theta << "\t" << r_ellipsoid << "\t" << r_atomistic << "\t" << element_type->r_diff[t*phi_res+p] << "\n";
#endif
			}
		}
		for(unsigned int i=0; i<2; i++){
			double cost,sint;
			if(i==0){
				cost=1.0;
				sint=0.0;
				phi=0.0;
			} else{
				cost=-1.0;
				sint=0.0;
				phi=0.0;
			}
			direction.vec[0]=sint*cos(phi);
			direction.vec[1]=sint*sin(phi);
			direction.vec[2]=cost;
			direction=element->rot*direction; // get into ellipsoid coordinate system
			double r_ellipsoid=touch_sphere_sigma(direction,element_type->saxes,element->rot,Rtest);
			r.vec[0]=0.0;
			r.vec[1]=r_ellipsoid*2.0;
			success&=VLJatomisticZero(element->center,combine_elements,count,direction,r,Rtest);
			double r_atomistic=r.vec[2];
			element_type->r_diff[theta_res*phi_res+i]=r_atomistic-r_ellipsoid;
			
			r.vec[0]=r_atomistic;
			r.vec[1]=r_atomistic*2.0;
			success&=dVLJatomisticZero(element->center,combine_elements,count,direction,r,Rtest);
			w=(r.vec[2]-r_atomistic)/(pow(2.0,1.0/(parameters.LJexp[1]-parameters.LJexp[0]))-1.0)-Rtest;
			avg_w += w;
			avg_w2 += w*w;
			element_type->width_diff[theta_res*phi_res+i]=w;
		}
		cout << "\n";
		avg_w/=theta_res*phi_res+2; // avg_w = (sigma+rT)*2^1/6 - (sigma+rT) => sigma = 1.0/(2^1/6-1)*avg_w - rT
		avg_w2/=theta_res*phi_res+2;
		if(avg_w2-avg_w*avg_w<EPS) avg_w2=avg_w*avg_w;
		element_type->avg_width=avg_w;
		cout << "\t\t\t-> Average width parameter for ellipsoid: " << element_type->avg_width << " +/- " << sqrt(avg_w2-avg_w*avg_w) << " Angström\n";
#if DEBUG_LEVEL>1
		cout << "\t\t\t-> Determined ellipsoid shape difference to fully-atomistic model (test sphere radius: " << Rtest << " Angström)\n";
#endif
	}
#ifdef RT_TEST_LOOP
		datout << Rtest << "\t" << element_type->saxes.V3Str('\t') << "\t" << newV << "\t" << element_type->avg_width << "\t" << element_type->sphericity << "\n";
		datout.flush();
		if(Rtest<9.99){
			Rtest+=0.1;
		} else{
			if(Rtest<19.99){
				Rtest+=1.0;
			} else{
				if(Rtest<49.99){
					Rtest+=2.0;
				} else{
					if(Rtest<99.99){
						Rtest+=5.0;
					} else{
						if(Rtest<199.99){
							Rtest+=10.0;
						} else{
							if(Rtest<499.99){
								Rtest+=20.0;
							} else{
								if(Rtest<999.99){
									Rtest+=50.0;
								} else{
									if(Rtest<1999.99){
										Rtest+=100.0;
									} else{
										if(Rtest<4999.99){
											Rtest+=200.0;
										} else{
											Rtest+=500.0;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	} while(Rtest<10000.1);
	datout.close();
	exit(42);
#endif
#if DEBUG_LEVEL>1
	cout << "\t\t-> Use LOD ellipsoid to determine Lennard-Jones epsilon through weighted RSS ...\n";
#endif
	// Now that the LOD ellipsoid is known (apart from a constant factor), we can calculate the average Lennard-Jones epsilon as a weighted average
	// using <eps> = sum(w_i*eps_i)/sum(w_i)
	// with weights (arbitrarly defined as) w_i = (|x_i - <x>|+r_min^i)/r_min^LOD (was: w_i = |x_i - <x>|+Rg) <- idea for this formula is to give a bit more weigth to LJ spheres closer to ellipsoid surface
	element_type->Vvdw=0.0;
	wsum=0.0;
	double avgvdw=0.0;
	double maxvdw=0.0;
	element_type->Vvdw_coefficients=NULL;
	if(newgroup->Type->LOD->match_epsilon && (!((fabs(element_type->saxes.vec[0]-element_type->saxes.vec[1])<EPS) && (fabs(element_type->saxes.vec[1]-element_type->saxes.vec[2])<EPS)))){
		element_type->nr_Vvdw_coefficients=8;
		element_type->Vvdw_coefficients=new double[8];
		element_type->Vvdw_coefficients[0]=element_type->rT;
	} else element_type->nr_Vvdw_coefficients=0;
	for(unsigned int i=1; i<element_type->nr_Vvdw_coefficients; i++) element_type->Vvdw_coefficients[i]=0.0;
	for(unsigned int i=0; i<count; i++){
		Element* curr=&parameters.group_elements[combine_elements[i]];
		double w=1.0;
		Vec3 dir=curr->center-element->center;
		double dist=dir.V3Norm();
		if(dist>EPS) w=(dist+EllipsoidRmin(dir,element->MyType->saxes,element->rot))/EllipsoidRmin(dir,element_type->saxes,element->rot);
		wsum+=w;
		// RSS sqrt(LJ epsilon) = sqrt(sum sqrt(LJ epsilon)^2) => LJ eps. = sum LJ eps.
		// -- in other words: RSS = sqrt(N)*RMS
		avgvdw+=curr->MyType->Vvdw;
		// in analogy to the above, weighted RSS = sqrt(N)*weighted RMS
		// weighted RSS sqrt(LJ epsilon) = sqrt(N)*sqrt(sum w*sqrt(LJ epsilon)^2*/sum w) => LJ eps. = (sum w*LJ eps.)*N/sum w
		element_type->Vvdw+=w*curr->MyType->Vvdw;
		double sqrt_eps_i=sqrt(curr->MyType->Vvdw);
		if(sqrt_eps_i>maxvdw) maxvdw=sqrt_eps_i;
		if(element_type->nr_Vvdw_coefficients==8){
			double s_i=average(curr->MyType->saxes.vec,3);
			element_type->Vvdw_coefficients[7]+=sqrt_eps_i;
			element_type->Vvdw_coefficients[6]+=sqrt_eps_i*6.0*s_i;
			element_type->Vvdw_coefficients[5]+=sqrt_eps_i*15.0*qqpwr(s_i,2);
			element_type->Vvdw_coefficients[4]+=sqrt_eps_i*20.0*qqpwr(s_i,3);
			element_type->Vvdw_coefficients[3]+=sqrt_eps_i*15.0*qqpwr(s_i,4);
			element_type->Vvdw_coefficients[2]+=sqrt_eps_i*6.0*qqpwr(s_i,5);
			element_type->Vvdw_coefficients[1]+=sqrt_eps_i*qqpwr(s_i,6);
		}
	}
	element_type->Vvdw*=count/wsum;
#if DEBUG_LEVEL>1
	cout << "\t\t\t-> Lennard-Jones epsilon is " << element_type->Vvdw << " (simple RSS: " << avgvdw;
	if((element_type->nr_Vvdw_coefficients>0) && parameters.fit2lod){
		cout << "; coefficient guesses: ";
		for(unsigned int i=0; i<element_type->nr_Vvdw_coefficients; i++){
			if(i>0) cout << ", ";
			cout << element_type->Vvdw_coefficients[i];
		}
	}
	cout << ")\n";
#endif
	if(!parameters.fit2lod){
		if(newgroup->Type->LOD->match_epsilon && (newgroup->levelofdetail>=newgroup->Type->LOD->start_level)){
			cout << "\t\t-> Adjusting Lennard-Jones epsilon to match fully-atomistic potential energies.\n";
			if(origEps<0.0){
				cout << "\n****************************************************************\n* Please run \"fit2lod <configuration file>\" in order to set up *\n* the LOD ellipsoid parameters for this simulation.            *\n****************************************************************\n\n";
				exit(8);
			}
			element_type->Vvdw=origEps;
			cout << "\t\t\t-> Lennard-Jones epsilon: " << element_type->Vvdw << "\n";
		}
	}
	cout << "\t<- Done.\n";
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void MC_Config::SetEllipsoidLinks(Element_Group* group, Element_Group* newgroup, int* combine_elements, unsigned int* counters)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	cout << "\t-> Creating links for <" << newgroup->Type->name << ">\n";
	
	int* elements_ellipses=new int[group->nr_elements];
	memset(elements_ellipses,0xFF,group->nr_elements*sizeof(int)); // puts -1 in each field -- for finding elements which don't belong to an ellipse
	for(unsigned int i=0; i<newgroup->nr_elements; i++){
		for(unsigned int j=0; j<(unsigned int)combine_elements[counters[i]]; j++){
			if(elements_ellipses[combine_elements[counters[i]+j+1]-1]<0){
				elements_ellipses[combine_elements[counters[i]+j+1]-1]=i;
			} else{ // Need to catch elements which are in more than one ellipsoid (otherwise the ring-finding algorithm will go into an infinite loop)
				cout << "ERROR: Group element " << combine_elements[counters[i]+j+1] << " appears in LOD ellipsoid " << elements_ellipses[combine_elements[counters[i]+j+1]-1]+1 << " and " << i+1 << ".\n";
				exit(4);
			}
		}
	} // we now have an array telling what ellipse each element is in (or if it's in none)
	
	for(unsigned int i=0; i<group->nr_elements; i++){
		if(elements_ellipses[i]<0){ // element is in no LOD ellipse - bail out and let user decide
			cout << "\t\t-> Individual element " << i+1 << " of <" << group->Type->name << "> is not assigned. Please make it part of a combined ellipsoid.\n";
			exit(8);
		}
	}
	
	unsigned int* multilinks=new unsigned int[newgroup->nr_elements*newgroup->nr_elements];
	memset(multilinks,0,newgroup->nr_elements*newgroup->nr_elements*sizeof(unsigned int));
	// Finally assign links
	for(unsigned int m=0; m<newgroup->nr_elements; m++){
		for(unsigned int n=0; n<(unsigned int)combine_elements[counters[m]]; n++){
			unsigned int i=combine_elements[counters[m]+n+1]-1;
			Element* element=&parameters.group_elements[group->elements[i]];
			Element* ellipse=&parameters.group_elements[newgroup->elements[elements_ellipses[i]]];
			for(unsigned int j=0; j<element->nr_interactions; j++){
				if(elements_ellipses[element->interactions[j].partner]!=elements_ellipses[i]){ // found link outside current ellipse
					multilinks[elements_ellipses[i]*newgroup->nr_elements+elements_ellipses[element->interactions[j].partner]]++;
					bool can_rotate=!(element->interactions[j].bond_order>1.0) && !((int)element->interactions[j].bond_order==-1);
#if DEBUG_LEVEL>1
					cout << "\t\t-> Creating";
					if(element->interactions[j].fixed) cout << " static";
					cout << " link between ellipsoids " << elements_ellipses[i]+1 << " and " << elements_ellipses[element->interactions[j].partner]+1 << " (elements " << i+1 << " and " << element->interactions[j].partner+1 << ", bond order: ";
					if(element->interactions[j].bond_order>0.0) cout  << element->interactions[j].bond_order; else cout << bondorder_types[(unsigned int)-element->interactions[j].bond_order];
					cout << ", ";
					if(!can_rotate){
						cout << "non-";
					}
					cout << "rotable)\n";
#endif
					if(!can_rotate && (element->interactions[j].nr_potentials==0)){
						cout << "\t\t\t-> Fixing non-rotable link.\n";
						element->interactions[j].fixed=true;
					}
					ellipse->nr_interactions++;
					ellipse->interactions=(Interaction*)realloc(ellipse->interactions,ellipse->nr_interactions*sizeof(Interaction));
					if(!ellipse->interactions){
						cout << "Not enough memory to add link.\n";
						exit(3);
					}
					Interaction* interaction=&ellipse->interactions[ellipse->nr_interactions-1];
					interaction->fixed=element->interactions[j].fixed;
					interaction->allow_bond_stretch=false; // for the moment, since no stretch potential has been assigned, connections are fixed in distance
					interaction->nr_potentials=0; // for the moment, potentials pointers have not been read
					interaction->potentials=NULL;
					interaction->partner=elements_ellipses[element->interactions[j].partner]; // which element we're connecting to
					
					interaction->initial_location=new Vec3;
					// for a rotation matrix R^T = R^(-1)
					*interaction->initial_location=ellipse->rot.M3Transpose()*(element->center+(element->rot*(*element->interactions[j].initial_location))-ellipse->center);
#if DEBUG_LEVEL>2
					cout << "\t\t\t-> Link location: (" << (ellipse->center+ellipse->rot*(*interaction->initial_location)).V3Str(',') << ")\n";
#endif
					if(*interaction->initial_location!=Vec3(0.0)){
						interaction->location = new Vec3;
						*interaction->location=ellipse->rot*(*interaction->initial_location);
					} else interaction->location=NULL;
					
					interaction->initial_normal=new Vec3;
					*interaction->initial_normal=ellipse->rot.M3Transpose()*(*element->interactions[j].initial_normal);
					if(*interaction->initial_normal!=Vec3(0.0)){
						interaction->normal = new Vec3;
						*interaction->normal=ellipse->rot*(*interaction->initial_normal);
					} else interaction->normal=NULL;
					
					interaction->initial_tangent=new Vec3;
					*interaction->initial_tangent=ellipse->rot.M3Transpose()*(*element->interactions[j].initial_tangent);
					if(*interaction->initial_tangent!=Vec3(0.0)){
						interaction->tangent = new Vec3;
						*interaction->tangent=ellipse->rot*(*interaction->initial_tangent);
					} else interaction->tangent=NULL;
					
					if(m>interaction->partner){ // can only determine back link if link has already been created ...
#if DEBUG_LEVEL>2
						cout << "\t\t\t-> Registering existing mutual back links ";
#endif
						bool found=false;
						Element* partner=&parameters.group_elements[newgroup->elements[interaction->partner]];
						unsigned int k=0;
						unsigned int linkcount=0;
						while(k<partner->nr_interactions){
							if(partner->interactions[k].partner==(unsigned int)elements_ellipses[i]){
								linkcount++;
								// want to find n'th link on other side if from this side n have been created
								// remember: backlinks are only set if links to current element have already been created
								if(linkcount==multilinks[elements_ellipses[i]*newgroup->nr_elements+elements_ellipses[element->interactions[j].partner]]){
									found=true;
									break;
								}
							}
							k++;
						}
						if(found) interaction->back_link=k; else{ // should not happen
							cout << "Could not find link back from " << element->interactions[j].partner+1 << " to " << i+1 << ".\n";
							exit(3);
						}
						partner->interactions[k].back_link=ellipse->nr_interactions-1; // don't forgot to set way back on other end
#if DEBUG_LEVEL>2
						cout << "(" << interaction->partner+1 << " and " << partner->interactions[interaction->back_link].partner+1 << ")\n";
#endif
					}
				}
			}
		}
	}
	delete[] multilinks;
	cout << "\t<- Done.\n";
	newgroup->Type->LOD->element_in_ellipsoid[newgroup->levelofdetail-1]=elements_ellipses;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void MC_Config::SetEllipsoidCharges(Element_Group* group, Element_Group* newgroup)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	cout << "\t-> Creating charge/dipole distribution for <" << newgroup->Type->name << ">\n";
	double* qs = new double[newgroup->nr_elements];
	double* wcsum = new double[newgroup->nr_elements];
	
	for(unsigned int i=0; i<newgroup->nr_elements; i++){
		double reduce_es=newgroup->Type->LOD->reduce_electrostatics[newgroup->levelofdetail-1][i];
		qs[i]=0.0;
		wcsum[i]=0.0;
		Element* ellipse=&parameters.group_elements[newgroup->elements[i]];
		ellipse->MyType->initial_dipole=Vec3(0.0);
		if((fabs(reduce_es-1.0)<EPS) || (fabs(reduce_es-2.0)<EPS)){
			if((ellipse->MyType->nr_charges==0) && (fabs(ellipse->MyType->dipole)<EPS)){ // otherwise we got a copied element here
				ellipse->MyType->nr_charges=1;
				ellipse->MyType->q=(double*)realloc(ellipse->MyType->q,sizeof(double)*ellipse->MyType->nr_charges);
				if(!ellipse->MyType->q){
					cout << "Not enough memory for charges.\n";
					exit(3);
				}
				ellipse->MyType->q[0]=0.0;
				ellipse->MyType->q_pos=(Vec3*)realloc(ellipse->MyType->q_pos,sizeof(Vec3)*ellipse->MyType->nr_charges);
				if(!ellipse->MyType->q_pos){
					cout << "Not enough memory for charge positions.\n";
					exit(3);
				}
				ellipse->MyType->q_pos[0]=Vec3(0.0);
			}
		}
	}
	Vec3 original_dipole=Vec3(0.0);
	for(unsigned int i=0; i<group->nr_elements; i++){
		Element* element=&parameters.group_elements[group->elements[i]];
		original_dipole+=element->rot*element->MyType->initial_dipole;
		for(unsigned int j=0; j<element->MyType->nr_charges; j++){
			int el=newgroup->Type->LOD->element_in_ellipsoid[newgroup->levelofdetail-1][i];
			wcsum[el]+=fabs(element->MyType->q[j]);
			qs[el]+=element->MyType->q[j];
			original_dipole+=(element->center+element->rot*element->MyType->q_pos[j])*element->MyType->q[j];
		}
	}
	Vec3 r;
	double factor=1.0;
	bool neutralize=false;
	Mat33* quads=new Mat33[newgroup->nr_elements];
	unsigned int* nr_charges_per_ellipsoid=new unsigned int[newgroup->nr_elements];
	for(unsigned int i=0; i<newgroup->nr_elements; i++){
		quads[i].M3Zeros();
		nr_charges_per_ellipsoid[i]=0;
	}
	for(unsigned int i=0; i<group->nr_elements; i++){
		Element* element=&parameters.group_elements[group->elements[i]];
		unsigned int newnr=newgroup->Type->LOD->element_in_ellipsoid[newgroup->levelofdetail-1][i];
		Element* ellipse=&parameters.group_elements[newgroup->elements[newnr]];
		double reduce_es=newgroup->Type->LOD->reduce_electrostatics[newgroup->levelofdetail-1][newgroup->Type->LOD->element_in_ellipsoid[newgroup->levelofdetail-1][i]];
		if(element->mytype!=ellipse->mytype){
			ellipse->MyType->initial_dipole+=element->MyType->initial_dipole;
			for(unsigned int j=0; j<element->MyType->nr_charges; j++){ // assign new charges/dipoles to each ellipsoid
				if(fabs(element->MyType->q[j])>EPS) nr_charges_per_ellipsoid[newnr]++;
				switch((int)reduce_es){
					default:
					case -1:
						if(reduce_es<0.0) factor=fabs(reduce_es);
						if(factor-1.0>EPS){
							if(newgroup->nr_elements>1) neutralize=true;
						}
						ellipse->MyType->nr_charges++;
						ellipse->MyType->q=(double*)realloc(ellipse->MyType->q,sizeof(double)*ellipse->MyType->nr_charges);
						if(!ellipse->MyType->q){
							cout << "Not enough memory for charges.\n";
							exit(3);
						}
						if(neutralize){
							ellipse->MyType->q[ellipse->MyType->nr_charges-1]=(element->MyType->q[j]-fabs(element->MyType->q[j])/wcsum[newgroup->Type->LOD->element_in_ellipsoid[newgroup->levelofdetail-1][i]]*qs[newgroup->Type->LOD->element_in_ellipsoid[newgroup->levelofdetail-1][i]])*factor;
						} else ellipse->MyType->q[ellipse->MyType->nr_charges-1]=element->MyType->q[j]*factor;
						ellipse->MyType->q_pos=(Vec3*)realloc(ellipse->MyType->q_pos,sizeof(Vec3)*ellipse->MyType->nr_charges);
						if(!ellipse->MyType->q_pos){
							cout << "Not enough memory for charge positions.\n";
							exit(3);
						}
						ellipse->MyType->q_pos[ellipse->MyType->nr_charges-1]=(ellipse->rot.M3Transpose()*(element->center+(element->rot*element->MyType->q_pos[j])-ellipse->center))/factor;
						break;
					case 0: // only dipole, ignore left-over charges
						neutralize=true; // need to make sure result is neutral in this case as well
						r=ellipse->rot.M3Transpose()*(element->center+(element->rot*element->MyType->q_pos[j])-ellipse->center);
						ellipse->MyType->initial_dipole+=r*(element->MyType->q[j]-fabs(element->MyType->q[j])/wcsum[newgroup->Type->LOD->element_in_ellipsoid[newgroup->levelofdetail-1][i]]*qs[newgroup->Type->LOD->element_in_ellipsoid[newgroup->levelofdetail-1][i]]);
//						ellipse->MyType->initial_dipole+=r*element->MyType->q[j];
						break;
					case 1:
					case 2:
						ellipse->MyType->q[0]+=element->MyType->q[j];
						r=ellipse->rot.M3Transpose()*(element->center+(element->rot*element->MyType->q_pos[j])-ellipse->center);
						ellipse->MyType->initial_dipole+=r*element->MyType->q[j];
						quads[newnr].mat[0][0]+=element->MyType->q[j]*(3.0*r.vec[0]*r.vec[0]-(r*r));
						quads[newnr].mat[0][1]+=element->MyType->q[j]*(3.0*r.vec[0]*r.vec[1]);
						quads[newnr].mat[0][2]+=element->MyType->q[j]*(3.0*r.vec[0]*r.vec[2]);
						quads[newnr].mat[1][1]+=element->MyType->q[j]*(3.0*r.vec[1]*r.vec[1]-(r*r));
						quads[newnr].mat[1][2]+=element->MyType->q[j]*(3.0*r.vec[1]*r.vec[2]);
						quads[newnr].mat[2][2]+=element->MyType->q[j]*(3.0*r.vec[2]*r.vec[2]-(r*r));
						break;
				}
			}
		}
	}
	Vec3 missing_dipole=Vec3(0.0);
	double wsum=0.0;
	if(neutralize){
		// calculate missing dipole moment
		missing_dipole=original_dipole;
		for(unsigned int i=0; i<newgroup->nr_elements; i++){
			Element* ellipse=&parameters.group_elements[newgroup->elements[i]];
			missing_dipole-=ellipse->rot*ellipse->MyType->initial_dipole;
			for(unsigned int j=0; j<ellipse->MyType->nr_charges; j++){
				missing_dipole-=(ellipse->center+ellipse->rot*ellipse->MyType->q_pos[j])*ellipse->MyType->q[j];
			}
		}
		for(unsigned int i=0; i<newgroup->nr_elements; i++){
			Element* ellipse=&parameters.group_elements[newgroup->elements[i]];
			Vec3 dipole=(ellipse->rot*ellipse->MyType->initial_dipole);
			for(unsigned int j=0; j<ellipse->MyType->nr_charges; j++){
				dipole+=(ellipse->center+ellipse->rot*ellipse->MyType->q_pos[j])*ellipse->MyType->q[j];
			}
			wsum+=fabs(dipole*missing_dipole);
		}
	}
#if DEBUG_LEVEL>2
	cout << "(" << original_dipole.V3Str(',') << "), (" << missing_dipole.V3Str(',') << ")\n";
#endif
	Vec3 new_dipole=Vec3(0.0);
	for(unsigned int i=0; i<newgroup->nr_elements; i++){
		Element* ellipse=&parameters.group_elements[newgroup->elements[i]];
		new_dipole+=ellipse->rot*ellipse->MyType->initial_dipole;
		for(unsigned int j=0; j<ellipse->MyType->nr_charges; j++) new_dipole+=(ellipse->center+ellipse->rot*ellipse->MyType->q_pos[j])*ellipse->MyType->q[j];
	}
	double prev_reduce=newgroup->Type->LOD->reduce_electrostatics[newgroup->levelofdetail-1][0];
	bool info_shown=false;
	for(unsigned int i=0; i<newgroup->nr_elements; i++){
		bool calculate_quadrupole=false;
		double reduce_es=newgroup->Type->LOD->reduce_electrostatics[newgroup->levelofdetail-1][i];
		Element* ellipse=&parameters.group_elements[newgroup->elements[i]];
#if DEBUG_LEVEL>2
		cout << "(" << (ellipse->rot*ellipse->MyType->initial_dipole).V3Str(',') << ") -> (";
#endif
		if(neutralize){
			Vec3 dipole=(ellipse->rot*ellipse->MyType->initial_dipole);
			for(unsigned int j=0; j<ellipse->MyType->nr_charges; j++){
				dipole+=(ellipse->center+ellipse->rot*ellipse->MyType->q_pos[j])*ellipse->MyType->q[j];
			}
			ellipse->MyType->initial_dipole+=ellipse->rot.M3Transpose()*missing_dipole*(fabs(dipole*missing_dipole)/wsum);
		}
#if DEBUG_LEVEL>2
		cout << (ellipse->rot*ellipse->MyType->initial_dipole).V3Str(',') << ")\n";
#endif
		if(ellipse->MyType->initial_dipole*ellipse->MyType->initial_dipole<EPS){ // if dipole is too small ignore it
			ellipse->MyType->initial_dipole=Vec3(0.0);
		}
		ellipse->dipole=ellipse->rot*ellipse->MyType->initial_dipole;
		ellipse->MyType->dipole=ellipse->MyType->initial_dipole.V3Norm();
		ellipse->MyType->hasmu=(fabs(ellipse->MyType->dipole)>EPS);
		if(ellipse->MyType->nr_charges==1){
			if(fabs(ellipse->MyType->q[0]/e_in_esu)<EPS){ // no charge after all
				ellipse->MyType->nr_charges=0;
				free(ellipse->MyType->q);
				ellipse->MyType->q=NULL;
				free(ellipse->MyType->q_pos);
				ellipse->MyType->q_pos=NULL;
			}
		}
		switch((int)reduce_es){
			default:
			case -1:
				if(reduce_es<0.0) factor=fabs(reduce_es);
				if(factor-1.0>EPS){
					if((fabs(prev_reduce-reduce_es)>EPS) || !info_shown){
						cout << "\t\t-> Scaling original charge distances by a factor of " << factor << "\n";
					}
					if(newgroup->nr_elements>1){
						neutralize=true;
						if((fabs(prev_reduce-reduce_es)>EPS) || !info_shown){
							cout << "\t\t\t-> Need to neutralize ellipsoid charges to maintain overall dipole moment.\n\t\t\t   This also may create a an additional dipole per ellipsoid.\n";
						}
					}
				} else{
					if((fabs(prev_reduce-reduce_es)>EPS) || !info_shown){
						cout << "\t\t-> Keeping original charges\n";
					}
				}
				if(!newgroup->Type->group_dipole){
					if((fabs(prev_reduce-reduce_es)>EPS) || !info_shown){
						cout << "\t\t\t-> Switching to dipole calculation relative to group.\n";
					}
					newgroup->Type->group_dipole=true;
				}
				break;
			case 0: // only dipole, ignore left-over charges
				if((fabs(prev_reduce-reduce_es)>EPS) || !info_shown){
					cout << "\t\t-> Reducing charges to a center dipole (adjusted to match overall group dipole).\n";
				}
				break;
			case 1:
				if((fabs(prev_reduce-reduce_es)>EPS) || !info_shown){
					cout << "\t\t-> Reducing charges to a center dipole and a left-over charge (for non-neutral ellipsoids).\n";
				}
				break;
			case 2:
				calculate_quadrupole=true;
				if((fabs(prev_reduce-reduce_es)>EPS) || !info_shown){
					cout << "\t\t-> Reducing charges to a center dipole, a quadrupole represented by charges, and a left-over charge (for non-neutral ellipsoids).\n";
				}
				break;
		}
		if(calculate_quadrupole){
			if(nr_charges_per_ellipsoid[i]>4){
#if DEBUG_LEVEL>2
				cout << "\t-> Ellipsoid " << i+1 << ":\n";
#endif
				quads[i].mat[1][0]=quads[i].mat[0][1];
				quads[i].mat[2][0]=quads[i].mat[0][2];
				quads[i].mat[2][1]=quads[i].mat[1][2];
				quads[i]*=0.5;
#if DEBUG_LEVEL>2
				cout << "\t\t-> Tensor:\n" << quads[i].M3Str() << "\n";
#endif
				CVec3 cew=quads[i].Eigenvalues();
				if(cew.Im()*cew.Im()>EPS){ // Quadrupol tensor eigenvalues need to be real
					cout << "ERROR: Quadropol tensor eigenvalues need to be real.\n";
					cout << "Quadrupol tensor:\n";
					cout << quads[i].M3Str() << "\n";
					cout << "Complex eigenvalues: " << cew.CV3Str(',') << "\n";
					exit(4);
				}
				Vec3 ew=cew.Re();
				Mat33 qrot=RotFromEigenvectors(quads[i].Eigenvectors(ew));
#if DEBUG_LEVEL>2
				cout << "\t\t-> Eigenvalues: " << ew.V3Str(',') << "\n";
				cout << "\t\t-> Rotation matrix:\n";
				cout << qrot.M3Str() << "\n";
#endif
				// Determine charges and charge locations
				//     ( 2q_x*a²-q_yb²-q_zc²   0             0           )
				// Q = (        0      2q_y*b²-q_xa²-q_zc²   0           )
				//     (        0                0   2q_z*c²-q_xa²-q_yb² )
				// q_x + q_y + q_z = 0
				// Let eigenvalues and maximum eigenvalue magnitude determine relative charge magnitudes:
				// q_i = lambda_i/|lambda_ref|*q_exp <- q_exp is expansion charge used for quadrupol
				//
				// Now need to solve for a², b², c²:
				// 2q_x*a²-q_yb²-q_zc² = lambda_x
				// 2q_y*b²-q_xa²-q_zc² = lambda_y
				// 2q_z*c²-q_xa²-q_yb² = lambda_z
				//
				//    ( 2 -1 -1) (a²*q_x) (lambda_x)
				// => (-1  2 -1)*(b²*q_y)=(lambda_y)
				//    (-1 -1  2) (c²*q_z) (lambda_z)
				Mat33 dists;
				dists.mat[0][0]=2.0;
				dists.mat[1][0]=-1.0;
				dists.mat[2][0]=-1.0;
				
				dists.mat[0][1]=-1.0;
				dists.mat[1][1]=2.0;
				dists.mat[2][1]=-1.0;
				
				dists.mat[0][2]=-1.0;
				dists.mat[1][2]=-1.0;
				dists.mat[2][2]=2.0;
				
				Vec3 solution=ew;
				dists.M3GEPP(solution);
				unsigned int l_max=0;
				for(unsigned int j=1; j<3; j++){
					if(fabs(ew.vec[j])>fabs(ew.vec[l_max])) l_max=j;
				}
				BackSubstitute(dists,solution,l_max);
				// now go over solution and determine where charges are located
				double qsum=0.0;
				double q, pos;
				Vec3 position;
				cout << "\t\t\t-> Ellipsoid " << i+1 << ": Placing quadrupole charges (quadrupole diagonal elements: " << (ew/e_in_esu).V3Str(',') << " eAng^2)\n";
				if(ellipse->MyType->nr_charges==0){
					ellipse->MyType->nr_charges=1;
					ellipse->MyType->q=(double*)realloc(ellipse->MyType->q,sizeof(double)*ellipse->MyType->nr_charges);
					if(!ellipse->MyType->q){
						cout << "Not enough memory for charges.\n";
						exit(3);
					}
					ellipse->MyType->q[0]=0.0;
					ellipse->MyType->q_pos=(Vec3*)realloc(ellipse->MyType->q_pos,sizeof(Vec3)*ellipse->MyType->nr_charges);
					if(!ellipse->MyType->q_pos){
						cout << "Not enough memory for charge positions.\n";
						exit(3);
					}
					ellipse->MyType->q_pos[0]=Vec3(0.0);
				}
				Mat33 newquad(0.0);
				for(unsigned int j=0; j<3; j++){
					if(fabs(solution.vec[j])>EPS){
						q=parameters.quadrupole_charge_factor*fabs(ew.vec[j]/ew.vec[l_max])*solution.vec[j]/fabs(solution.vec[j]); // sign of q is determined by solution
						pos=sqrt(solution.vec[j]/q);
						if(pos*pos>1E-4){ // only interested in contributions to quadrupol tensor that are large enough
							ellipse->MyType->nr_charges+=2;
							ellipse->MyType->q=(double*)realloc(ellipse->MyType->q,sizeof(double)*ellipse->MyType->nr_charges);
							if(!ellipse->MyType->q){
								cout << "Not enough memory for charges.\n";
								exit(3);
							}
							ellipse->MyType->q[ellipse->MyType->nr_charges-2]=q;
							ellipse->MyType->q[ellipse->MyType->nr_charges-1]=q;
							ellipse->MyType->q_pos=(Vec3*)realloc(ellipse->MyType->q_pos,sizeof(Vec3)*ellipse->MyType->nr_charges);
							if(!ellipse->MyType->q_pos){
								cout << "Not enough memory for charge positions.\n";
								exit(3);
							}
							qsum+=2.0*q;
							cout << "\t\t\t\t-> Two charges (" << q/e_in_esu << " e) at +/- " << pos << " Angström in ";
							switch(j){
								case 0: cout << "x";
									position=Vec3(1.0,0.0,0.0)*pos;
									break;
								case 1: cout << "y";
									position=Vec3(0.0,1.0,0.0)*pos;
									break;
								case 2: cout << "z";
									position=Vec3(0.0,0.0,1.0)*pos;
									break;
							}
							position=qrot*position;
							ellipse->MyType->q_pos[ellipse->MyType->nr_charges-2]=position;
							ellipse->MyType->q_pos[ellipse->MyType->nr_charges-1]=position*(-1.0);
							cout << "-direction (in ellipsoid frame: " << position.V3Str(',') << ").\n";
							if(position.V3Norm()>0.8*EllipsoidRmin(position,ellipse->MyType->saxes,Mat33())){
								cout << "ERROR: Charges are located too close to the ellipsoid surface.\nIncrease \"quadrupole_charge_factor\" in [Simulation Parameters] to compensate.\n";
								exit(8);
							}
							newquad.mat[0][0]+=2.0*q*(3.0*position.vec[0]*position.vec[0]-(position*position));
							newquad.mat[0][1]+=2.0*q*(3.0*position.vec[0]*position.vec[1]);
							newquad.mat[0][2]+=2.0*q*(3.0*position.vec[0]*position.vec[2]);
							newquad.mat[1][1]+=2.0*q*(3.0*position.vec[1]*position.vec[1]-(position*position));
							newquad.mat[1][2]+=2.0*q*(3.0*position.vec[1]*position.vec[2]);
							newquad.mat[2][2]+=2.0*q*(3.0*position.vec[2]*position.vec[2]-(position*position));
						}
					}
				}
				// center charge does not contribute to quadrupole
				ellipse->MyType->q[0]+=-qsum;
				cout << "\t\t\t\t-> Charge (" << -qsum/e_in_esu << " e) at center.\n";
				newquad.mat[1][0]=newquad.mat[0][1];
				newquad.mat[2][0]=newquad.mat[0][2];
				newquad.mat[2][1]=newquad.mat[1][2];
				newquad*=0.5;
#if DEBUG_LEVEL>2
				cout << "New quadrupol tensor:\n" << newquad.M3Str() << "\n";
#endif
				if(!M3equal(newquad,quads[i],2*nr_charges_per_ellipsoid[i]*parameters.quadrupole_charge_factor*sqrt(EPS))){
					cout << "ERROR: Charges do not represent quadrupol moment properly within error of " << 2*nr_charges_per_ellipsoid[i]*parameters.quadrupole_charge_factor*sqrt(EPS) << ".\n";
					cout << "Original quadrupol tensor:\n" << quads[i].M3Str() << "\nDifference tensor:\n" << (quads[i]-newquad).M3Str() << "\n";
					exit(3);
				}
			} else{
				cout << "\t\t\t-> Ellipsoid " << i+1 << ":\n\t\t\t\t-> Not enough underlying charges to warrant a quadrupole.\n";
			}
		}
		info_shown=true;
		prev_reduce=reduce_es;
		cout << "\t\t\t-> Ellipsoid " << i+1 << " has " << ellipse->MyType->nr_charges;
		if(ellipse->MyType->nr_charges==1) cout << " charge " << "(" << ellipse->MyType->q[0]/e_in_esu << " e) "; else cout << " charges ";
		if(ellipse->MyType->hasmu) cout << "and a dipole of " << ellipse->MyType->dipole << " Debye (" << ellipse->dipole.V3Str(',') << ").\n"; else cout << "and no dipole.\n";
	}
	delete[] quads;
	delete[] nr_charges_per_ellipsoid;
	delete[] qs;
	delete[] wcsum;
	cout << "\t<- Done.\n";
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void MC_Config::LoadLevelofDetail(char* confbase, Element_Group* group, int nr)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	string storage_filename=parameters.fileout+".conf.dat";
	readfile input;
	input.filename=storage_filename;
	input.directory="";
	input.file.open(storage_filename.c_str(), fstream::in);
	char* conf=confbase;
	if(!input.file.fail()){
		unsigned int pos=0;
		string subname=group->Type->name;
		char* section=GetSection(&input,"Detail",pos,subname);
		if(section){
			conf=new char[strlen(section)+strlen(confbase)+1];
			conf[strlen(section)+strlen(confbase)]='\0';
			strncpy(conf,section,strlen(section));
			char* temp=conf+strlen(section);
			strncpy(temp,confbase,strlen(confbase));
			delete[] section;
		}
		input.file.close();
	}
	cout << "-> Check for user-specified components:\n";
	double* components=get_flex_tupel("components",conf);
	unsigned int count=0;
	if(components){
		parameters.lods[nr]->nr_components=0;
		while(components[count]>0){
			count+=(unsigned int)components[count]+1;
			parameters.lods[nr]->nr_components++;
		}
		if(parameters.lods[nr]->nr_components==0){
			cout << "No components specified after keyword.\n";
			exit(4);
		}
		if(parameters.lods[nr]->nr_components>max_nr_components) max_nr_components=parameters.lods[nr]->nr_components;
		parameters.lods[nr]->nr_elements=new unsigned int[parameters.lods[nr]->nr_components];
		parameters.lods[nr]->component_elements=new unsigned int*[parameters.lods[nr]->nr_components];
		parameters.lods[nr]->component_show_dipole=new bool[parameters.lods[nr]->nr_components];
		parameters.lods[nr]->component_color=new Vec3[parameters.lods[nr]->nr_components];
		parameters.lods[nr]->component_transparency=new double[parameters.lods[nr]->nr_components];
		parameters.lods[nr]->component_texture=new string[parameters.lods[nr]->nr_components];
		parameters.lods[nr]->component_texture_is_dipole=new bool[parameters.lods[nr]->nr_components];
		parameters.lods[nr]->element_in_component=new int[group->nr_elements];
		memset(parameters.lods[nr]->element_in_component,0xFF,group->nr_elements*sizeof(int)); // puts -1 in each field
		count=0;
		unsigned int current=0;
		while(components[count]>0){
			SetParam(parameters.lods[nr]->component_show_dipole[current],("component_show_dipole."+int2str(current+1)).c_str(),conf,false);
			SetParam(parameters.lods[nr]->component_color[current].vec,("component_color."+int2str(current+1)).c_str(),conf,Vec3(1.0));
			SetParam(parameters.lods[nr]->component_transparency[current],("component_transparency."+int2str(current+1)).c_str(),conf,-1.0);
			SetParam(parameters.lods[nr]->component_texture[current],("component_texture."+int2str(current+1)).c_str(),conf,"");
			SetParam(parameters.lods[nr]->component_texture_is_dipole[current],("component_texture_is_dipole."+int2str(current+1)).c_str(),conf,parameters.lods[nr]->component_show_dipole[current]);
			parameters.lods[nr]->nr_elements[current]=(unsigned int)components[count];
			parameters.lods[nr]->component_elements[current]=new unsigned int[parameters.lods[nr]->nr_elements[current]];
			for(unsigned int i=0; i<parameters.lods[nr]->nr_elements[current]; i++){
				parameters.lods[nr]->component_elements[current][i]=(unsigned int)components[count+i+1]-1;
				int curr=(int)components[count+i+1]-1;
				if((curr<0) || (curr>(int)group->nr_elements)){
					cout << "Specified elements need to be between 1 and the total number of elements in the group.\n";
					exit(2);
				}
				parameters.lods[nr]->element_in_component[(unsigned int)curr]=current;
			}
			current++;
			count+=(unsigned int)components[count]+1;
		}
		group->Type->LOD=parameters.lods[nr]; // group only gets an LOD associated if user specifies components ...
		delete[] components;
		// If components are specified we can also use them to specify which dipole to use for order calculation
		SetParam(parameters.lods[nr]->order_dipole_component,"order_dipole_component",conf,0); // default is 0 (whole group)
	} else{
		parameters.lods[nr]->order_dipole_component=0;
		parameters.lods[nr]->nr_components=0;
		parameters.lods[nr]->nr_elements=NULL;
		parameters.lods[nr]->component_elements=NULL;
		parameters.lods[nr]->component_show_dipole=NULL;
		parameters.lods[nr]->component_color=NULL;
		parameters.lods[nr]->component_transparency=NULL;
		parameters.lods[nr]->component_texture=NULL;
	}
	cout << "\t-> " << parameters.lods[nr]->nr_components << " component(s) defined.\n";
	SetParam(parameters.lods[nr]->levels,"levels",conf,0);
	if(parameters.lods[nr]->levels>max_levelofdetail) max_levelofdetail=parameters.lods[nr]->levels;
	if(parameters.lods[nr]->levels==0){
		parameters.lods[nr]->groups=new Element_Group*;
		parameters.lods[nr]->groups[0]=group;
		cout << "-> No levels of detail specified.\n";
	} else{
		group->Type->LOD=parameters.lods[nr]; // ... or level of detail groups
		double rT=parameters.rT;
		SetParam(parameters.lods[nr]->start_level,"start_level",conf,0);
		SetParam(parameters.lods[nr]->end_level,"end_level",conf,parameters.lods[nr]->levels);
		SetParam(parameters.lods[nr]->gyration_sphere,"gyration_sphere",conf,false);
		SetParam(parameters.lods[nr]->symmetry_center,"symmetry_center",conf,false);
		SetParam(parameters.lods[nr]->symmetry_axes,"symmetry_axes",conf,true);
		SetParam(parameters.lods[nr]->match_volume,"match_volume",conf,true);
		SetParam(parameters.lods[nr]->match_epsilon,"match_epsilon",conf,true);
		SetParam(parameters.lods[nr]->scale_original_vdw,"scale_original_vdw",conf,1.0);
		SetParam(parameters.lods[nr]->show_internal_bonds,"show_internal_bonds",conf,false);
		SetParam(parameters.lods[nr]->internal_charge_colors,"internal_charge_colors",conf,false);
		SetParam(rT,"rT",conf,rT);
		parameters.lods[nr]->reduce_electrostatics=new double*[parameters.lods[nr]->levels];
		
		parameters.lods[nr]->distances = new double[parameters.lods[nr]->levels+1];
		parameters.lods[nr]->nr_ellipsoids = new unsigned int[parameters.lods[nr]->levels];
		parameters.lods[nr]->ellipsoid_counts = new unsigned int*[parameters.lods[nr]->levels];
		parameters.lods[nr]->element_groups = new int*[parameters.lods[nr]->levels];
		parameters.lods[nr]->element_in_ellipsoid = new int*[parameters.lods[nr]->levels];
		parameters.lods[nr]->volumes = new double*[parameters.lods[nr]->levels];
		parameters.lods[nr]->epsilons = new double*[parameters.lods[nr]->levels];
		parameters.lods[nr]->placement_deltas = new Vec3[parameters.lods[nr]->levels];
		parameters.lods[nr]->groups = new Element_Group*[parameters.lods[nr]->levels+1];
		
		parameters.lods[nr]->visual_original = new bool[parameters.lods[nr]->levels];
		parameters.lods[nr]->keep_original_potentials = new bool[parameters.lods[nr]->levels];
		parameters.lods[nr]->use_epsilon_texture = new bool[parameters.lods[nr]->levels];
		parameters.lods[nr]->transparency = new double[parameters.lods[nr]->levels];
		parameters.lods[nr]->inside_transparency = new double[parameters.lods[nr]->levels];
		
		parameters.lods[nr]->groups[0]=group;
		for(unsigned int i=0; i<parameters.lods[nr]->levels; i++){
			cout << "-> Reading settings for detail level " << i+1 << "\n";
			SetParam(parameters.lods[nr]->distances[i],("lod_distance."+int2str(i+1)).c_str(),conf,0.0);
			SetParam(parameters.lods[nr]->visual_original[i],("visual_original."+int2str(i+1)).c_str(),conf,true);
			SetParam(parameters.lods[nr]->keep_original_potentials[i],("keep_original_potentials."+int2str(i+1)).c_str(),conf,false);
			SetParam(parameters.lods[nr]->use_epsilon_texture[i],("use_epsilon_texture."+int2str(i+1)).c_str(),conf,false);
			SetParam(parameters.lods[nr]->transparency[i],("transparency."+int2str(i+1)).c_str(),conf,group->Type->transparency);
			SetParam(parameters.lods[nr]->inside_transparency[i],("inside_transparency."+int2str(i+1)).c_str(),conf,0.0);
			
			double* element_groups=get_flex_tupel(("lod_partitions."+int2str(i+1)).c_str(),conf);
			if(!element_groups){
				cout << "<lod_partitions."+int2str(i+1) << "> not defined.\n";
				exit(2);
			}
			count=0;
			unsigned nr_ellipsoids=0;
			while(element_groups[count]>0){
				count+=(unsigned int)element_groups[count]+1;
				nr_ellipsoids++;
			}
			parameters.lods[nr]->element_groups[i]=new int[count+1];
			for(unsigned int j=0; j<=count; j++) parameters.lods[nr]->element_groups[i][j]=(int)element_groups[j];
			delete[] element_groups;
			parameters.lods[nr]->reduce_electrostatics[i]=new double[nr_ellipsoids];
			for(unsigned int j=0; j<nr_ellipsoids; j++) parameters.lods[nr]->reduce_electrostatics[i][j]=-1.0;
			double* ellipsoid_reduce_es=get_flex_tupel(("reduce_electrostatics."+int2str(i+1)).c_str(),conf);
			if(ellipsoid_reduce_es){
				unsigned int curr_e=0;
				count=0;
				cout << "\t-> Using the following charge reduction parameters: ";
				while((ellipsoid_reduce_es[count]>0) && (curr_e<nr_ellipsoids)){
					for(unsigned int j=0; j<ellipsoid_reduce_es[count]; j++){
						if(curr_e<nr_ellipsoids){
							parameters.lods[nr]->reduce_electrostatics[i][curr_e]=ellipsoid_reduce_es[count+j+1];
							if(curr_e>0) cout << ", ";
							cout << parameters.lods[nr]->reduce_electrostatics[i][curr_e];
						}
						curr_e++;
					}
					count+=(unsigned int)ellipsoid_reduce_es[count]+1;
				}
				cout << "\n";
				if(curr_e==1){ // Apply value to all ellipsoids
					cout << "NOTE: Defined charge reduction parameter is applied to all ellipsoids.\n";
					for(unsigned int j=1; j<nr_ellipsoids; j++) parameters.lods[nr]->reduce_electrostatics[i][j]=parameters.lods[nr]->reduce_electrostatics[i][0];
				} else{
					if(curr_e<nr_ellipsoids) cout << "WARNING: Did only specify " << curr_e << " charge reduction parameters but " << nr_ellipsoids << " ellipsoids exist. Using -1 as default for leftover ellipsoids.\n";
					if(curr_e>nr_ellipsoids) cout << "WARNING: Too many charge reduction parameters specified for " << nr_ellipsoids << " ellipsoids. Ignoring additionals.\n";
				}
				delete[] ellipsoid_reduce_es;
			} else cout << "NOTE: No charge reduction parameter specified, default -1 is used for all ellipsoids.\n";
			parameters.lods[nr]->volumes[i]=new double[nr_ellipsoids];
			for(unsigned int j=0; j<nr_ellipsoids; j++) parameters.lods[nr]->volumes[i][j]=-1.0;
			double* ellipsoid_volumes=get_MxN_tupel(("lod_volumes."+int2str(i+1)).c_str(),conf,nr_ellipsoids,1,NULL);
			if(ellipsoid_volumes){
				for(unsigned int j=0; j<nr_ellipsoids; j++) parameters.lods[nr]->volumes[i][j]=ellipsoid_volumes[j];
				delete[] ellipsoid_volumes;
			}
			parameters.lods[nr]->epsilons[i]=new double[nr_ellipsoids];
			for(unsigned int j=0; j<nr_ellipsoids; j++) parameters.lods[nr]->epsilons[i][j]=-1.0;
			double* ellipsoid_epsilons=get_MxN_tupel(("lod_epsilons."+int2str(i+1)).c_str(),conf,nr_ellipsoids,1,NULL);
			if(ellipsoid_epsilons){
				for(unsigned int j=0; j<nr_ellipsoids; j++) parameters.lods[nr]->epsilons[i][j]=ellipsoid_epsilons[j];
				delete[] ellipsoid_epsilons;
			}
			parameters.num_groups++;
			parameters.groups = (Element_Group**)realloc(parameters.groups,parameters.num_groups*sizeof(Element_Group*));
			if(parameters.groups==NULL){ // Danger, danger
				cout << "Could not find a memory block large enough to extend groups for level of detail.\n";
				exit(1);
			}
			Element_Group* newgroup=new Element_Group;
			newgroup->type=parameters.num_groups-1;
			newgroup->Type=new Element_Group_Type;
			newgroup->Type->name=group->Type->name+".LOD"+int2str(i+1);
			newgroup->Type->dof=group->Type->dof;
			newgroup->Type->allow_bond_bend=group->Type->allow_bond_bend;
			newgroup->Type->group_dipole=group->Type->group_dipole;
			newgroup->Type->rand_independent=group->Type->rand_independent;
			newgroup->Type->show_just_dipoles=group->Type->show_just_dipoles;
			newgroup->Type->group_dipole_color=group->Type->group_dipole_color;
			newgroup->Type->visual_PBCs=group->Type->visual_PBCs;
			newgroup->Type->transparency=parameters.lods[nr]->transparency[i];
			newgroup->Type->label_elements=group->Type->label_elements;
			newgroup->Type->density=group->Type->density;
			newgroup->Type->packing_density=group->Type->packing_density;
			newgroup->Type->bond_range=group->Type->bond_range;
			newgroup->Type->calculate_order=group->Type->calculate_order;
			SetParam(newgroup->Type->adjust_overlap,("adjust_overlap."+int2str(i+1)).c_str(),conf,0.5);
			for(unsigned int j=0; j<SPECIAL_DISTANCES; j++){
				newgroup->Type->range_factors[j].distance=group->Type->range_factors[j].distance;
				newgroup->Type->range_factors[j].factor=group->Type->range_factors[j].factor;
			}
			newgroup->Type->range_field=NULL;
			if(parameters.lods[nr]->start_level==i+1){
				newgroup->number=group->number;
				group->number=0;
			} else newgroup->number=0;
			parameters.lods[nr]->groups[i+1]=newgroup;
			newgroup->levelofdetail=i+1;
			newgroup->Type->LOD=parameters.lods[nr];
			parameters.groups[newgroup->type]=newgroup;
			newgroup->nr_elements=0;
			newgroup->elements=NULL;
			newgroup->nr_potentials=0;
			newgroup->potentials=NULL;
			count=0;
			parameters.lods[nr]->nr_ellipsoids[i]=0;
			parameters.lods[nr]->ellipsoid_counts[i]=(unsigned int*)malloc(sizeof(unsigned int));
			parameters.lods[nr]->ellipsoid_counts[i][0]=count;
			while(parameters.lods[nr]->element_groups[i][count]>0){
				unsigned int* combine_elements = new unsigned int[parameters.lods[nr]->element_groups[i][count]];
				for(unsigned j=0; j<(unsigned int)parameters.lods[nr]->element_groups[i][count]; j++){
					combine_elements[j]=parameters.lods[nr]->element_groups[i][count+j+1];
					if((combine_elements[j]>0) && (combine_elements[j]<=group->nr_elements)){
						combine_elements[j]=group->elements[combine_elements[j]-1];
					} else{
						cout << "Specified elements need to be between 1 and the total number of elements in the group.\n";
						exit(2);
					}
				}
				string specifier=":"+int2str(parameters.lods[nr]->nr_ellipsoids[i]+1);
				CombineElements(group,newgroup,combine_elements,parameters.lods[nr]->element_groups[i][count],specifier,parameters.lods[nr]->volumes[i][parameters.lods[nr]->nr_ellipsoids[i]],parameters.lods[nr]->epsilons[i][parameters.lods[nr]->nr_ellipsoids[i]],rT);
				delete[] combine_elements;
				count+=parameters.lods[nr]->element_groups[i][count]+1;
				parameters.lods[nr]->nr_ellipsoids[i]++;
				parameters.lods[nr]->ellipsoid_counts[i]=(unsigned int*)realloc(parameters.lods[nr]->ellipsoid_counts[i],(parameters.lods[nr]->nr_ellipsoids[i]+1)*sizeof(unsigned int));
				if(!parameters.lods[nr]->ellipsoid_counts[i]){
					cout << "Not enough memory for ellipse storage.\n";
					exit(2);
				}
				parameters.lods[nr]->ellipsoid_counts[i][parameters.lods[nr]->nr_ellipsoids[i]]=count;
			}
			SetEllipsoidLinks(group,newgroup,parameters.lods[nr]->element_groups[i],parameters.lods[nr]->ellipsoid_counts[i]);
			SetEllipsoidCharges(group,newgroup);
			if(newgroup->nr_elements>max_LOD_elements) max_LOD_elements=newgroup->nr_elements;
			SetParam(newgroup->Type->allow_bond_bend,("allow_bond_bend."+int2str(i+1)).c_str(),conf,group->Type->allow_bond_bend);
			AssignFixedElements(newgroup);
			SetParam(newgroup->Type->rand_elements,("rand_elements."+int2str(i+1)).c_str(),conf,(unsigned int)round((double)group->Type->rand_elements/(double)group->Type->nr_movable*newgroup->Type->nr_movable));
#if DEBUG_LEVEL>1
			cout << "-> Number of elements randomized for detail level " << i+1 << " group is: " << newgroup->Type->rand_elements << "\n";
#endif
			SetParam(newgroup->Type->nr_rand_per_cycle,("nr_rand_per_cycle."+int2str(i+1)).c_str(),conf,sqrt(newgroup->Type->nr_movable));
#if DEBUG_LEVEL>1
			cout << "-> Number of elements randomized on average per cycle for detail level " << i+1 << " group is: " << newgroup->Type->nr_rand_per_cycle << "\n";
#endif
			if(newgroup->Type->nr_movable>0){
				newgroup->Type->inv_sqrt_nrpc = sqrt(1.0/newgroup->Type->nr_rand_per_cycle);
				if(newgroup->Type->inv_sqrt_nrpc<EPS) newgroup->Type->inv_sqrt_nrpc=1.0;
			} else newgroup->Type->inv_sqrt_nrpc=1.0;
			if(!parameters.fit2lod && newgroup->Type->LOD->match_epsilon){
				// read in Lennard-Jones epsilon coefficients
				double* eps_coeff=get_flex_tupel(("lod_epsilon_of_rT."+int2str(i+1)).c_str(),conf);
				if(eps_coeff){
					cout << "-> Setting Lennard Jones epsilon coefficients for newly created LOD ellipsoids.\n";
					count=0;
					unsigned int j=0;
					while(eps_coeff[count]>=0){
						unsigned int e_nr=(unsigned int)eps_coeff[count];
						if(e_nr==parameters.group_elements[newgroup->elements[j]].MyType->nr_Vvdw_coefficients){
							if(j<newgroup->nr_elements){
								if(parameters.group_elements[newgroup->elements[j]].MyType->Vvdw_coefficients!=NULL) delete[] (parameters.group_elements[newgroup->elements[j]].MyType->Vvdw_coefficients);
								parameters.group_elements[newgroup->elements[j]].MyType->Vvdw_coefficients=new double[e_nr];
								cout << "\t-> Ellipsoid <" << parameters.group_elements[newgroup->elements[j]].MyType->name << ">: ";
								for(unsigned int l=0; l<e_nr; l++){
									parameters.group_elements[newgroup->elements[j]].MyType->Vvdw_coefficients[l]=eps_coeff[count+1+l];
									if(l>0) cout << ", ";
									cout << parameters.group_elements[newgroup->elements[j]].MyType->Vvdw_coefficients[l];
								}
								cout << "\n";
								count+=e_nr+1;
							} else{
								cout << "WARNING: Different number of LOD ellipoids in configuration.\nPlease make sure you run \"fit2lod <configuration file>\".\n";
								break;
							}
						} else{
							if(e_nr>0){ // Single spheres will not get coefficients which is OK
								cout << "WARNING: Number of epsilon coefficients is not in conformance with this code.\nPlease make sure you run \"fit2lod <configuration file>\".\n";
								break;
							}
						}
						j++;
					}
					delete[] eps_coeff;
				}
				double* IA_coeff=get_flex_tupel(("lod_IA."+int2str(i+1)).c_str(),conf);
				if(IA_coeff){
					cout << "-> Setting Interaction Area adjustment coefficients for newly created LOD ellipsoids.\n";
					count=0;
					unsigned int j=0;
					while(IA_coeff[count]>0){
						unsigned int e_nr=(unsigned int)IA_coeff[count];
						if(e_nr==7){
							if(j<newgroup->nr_elements){
								if(parameters.group_elements[newgroup->elements[j]].MyType->IA_coefficients!=NULL) delete[] (parameters.group_elements[newgroup->elements[j]].MyType->IA_coefficients);
								parameters.group_elements[newgroup->elements[j]].MyType->IA_coefficients=new double[e_nr];
								cout << "\t-> Ellipsoid <" << parameters.group_elements[newgroup->elements[j]].MyType->name << ">: ";
								for(unsigned int l=0; l<e_nr; l++){
									parameters.group_elements[newgroup->elements[j]].MyType->IA_coefficients[l]=IA_coeff[count+1+l];
									if(l>0) cout << ", ";
									cout << parameters.group_elements[newgroup->elements[j]].MyType->IA_coefficients[l];
								}
								cout << "\n";
								count+=e_nr+1;
							} else{
								cout << "WARNING: Different number of LOD ellipoids in configuration.\nPlease make sure you run \"fit2lod <configuration file>\".\n";
								break;
							}
						} else{
							cout << "WARNING: Number of IA coefficients is not in conformance with this code.\nPlease make sure you run \"fit2lod <configuration file>\".\n";
							break;
						}
						j++;
					}
					delete[] IA_coeff;
				}
				if(parameters.lods[nr]->use_epsilon_texture[i]){
					// set textures if there are any
					cout << "-> Setting epsilon textures for newly created LOD ellipsoids (may take a while) ...\n";
					double* textures=get_flex_tupel(("lod_textures."+int2str(i+1)).c_str(),conf);
					if(textures){
						count=0;
						unsigned int j=0;
						while(textures[count]>0){
							unsigned int e_nr=(unsigned int)textures[count];
							if(e_nr==phi_res*theta_res+4){
								if(j<newgroup->nr_elements){
									cout << "\t-> Setting texture for ellipsoid <" << parameters.group_elements[newgroup->elements[j]].MyType->name << "> ... ";
									cout.flush();
									parameters.group_elements[newgroup->elements[j]].MyType->eps_texture=new double[e_nr];
									for(unsigned int l=0; l<e_nr; l++) parameters.group_elements[newgroup->elements[j]].MyType->eps_texture[l]=textures[count+1+l];
									cout << "Done\n";
									count+=e_nr+1;
								} else{
									cout << "WARNING: Different number of LOD ellipoids in configuration.\nPlease make sure you run \"fit2lod <configuration file>\".\n";
									break;
								}
							} else{
								cout << "WARNING: Number of pixels in epsilon texture is not in conformance with this code.\nPlease make sure you run \"fit2lod <configuration file>\".\n";
								break;
							}
							j++;
						}
						delete[] textures;
					}
				}
			}
			cout << "<- Done.\n";
			// Now is a good time for that ...
			newgroup->Type->rot_notrans=group->Type->rot_notrans;
			newgroup->Type->still=group->Type->still;
			newgroup->Type->mass=0.0;
			for(unsigned int j=0; j<newgroup->nr_elements; j++) newgroup->Type->mass+=parameters.group_elements[newgroup->elements[j]].MyType->mass;
			cout << "-> Mass of detail level " << i+1 << " group is " << newgroup->Type->mass << " g/mol\n";
			SetGroupRangeField(newgroup,parameters.group_elements);
			DetermineGroupRings(newgroup,parameters.group_elements);
			newgroup->Type->Volume=-1;
			
			// Calculate center of group
			parameters.group_centers=(Vec3*)realloc(parameters.group_centers,parameters.num_groups*sizeof(Vec3));
			if(!parameters.group_centers){
				cout << "Not enough memory for new group centers.\n";
				exit(2);
			}
			parameters.group_centers[newgroup->type] = GroupCenter(newgroup,parameters.group_elements);
#if DEBUG_LEVEL>1
			cout << "-> Center (cartesian average of element centers) of detail level " << i+1 << " group is: (" << parameters.group_centers[newgroup->type].vec[0] << ", " << parameters.group_centers[newgroup->type].vec[1] << ", " << parameters.group_centers[newgroup->type].vec[2] << ")\n";
#endif
			// Calculate external sphere
			parameters.group_radii=(double*)realloc(parameters.group_radii,parameters.num_groups*sizeof(double));
			if(!parameters.group_radii){
				cout << "Not enough memory for new group radii.\n";
				exit(2);
			}
			parameters.group_radii[newgroup->type] = GroupRadius(newgroup,parameters.group_elements,parameters.group_centers[newgroup->type]);
#if DEBUG_LEVEL>1
			cout << "-> Effective radius of bounding sphere able to contain detail level " << i+1 << " group is: " << parameters.group_radii[newgroup->type] << " Angström\n";
#endif
			parameters.lods[nr]->placement_deltas[i]=parameters.group_centers[newgroup->type]-parameters.group_centers[group->type];
#if DEBUG_LEVEL>1
			cout << "-> Placement difference to original is: (" << parameters.lods[nr]->placement_deltas[i].vec[0] << ", " << parameters.lods[nr]->placement_deltas[i].vec[1] << ", " << parameters.lods[nr]->placement_deltas[i].vec[2] << ")\n";
#endif
		}
	}
	if(conf!=confbase) delete[] conf;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void MC_Config::LoadElementTypes(char* conf, int nr, bool load_charges)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	parameters.element_types[nr]->eps_texture=NULL;
	parameters.element_types[nr]->r_diff=NULL;
	parameters.element_types[nr]->width_diff=NULL;
	SetParam(parameters.element_types[nr]->avg_width,"avg_width",conf,0.0);
	parameters.element_types[nr]->sphericity=1.0;
	parameters.element_types[nr]->IA_average=1.0;
	string archetype="";
	SetParam(archetype,"archetype",conf,"");
	parameters.element_types[nr]->archetype=-1;
	if(archetype!=""){
		bool found=false;
		for(parameters.element_types[nr]->archetype=0; parameters.element_types[nr]->archetype<nr-1; parameters.element_types[nr]->archetype++){
			if(compare_strings(parameters.element_types[parameters.element_types[nr]->archetype]->name.c_str(),archetype.c_str())){
				found=true;
				break;
			}
		}
		if(!found){
			cout << "Could not find archetype <" << archetype << ">.\n";
			exit(3);
		}
	}
	Mat33 z;
	z.M3Zeros();
	// parameters.element_types[nr]->archetype now contains type nr of archetype (if parameters.element_types[nr]->archetype>=0)
	if(parameters.element_types[nr]->archetype<0){
		SetParam(parameters.element_types[nr]->number,"number",conf,0);
		SetParam(parameters.element_types[nr]->dof,"dof",conf,5);
		SetParam(parameters.element_types[nr]->saxes.vec,"r0",conf,Vec3(0.0));
		SetParam(parameters.element_types[nr]->mass,"mass",conf,0.0); // default is unknown mass
		SetParam(parameters.element_types[nr]->MOI,"MOI",conf,z); // per default atoms can be seen as spheres with "point" mass at center (which means no moment of inertia)
		SetParam(parameters.element_types[nr]->mu_pos.vec,"m0",conf,Vec3(0.0));
		if(fabs(parameters.element_types[nr]->mu_pos.V3Norm())>EPS) parameters.offctrmu=true;
		SetParam(parameters.element_types[nr]->dipole,"muz",conf,0.0);
		SetParam(parameters.element_types[nr]->initial_dipole.vec,"dipole",conf,Vec3(0.0,0.0,parameters.element_types[nr]->dipole));
		if(parameters.element_types[nr]->initial_dipole!=Vec3(0.0,0.0,parameters.element_types[nr]->dipole)){
			cout << "\t -> Using provided dipole vector.\n";
			parameters.element_types[nr]->dipole=parameters.element_types[nr]->initial_dipole.V3Norm();
		}
		SetParam(parameters.element_types[nr]->Vvdw,"vdw",conf);
		SetParam(parameters.element_types[nr]->rot_notrans,"rot_notrans",conf,false);
		SetParam(parameters.element_types[nr]->still,"stillmembers",conf,false);
		SetParam(parameters.element_types[nr]->calculate_order,"calculate_order",conf,false);
		
		SetParam(parameters.element_types[nr]->color.vec,"color",conf,Vec3(1.0)); // default color is white
		SetParam(parameters.element_types[nr]->transparency,"transparency",conf,0.0); // default is no transparency
		SetParam(parameters.element_types[nr]->label_elements,"label_elements",conf,true); // default is label elements when transparent
		SetParam(parameters.element_types[nr]->texture,"texture",conf,""); // no texture by default
		SetParam(parameters.element_types[nr]->texture_is_dipole,"texture_is_dipole",conf,true); // texture is treated to point in dipole direction
	} else{
		SetParam(parameters.element_types[nr]->number,"number",conf,parameters.element_types[parameters.element_types[nr]->archetype]->number);
		SetParam(parameters.element_types[nr]->dof,"dof",conf,parameters.element_types[parameters.element_types[nr]->archetype]->dof);
		SetParam(parameters.element_types[nr]->saxes.vec,"r0",conf,parameters.element_types[parameters.element_types[nr]->archetype]->saxes);
		SetParam(parameters.element_types[nr]->mass,"mass",conf,parameters.element_types[parameters.element_types[nr]->archetype]->mass);
		SetParam(parameters.element_types[nr]->MOI,"MOI",conf,parameters.element_types[parameters.element_types[nr]->archetype]->MOI);
		SetParam(parameters.element_types[nr]->mu_pos.vec,"m0",conf,parameters.element_types[parameters.element_types[nr]->archetype]->mu_pos);
		if(fabs(parameters.element_types[nr]->mu_pos.V3Norm())>EPS) parameters.offctrmu=true;
		SetParam(parameters.element_types[nr]->dipole,"muz",conf,parameters.element_types[parameters.element_types[nr]->archetype]->dipole);
		SetParam(parameters.element_types[nr]->initial_dipole.vec,"dipole",conf,parameters.element_types[parameters.element_types[nr]->archetype]->initial_dipole);
		if(parameters.element_types[nr]->initial_dipole!=parameters.element_types[parameters.element_types[nr]->archetype]->initial_dipole){
			cout << "\t -> Using provided dipole vector.\n";
			parameters.element_types[nr]->dipole=parameters.element_types[nr]->initial_dipole.V3Norm();
		}
		SetParam(parameters.element_types[nr]->Vvdw,"vdw",conf,parameters.element_types[parameters.element_types[nr]->archetype]->Vvdw);
		SetParam(parameters.element_types[nr]->rot_notrans,"rot_notrans",conf,parameters.element_types[parameters.element_types[nr]->archetype]->rot_notrans);
		SetParam(parameters.element_types[nr]->still,"stillmembers",conf,parameters.element_types[parameters.element_types[nr]->archetype]->still);
		SetParam(parameters.element_types[nr]->calculate_order,"calculate_order",conf,parameters.element_types[parameters.element_types[nr]->archetype]->calculate_order);
		SetParam(parameters.element_types[nr]->color.vec,"color",conf,parameters.element_types[parameters.element_types[nr]->archetype]->color);
		SetParam(parameters.element_types[nr]->transparency,"transparency",conf,parameters.element_types[parameters.element_types[nr]->archetype]->transparency);
		SetParam(parameters.element_types[nr]->texture,"texture",conf,parameters.element_types[parameters.element_types[nr]->archetype]->texture);
		SetParam(parameters.element_types[nr]->texture_is_dipole,"texture_is_dipole",conf,parameters.element_types[parameters.element_types[nr]->archetype]->texture_is_dipole);
	}
	
	// set some values from just read parameters
	parameters.element_types[nr]->thistype=nr;
	parameters.element_types[nr]->hasmu=(fabs(parameters.element_types[nr]->dipole)>EPS);
	
	// create charges
	parameters.element_types[nr]->nr_charges=0;
	parameters.element_types[nr]->q_pos=NULL;
	parameters.element_types[nr]->q=NULL;
	if(load_charges){
#if DEBUG_LEVEL>0
		cout << "-> Creating charges for element.\n";
#endif
		double* charges = get_flex_tupel("charges",conf);
		if(((parameters.element_types[nr]->archetype>=0) && (parameters.element_types[parameters.element_types[nr]->archetype]->nr_charges>0)) && (charges==NULL)){
			parameters.element_types[nr]->nr_charges=parameters.element_types[parameters.element_types[nr]->archetype]->nr_charges;
			parameters.element_types[nr]->q_pos=parameters.element_types[parameters.element_types[nr]->archetype]->q_pos;
			parameters.element_types[nr]->q=parameters.element_types[parameters.element_types[nr]->archetype]->q;
		}
		if(charges!=NULL){
			string item;
			int i=0;
			while((int)charges[i]!=-1){ // go through all charges (no matter how they are divided)
				for(int j=0; j<charges[i]; j++){
					// assign charge if not zero
					if(fabs(charges[i+j+1])>EPS){
						parameters.element_types[nr]->nr_charges++;
						parameters.element_types[nr]->q=(double*)realloc(parameters.element_types[nr]->q,sizeof(double)*parameters.element_types[nr]->nr_charges);
						if(parameters.element_types[nr]->q==NULL){
							cout << "Could not get enough memory to hold " << parameters.element_types[nr]->nr_charges << " partial charges.\n";
							exit(3);
						}
						parameters.element_types[nr]->q[parameters.element_types[nr]->nr_charges-1]=e_in_esu*charges[i+j+1]; // partial charges in e are now in esu
						
						// assign charge location
						item="charge_location."+to_string(parameters.element_types[nr]->nr_charges);
						parameters.element_types[nr]->q_pos=(Vec3*)realloc(parameters.element_types[nr]->q_pos,sizeof(Vec3)*parameters.element_types[nr]->nr_charges);
						if(parameters.element_types[nr]->q_pos==NULL){
							cout << "Could not get enough memory to hold " << parameters.element_types[nr]->nr_charges << " charge locations.\n";
							exit(3);
						}
						SetParam(parameters.element_types[nr]->q_pos[parameters.element_types[nr]->nr_charges-1].vec,item.c_str(),conf,Vec3(0.0));
					}
				}
				i+=(int)charges[i]+1;
			}
#if DEBUG_LEVEL>0
			cout << "<- Created " << parameters.element_types[nr]->nr_charges;
			if(parameters.element_types[nr]->nr_charges>1) cout << " charges.\n"; else cout << " charge.\n";
#endif
			delete[] charges;
		} else{
#if DEBUG_LEVEL>0
			cout << "<- No charges specified.\n";
#endif
		}
	}
	parameters.element_types[nr]->nr_Vvdw_coefficients=0;
	parameters.element_types[nr]->Vvdw_coefficients=NULL;
	parameters.element_types[nr]->IA_coefficients=NULL;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
} // element types are set up now - group creation can follow ...

char* MC_Config::Mol2Convert(char* content, char* conf, string groupname)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	string conversion="[Mol2 conversion]\n";
	// first, get element definition section
	char* pos=find_string(content,"@<TRIPOS>ATOM")+strlen("@<TRIPOS>ATOM")+1; // pos is on first line
	char* end=find_string(pos,"@");
	if(end==NULL) end=content+strlen(content);
	// now go through it line by line
	char* temp;
	string elements="elements = ";
	string charges="import_charges = {";
	string item;
	Vec3 v;
	bool getposition;
	
	bool first=true;
	int count=0;
	while((pos<end) && (pos!=NULL)){
		temp=strchr(pos,'\n');
		if(temp==NULL) temp=end; // last line
		if(temp>pos){ // interpret current line
			count++;
			item="position."+to_string(count);
			conversion+="\t"+item+" = (";
			SetParam(v.vec,item.c_str(),conf,Vec3(0.0)); // read user-defined position
			getposition=false;
			if(v*v<1E-10) getposition=true;
			while((*pos==' ') || (*pos=='\t')) pos++; // whitespace removal
			string atomid="";
			while((*pos!=' ') && (*pos!='\t')){ // first column: atom id
				atomid+=*pos; // copy atom id
				pos++;
			}
			while((*pos==' ') || (*pos=='\t')) pos++; // more whitespace removal
			while((*pos!=' ') && (*pos!='\t')) pos++; // second column: atom name (not used b/c mostly arbitrary)
			while((*pos==' ') || (*pos=='\t')) pos++; // whitespace removal to the next entry
			while((*pos!=' ') && (*pos!='\t')){ // third column: position x
				if(getposition) conversion+=*pos; // copy x
				pos++;
			}
			conversion+=",";
			while((*pos==' ') || (*pos=='\t')) pos++; // whitespace removal to the next entry
			while((*pos!=' ') && (*pos!='\t')){ // fourth column: position y
				if(getposition) conversion+=*pos; // copy y
				pos++;
			}
			conversion+=",";
			while((*pos==' ') || (*pos=='\t')) pos++; // whitespace removal to the next entry
			while((*pos!=' ') && (*pos!='\t')){ // fifth column: position z
				if(getposition) conversion+=*pos; // copy z
				pos++;
			}
			if(!getposition) conversion+=to_string(v.vec[0])+","+to_string(v.vec[1])+","+to_string(v.vec[2]);
			conversion+=")\n";
			
			while((*pos==' ') || (*pos=='\t')) pos++; // whitespace removal to the next entry
			if(!first) elements+=","; else first=false; // first element does not start with a comma
			string element_name=groupname+" ";
			string type="";
			while((*pos!=' ') && (*pos!='\t')){ // sixth column: atom type
				type+=*pos; // copy element name
				pos++;
			}
			element_name+=type+" "+atomid;
			// does such an element currently exist (only if user defined it)
			bool element_exists=false;
			for(unsigned int i=0; i<parameters.num_element_types; i++){
				if(compare_strings(parameters.element_types[i]->name.c_str(),element_name.c_str())){
					element_exists=true;
					break;
				}
			}
			if(!element_exists){
				parameters.num_element_types++;
				parameters.element_types=(Element_Type**)realloc(parameters.element_types,parameters.num_element_types*sizeof(Element_Type*));
				if(!parameters.element_types){
					cout << "Not enough memory to create an additional element type.\n";
					exit(1);
				}
				parameters.element_types[parameters.num_element_types-1]=new Element_Type;
				parameters.element_types[parameters.num_element_types-1]->name=element_name;
				string newelement="archetype = "+type+"\n";
				char* element=new char[newelement.length()+1];
				strcpy(element,newelement.c_str());
				LoadElementTypes(element,parameters.num_element_types-1,false); // want no charges because we'll override (and create) them later
				delete[] element;
			}
			elements+=element_name;
			while((*pos==' ') || (*pos=='\t')) pos++; // whitespace removal to the next entry
			while((*pos!=' ') && (*pos!='\t')) pos++; // seventh column: not used
			while((*pos==' ') || (*pos=='\t')) pos++; // whitespace removal to the next entry
			while((*pos!=' ') && (*pos!='\t')) pos++; // eight column: not used
			
			if(count>1) charges+=",";
			while((*pos==' ') || (*pos=='\t')) pos++; // whitespace removal to the next entry
			while(((*pos!=' ') && (*pos!='\t')) && ((*pos!='\r') && (*pos!='\n'))){ // nineth column: charge
				charges+=*pos;
				pos++;
			}
		}
		pos=temp+1;
	}
	charges+="}\n";
	conversion+="\t"+charges;
	char* value = GetItemValue(conf,"elements",true);
	temp=value;
	first=true;
	int useradd=0;
	while(temp!=NULL){
		if(first){
			useradd++;
			elements+=","; // user defined elements are added with a comma (after mol2 elements)
			first=false;
		}
		elements+=*temp;
		if(*temp==',') useradd++; // count added elements
		temp++;
		if(*temp=='\0') temp=NULL;
	}
	if(value) delete[] value;
	for(int i=0; i<useradd; i++){
		item="position."+to_string(count+i+1);
		SetParam(v.vec,item.c_str(),conf,Vec3(0.0));
		// only write to configuration if specified
		if(v*v>1E-10) conversion+="\t"+item+" = ("+to_string(v.vec[0])+","+to_string(v.vec[1])+","+to_string(v.vec[2])+")\n";
	}
	conversion+="\t"+elements+"\n";
	
	pos=find_string(content,"@<TRIPOS>BOND")+strlen("@<TRIPOS>ATOM")+1;
	end=find_string(pos,"@");
	if(end==NULL) end=content+strlen(content);
	
	string* connections = new string[count];
	string* bond_orders = new string[count];
	int current;
	while((pos<end) && (pos!=NULL)){
		temp=strchr(pos,'\n');
		if(temp==NULL) temp=end; // last line
		if(temp>pos){ // interpret current line
			while((*pos==' ') || (*pos=='\t')) pos++; // whitespace removal
			while((*pos!=' ') && (*pos!='\t')) pos++; // don't want first column (running number)
			while((*pos==' ') || (*pos=='\t')) pos++; // more whitespace removal
			item="";
			while((*pos!=' ') && (*pos!='\t')){ // second column: element type
				item+=*pos; // copy item
				pos++;
			}
			if(from_string(current,item)){
				cout << "Problem with connectivity: Can not read element number.\n";
				exit(3);
			}
			if(current>count){
				cout << "Problem with connectivity: Element number larger than existing elements.\n";
				exit(3);
			}
			if(connections[current-1]!=""){
				connections[current-1]+=",";
				bond_orders[current-1]+=",";
			}
			while((*pos==' ') || (*pos=='\t')) pos++; // more whitespace removal
			item="";
			while((*pos!=' ') && (*pos!='\t')){ // third column: element type
				item+=*pos; // copy item
				pos++;
			}
			connections[current-1]+=item;
			while((*pos==' ') || (*pos=='\t')) pos++; // more whitespace removal
			item="";
			while((*pos!=' ') && (*pos!='\t') && (*pos!='\n')){ // fourth column: bond order
				item+=*pos; // copy item
				pos++;
			}
			bond_orders[current-1]+=item;
		}
		pos=temp+1;
	}
	conversion+="\tconnectivity = {";
	string bo="\tbond_order = {";
	first=true;
	for(int i=0; i<count; i++){
		if(!first){
			conversion+="|";
			bo+="|";
		} else first=false;
		conversion+=connections[i];
		bo+=bond_orders[i];
	}
	delete[] connections; // done with it, may as well clean up ...
	delete[] bond_orders;
	
	value = GetItemValue(conf,"connectivity",true);
	temp = value;
	first=true;
	while(temp!=NULL){
		if(first){
			temp++; // get rid of opening bracket
			conversion+="|"; // user defined elements are added with a comma
			first=false;
		}
		conversion+=*temp;
		temp++;
		if(*temp=='\0') temp=NULL;
	}
	delete[] value;
	value = GetItemValue(conf,"bond_order",true);
	temp = value;
	first=true;
	while(temp!=NULL){
		if(first){
			temp++; // get rid of opening bracket
			bo+="|"; // user defined elements are added with a comma
			first=false;
		}
		bo+=*temp;
		temp++;
		if(*temp=='\0') temp=NULL;
	}
	delete[] value;
	count+=useradd; // just for the record
	if(first){
		conversion+="}\n";
		bo+="}\n";
	} else{
		conversion+="\n";
		bo+="\n";
	}
	conversion+=bo;
	
	char* result = new char[conversion.length()+1];
	strcpy(result,conversion.c_str());
	
	delete[] content;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return result;
}

char* MC_Config::PDBConvert(char* content, char* conf, string groupname)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	cout << "-> Patience young padawan. Feature you are looking for not yet written it is.\n";
	exit(42);
	
	char* result=NULL;
	
	delete[] content;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return result;
}

inline char* ReadFile(string filename)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	unsigned int size;
	ifstream file(filename.c_str(),ifstream::in);
	if (file.fail()==true){
		cout << "Could not open file \"" << filename << "\".\n";
		exit(1);
	}
	// Get file size
	file.seekg(0,ifstream::end);
	size=file.tellg();
	file.seekg(0);
	// Sanity checks
	if (size==0){
		cout << "File \"" << filename << "\" is empty.\n";
		exit(1);
	}
	if (size>MAXCONFSIZE){ // config files should not be bigger than MAXCONFSIZE/(1024*1024) MB
		cout << "File \"" << filename << "\" is too big (>" << MAXCONFSIZE/(1024*1024) << " MB).\n";
		exit(1);
	}
	// Read content
	char* content = new char[size+1];
	file.read(content,size);
	content[size]='\0';
	file.close();
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return content;
}

void MC_Config::GetSpecialDistances(double* special, int nr)
{
	unsigned int count=0;
	unsigned int total=0;
	while((int)special[count]>0){
		if(total>=SPECIAL_DISTANCES){
			cout << "Too many bond distance factors, only up to " << SPECIAL_DISTANCES << " are possible.\n";
			exit(2);
		}
		if((int)special[count]!=2){
			cout << "Need exactly two parameters for each entry in <bond_distance_factors> (e.g. {3,0.5|4,0.75}).\n";
			exit(2);
		}
		double d=special[3*total+1];
		int i=(int)qround(d);
		if(fabs(i-d)>EPS) cout << "WARNING: First parameter in each key of <bond_distance_factors> needs to be an integer (rounding to get one).\n";
		if(i<1){
			cout << "First parameter in each key of <bond_distance_factors> needs to be an integer greater or equal to 1.\n";
			exit(2);
		}
		parameters.groups[nr]->Type->range_factors[total].distance=(unsigned int)i;
		parameters.groups[nr]->Type->range_factors[total].factor=(float)special[3*total+2];
		total++;
		count+=3;
	}
	delete[] special;
}

bool MC_Config::LoadGroup(char* conf, int nr, bool rerun)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	Element_Group* group=parameters.groups[nr];
	Element_Group* lodgroup=NULL;
	
	unsigned int i,j,k,count;
	bool finished=true;
	bool fileimport=false;
	char* content=NULL;
	if(!rerun){
		SetParam(group->number,"number",conf,0);
		group->levelofdetail=0;
		group->Type->LOD=NULL; // reasonable default
		SetParam(group->Type->dof,"dof",conf,6);
		SetParam(group->Type->allow_bond_bend,"allow_bond_bend",conf,false);
		SetParam(group->Type->group_dipole,"group_dipole",conf,false);
		SetParam(group->Type->rand_independent,"rand_independent",conf,false);
		SetParam(group->Type->rot_notrans,"rot_notrans",conf,false);
		SetParam(group->Type->still,"stillmembers",conf,false);
		SetParam(group->Type->add_start_LOD,"add_start_LOD",conf,true);
		SetParam(group->Type->visual_PBCs,"visual_PBCs",conf,true);
		SetParam(group->Type->transparency,"transparency",conf,-1.0); // default is no special group transparency
		SetParam(group->Type->show_just_dipoles,"show_just_dipoles",conf,false);
		SetParam(group->Type->density,"density",conf,0.0); // default is unknown density
		SetParam(group->Type->label_elements,"label_elements",conf,true); // default is label elements when transparent
		SetParam(group->Type->calculate_order,"calculate_order",conf,false);
		SetParam(group->Type->adjust_overlap,"adjust_overlap",conf,-1.0);
		double* special_distances = get_flex_tupel("bond_distance_factors",conf);
		for(i=0; i<SPECIAL_DISTANCES; i++){ // set reasonable default
			group->Type->range_factors[i].distance=0;
			group->Type->range_factors[i].factor=0.0;
		}
		if(special_distances) GetSpecialDistances(special_distances,nr);
		string filename;
		SetParam(filename,"structure_mol2",conf,""); // get filename from configuration file, mol2 format first (default is no filename)
		if(filename!=""){
			fileimport=true;
#if DEBUG_LEVEL>0
			cout << "-> Reading group elements, positions, and connectivity from Mol2 file <" << filename << ">.\n";
#endif
			content = Mol2Convert(ReadFile(parameters.configdir+filename),conf,group->Type->name);
#if DEBUG_LEVEL>2
			cout << content << "<- EOC.\n";
#endif
		} else{
			SetParam(filename,"structure_pdb",conf,""); // get filename from configuration file, mol2 format first (default is no filename)
			if(filename!=""){
				fileimport=true;
#if DEBUG_LEVEL>0
				cout << "-> Reading group elements, positions, and connectivity from PDB file <" << filename << ">.\n";
#endif
				content = PDBConvert(ReadFile(parameters.configdir+filename),conf,group->Type->name);
#if DEBUG_LEVEL>2
				cout << content << "<- EOC.\n";
#endif
			}
		}
	}
	string element_list, item;
	// Create elements belonging to this particular group
	if(fileimport) SetParam(element_list,"elements",content); else SetParam(element_list,"elements",conf);
	// Split list and assign group elements
	element_list+=','; // cheat a bit and keep code simple ;-)
	if(!rerun){
		group->nr_elements=0; // initially, there are no elements ...
		group->elements=NULL; // needed for realloc to work
	}
	unsigned int rerun_nr_elements=group->nr_elements;
	unsigned int* rerun_elements=NULL;
	unsigned int delete_nr_elements=0;
	unsigned int* delete_elements=NULL; // to keep track of group elements for deletion (dynamic array)
	unsigned int* group_startnr=NULL; 
	unsigned int* lodgroup_startnr=NULL; 
	bool copy_lod=false;
	Element_Group** addgroups=NULL;
	int group_elements_start=nr_group_elements; // how many associated group elements had been created before this group
	
	bool element=false;
	bool found, found_group;
	unsigned int added_groups=0;
	string element_name = "";
	
	int* group_connections=NULL;
	if(rerun){
		if(group->nr_potentials>0){
			free(group->potentials);
			group->nr_potentials=0;
			group->potentials=NULL;
		}
		rerun_elements=new unsigned int[rerun_nr_elements];
		group_startnr=new unsigned int[rerun_nr_elements];
		lodgroup_startnr=new unsigned int[rerun_nr_elements];
		memcpy(rerun_elements,group->elements,sizeof(unsigned int)*rerun_nr_elements);
		memset(group_startnr,0,sizeof(unsigned int)*rerun_nr_elements);
		memset(lodgroup_startnr,0,sizeof(unsigned int)*rerun_nr_elements);
		addgroups=new Element_Group*[rerun_nr_elements];
		memset(addgroups,0,sizeof(Element_Group*)*rerun_nr_elements);
		// read in super group connection data
		group_connections=new int[rerun_nr_elements*rerun_nr_elements];
		memset(group_connections,0,rerun_nr_elements*rerun_nr_elements*sizeof(int)); // puts 0 in each field
		for(i=0; i<rerun_nr_elements; i++){
			for(j=0; j<rerun_nr_elements; j++){
				if(i!=j){
					item="group_connection."+int2str(i+1)+"-"+int2str(j+1);
					double* connection=get_MxN_tupel(item.c_str(),conf,2,1,NULL,false);
					if(connection){
						if(abs(group_connections[i*rerun_nr_elements+j])+abs(group_connections[j*rerun_nr_elements+i])>0){
							cout << "ERROR: There already is a connection between entities " << i+1 << " and " << j+1 << ". Only one definition allowed between any two entities.\n";
							exit(4);
						}
						if(((fabs(connection[0]-(int)connection[0])>EPS) || (fabs(connection[1]-(int)connection[1])>EPS)) || (((int)connection[0]==0) || ((int)connection[1]==0))){ // sanity check input to be integers
							cout << "ERROR in \"" << item << "\": Please use only non-zero integers numbered with respect to the corresponding group.\n-> Positive numbers denote individual elements, negative numbers denote connection sites.\n";
							exit(5);
						}
						group_connections[i*rerun_nr_elements+j]=(int)connection[0];
						group_connections[j*rerun_nr_elements+i]=(int)connection[1];
						delete[] connection;
					}
				}
			}
		}
	}
	count=0;
	for(i=0; i<element_list.length(); i++){
		if((element_list[i]==',') || (element_list[i]==';')){ // element_name finished
			// trailing white space removal
			j=element_name.length()-1;
			while((element_name[j]==' ') || (element_name[j]=='\t')){ j--; }
			j++;
			element_name.erase(j,element_name.length()-j); // element_name now contains name of element
			count++;
			if(!rerun){
				// Create element for group and link to it in the group
#if DEBUG_LEVEL>0
				cout << "Creating " << group->Type->name << " entity: " << element_name << "\n";
#endif
				nr_group_elements++; // Adding one more element, increasing running total
				parameters.group_elements = (Element*)realloc(parameters.group_elements,nr_group_elements*sizeof(Element));
				if(parameters.group_elements==NULL){ // Danger, danger
					cout << "ERROR: Could not find a memory block large enough to hold a total of " << nr_group_elements << " group associated elements.\n";
					exit(1);
				}
				group->nr_elements++;
				group->elements = (unsigned int*)realloc(group->elements,group->nr_elements*sizeof(unsigned int));
				if(group->elements==NULL){ // More danger ...
					cout << "ERROR: Not enough memory to hold elements numbers of group " << group->Type->name << "\n";
					exit(1);
				}
				group->elements[group->nr_elements-1]=nr_group_elements-1;
			}
			// find element type and link just created element to it
			found=false;
			for(j=0; j<parameters.num_element_types; j++){
				if(compare_strings(parameters.element_types[j]->name.c_str(),element_name.c_str())){
					found=true;
					break;
				}
			}
			found_group=false;
			if(!found){ // Two possibilities here, the name may be a group or it really does not exist
				for(j=0; j<parameters.num_groups; j++){
					if(compare_strings(parameters.groups[j]->Type->name.c_str(),element_name.c_str())){
						found_group=true;
						break;
					}
				}
				if(!found_group){
					if(rerun){ // if this is a rerun then the element/group does not exist
						cout << "ERROR: Entity '" << element_name << "' does not exist.\n";
						exit(1);
					} else{ // otherwise, come back to it later
						cout << "-> Could not (yet) find entity '" << element_name << "'.\n";
					}
				}
				if(!rerun) finished=false;
			}
			if(found_group && rerun){ // now things get interesting: add a group
				Element_Group* add_group=parameters.groups[j];
				Element_Group* base_add_group=add_group;
				if(add_group->Type->LOD){
					// if group entity has LOD version and is the base group: incorporate LOD start level
					if(add_group->levelofdetail==0 && group->Type->add_start_LOD) add_group=add_group->Type->LOD->groups[add_group->Type->LOD->start_level];
					base_add_group=add_group->Type->LOD->groups[0];
					if(add_group->levelofdetail>0){
						if(!group->Type->LOD){
							cout << "-> Adding Level-of-detail feature\n";
							parameters.num_levelofdetail++;
							parameters.lods=(Level_of_Detail**)realloc(parameters.lods,parameters.num_levelofdetail*sizeof(Level_of_Detail*));
							if(!parameters.lods){
								cout << "ERROR: Out of memory adding an LOD model.\n";
								exit(11);
							}
							parameters.lods[parameters.num_levelofdetail-1]=new Level_of_Detail;
							group->Type->LOD=parameters.lods[parameters.num_levelofdetail-1];
							group->Type->LOD->groups=new Element_Group*[2];
							group->Type->LOD->groups[0]=group;
							group->Type->LOD->levels=1;
							group->Type->LOD->start_level=1;
							group->Type->LOD->end_level=1;
							SetParam(group->Type->LOD->symmetry_center,"symmetry_axes",conf,false);
							SetParam(group->Type->LOD->match_volume,"match_volume",conf,false);
							SetParam(group->Type->LOD->symmetry_axes,"symmetry_axes",conf,true);
							SetParam(group->Type->LOD->match_epsilon,"match_epsilon",conf,false);
							// set some flags (all start out zero, but added group's value are combined using OR)
							group->Type->LOD->keep_original_potentials=new bool;
							*group->Type->LOD->keep_original_potentials=true;
							group->Type->LOD->use_epsilon_texture=new bool;
							SetParam(*group->Type->LOD->use_epsilon_texture,"use_epsilon_texture",conf,false);
							group->Type->LOD->visual_original=new bool;
							*group->Type->LOD->visual_original=false;
							group->Type->LOD->show_internal_bonds=false;
							group->Type->LOD->internal_charge_colors=false;
							// visual aspects are averaged while building
							group->Type->LOD->transparency=new double;
							*group->Type->LOD->transparency=0.0;
							group->Type->LOD->inside_transparency=new double;
							*group->Type->LOD->inside_transparency=0.0;
							group->Type->LOD->scale_original_vdw=added_groups;
							// reasonable defaults for other stuff
							group->Type->LOD->reduce_electrostatics=new double*; // electrostatics are already reduced based on chosen LOD model
							group->Type->LOD->reduce_electrostatics[0]=NULL;
							group->Type->LOD->distances=NULL; // not used currently (TODO)
							group->Type->LOD->volumes = new double*[parameters.lods[nr]->levels];
							group->Type->LOD->epsilons = new double*[parameters.lods[nr]->levels];
							group->Type->LOD->placement_deltas=(Vec3*)malloc(group->Type->LOD->levels*sizeof(Vec3));
							// no components at the moment ...
							group->Type->LOD->order_dipole_component=0;
							group->Type->LOD->nr_components=0;
							group->Type->LOD->nr_elements=NULL;
							group->Type->LOD->component_elements=NULL;
							group->Type->LOD->element_in_component=NULL;
							group->Type->LOD->component_show_dipole=NULL;
							group->Type->LOD->component_color=NULL;
							group->Type->LOD->component_transparency=NULL;
							group->Type->LOD->component_texture=NULL;
							// the remaining LOD parameters are set by "AddGroup2Group" but need to be initialized
							group->Type->LOD->nr_ellipsoids=new unsigned int[group->Type->LOD->levels]; // "group->Type->LOD->levels" is a fancy way of writing "1" (see above), but I may want to extend that later, so ...
							group->Type->LOD->nr_ellipsoids[0]=0;
							group->Type->LOD->ellipsoid_counts=new unsigned int*[group->Type->LOD->levels];
							group->Type->LOD->ellipsoid_counts[0]=NULL;
							group->Type->LOD->element_groups=new int*[group->Type->LOD->levels];
							group->Type->LOD->element_groups[0]=NULL;
							group->Type->LOD->element_in_ellipsoid=new int*[group->Type->LOD->levels];
							group->Type->LOD->element_in_ellipsoid[0]=NULL;
							// create new LOD group
							parameters.num_groups++;
							parameters.groups = (Element_Group**)realloc(parameters.groups,parameters.num_groups*sizeof(Element_Group*));
							if(parameters.groups==NULL){ // Danger, danger
								cout << "Could not find a memory block large enough to extend groups for level of detail.\n";
								exit(1);
							}
							lodgroup=new Element_Group;
							lodgroup->type=parameters.num_groups-1;
							lodgroup->Type=new Element_Group_Type;
							lodgroup->Type->name=group->Type->name+".LOD1";
							lodgroup->Type->dof=group->Type->dof;
							lodgroup->Type->rand_elements=0;
							lodgroup->Type->allow_bond_bend=group->Type->allow_bond_bend;
							lodgroup->Type->group_dipole=group->Type->group_dipole;
							lodgroup->Type->show_just_dipoles=group->Type->show_just_dipoles;
							lodgroup->Type->rand_independent=group->Type->rand_independent;
							lodgroup->Type->still=group->Type->still;
							lodgroup->Type->rot_notrans=group->Type->rot_notrans;
							lodgroup->Type->add_start_LOD=group->Type->add_start_LOD;
							lodgroup->Type->visual_PBCs=group->Type->visual_PBCs;
							lodgroup->Type->transparency=group->Type->LOD->transparency[0];
							lodgroup->Type->label_elements=group->Type->label_elements;
							lodgroup->Type->mass=0.0;
							lodgroup->Type->density=group->Type->density;
							lodgroup->Type->packing_density=group->Type->packing_density;
							lodgroup->Type->bond_range=group->Type->bond_range;
							lodgroup->Type->calculate_order=group->Type->calculate_order;
							SetParam(lodgroup->Type->adjust_overlap,"adjust_overlap",conf,0.5);
							for(j=0; j<SPECIAL_DISTANCES; j++){
								lodgroup->Type->range_factors[j].distance=group->Type->range_factors[j].distance;
								lodgroup->Type->range_factors[j].factor=group->Type->range_factors[j].factor;
							}
							lodgroup->Type->elements_in_ring=NULL;
							lodgroup->Type->nr_rings=0;
							lodgroup->Type->rings=NULL;
							lodgroup->Type->nr_connection_sites=0;
							lodgroup->Type->connection_sites=NULL;
							lodgroup->Type->range_field=NULL;
							lodgroup->number=group->number;
							group->number=0;
							group->Type->LOD->groups[1]=lodgroup;
							lodgroup->levelofdetail=1;
							lodgroup->Type->LOD=group->Type->LOD;
							parameters.groups[lodgroup->type]=lodgroup;
							lodgroup->nr_elements=0;
							lodgroup->elements=NULL;
							lodgroup->nr_potentials=0;
							lodgroup->potentials=NULL;
							copy_lod=true;
							if(count>1){ // Need to take care of (aka add) non-LOD'ed stuff added before
								AddGroup2Group(lodgroup,group,0);
							}
						}
					}
				}
#if DEBUG_LEVEL>0
				cout << "-> Adding group <" << add_group->Type->name << ">\n";
#endif
#if DEBUG_LEVEL>2
				for(k=0; k<base_add_group->Type->nr_connection_sites; k++) cout << base_add_group->Type->connection_sites[k] << "\n";
#endif
				added_groups++;
				group->Type->bond_range+=add_group->Type->bond_range;
				// store placeholder element number to delete them later
				delete_nr_elements++;
				delete_elements=(unsigned int*)realloc(delete_elements,delete_nr_elements*sizeof(unsigned int));
				if(!delete_elements){
					cout << "ERROR: O memory, where art thou?\n";
					exit(4);
				}
				delete_elements[delete_nr_elements-1]=count-1;
				// add group elements of base group
				addgroups[count-1]=add_group; // only one pointer needed here b/c LOD groups contain link back to their base group anyway
				group_startnr[count-1]=group->nr_elements;
				AddGroup2Group(group,base_add_group,0);
				// add LOD group elements
				if(group->Type->LOD){
					lodgroup_startnr[count-1]=lodgroup->nr_elements;
					AddGroup2Group(lodgroup,add_group,group_startnr[count-1]);
					if(add_group->levelofdetail>0){
						*group->Type->LOD->transparency+=add_group->Type->LOD->transparency[add_group->levelofdetail-1];
						group->Type->LOD->scale_original_vdw+=add_group->Type->LOD->scale_original_vdw;
						*group->Type->LOD->visual_original|=add_group->Type->LOD->visual_original[add_group->levelofdetail-1];
						group->Type->LOD->show_internal_bonds|=add_group->Type->LOD->show_internal_bonds;
						group->Type->LOD->internal_charge_colors|=add_group->Type->LOD->internal_charge_colors;
					}
				} else lodgroup_startnr[count-1]=group_startnr[count-1]; // when LOD feature is switched on (b/c LOD group is added at some point) previous elements are copied 1:1
				// if group is connected to something that has already been added connect and rotate ...
				for(k=0; k<parameters.group_elements[rerun_elements[count-1]].nr_interactions; k++){
					unsigned int partnernr=parameters.group_elements[rerun_elements[count-1]].interactions[k].partner;
					unsigned int nr_potentials=parameters.group_elements[rerun_elements[count-1]].interactions[k].nr_potentials;
					int curr_connection=group_connections[(count-1)*rerun_nr_elements+partnernr];
					int partner_connection=group_connections[rerun_nr_elements*partnernr+(count-1)];
#if DEBUG_LEVEL>2
					cout << "-> " << partnernr+1 << " (" << curr_connection << "<->" << partner_connection << "), (potentials: " << nr_potentials << ")\n";
#endif
					if(partnernr<count){ // ... this works b/c both connections are stored on read-in
						Element_Group* partnergroup=addgroups[partnernr];
						Element_Group* basepartner=partnergroup;
						if(partnergroup->Type->LOD) basepartner=partnergroup->Type->LOD->groups[0];
#if DEBUG_LEVEL>0
						cout << "\t-> Setting up connection with ";
						if(partnergroup) cout << "<" << partnergroup->Type->name << ">\n"; else cout << "element " << partnernr+1 << "\n";
#endif
						// sanity-check connection points
						if((curr_connection>0) && (curr_connection>(int)base_add_group->nr_elements)){
							cout << "ERROR: <" << base_add_group->Type->name << "> has only " << base_add_group->nr_elements << " elements.\n";
							exit(3);
						}
						if((curr_connection<0) && (-curr_connection>(int)base_add_group->Type->nr_connection_sites)){
							cout << "ERROR: <" << base_add_group->Type->name << "> has only " << base_add_group->Type->nr_connection_sites << " connection sites.\n";
							exit(3);
						}
						if((partner_connection>0) && (partner_connection>(int)basepartner->nr_elements)){
							cout << "ERROR: <" << basepartner->Type->name << "> has only " << basepartner->nr_elements << " elements.\n";
							exit(3);
						}
						if((partner_connection<0) && (-partner_connection>(int)basepartner->Type->nr_connection_sites)){
							cout << "ERROR: <" << basepartner->Type->name << "> has only " << basepartner->Type->nr_connection_sites << " connection sites.\n";
							exit(3);
						}
						if((curr_connection==0) || (partner_connection==0)){
							cout << "ERROR: No connection site specified for connecting groups " << count << " and " << partnernr+1 << ".\nPlease specify using \"group_connection." << count << "-" << partnernr+1 << "\"" << " or \"group_connection." << partnernr+1 <<"-" <<  count << "\".\n";
							exit(3);
						}
						// take care of fully-atomistic model
						unsigned int connection_A, connection_B; // A is element on previous group, B current
						if(curr_connection>0) connection_B=curr_connection-1; else connection_B=base_add_group->Type->connection_sites[(unsigned int)(-curr_connection-1)]-1;
						if(partner_connection>0) connection_A=partner_connection-1; else connection_A=basepartner->Type->connection_sites[(unsigned int)(-partner_connection-1)]-1;
						connection_B+=group_startnr[count-1];
						connection_A+=group_startnr[partnernr];
						Element* A=&parameters.group_elements[group->elements[connection_A]];
						Element* B=&parameters.group_elements[group->elements[connection_B]];
						cout << "\t\t-> Connection sites: " << A->MyType->name << " <-> " << B->MyType->name << "\n";
						// determine how many bonds each partner has and calculate average bond vector
						bool A_onebond=(A->nr_interactions==1);
						bool B_onebond=(B->nr_interactions==1);
						if(!A_onebond && !B_onebond){ // enforce that both A and B are same type (otherwise this combination does not work)
							if(A->mytype!=B->mytype){
								cout << "ERROR: Can only connect two multibond sites if their connection sites have the same type.\n";
								exit(7);
							}
						}
						if((A->nr_interactions==0) || (B->nr_interactions==0)){ // problem
							cout << "ERROR: Connection sites are not connected to group.\n";
							exit(3);
						}
						Vec3 Avec(0.0);
						Vec3 CAloc=A->center;
						double CAfrac=0.0;
						unsigned int CAnr=connection_A;
						unsigned int CAlink=0;
						Element* CA=A;
						for(unsigned int l=0; l<A->nr_interactions; l++){
							Element* C=&parameters.group_elements[group->elements[A->interactions[l].partner]];
							Avec+=A->center-C->center+A->rot*(*A->interactions[l].initial_location)-C->rot*(*C->interactions[A->interactions[l].back_link].initial_location);
							if(A_onebond){
								CA=C;
								CAnr=A->interactions[l].partner;
								CAlink=A->interactions[l].back_link;
								CAloc=CA->center+CA->rot*(*CA->interactions[l].initial_location);
								CAfrac=C->MyType->saxes.vec[0]*C->MyType->saxes.vec[1]*C->MyType->saxes.vec[2]/(A->MyType->saxes.vec[0]*A->MyType->saxes.vec[1]*A->MyType->saxes.vec[2]+C->MyType->saxes.vec[0]*C->MyType->saxes.vec[1]*C->MyType->saxes.vec[2]);
							}
						}
						Avec/=A->nr_interactions;
						double lCA=CAfrac*Avec.V3Norm();
						if(CAfrac<=EPS){ // for multiple bonds on the connection site use half the semiaxis as the bond length contribution
							lCA=average(A->MyType->saxes.vec,3)/2;
						}
#if DEBUG_LEVEL>2
						cout << Avec.V3Str() << " (" << Avec.V3Norm() << ", " << lCA << ")\n";
#endif
						Avec/=-Avec.V3Norm(); // unit vector pointing away from connection on A (aka towards connection from B site)
						Vec3 Bvec(0.0);
						Vec3 CBloc=B->center;
						double CBfrac=0.0;
						unsigned int CBnr=connection_B;
						unsigned int CBlink=0;
						Element* CB=B;
						for(unsigned int l=0; l<B->nr_interactions; l++){
							Element* C=&parameters.group_elements[group->elements[B->interactions[l].partner]];
							Bvec+=B->center-C->center+B->rot*(*B->interactions[l].initial_location)-C->rot*(*C->interactions[B->interactions[l].back_link].initial_location);
							if(B_onebond){
								CB=C;
								CBnr=B->interactions[l].partner;
								CBlink=B->interactions[l].back_link;
								CBloc=CB->center+CB->rot*(*CB->interactions[l].initial_location);
								CBfrac=C->MyType->saxes.vec[0]*C->MyType->saxes.vec[1]*C->MyType->saxes.vec[2]/(B->MyType->saxes.vec[0]*B->MyType->saxes.vec[1]*B->MyType->saxes.vec[2]+C->MyType->saxes.vec[0]*C->MyType->saxes.vec[1]*C->MyType->saxes.vec[2]);
							}
						}
						Bvec/=B->nr_interactions;
						double lCB=CBfrac*Bvec.V3Norm();
						if(CBfrac<=EPS){
							lCB=average(B->MyType->saxes.vec,3)/2;
						}
#if DEBUG_LEVEL>2
						cout << Bvec.V3Str() << " (" << Bvec.V3Norm() << ", " << lCB << ")\n";
#endif
						Bvec/=Bvec.V3Norm(); // unit vector pointing towards connection on B
						double bl=lCA+lCB;
						if(!A_onebond && !B_onebond){
							cout  << "\t\t\t-> Bond between two similar typed multibonded sites.\n"; // no bond length changes here
						} else cout << "\t\t\t-> Bond between " << CA->MyType->name << " and " << CB->MyType->name << " (bond length " << bl << " Angström)\n";
						// check if potentials from group definition need to be preserved and bondlength needs adjustment
						if(nr_potentials>0){
							if(nr_potentials>1){ // only allow one potential
								cout << "ERROR: Only one potential allowed per bond.\n";
								exit(2);
							}
							if(parameters.group_elements[rerun_elements[count-1]].interactions[k].potentials[0]->n!=2){ // artificially limit ourselves to two partners (otherwise things get more complex than I want atm)
								cout << "ERROR: Only two-element potentials are allowed at the moment for super groups.\n";
								exit(3);
							}
							// for all two-element potentials the second parameter is the optimal bond length
							bl=parameters.group_elements[rerun_elements[count-1]].interactions[k].potentials[0]->parameters[1];
							cout << "\t\t\t\t -> Bond length adjusted based on super group potential: " << bl << " Angström\n";
						}
						Vec3 transGroup=CAloc-Avec*bl-CBloc; // remember that Avec points toward A from B
						Mat33 rotGroup=RotAtoB(Bvec,Avec); // rotation matrix to rotate B onto A
#if DEBUG_LEVEL>2
						cout << "Translate:\n(" << transGroup.V3Str(',') << ")\nRotate:\n" << rotGroup.M3Str() << "\n";
#endif
						for(unsigned int l=group_startnr[count-1]; l<group->nr_elements; l++){
							Element* curr_element=&parameters.group_elements[group->elements[l]];
							// rotate around CB's bond location
							RotateElement(curr_element,CBloc,rotGroup);
							// translate so bond location is where it should be
							curr_element->center+=transGroup;
						}
						// do same thing for LOD group (if existent)
						if(group->Type->LOD){
							for(unsigned int l=lodgroup_startnr[count-1]; l<lodgroup->nr_elements; l++){
								Element* curr_element=&parameters.group_elements[lodgroup->elements[l]];
								// rotate around CB's bond location
								RotateElement(curr_element,CBloc,rotGroup);
								// translate so bond location is where it should be
								curr_element->center+=transGroup;
							}
						}
						// now get rid of elements not needed anymore and take care of connections
						if(A_onebond){ // single bonded A gets removed
							delete_nr_elements++;
							delete_elements=(unsigned int*)realloc(delete_elements,delete_nr_elements*sizeof(unsigned int));
							if(!delete_elements){
								cout << "ERROR: O memory, where art thou?\n";
								exit(4);
							}
							delete_elements[delete_nr_elements-1]=connection_A;
							// take care of charge of deleted element - put on opposing end (which takes its place)
							CB->MyType=NewType(CB->MyType);
							CB->mytype=CB->MyType->thistype;
							Vec3 residual_dipole=A->rot*A->MyType->initial_dipole;
							double residual_q=0.0;
							for(unsigned int l=0; l<A->MyType->nr_charges; l++){
								residual_dipole+=(A->rot*A->MyType->q_pos[l])*A->MyType->q[l];
								residual_q+=A->MyType->q[l];
							}
							CB->MyType->initial_dipole+=CB->rot.M3Transpose()*residual_dipole;
							bool found_bond=false;
							for(unsigned int l=0; l<CB->MyType->nr_charges; l++){
								if(CB->MyType->q_pos[l].V3Norm()<EPS){
									found_bond=true;
									CB->MyType->q[l]+=residual_q;
									break;
								}
							}
							if(!found_bond && (residual_q>EPS)){
								CB->MyType->nr_charges++;
								CB->MyType->q_pos=(Vec3*)realloc(CB->MyType->q_pos,CB->MyType->nr_charges*sizeof(Vec3));
								CB->MyType->q=(double*)realloc(CB->MyType->q,CB->MyType->nr_charges*sizeof(double));
								if(!CB->MyType->q_pos || !CB->MyType->q_pos){
									cout << "ERROR: Not enough memory to add charge.\n";
									exit(11);
								}
								CB->MyType->q_pos[CB->MyType->nr_charges-1]=Vec3(0.0);
								CB->MyType->q[CB->MyType->nr_charges-1]=residual_q;
							}
						}
						if(B_onebond){ // single bonded B gets removed
							delete_nr_elements++;
							delete_elements=(unsigned int*)realloc(delete_elements,delete_nr_elements*sizeof(unsigned int));
							if(!delete_elements){
								cout << "ERROR: O memory, where art thou?\n";
								exit(4);
							}
							delete_elements[delete_nr_elements-1]=connection_B;
							// take care of charge of deleted element - put on opposing end (which takes its place)
							CA->MyType=NewType(CA->MyType);
							CA->mytype=CA->MyType->thistype;
							Vec3 residual_dipole=B->rot*B->MyType->initial_dipole;
							double residual_q=0.0;
							for(unsigned int l=0; l<B->MyType->nr_charges; l++){
								residual_dipole+=(B->rot*B->MyType->q_pos[l])*B->MyType->q[l];
								residual_q+=B->MyType->q[l];
							}
							CA->MyType->initial_dipole+=CA->rot.M3Transpose()*residual_dipole;
							bool found_bond=false;
							for(unsigned int l=0; l<CA->MyType->nr_charges; l++){
								if(CA->MyType->q_pos[l].V3Norm()<EPS){
									found_bond=true;
									CA->MyType->q[l]+=residual_q;
									break;
								}
							}
							if(!found_bond && (residual_q>EPS)){
								CA->MyType->nr_charges++;
								CA->MyType->q_pos=(Vec3*)realloc(CA->MyType->q_pos,CA->MyType->nr_charges*sizeof(Vec3));
								CA->MyType->q=(double*)realloc(CA->MyType->q,CA->MyType->nr_charges*sizeof(double));
								if(!CA->MyType->q_pos || !CA->MyType->q_pos){
									cout << "ERROR: Not enough memory to add charge.\n";
									exit(11);
								}
								CA->MyType->q_pos[CA->MyType->nr_charges-1]=Vec3(0.0);
								CA->MyType->q[CA->MyType->nr_charges-1]=residual_q;
							}
						}
						if(!A_onebond && !B_onebond){ // in case both are multibonded (and A=B) B can be removed
							delete_nr_elements++;
							delete_elements=(unsigned int*)realloc(delete_elements,delete_nr_elements*sizeof(unsigned int));
							if(!delete_elements){
								cout << "ERROR: O memory, where art thou?\n";
								exit(4);
							}
							delete_elements[delete_nr_elements-1]=connection_B;
							// take care of connections and charges
							// <TODO>
							// at the moment bail out b/c implementation is quite complex for LOD systems
							cout << "ERROR: Fusion feature not implemented yet.\n";
							exit(8);
						} else{
							// take care of connections
							if(!A_onebond){ // link needs to be added on top of existing ones
								CAlink=CA->nr_interactions;
								CA->nr_interactions++;
								CA->interactions=(Interaction*)realloc(CA->interactions,CA->nr_interactions*sizeof(Interaction));
								if(!CA->interactions){
									cout << "ERROR: Not enough memory to create new link.\n";
									exit(9);
								}
								CA->interactions[CAlink].fixed=false;
								CA->interactions[CAlink].allow_bond_stretch=false;
								CA->interactions[CAlink].allow_bond_bend=false;
								CA->interactions[CAlink].bond_order=1.0;
								CA->interactions[CAlink].nr_potentials=0;
								CA->interactions[CAlink].potentials=NULL;
								CA->interactions[CAlink].fixed=false;
								CA->interactions[CAlink].location=NULL;
								CA->interactions[CAlink].normal=NULL;
								CA->interactions[CAlink].tangent=NULL;
								CA->interactions[CAlink].initial_location=new Vec3;
								*CA->interactions[CAlink].initial_location=Vec3(0.0);
								CA->interactions[CAlink].initial_normal=new Vec3;
								*CA->interactions[CAlink].initial_normal=Vec3(0.0);
								CA->interactions[CAlink].initial_tangent=new Vec3;
								*CA->interactions[CAlink].initial_tangent=Vec3(0.0);
							}
							if(!B_onebond){ // link needs to be added on top of existing ones
								CBlink=CB->nr_interactions;
								CB->nr_interactions++;
								CB->interactions=(Interaction*)realloc(CB->interactions,CB->nr_interactions*sizeof(Interaction));
								if(!CB->interactions){
									cout << "ERROR: Not enough memory to create new link.\n";
									exit(9);
								}
								CB->interactions[CBlink].fixed=false;
								CB->interactions[CBlink].allow_bond_stretch=false;
								CB->interactions[CBlink].allow_bond_bend=false;
								CB->interactions[CBlink].bond_order=1.0;
								CB->interactions[CBlink].nr_potentials=0;
								CB->interactions[CBlink].potentials=NULL;
								CB->interactions[CBlink].fixed=false;
								CB->interactions[CBlink].location=NULL;
								CB->interactions[CBlink].normal=NULL;
								CB->interactions[CBlink].tangent=NULL;
								CB->interactions[CBlink].initial_location=new Vec3;
								*CB->interactions[CBlink].initial_location=Vec3(0.0);
								CB->interactions[CBlink].initial_normal=new Vec3;
								*CB->interactions[CBlink].initial_normal=Vec3(0.0);
								CB->interactions[CBlink].initial_tangent=new Vec3;
								*CB->interactions[CBlink].initial_tangent=Vec3(0.0);
							}
							CA->interactions[CAlink].partner=CBnr;
							CA->interactions[CAlink].back_link=CBlink;
							CB->interactions[CBlink].partner=CAnr;
							CB->interactions[CBlink].back_link=CAlink;
						}
						if(nr_potentials>0){
							Interaction_Potential* orig_potential=parameters.group_elements[rerun_elements[count-1]].interactions[k].potentials[0];
							// preserve partner order as well (i.e. important for "chainspring" potential)
							unsigned int shift=1;
							if(orig_potential->partners[1]==count-1) shift=0;
#if DEBUG_LEVEL>2
							cout << orig_potential->partners[0] << ", " << orig_potential->partners[1] << " -> " << shift << "\n";
#endif
							// create new potential
							parameters.num_potentials++;
							parameters.potentials=(Interaction_Potential**)realloc(parameters.potentials,parameters.num_potentials*sizeof(Interaction_Potential*));
							if(!parameters.potentials){
								cout << "ERROR: Not enough memory to create additional potential.\n";
								exit(11);
							}
							Interaction_Potential* new_potential=new Interaction_Potential;
							new_potential->name=orig_potential->name;
							new_potential->type=orig_potential->type;
							new_potential->n=orig_potential->n;
							new_potential->nr_parameters=orig_potential->nr_parameters;
							new_potential->parameters=orig_potential->parameters; // copy just pointer here (parameters don't change)
							new_potential->partners=new unsigned int[2];
							new_potential->partners[shift]=CAnr;
							new_potential->partners[(1+shift)%2]=CBnr;
#if DEBUG_LEVEL>2
							cout << CAnr << " (" << CAlink << ") : " << CBnr << " (" << CBlink  << ")\n";
#endif
							new_potential->links=new unsigned int[4];
							new_potential->links[0]=0; new_potential->links[3]=0;
							new_potential->links[shift*2+(1+shift)%2]=CAlink; // link nr from CA->CB
							new_potential->links[((1+shift)%2)*2+shift]=CBlink; // link nr from CB->CA
							// Store new potential
							parameters.potentials[parameters.num_potentials-1]=new_potential;
							// Attach to group
							group->nr_potentials++;
							group->potentials=(Interaction_Potential**)realloc(group->potentials,group->nr_potentials*sizeof(Interaction_Potential*));
							if(!group->potentials){
								cout << "ERROR: Not enough memory to create additional potential.\n";
								exit(11);
							}
							group->potentials[group->nr_potentials-1]=new_potential;
							// Attach to both links (forward and backward)
							if(CA->interactions[CAlink].nr_potentials==0) CA->interactions[CAlink].potentials=NULL;
							CA->interactions[CAlink].nr_potentials++;
							CA->interactions[CAlink].potentials=(Interaction_Potential**)realloc(CA->interactions[CAlink].potentials,CA->interactions[CAlink].nr_potentials*sizeof(Interaction_Potential*));
							if(!CA->interactions[CAlink].potentials){
								cout << "ERROR: Not enough memory to create additional potential.\n";
								exit(11);
							}
							CA->interactions[CAlink].potentials[CA->interactions[CAlink].nr_potentials-1]=new_potential;
							if(CB->interactions[CBlink].nr_potentials==0) CB->interactions[CBlink].potentials=NULL;
							CB->interactions[CBlink].nr_potentials++;
							CB->interactions[CBlink].potentials=(Interaction_Potential**)realloc(CB->interactions[CBlink].potentials,CB->interactions[CBlink].nr_potentials*sizeof(Interaction_Potential*));
							if(!CB->interactions[CBlink].potentials){
								cout << "ERROR: Not enough memory to create additional potential.\n";
								exit(11);
							}
							CB->interactions[CBlink].potentials[CB->interactions[CBlink].nr_potentials-1]=new_potential;
						}
						// rotate just added group around bond in order to minimize energy
						Vec3 location=CB->center+CB->rot*(*CB->interactions[CBlink].initial_location);
						double mintheta=0.0;
						bool* connection_sites=new bool[group->nr_elements-group_startnr[count-1]+1];
						for(unsigned int l=0; l<group->nr_elements-group_startnr[count-1]+1; l++) connection_sites[l]=false;
						for(unsigned int l=0; l<rerun_nr_elements; l++){
							int connection=group_connections[(count-1)*rerun_nr_elements+l];
							if(connection!=0){
								if((connection>0) && (curr_connection<=(int)base_add_group->nr_elements)) connection_sites[curr_connection-1]=true;
								if((connection<0) && (-curr_connection<=(int)base_add_group->Type->nr_connection_sites)) connection_sites[base_add_group->Type->connection_sites[(unsigned int)(-curr_connection-1)]]=true;
							}
						}
						double minE=GroupLJdisp(group,delete_elements,delete_nr_elements,connection_sites,group_startnr[0],group_startnr[count-1],mintheta,Avec,location);
						for(unsigned int l=1; l<SUPERGROUP_ROTATE_STEPS; l++){
							double theta=(double)l/((double)SUPERGROUP_ROTATE_STEPS)*2.0*pi;
							double E=GroupLJdisp(group,delete_elements,delete_nr_elements,connection_sites,group_startnr[0],group_startnr[count-1],theta,Avec,location);
#if DEBUG_LEVEL>2
							cout << l << ", " << theta << ": " << E << "\n";
#endif
							if(E<minE){
								mintheta=theta;
								minE=E;
							}
						}
						delete[] connection_sites;
						// rotate around bond to minimum dispersive energy
						cout << "\t\t\t-> Minimum dispersive energy found at " << mintheta/pi << "*pi rad rotated around bond.\n";
						rotGroup=AxisAngle2Rot(Avec,mintheta);
						for(unsigned int l=group_startnr[count-1]; l<group->nr_elements; l++){
							// rotate around CB's bond location
							RotateElement(&parameters.group_elements[group->elements[l]],location,rotGroup);
						}
						// again apply for LOD group as well
						if(group->Type->LOD){
							for(unsigned int l=lodgroup_startnr[count-1]; l<lodgroup->nr_elements; l++){
								// rotate around CB's bond location
								RotateElement(&parameters.group_elements[lodgroup->elements[l]],location,rotGroup);
							}
						}
						// add LOD model to base model
						if(group->Type->LOD){
							// reuse variables here
							unsigned int CA_LOD_nr=lodgroup->Type->LOD->element_in_ellipsoid[lodgroup->levelofdetail-1][CAnr];
							Element* CA_LOD=&parameters.group_elements[lodgroup->elements[CA_LOD_nr]];
							unsigned int CB_LOD_nr=lodgroup->Type->LOD->element_in_ellipsoid[lodgroup->levelofdetail-1][CBnr];
							Element* CB_LOD=&parameters.group_elements[lodgroup->elements[CB_LOD_nr]];
							unsigned int CA_LOD_link=0;
							unsigned int CB_LOD_link=0;
							bool A_newlink=false;
							bool B_newlink=false;
#if DEBUG_LEVEL>0
							cout << "\t\t\t-> LOD connection: " << CA_LOD->MyType->name << " <-> " << CB_LOD->MyType->name << "\n";
#endif
							if(A_onebond){ // existing group element is taken out
								// find ellipsoid the connection site element is from
								int prev_count;
								bool found_prev=false;
								for(prev_count=count-1; prev_count>=0; prev_count--){
									if(addgroups[prev_count]!=NULL){
										if(connection_A>=group_startnr[(unsigned int)prev_count]){
											found_prev=true;
											break;
										}
									}
								}
								int lod_delete_nr=-(connection_A+1);
								if(found_prev){
									Element_Group* prev_group=addgroups[prev_count];
									if(prev_group->levelofdetail>0){
										unsigned int recalc_ellipsoid=prev_group->Type->LOD->element_in_ellipsoid[prev_group->levelofdetail-1][connection_A-group_startnr[prev_count]];
										unsigned int eg_index=prev_group->Type->LOD->ellipsoid_counts[prev_group->levelofdetail-1][recalc_ellipsoid];
										unsigned int nr_ellipsoid_elements=prev_group->Type->LOD->element_groups[prev_group->levelofdetail-1][eg_index];
										if(nr_ellipsoid_elements==1){ // easy case: the deleted connection site also deletes the LOD ellipsoid
											lod_delete_nr=-(lodgroup_startnr[prev_count]+recalc_ellipsoid+1); // negative number indicates deletion
										} else{
											lod_delete_nr=lodgroup_startnr[prev_count]+recalc_ellipsoid+1; // positive number indicates recalculation
#if DEBUG_LEVEL>2
											cout << recalc_ellipsoid+1 << ", " << nr_ellipsoid_elements << "\n";
#endif
										}
									} else lod_delete_nr=-(lodgroup_startnr[prev_count]+connection_A-group_startnr[prev_count]+1);
								}
								if((lod_delete_nr<0) && (CA_LOD->nr_interactions>0)){ // deleting an element means there must've been a link to it one can find
									bool found_link=false;
									for(unsigned int l=0; l<CA_LOD->nr_interactions; l++){
										if((int)CA_LOD->interactions[l].partner+1==-lod_delete_nr){ // found it
											found_link=true;
											CA_LOD_link=l;
										}
									}
									if(!found_link) A_newlink=true; // shouldn't happen, but just to be safe ...
								} else A_newlink=true;
							} else A_newlink=true;
							if(A_newlink){
								CA_LOD_link=CA_LOD->nr_interactions;
								CA_LOD->nr_interactions++;
								CA_LOD->interactions=(Interaction*)realloc(CA_LOD->interactions,CA_LOD->nr_interactions*sizeof(Interaction));
								if(!CA_LOD->interactions){
									cout << "ERROR: Not enough memory to create new link.\n";
									exit(9);
								}
								CA_LOD->interactions[CA_LOD_link].fixed=false;
								CA_LOD->interactions[CA_LOD_link].allow_bond_stretch=false;
								CA_LOD->interactions[CA_LOD_link].allow_bond_bend=false;
								CA_LOD->interactions[CA_LOD_link].bond_order=1.0;
								CA_LOD->interactions[CA_LOD_link].nr_potentials=0;
								CA_LOD->interactions[CA_LOD_link].potentials=NULL;
								CA_LOD->interactions[CA_LOD_link].fixed=false;
							}
							if(B_onebond){ // group to be added is LOD group and one element is taken out
								int lod_delete_nr=-(lodgroup_startnr[count-1]+connection_B-group_startnr[count-1]+1);
								if(add_group->levelofdetail>0){
									// find ellipsoid the connection site element is from
									unsigned int recalc_ellipsoid=add_group->Type->LOD->element_in_ellipsoid[add_group->levelofdetail-1][connection_B-group_startnr[count-1]];
									unsigned int eg_index=add_group->Type->LOD->ellipsoid_counts[add_group->levelofdetail-1][recalc_ellipsoid];
									unsigned int nr_ellipsoid_elements=add_group->Type->LOD->element_groups[add_group->levelofdetail-1][eg_index];
									if(nr_ellipsoid_elements==1){ // easy case: the deleted connection site also deletes the LOD ellipsoid
										lod_delete_nr=-(lodgroup_startnr[count-1]+recalc_ellipsoid+1); // negative number indicates deletion
									} else{
										lod_delete_nr=lodgroup_startnr[count-1]+recalc_ellipsoid+1; // positive number indicates recalculation
#if DEBUG_LEVEL>2
										cout << recalc_ellipsoid+1 << ", " << nr_ellipsoid_elements << "\n";
#endif
									}
								}
								if((lod_delete_nr<0) && (CB_LOD->nr_interactions>0)){ // deleting an element means there must've been a link to it one can find
									bool found_link=false;
									for(unsigned int l=0; l<CB_LOD->nr_interactions; l++){
										if((int)CB_LOD->interactions[l].partner+1==-lod_delete_nr){ // found it
											found_link=true;
											CB_LOD_link=l;
										}
									}
									if(!found_link) B_newlink=true; // shouldn't happen, but just to be safe ...
								} else B_newlink=true;
							} else B_newlink=true;
							if(B_newlink){
								CB_LOD_link=CB_LOD->nr_interactions;
								CB_LOD->nr_interactions++;
								CB_LOD->interactions=(Interaction*)realloc(CB_LOD->interactions,CB_LOD->nr_interactions*sizeof(Interaction));
								if(!CB_LOD->interactions){
									cout << "ERROR: Not enough memory to create new link.\n";
									exit(9);
								}
								CB_LOD->interactions[CB_LOD_link].fixed=false;
								CB_LOD->interactions[CB_LOD_link].allow_bond_stretch=false;
								CB_LOD->interactions[CB_LOD_link].allow_bond_bend=false;
								CB_LOD->interactions[CB_LOD_link].bond_order=1.0;
								CB_LOD->interactions[CB_LOD_link].nr_potentials=0;
								CB_LOD->interactions[CB_LOD_link].potentials=NULL;
								CB_LOD->interactions[CB_LOD_link].fixed=false;
							}
							// take care of connections
							if(!A_onebond && !B_onebond){ // not implemented yet and error thrown earlier already
								// <TODO>
							} else{
								// CAloc-Avec*bl
								CA_LOD->interactions[CA_LOD_link].initial_location=new Vec3;
								*CA_LOD->interactions[CA_LOD_link].initial_location=(CA_LOD->rot.M3Transpose()*(CAloc-CA_LOD->center));
								CA_LOD->interactions[CA_LOD_link].initial_normal=new Vec3;
								*CA_LOD->interactions[CA_LOD_link].initial_normal=Vec3(0.0);
								CA_LOD->interactions[CA_LOD_link].initial_tangent=new Vec3;
								*CA_LOD->interactions[CA_LOD_link].initial_tangent=Vec3(0.0);
								CA_LOD->interactions[CA_LOD_link].location=new Vec3;
								*CA_LOD->interactions[CA_LOD_link].location=(CA_LOD->rot*(*CA_LOD->interactions[CA_LOD_link].initial_location));
								CA_LOD->interactions[CA_LOD_link].normal=NULL;
								CA_LOD->interactions[CA_LOD_link].tangent=NULL;
								
								CB_LOD->interactions[CB_LOD_link].initial_location=new Vec3;
								*CB_LOD->interactions[CB_LOD_link].initial_location=(CB_LOD->rot.M3Transpose()*(CAloc-Avec*bl-CB_LOD->center));
								CB_LOD->interactions[CB_LOD_link].initial_normal=new Vec3;
								*CB_LOD->interactions[CB_LOD_link].initial_normal=Vec3(0.0);
								CB_LOD->interactions[CB_LOD_link].initial_tangent=new Vec3;
								*CB_LOD->interactions[CB_LOD_link].initial_tangent=Vec3(0.0);
								CB_LOD->interactions[CB_LOD_link].location=new Vec3;
								*CB_LOD->interactions[CB_LOD_link].location=(CB_LOD->rot*(*CB_LOD->interactions[CB_LOD_link].initial_location));
								CB_LOD->interactions[CB_LOD_link].normal=NULL;
								CB_LOD->interactions[CB_LOD_link].tangent=NULL;
								
								CA_LOD->interactions[CA_LOD_link].partner=CB_LOD_nr;
								CA_LOD->interactions[CA_LOD_link].back_link=CB_LOD_link;
								CB_LOD->interactions[CB_LOD_link].partner=CA_LOD_nr;
								CB_LOD->interactions[CB_LOD_link].back_link=CA_LOD_link;
#if DEBUG_LEVEL>2
								cout << (CB_LOD->center-CA_LOD->center+CB_LOD->rot*(*CB_LOD->interactions[CB_LOD_link].initial_location)-CA_LOD->rot*(*CA_LOD->interactions[CA_LOD_link].initial_location)).V3Norm() << "\n";
#endif
							}
						}
					}
				}
			}
			// need to copy non-LOD'ed stuff over to the associated LOD group (if there is one) as well (basegroup elements are already in there from earlier pass)
			if(copy_lod && !found_group){ // this is really all the logic needed here since we already made sure that if found_group is false found will be true
#if DEBUG_LEVEL>0
				cout << "-> Adding element " << parameters.group_elements[group->elements[count-1]].MyType->name << " to LOD level\n";
#endif
				lodgroup->nr_elements++;
				lodgroup->elements = (unsigned int*)realloc(group->elements,group->nr_elements*sizeof(unsigned int));
				if(!lodgroup->elements){
					cout << "ERROR: Not enough memory to hold elements numbers of group " << group->Type->name << "\n";
					exit(1);
				}
				nr_group_elements++;
				parameters.group_elements = (Element*)realloc(parameters.group_elements,nr_group_elements*sizeof(Element));
				if(!parameters.group_elements){ // Danger, danger
					cout << "ERROR: Could not find a memory block large enough to hold a total of " << nr_group_elements << " group associated elements.\n";
					exit(1);
				}
				lodgroup->elements[lodgroup->nr_elements-1]=nr_group_elements-1;
				CopyElement(&parameters.group_elements[group->elements[count-1]],&parameters.group_elements[lodgroup->elements[lodgroup->nr_elements-1]],lodgroup->nr_elements-count);
				parameters.group_elements[lodgroup->elements[lodgroup->nr_elements-1]].group=lodgroup;
			}
			if(!rerun){ // add one element (or fake element in case of super group)
				// do accounting
				parameters.group_elements[nr_group_elements-1].nr_interactions=0; // none read yet
				parameters.group_elements[nr_group_elements-1].interactions=NULL; // same realloc-needing-NULL-thing
				parameters.group_elements[nr_group_elements-1].group=group;
				if(!found || found_group){ // Entity is either not existent yet or a group: We have a super group
					parameters.group_elements[nr_group_elements-1].mytype=0;
					parameters.group_elements[nr_group_elements-1].MyType=NULL;
				} else{
					if(!finished){ // do not allow to add element after a group has been added (too complicated at the moment ...)
						cout << "ERROR: Adding an element to a super group once the first group has been added is not allowed.\n";
						exit(5);
					}
					parameters.group_elements[nr_group_elements-1].mytype=parameters.element_types[j]->thistype;
					parameters.group_elements[nr_group_elements-1].MyType=parameters.element_types[j];
					
					// read user definable stuff
					
					// position of element in group
					item="position."+to_string(group->nr_elements);
					if(fileimport) SetParam(parameters.group_elements[nr_group_elements-1].center.vec,item.c_str(),content,Vec3(0.0)); else SetParam(parameters.group_elements[nr_group_elements-1].center.vec,item.c_str(),conf,Vec3(0.0));
					
					// rotation angles
					Vec4 rot_angles; // axis, angle
					item="rotation."+to_string(group->nr_elements);
					SetParam(rot_angles.vec,4,item.c_str(),conf,Vec4(1.0,0.0,0.0,0.0).vec);
					parameters.group_elements[nr_group_elements-1].rot=AxisAngle2Rot(rot_angles);
					
					parameters.group_elements[nr_group_elements-1].dipole=parameters.group_elements[nr_group_elements-1].rot*parameters.group_elements[nr_group_elements-1].MyType->initial_dipole;
				}
			}
			element_name="";
			element=false;
		} else{
			if((element) || ((element_list[i]!=' ') && (element_list[i]!='\t'))){
				element=true;
				element_name+=element_list[i];
			}
		}
	} // finished creating group elements
	if(group_connections) delete[] group_connections;
	
	// Now create connectivity
	if(rerun){
		unsigned int lgnr=group->nr_elements; // delete_from_array changes that number so we need a copy
		// now sort by element nr
#if DEBUG_LEVEL>1
		cout << "-> Deleting placeholder elements and removed connection site elements ... ";
		cout.flush();
#endif
		for(i=0; i<delete_nr_elements-1; i++){
			for(j=0; j<delete_nr_elements-i-1; j++){
				if(delete_elements[j]>delete_elements[j+1]){ // current element number is larger then next, exchange
					delete_elements[j]^=delete_elements[j+1]; // a xor b
					delete_elements[j+1]^=delete_elements[j]; // b xor (a xor b) = a
					delete_elements[j]^=delete_elements[j+1]; // (a xor b) xor a = b
				}
			} // last element is biggest of them all (hence no need to go there again)
		}
		// delete elements that are not in use anymore (start from last, so shifting things preserves array numbers
		for(int di=delete_nr_elements-1; di>=0; di--){
			// all elements above the elements to be deleted will be shifted down one place ...
			for(j=group->nr_elements-1; j>delete_elements[di]; j--){
				Element* curr_element=&parameters.group_elements[group->elements[j]];
				for(k=0; k<curr_element->nr_interactions; k++){ // ... go through all interactions on these elements ...
					Interaction* inter=&curr_element->interactions[k];
					// ... adjust partner numbers down if above the element to be deleted
					if(inter->partner>delete_elements[di]){
						inter->partner--;
					} else{ // do not forget that partner nr changes on other side if it is below the element to be deleted
						parameters.group_elements[group->elements[inter->partner]].interactions[inter->back_link].partner--;
					}
				}
			}
			for(j=0; j<group->nr_potentials; j++){
				for(unsigned int l=0; l<group->potentials[j]->n; l++){
					if(group->potentials[j]->partners[l]>delete_elements[di]) group->potentials[j]->partners[l]--;
				}
			}
			// finally delete element (and shift down elements)
			delete_from_array(&(group->elements),group->nr_elements,delete_elements[di]);
		}
#if DEBUG_LEVEL>1
		cout << "Done.\n";
#endif
		if(group->Type->LOD){
			unsigned int* recalc_ellipsoid=(unsigned int*)malloc(lodgroup->nr_elements*sizeof(unsigned int));
			if(!recalc_ellipsoid){
				cout << "ERROR: Not enough memory to store which ellipsoids to recalculate.\n";
				exit(11);
			}
			memset(recalc_ellipsoid,0,lodgroup->nr_elements*sizeof(unsigned int));
			*group->Type->LOD->transparency/=added_groups;
			lodgroup->Type->transparency=*group->Type->LOD->transparency;
			group->Type->LOD->scale_original_vdw/=added_groups;
#if DEBUG_LEVEL>1
			cout << "-> Deleting placeholder LOD elements and recalculating ellipsoids if necessary\n";
#endif
			// delete elements first
			Level_of_Detail* LOD=lodgroup->Type->LOD;
			unsigned int egnr=0;
			for(i=0; i<LOD->nr_ellipsoids[lodgroup->levelofdetail-1]; i++) egnr+=LOD->element_groups[lodgroup->levelofdetail-1][LOD->ellipsoid_counts[lodgroup->levelofdetail-1][i]]+1;
			for(int di=delete_nr_elements-1; di>=0; di--){
				unsigned int del_nr=delete_elements[di]+1;
				// adjust LOD meta-information
				for(int l=LOD->nr_ellipsoids[lodgroup->levelofdetail-1]-1; l>=0; l--){
					for(int m=LOD->element_groups[lodgroup->levelofdetail-1][LOD->ellipsoid_counts[lodgroup->levelofdetail-1][l]]-1; m>=0; m--){ // need an int here b/c value can get negative
						if(LOD->element_groups[lodgroup->levelofdetail-1][LOD->ellipsoid_counts[lodgroup->levelofdetail-1][l]+m+1]>(int)del_nr){
							LOD->element_groups[lodgroup->levelofdetail-1][LOD->ellipsoid_counts[lodgroup->levelofdetail-1][l]+m+1]--;
						} else{
							if(LOD->element_groups[lodgroup->levelofdetail-1][LOD->ellipsoid_counts[lodgroup->levelofdetail-1][l]+m+1]==(int)del_nr){
								delete_from_array(&(LOD->element_groups[lodgroup->levelofdetail-1]),egnr,LOD->ellipsoid_counts[lodgroup->levelofdetail-1][l]+m+1);
								LOD->element_groups[lodgroup->levelofdetail-1][LOD->ellipsoid_counts[lodgroup->levelofdetail-1][l]]--;
								// get rid of ellipsoids that don't exist anymore
								bool delete_ellipsoid=false;
								if(LOD->element_groups[lodgroup->levelofdetail-1][LOD->ellipsoid_counts[lodgroup->levelofdetail-1][l]]==0){
#if DEBUG_LEVEL>2
									cout << "\t-> Deleting LOD element nr " << l+1 << " ... ";
									cout.flush();
#endif
									// all elements above the elements to be deleted will be shifted down one place ...
									for(j=lodgroup->nr_elements-1; (int)j>l; j--){
										Element* curr_element=&parameters.group_elements[lodgroup->elements[j]];
										for(k=0; k<curr_element->nr_interactions; k++){ // ... go through all interactions on these elements ...
											Interaction* inter=&curr_element->interactions[k];
											// ... adjust partner numbers down if above the element to be deleted
											if((int)inter->partner>l){
												inter->partner--;
											} else{ // do not forget that partner nr changes on other side if it is below the element to be deleted
												parameters.group_elements[lodgroup->elements[inter->partner]].interactions[inter->back_link].partner--;
											}
										}
									}
									for(j=0; j<lodgroup->nr_potentials; j++){
										for(k=0; k<lodgroup->potentials[j]->n; k++){
											if((int)lodgroup->potentials[j]->partners[k]>l) lodgroup->potentials[j]->partners[k]--;
										}
									}
									unsigned int lnr=lodgroup->nr_elements; // need a copy b/c delete_from_array changes it
									delete_from_array(&recalc_ellipsoid,lnr,(unsigned int)l);
									lnr=lodgroup->nr_elements; // need another copy (see above)
									delete_from_array(&(LOD->reduce_electrostatics[lodgroup->levelofdetail-1]),lnr,(unsigned int)l);
									delete_from_array(&(lodgroup->elements),lodgroup->nr_elements,(unsigned int)l);
									delete_from_array(&(LOD->element_groups[lodgroup->levelofdetail-1]),egnr,LOD->ellipsoid_counts[lodgroup->levelofdetail-1][l]);
									delete_from_array(&(LOD->ellipsoid_counts[lodgroup->levelofdetail-1]),LOD->nr_ellipsoids[lodgroup->levelofdetail-1],l);
									delete_ellipsoid=true;
									// need to adjust element_in_ellipsoid
									for(j=0; j<lgnr; j++) if(LOD->element_in_ellipsoid[lodgroup->levelofdetail-1][j]>l) LOD->element_in_ellipsoid[lodgroup->levelofdetail-1][j]--;
#if DEBUG_LEVEL>2
									cout << "Done.\n";
#endif
								} else{ // LOD ellipsoid needs to be recalculated later
									recalc_ellipsoid[l]=1;
								}
								for(k=l+(!delete_ellipsoid); k<LOD->nr_ellipsoids[lodgroup->levelofdetail-1]; k++){ // adjust ellipsoid_counts
									LOD->ellipsoid_counts[lodgroup->levelofdetail-1][k]-=1+(delete_ellipsoid);
								}
							}
						}
					}
				}
				// delete what can be deleted
				delete_from_array(&(LOD->element_in_ellipsoid[lodgroup->levelofdetail-1]),lgnr,del_nr-1);
			}
			string storage_filename=parameters.fileout+".conf.dat";
			readfile input;
			input.filename=storage_filename;
			input.directory="";
			input.file.open(storage_filename.c_str(), fstream::in);
			char* confext=conf;
			if(!input.file.fail()){
				unsigned int pos=0;
				string subname=group->Type->name;
				char* section=GetSection(&input,"Detail",pos,subname);
				if(section){
					confext=new char[strlen(section)+strlen(conf)+1];
					confext[strlen(section)+strlen(conf)]='\0';
					strncpy(confext,section,strlen(section));
					char* temp=confext+strlen(section);
					strncpy(temp,conf,strlen(conf));
					delete[] section;
				}
				input.file.close();
			}
			LOD->volumes[lodgroup->levelofdetail-1]=new double[lodgroup->nr_elements];
			for(j=0; j<lodgroup->nr_elements; j++) LOD->volumes[lodgroup->levelofdetail-1][j]=-1.0;
			LOD->epsilons[lodgroup->levelofdetail-1]=new double[lodgroup->nr_elements];
			for(j=0; j<lodgroup->nr_elements; j++) LOD->epsilons[lodgroup->levelofdetail-1][j]=-1.0;
			double* ellipsoid_volumes=get_MxN_tupel(("lod_volumes."+int2str(lodgroup->levelofdetail)).c_str(),confext,lodgroup->nr_elements,1,NULL);
			double* ellipsoid_epsilons=get_MxN_tupel(("lod_epsilons."+int2str(lodgroup->levelofdetail)).c_str(),confext,lodgroup->nr_elements,1,NULL);
			SetParam(LOD->gyration_sphere,"gyration_sphere",conf,false);
			// Calculate center of group
			parameters.group_centers[nr] = GroupCenter(group,parameters.group_elements);
			// Determine gyration tensor because we need it sooner than normal ...
			Vec3 points[6];
			Mat33 S;
			S.M3Zeros();
			double gcount=0.0;
			for(i=0; i<group->nr_elements; i++){
				Element* curr=&parameters.group_elements[group->elements[i]];
				Vec3 center=curr->center-parameters.group_centers[nr];
				for(unsigned int l=0; l<3; l++){
					points[l<<1].vec[0]=(center.vec[0]+curr->rot.mat[l][0]*curr->MyType->saxes.vec[0]);
					points[(l<<1)+1].vec[0]=(center.vec[0]-curr->rot.mat[l][0]*curr->MyType->saxes.vec[0]);
					
					points[l<<1].vec[1]=(center.vec[1]+curr->rot.mat[l][1]*curr->MyType->saxes.vec[1]);
					points[(l<<1)+1].vec[1]=(center.vec[1]-curr->rot.mat[l][1]*curr->MyType->saxes.vec[1]);
					
					points[l<<1].vec[2]=(center.vec[2]+curr->rot.mat[l][2]*curr->MyType->saxes.vec[2]);
					points[(l<<1)+1].vec[2]=(center.vec[2]-curr->rot.mat[l][2]*curr->MyType->saxes.vec[2]);
				}
				for(k=0; k<6; k++){
					for(unsigned int m=0; m<3; m++){
						for(unsigned int n=0; n<3; n++){
							// S_mn = sum_i r_m(i)*r_n(i) (corrected for center
							S.mat[m][n]+=points[k].vec[m]*points[k].vec[n];
						}
					}
					gcount+=1.0;
				}
			}
			S*=3.0/gcount; // since this gets square rooted later -> *3^(1/2) gives effective radius and principle axes (is needed to pass sanity check of single ellipsoid points)
			// get eigenvalues for ellipse semiaxis
			CVec3 eigenvalues=S.Eigenvalues();
			Vec3 semiaxis=eigenvalues.Re(); // lambda^2's
			group->rT=sqrt(EllipsoidSurface(sqrt(semiaxis.vec[0]),sqrt(semiaxis.vec[1]),sqrt(semiaxis.vec[2]))/(4.0*pi));
			double rT=parameters.rT;
			SetParam(rT,"rT",conf,rT);
			for(i=0; i<lodgroup->nr_elements; i++){
				Element* source=&parameters.group_elements[lodgroup->elements[i]];
				if(lodgroup->Type->LOD->reduce_electrostatics){
					source->MyType->nr_charges=0;
					source->MyType->dipole=0.0;
				}
				if(ellipsoid_volumes) LOD->volumes[lodgroup->levelofdetail-1][i]=ellipsoid_volumes[i];
				if(ellipsoid_epsilons) LOD->epsilons[lodgroup->levelofdetail-1][i]=ellipsoid_epsilons[i];
				if(recalc_ellipsoid[i]>0){
					unsigned int nrce=(unsigned int)LOD->element_groups[lodgroup->levelofdetail-1][LOD->ellipsoid_counts[lodgroup->levelofdetail-1][i]];
					unsigned int* combine_elements = new unsigned int[nrce];
					for(j=0; j<nrce; j++){
						combine_elements[j]=LOD->element_groups[lodgroup->levelofdetail-1][LOD->ellipsoid_counts[lodgroup->levelofdetail-1][i]+j+1];
						if((combine_elements[j]>0) && (combine_elements[j]<=group->nr_elements)){
							combine_elements[j]=group->elements[combine_elements[j]-1];
						} else{ // Should not happen here
							cout << "Specified elements need to be between 1 and the total number of elements in the group.\n";
							exit(2);
						}
					}
					string specifier=":"+int2str(i+1);
					CombineElements(group,lodgroup,combine_elements,nrce,specifier,LOD->volumes[lodgroup->levelofdetail-1][i],LOD->epsilons[lodgroup->levelofdetail-1][i],rT);
					source=&parameters.group_elements[lodgroup->elements[i]]; // address may have changed
					delete[] combine_elements;
					// Recalculated element is new element at end of group which we now need to properly embed into existing LOD group (including links, etc.)
					Element* dest=&parameters.group_elements[nr_group_elements-1];
					CopyElement(source,dest,0,false);
					// need to adjust bond locations to new rotation matrix
					for(j=0; j<dest->nr_interactions; j++){
						Interaction* inter=&dest->interactions[j];
						if(inter->location) *inter->initial_location=dest->rot.M3Transpose()*(source->center+source->rot*(*inter->initial_location)-dest->center);
						if(inter->normal) *inter->initial_normal*=dest->rot.M3Transpose()*source->rot;
						if(inter->tangent) *inter->initial_tangent*=dest->rot.M3Transpose()*source->rot;
					}
					// replace previous ellipsoid in LOD group
					lodgroup->elements[i]=nr_group_elements-1;
					// delete last element again
					lodgroup->nr_elements--;
					lodgroup->elements=(unsigned int*)realloc(lodgroup->elements,lodgroup->nr_elements*sizeof(unsigned int));
					if(!lodgroup->elements){
						cout << "ERROR: Not enough memory to shrink (!) LOD element block.\n";
						exit(11);
					}
				}
			}
			if(!parameters.fit2lod && lodgroup->Type->LOD->match_epsilon){
				// read in Lennard-Jones epsilons
				double* epsilons=get_flex_tupel(("lod_epsilons."+int2str(lodgroup->levelofdetail)).c_str(),confext);
				if(epsilons){
					cout << "-> Setting Lennard Jones epsilons for LOD ellipsoids.\n";
					if((unsigned int)epsilons[0]==lodgroup->nr_elements){
						for(j=0; j<(unsigned int)epsilons[0]; j++){
							parameters.group_elements[lodgroup->elements[j]].MyType->Vvdw=epsilons[j+1];
							cout << "\t-> Ellipsoid <" << parameters.group_elements[lodgroup->elements[j]].MyType->name << ">: " << parameters.group_elements[lodgroup->elements[j]].MyType->Vvdw << "\n";;
						}
					} else{
						cout << "ERROR: Not enough LJ epsilon parameters in *.conf.dat file. Please rerun fit2lod.\n";
						exit(3);
					}
					delete[] epsilons;
				}
				// read in Lennard-Jones epsilon coefficients
				double* eps_coeff=get_flex_tupel(("lod_epsilon_of_rT."+int2str(lodgroup->levelofdetail)).c_str(),confext);
				if(eps_coeff){
					cout << "-> Setting Lennard Jones epsilon coefficients for LOD ellipsoids.\n";
					count=0;
					j=0;
					while(eps_coeff[count]>0){
						unsigned int e_nr=(unsigned int)eps_coeff[count];
						if(e_nr==8){
							if(j<lodgroup->nr_elements){
								if(parameters.group_elements[lodgroup->elements[j]].MyType->Vvdw_coefficients!=NULL) delete[] (parameters.group_elements[lodgroup->elements[j]].MyType->Vvdw_coefficients);
								parameters.group_elements[lodgroup->elements[j]].MyType->Vvdw_coefficients=new double[e_nr];
								parameters.group_elements[lodgroup->elements[j]].MyType->nr_Vvdw_coefficients=e_nr;
								cout << "\t-> Ellipsoid <" << parameters.group_elements[lodgroup->elements[j]].MyType->name << ">: ";
								for(unsigned int l=0; l<e_nr; l++){
									parameters.group_elements[lodgroup->elements[j]].MyType->Vvdw_coefficients[l]=eps_coeff[count+1+l];
									if(l>0) cout << ", ";
									cout << parameters.group_elements[lodgroup->elements[j]].MyType->Vvdw_coefficients[l];
								}
								cout << "\n";
								count+=e_nr+1;
							} else{
								cout << "WARNING: Different number of LOD ellipoids in configuration.\nPlease make sure you run \"fit2lod <configuration file>\".\n";
								break;
							}
						} else{
							if(e_nr>0){ // Single spheres will not get coefficients which is OK
								cout << "WARNING: Number of epsilon coefficients is not in conformance with this code.\nPlease make sure you run \"fit2lod <configuration file>\".\n";
								break;
							}
						}
						j++;
					}
					delete[] eps_coeff;
				}
				double* IA_coeff=get_flex_tupel(("lod_epsilon_of_rT."+int2str(lodgroup->levelofdetail)).c_str(),confext);
				if(IA_coeff){
					cout << "-> Setting Interaction Area coefficients for LOD ellipsoids.\n";
					count=0;
					j=0;
					while(IA_coeff[count]>0){
						unsigned int e_nr=(unsigned int)IA_coeff[count];
						if(e_nr==7){
							if(j<lodgroup->nr_elements){
								if(parameters.group_elements[lodgroup->elements[j]].MyType->IA_coefficients!=NULL) delete[] (parameters.group_elements[lodgroup->elements[j]].MyType->IA_coefficients);
								parameters.group_elements[lodgroup->elements[j]].MyType->IA_coefficients=new double[e_nr];
								cout << "\t-> Ellipsoid <" << parameters.group_elements[lodgroup->elements[j]].MyType->name << ">: ";
								for(unsigned int l=0; l<e_nr; l++){
									parameters.group_elements[lodgroup->elements[j]].MyType->IA_coefficients[l]=IA_coeff[count+1+l];
									if(l>0) cout << ", ";
									cout << parameters.group_elements[lodgroup->elements[j]].MyType->IA_coefficients[l];
								}
								cout << "\n";
								count+=e_nr+1;
							} else{
								cout << "WARNING: Different number of LOD ellipoids in configuration.\nPlease make sure you run \"fit2lod <configuration file>\".\n";
								break;
							}
						} else{
							cout << "WARNING: Number of IA coefficients is not in conformance with this code.\nPlease make sure you run \"fit2lod <configuration file>\".\n";
							break;
						}
						j++;
					}
					delete[] IA_coeff;
				}
				if(LOD->use_epsilon_texture[lodgroup->levelofdetail-1]){
					// set textures if there are any
					cout << "-> Setting epsilon textures for LOD ellipsoids (may take a while) ...\n";
					double* textures=get_flex_tupel(("lod_textures."+int2str(lodgroup->levelofdetail)).c_str(),confext);
					if(textures){
						count=0;
						j=0;
						while(textures[count]>0){
							unsigned int e_nr=(unsigned int)textures[count];
							if(e_nr==phi_res*theta_res+4){
								if(j<lodgroup->nr_elements){
									cout << "\t-> Setting texture for ellipsoid <" << parameters.group_elements[lodgroup->elements[j]].MyType->name << "> ... ";
									cout.flush();
									parameters.group_elements[lodgroup->elements[j]].MyType->eps_texture=new double[e_nr];
									for(unsigned int l=0; l<e_nr; l++) parameters.group_elements[lodgroup->elements[j]].MyType->eps_texture[l]=textures[count+1+l];
									cout << "Done\n";
									count+=e_nr+1;
								} else{
									cout << "WARNING: Different number of LOD ellipoids in configuration.\nPlease make sure you run \"fit2lod <configuration file>\".\n";
									break;
								}
							} else{
								cout << "WARNING: Number of pixels in epsilon texture is not in conformance with this code.\nPlease make sure you run \"fit2lod <configuration file>\".\n";
								break;
							}
							j++;
						}
						delete[] textures;
					}
				}
			}
			if(confext!=conf) delete[] confext;
			if(ellipsoid_volumes) delete[] ellipsoid_volumes;
			if(ellipsoid_epsilons) delete[] ellipsoid_epsilons;
			if(lodgroup->Type->LOD->reduce_electrostatics) SetEllipsoidCharges(group,lodgroup);
			if(recalc_ellipsoid) free(recalc_ellipsoid);
		}
#if DEBUG_LEVEL>0
		cout << "-> Super group consists of " << group->nr_elements << " elements";
		if(group->Type->LOD) cout << " (LOD group: " << lodgroup->nr_elements << " elements)";
		cout << ".\n";
#endif
		group->Type->bond_range/=added_groups;
		SetParam(group->Type->bond_range,"bond_range",conf,group->Type->bond_range);
#if DEBUG_LEVEL>2
		for(i=0; i<group->nr_elements; i++){
			Element* element=&parameters.group_elements[group->elements[i]];
			if(element->MyType) cout << i+1 << "\t" << element->MyType->name << "\n"; else cout << i+1 << "\n";
			for(j=0; j<element->nr_interactions; j++) cout << "\t" << element->interactions[j].partner+1 << "\n";
		}
#endif
		finished=true;
	} else{
		double* connection_sites = get_flex_tupel("connection_sites",conf);
		if(connection_sites){
			if(!finished){
				cout << "ERROR: Connection sites can not be defined for super groups.\n";
				exit(3);
			}
			count=0;
			group->Type->nr_connection_sites=0;
			while(connection_sites[count]>0){
				group->Type->nr_connection_sites+=(unsigned int)connection_sites[count];
				count+=(unsigned int)connection_sites[count]+1;
			}
			if(group->Type->nr_connection_sites>0){
				group->Type->connection_sites=new unsigned int[group->Type->nr_connection_sites];
				count=0;
				i=0;
				j=0;
				while(connection_sites[count]>0){
					if(((unsigned int)connection_sites[count+i+1]>0) && ((unsigned int)connection_sites[count+i+1]<=group->nr_elements)){
						group->Type->connection_sites[j]=(unsigned int)connection_sites[count+i+1];
					} else{
						cout << "ERROR: Connection sites can only be chosen from the 1-" << group->nr_elements << " elements of the group.\n";
						exit(3);
					}
					i++;
					j++;
					if(i==(unsigned int)connection_sites[count]){
						i=0;
						count+=(unsigned int)connection_sites[count]+1;
					}
				}
			}
			delete[] connection_sites;
		}
		double* connectivity=NULL;
		double* bond_order=NULL;
		if(fileimport) connectivity = get_flex_tupel("connectivity",content); else connectivity = get_flex_tupel("connectivity",conf);
		if(connectivity==NULL){
			cout << "WARNING: Connectivity for group <" << group->Type->name << "> not specified.\n";
			if(!group->Type->rand_independent) cout << "Group elements will be randomized independently.\n";
			group->Type->rand_independent=true;
		} else{ // only look for bond_order in case connectivity is defined
			char** bondorder_typenames=new char*[nr_bondorder_types];
			unsigned int* max_type_nr= new unsigned int[nr_bondorder_types];
			for(i=0; i<nr_bondorder_types; i++){
				string name=bondorder_types[i];
				bondorder_typenames[i]=stringNcopy(name.c_str(),name.length());
				max_type_nr[i]=0;
			}
			if(fileimport) bond_order = get_flex_tupel("bond_order",content,bondorder_typenames,max_type_nr,nr_bondorder_types); else bond_order = get_flex_tupel("bond_order",conf,bondorder_typenames,max_type_nr,nr_bondorder_types);
			// clean up
			for(i=0; i<nr_bondorder_types; i++) delete[] bondorder_typenames[i];
			delete[] bondorder_typenames;
			delete[] max_type_nr;
		}
		unsigned int* connected = new unsigned int[group->nr_elements-1]; // hold elements to be connect to another element for sorting (there can only be a maximum of n-1 links)
		double* connection_order = new double[group->nr_elements-1];
		count = 0;
		unsigned int connections, additionals;
		Interaction* interaction;
		
		double* charge_override=NULL;
		charge_override = get_flex_tupel("charge_override",conf);
		if(fileimport && !(charge_override)) charge_override = get_flex_tupel("import_charges",content);
		
		double* static_link_conf = get_flex_tupel("static_links",conf);
		unsigned int nr_statics=0;
		unsigned int* static_links=NULL;
		if(static_link_conf){
			while(static_link_conf[count]>0){
				if((unsigned int)static_link_conf[count]==2){
					if(!((((unsigned int)static_link_conf[count+1]<1) || ((unsigned int)static_link_conf[count+2]<1))
					|| (((unsigned int)static_link_conf[count+1]>group->nr_elements) || ((unsigned int)static_link_conf[count+2]>group->nr_elements)))){
						nr_statics++;
					}
				}
				count+=(unsigned int)static_link_conf[count]+1;
			}
			count=0;
			static_links=new unsigned int[nr_statics*2];
			nr_statics=0;
			while(static_link_conf[count]>0){
				if((unsigned int)static_link_conf[count]==2){
					if((((unsigned int)static_link_conf[count+1]<1) || ((unsigned int)static_link_conf[count+2]<1))
					|| (((unsigned int)static_link_conf[count+1]>group->nr_elements) || ((unsigned int)static_link_conf[count+2]>group->nr_elements))){
						cout << "WARNING: This group has elements number 1-" << group->nr_elements << ", ignoring static link entry #" << nr_statics+1 << ".\n";
					} else{
						static_links[nr_statics*2]=(unsigned int)static_link_conf[count+1];
						static_links[nr_statics*2+1]=(unsigned int)static_link_conf[count+2];
						nr_statics++;
					}
				} else{
					cout << "WARNING: A static link is only definable between two elements, ignoring entry #" << nr_statics+1 << ".\n";
				}
				count+=(unsigned int)static_link_conf[count]+1;
			}
			delete[] static_link_conf;
		}
		
		count=0;
		bool is_in_statics;
		for(i=0; i<group->nr_elements; i++){
			is_in_statics=false;
			for(j=0; j<2*nr_statics; j++){
				if(i+1==static_links[j]){
					is_in_statics=true;
					break;
				}
			}
			if(charge_override){
				if(i<(unsigned int)charge_override[0]){
					if(parameters.group_elements[group_elements_start+i].MyType->q) free(parameters.group_elements[group_elements_start+i].MyType->q);
					if(parameters.group_elements[group_elements_start+i].MyType->q_pos) free(parameters.group_elements[group_elements_start+i].MyType->q_pos);
					parameters.group_elements[group_elements_start+i].MyType->nr_charges=0;
					if(fabs(charge_override[i+1])>EPS){
						parameters.group_elements[group_elements_start+i].MyType->q=(double*)malloc(sizeof(double));
						parameters.group_elements[group_elements_start+i].MyType->q_pos=(Vec3*)malloc(sizeof(Vec3));
						parameters.group_elements[group_elements_start+i].MyType->nr_charges=1;
						parameters.group_elements[group_elements_start+i].MyType->q[0]=charge_override[i+1]*e_in_esu;
						parameters.group_elements[group_elements_start+i].MyType->q_pos[0]=Vec3(0.0);
					}
				}
			}
			if(connectivity){
				connections=(unsigned int)connectivity[count];
				if(bond_order){ // make sure the definition for bond_order follows that for the connectivity
					if(fabs(bond_order[count]-connectivity[count])>EPS){
						cout << "ERROR: Definitions for <bond_order> needs to follow that of <connectivity>.\n";
						exit(7);
					}
				}
			} else connections=0;
			// get user specified connections
			k=0;
			additionals=0;
			for(j=0; j<connections; j++){
				count++;
				if(((int)connectivity[count]<1) || ((int)connectivity[count]>(int)group->nr_elements)){
					cout << "Trying to connect group element " << i+1 << " to element number " << (int)connectivity[count] << " which does not exist. This group has elements number 1-" << group->nr_elements << ".\n";
					exit(3);
				}
				connected[k]=(int)connectivity[count];
				if(bond_order) connection_order[k]=bond_order[count]; else connection_order[k]=1.0;
				// take care of user-specified doublettes
				found=false;
				for(unsigned int l=0; l<k; l++){
					if(connected[l]==connected[k]){ // found existing link to same element -- pesky users ;-)
						found=true;
						break;
					}
				}
				if(found){ // complain and fix
					cout << "-> WARNING: Ignore multiple links to group element " << connected[k] << " from element " << i+1 << "\n";
					additionals++;
				} else k++;
			} // user specified connections are now in array connected
			connections-=additionals;
			
			// make sure to incorporate already existing interactions
			additionals=0;
			unsigned int partner;
			double pbo;
			for(j=0; j<parameters.group_elements[group_elements_start+i].nr_interactions; j++){
				partner=parameters.group_elements[group_elements_start+i].interactions[j].partner+1;
				pbo=parameters.group_elements[group_elements_start+i].interactions[j].bond_order;
				// try to find existing link in user specified list
				// - the ones found here are temporaries created earlier in the loop, look below)
				// - all links in connected will be recreated
				found=false;
				for(k=0; k<connections; k++){
					if(connected[k]==partner){
						found=true;
						break;
					}
				}
				if(!found){ // in case the user did not specify a link to an element which links to this element, add it to the list and complain
					additionals++;
					connected[connections+additionals-1]=partner;
					connection_order[connections+additionals-1]=pbo;
#if DEBUG_LEVEL>1
					cout << "-> Added previously defined connectivity: group element " << i+1 << " is connected to element " << partner << ".\n";
#endif
				}
			}
			connections+=additionals;
			
			if(connections>0){
				parameters.group_elements[group_elements_start+i].nr_interactions=connections;
				parameters.group_elements[group_elements_start+i].interactions = (Interaction*)realloc(parameters.group_elements[group_elements_start+i].interactions,connections*sizeof(Interaction));
			}
			interaction = parameters.group_elements[group_elements_start+i].interactions;
			
			for(j=0; j<connections; j++){
#if DEBUG_LEVEL>0
				cout << "-> Connecting group elements " << i+1;
				if(parameters.group_elements[group_elements_start+i].MyType) cout << " (" << parameters.group_elements[group_elements_start+i].MyType->name << ")";
				cout << " and " << connected[j];
				if(parameters.group_elements[group_elements_start+connected[j]-1].MyType) cout << " (" << parameters.group_elements[group_elements_start+connected[j]-1].MyType->name << ")";
				cout << ".\n";
#endif
				// do housekeeping first
				interaction[j].fixed=false;
				interaction[j].allow_bond_stretch=false; // for the moment, since no stretch potential has been assigned, connections are fixed in distance
				interaction[j].bond_order=connection_order[j];
				interaction[j].nr_potentials=0; // for the moment, potentials pointers have not been read
				interaction[j].potentials=NULL;
				interaction[j].partner=connected[j]-1; // which element we're connecting to
				// find out if link is static
				if(is_in_statics){
					for(k=0; k<nr_statics; k++){
						if(((i+1==static_links[k*2]) && (connected[j]==static_links[k*2+1]))
						|| ((i+1==static_links[k*2+1]) && (connected[j]==static_links[k*2]))){
							interaction[j].fixed=true;
							break;
						}
					}
				}
				// create temporary interaction on the other end as well (will get extended once we get to it)
				if(connected[j]>i+1){ // haven't gotten to element yet
					parameters.group_elements[group_elements_start+interaction[j].partner].nr_interactions++;
					parameters.group_elements[group_elements_start+interaction[j].partner].interactions = (Interaction*)realloc(parameters.group_elements[group_elements_start+interaction[j].partner].interactions,parameters.group_elements[group_elements_start+interaction[j].partner].nr_interactions*sizeof(Interaction));
					if(parameters.group_elements[group_elements_start+interaction[j].partner].interactions==NULL){
						cout << "Could not allocated memory for connection.\n";
						exit(3);
					}
					// set back link on current connection
					interaction[j].back_link=parameters.group_elements[group_elements_start+interaction[j].partner].nr_interactions-1;
					// set current element as partner in interaction on other side
					parameters.group_elements[group_elements_start+interaction[j].partner].interactions[interaction[j].back_link].partner=i;
					// set back link on partnering link
					parameters.group_elements[group_elements_start+interaction[j].partner].interactions[interaction[j].back_link].back_link=j;
					// set bond_order
					parameters.group_elements[group_elements_start+interaction[j].partner].interactions[interaction[j].back_link].bond_order=interaction[j].bond_order;
				} else{
					// need to find out if interactions to the current element on the other end has been created
					found=false;
					for(k=0; k<parameters.group_elements[group_elements_start+interaction[j].partner].nr_interactions; k++){
						if(parameters.group_elements[group_elements_start+interaction[j].partner].interactions[k].partner==i){
							found=true; // back link has been found on partner element
							interaction[j].back_link=k; // set back link on current link
							parameters.group_elements[group_elements_start+interaction[j].partner].interactions[k].back_link=j; // back link on other side (should already be set but better be safe)
							parameters.group_elements[group_elements_start+interaction[j].partner].interactions[k].fixed=interaction[j].fixed;
							break;
						}
					}
					if(!found){ // not found -- need to create it ...
#if DEBUG_LEVEL>0
						cout << "\t-> Corrected connectivity: group element " << connected[j] << " was not connected to element " << i+1 << ".\n";
#endif
						parameters.group_elements[group_elements_start+interaction[j].partner].nr_interactions++;
						parameters.group_elements[group_elements_start+interaction[j].partner].interactions = (Interaction*)realloc(parameters.group_elements[group_elements_start+interaction[j].partner].interactions,parameters.group_elements[group_elements_start+interaction[j].partner].nr_interactions*sizeof(Interaction));
						if(parameters.group_elements[group_elements_start+interaction[j].partner].interactions==NULL){
							cout << "Failed to allocated memory to correct connectivity.\n";
							exit(3);
						}
						k=parameters.group_elements[group_elements_start+interaction[j].partner].nr_interactions-1;
						// set back link on current link first
						interaction[j].back_link=k;
						// set bond order
						parameters.group_elements[group_elements_start+interaction[j].partner].interactions[k].bond_order=interaction[j].bond_order;
						// switch to other end now
						bool is_fixed=interaction[j].fixed;
						interaction = parameters.group_elements[group_elements_start+interaction[j].partner].interactions;
						
						interaction[k].fixed=is_fixed;
						interaction[k].allow_bond_stretch=false;
						interaction[k].nr_potentials=0;
						interaction[k].potentials=NULL;
						interaction[k].partner=i;
						interaction[k].back_link=j;
						
						// assign what else can be assigned
						item="connection_point."+to_string(connected[j])+"-"+to_string(i+1);
						interaction[k].initial_location = new Vec3;
						SetParam(interaction[k].initial_location->vec,item.c_str(),conf,Vec3(0.0));
						if(*interaction[k].initial_location!=Vec3(0.0)){
							interaction[k].location = new Vec3;
							*interaction[k].location=*interaction[k].initial_location;
						} else interaction[k].location=NULL;
						
						item="connection_normal."+to_string(connected[j])+"-"+to_string(i+1);
						interaction[k].initial_normal = new Vec3;
						SetParam(interaction[k].initial_normal->vec,item.c_str(),conf,Vec3(0.0));
						if(*interaction[k].initial_normal!=Vec3(0.0)){
							interaction[k].normal = new Vec3;
							*interaction[k].normal=*interaction[k].initial_normal;
						} else interaction[k].normal=NULL;
						
						item="connection_tangent."+to_string(connected[j])+"-"+to_string(i+1);
						interaction[k].initial_tangent = new Vec3;
						SetParam(interaction[k].initial_tangent->vec,item.c_str(),conf,Vec3(0.0));
						if(*interaction[k].initial_tangent!=Vec3(0.0)){
							interaction[k].tangent = new Vec3;
							*interaction[k].tangent=*interaction[k].initial_tangent;
						} else interaction[k].tangent=NULL;
						
						if(interaction[k].location!=NULL) *interaction[k].location=parameters.group_elements[group_elements_start+interaction[j].partner].rot*(*interaction[j].location);
						if(interaction[k].normal!=NULL) *interaction[k].normal=parameters.group_elements[group_elements_start+interaction[j].partner].rot*(*interaction[j].normal);
						if(interaction[k].tangent!=NULL) *interaction[k].tangent=parameters.group_elements[group_elements_start+interaction[j].partner].rot*(*interaction[j].tangent);
						
						interaction = parameters.group_elements[group_elements_start+i].interactions;
					}
				}
				
				// assign what else can be assigned
				item="connection_point."+to_string(i+1)+"-"+to_string(connected[j]);
				interaction[j].initial_location = new Vec3;
				SetParam(interaction[j].initial_location->vec,item.c_str(),conf,Vec3(0.0));
				if(*interaction[j].initial_location!=Vec3(0.0)){
					interaction[j].location = new Vec3;
					*interaction[j].location=*interaction[j].initial_location;
				} else interaction[j].location=NULL;
				
				item="connection_normal."+to_string(i+1)+"-"+to_string(connected[j]);
				interaction[j].initial_normal = new Vec3;
				SetParam(interaction[j].initial_normal->vec,item.c_str(),conf,Vec3(0.0));
				if(*interaction[j].initial_normal!=Vec3(0.0)){
					interaction[j].normal = new Vec3;
					*interaction[j].normal=*interaction[j].initial_normal;
				} else interaction[j].normal=NULL;
				
				item="connection_tangent."+to_string(i+1)+"-"+to_string(connected[j]);
				interaction[j].initial_tangent = new Vec3;
				SetParam(interaction[j].initial_tangent->vec,item.c_str(),conf,Vec3(0.0));
				if(*interaction[j].initial_tangent!=Vec3(0.0)){
					interaction[j].tangent = new Vec3;
					*interaction[j].tangent=*interaction[j].initial_tangent;
				} else interaction[j].tangent=NULL;
				
				if(interaction[j].location!=NULL) *interaction[j].location=parameters.group_elements[group_elements_start+i].rot*(*interaction[j].location);
				if(interaction[j].normal!=NULL) *interaction[j].normal=parameters.group_elements[group_elements_start+i].rot*(*interaction[j].normal);
				if(interaction[j].tangent!=NULL) *interaction[j].tangent=parameters.group_elements[group_elements_start+i].rot*(*interaction[j].tangent);
			}
			count++;
		} // Interactions have been set up
		if(static_links) delete[] static_links;
		// Set bond_range
		if(connectivity) group->Type->bond_range=0; else group->Type->bond_range=1;
		SetParam(group->Type->bond_range,"bond_range",conf,group->Type->bond_range);
		// clean up
		if(content) delete[] content;
		if(charge_override) delete[] charge_override;
		delete[] connected;
		delete[] connection_order;
		delete[] connectivity;
		delete[] bond_order;
		// potentials are non-existent yet (also only set to non-existent during first run of potential supergroup run)
		group->nr_potentials=0;
		group->potentials=NULL;
	}
	if(rerun_elements) delete[] rerun_elements;
	if(delete_elements) free(delete_elements); // dynamic in nature
	if(group_startnr) delete[] group_startnr;
	if(lodgroup_startnr) delete[] lodgroup_startnr;
	if(addgroups) delete[] addgroups; // no worries: just deleting pointers here
	
	if(finished){
		AssignFixedElements(group);
		SetParam(group->Type->rand_elements,"rand_elements",conf,(unsigned int)round(parameters.group_rand_frac*group->Type->nr_movable));
#if DEBUG_LEVEL>1
		cout << "-> Number of elements randomized for group type: " << group->Type->rand_elements << "\n";
#endif
		SetParam(group->Type->nr_rand_per_cycle,"nr_rand_per_cycle",conf,sqrt(group->Type->nr_movable));
#if DEBUG_LEVEL>1
		cout << "-> Number of elements randomized on average per cycle: " << group->Type->nr_rand_per_cycle << "\n";
#endif
		if(group->Type->nr_movable>0){
			group->Type->inv_sqrt_nrpc = sqrt(1.0/group->Type->nr_rand_per_cycle);
			if(group->Type->inv_sqrt_nrpc<EPS) group->Type->inv_sqrt_nrpc=1.0;
		} else group->Type->inv_sqrt_nrpc=1.0;
		SetParam(group->Type->ES_scale,"ES_scale",conf,1.0);
		if(fabs(group->Type->ES_scale-1.0)>EPS){
#if DEBUG_LEVEL>1
			cout << "-> Scaling group charges by " << group->Type->ES_scale << "\n";
#endif
		}
		// sanity check if overall group charge is zero
		double netcharge=0.0;
		double acsum=0.0;
		unsigned int nr_charges=0;
		group->Type->mass=0.0;
		for(i=0; i<group->nr_elements; i++){
			group->Type->mass+=parameters.group_elements[group->elements[i]].MyType->mass;
			for(j=0; j<parameters.group_elements[group->elements[i]].MyType->nr_charges; j++){
				if(fabs(group->Type->ES_scale-1.0)>EPS) parameters.group_elements[group->elements[i]].MyType->q[j]*=group->Type->ES_scale;
				netcharge+=parameters.group_elements[group->elements[i]].MyType->q[j];
				acsum+=fabs(parameters.group_elements[group->elements[i]].MyType->q[j]);
				nr_charges++;
			}
		}
		double newcharge=0.0;
		if(fabs(netcharge/e_in_esu-qround(netcharge/e_in_esu))>EPS){ // round charge to nearest integer
			cout << "Rounding charge of group <" << group->Type->name << "> to nearest integer charge of " << qround(netcharge/e_in_esu) << " e (charge currently is: " << netcharge/e_in_esu << " e)\n";
			double dq=netcharge-qround(netcharge/e_in_esu)*e_in_esu;
			cout << "-> Subtracting charges by weight from each charge.\n";
			for(i=0; i<group->nr_elements; i++){
				for(j=0; j<parameters.group_elements[group->elements[i]].MyType->nr_charges; j++){
					double aq=fabs(parameters.group_elements[group->elements[i]].MyType->q[j]);
					parameters.group_elements[group->elements[i]].MyType->q[j]-=(aq*dq)/acsum;
					newcharge+=parameters.group_elements[group->elements[i]].MyType->q[j];
				}
			}
		} else newcharge=netcharge;
		cout << "Charge of group: " << newcharge/e_in_esu << " e\n";
		Vec3 dipole=Vec3(0.0);
		Vec3 color(0.0);
		for(i=0; i<group->nr_elements; i++){
			color+=parameters.group_elements[group->elements[i]].MyType->color; // calculate average color while we're here
			dipole+=parameters.group_elements[group->elements[i]].dipole;
			for(j=0; j<parameters.group_elements[group->elements[i]].MyType->nr_charges; j++){
				dipole+=parameters.group_elements[group->elements[i]].center*parameters.group_elements[group->elements[i]].MyType->q[j];
			}
		}
		color/=group->nr_elements;
		SetParam(group->Type->group_dipole_color.vec,"group_dipole_color",conf,color);
#if DEBUG_LEVEL>1
		cout << "-> Group's dipole moment is (" << dipole.V3Str(',') << ") Debye (magnitude: " << dipole.V3Norm() << " D)\n";
		cout << "-> Group's molecular weight is " << group->Type->mass << " g/mol.\n";
#endif
		group->Type->Volume=-1;
		double default_packing=0.7; // default is 70% ellipse packing density
		if(((group->Type->density>0.0) && (group->Type->Volume>0.0)) && (group->Type->mass>0.0)) default_packing=group->Type->Volume/((group->Type->mass/group->Type->density)/NA*1E24);
		SetParam(group->Type->packing_density,"packing_density",conf,default_packing);
		if(fabs(group->Type->packing_density-default_packing)>EPS){
			if(group->Type->mass>0.0){
				if(group->Type->density>0.0){
					if(group->Type->Volume>0.0){ // problem
						cout << "Specified packing density is not consistent with specified group density, volume, or element molecular weights.\n";
						exit(4);
					} else group->Type->Volume=((group->Type->mass/group->Type->density)/NA*1E24)*group->Type->packing_density;
				} else{
					if(group->Type->Volume>0.0) group->Type->density=group->Type->Volume/(group->Type->mass/NA*1E24)/group->Type->packing_density;
				}
			}
		}
		// Now is a good time for that ...
		SetGroupRangeField(group,parameters.group_elements);
		DetermineGroupRings(group,parameters.group_elements);
		// Calculate center of group
		parameters.group_centers[nr] = GroupCenter(group,parameters.group_elements);
#if DEBUG_LEVEL>1
		cout << "-> Center (cartesian average of element centers) of group is: (" << parameters.group_centers[nr].vec[0] << ", " << parameters.group_centers[nr].vec[1] << ", " << parameters.group_centers[nr].vec[2] << ")\n";
#endif
		// Determine gyration tensor
		Vec3 points[6];
		Mat33 S;
		S.M3Zeros();
		double gcount=0.0;
		for(i=0; i<group->nr_elements; i++){
			Element* curr=&parameters.group_elements[group->elements[i]];
#if DEBUG_LEVEL>3
			cout << curr->MyType->name << "\n";
#endif
			Vec3 center=curr->center-parameters.group_centers[nr];
			for(unsigned int l=0; l<3; l++){
				points[l<<1].vec[0]=(center.vec[0]+curr->rot.mat[l][0]*curr->MyType->saxes.vec[0]);
				points[(l<<1)+1].vec[0]=(center.vec[0]-curr->rot.mat[l][0]*curr->MyType->saxes.vec[0]);
				
				points[l<<1].vec[1]=(center.vec[1]+curr->rot.mat[l][1]*curr->MyType->saxes.vec[1]);
				points[(l<<1)+1].vec[1]=(center.vec[1]-curr->rot.mat[l][1]*curr->MyType->saxes.vec[1]);
				
				points[l<<1].vec[2]=(center.vec[2]+curr->rot.mat[l][2]*curr->MyType->saxes.vec[2]);
				points[(l<<1)+1].vec[2]=(center.vec[2]-curr->rot.mat[l][2]*curr->MyType->saxes.vec[2]);
			}
			for(k=0; k<6; k++){
				for(unsigned int m=0; m<3; m++){
					for(unsigned int n=0; n<3; n++){
						// S_mn = sum_i r_m(i)*r_n(i) (corrected for center
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
		Vec3 semiaxis=eigenvalues.Re(); // lambda^2's
		// idea here is that at most one sees half the ellipsoid surface area - map to circle area (pi r^2 - not sphere which has additional factor 4)
		group->rT=sqrt(EllipsoidSurface(sqrt(semiaxis.vec[0]),sqrt(semiaxis.vec[1]),sqrt(semiaxis.vec[2]))/(4.0*pi));
#if DEBUG_LEVEL>1
		cout << "-> Group's test radius based on gyration tensor dimensions is " << group->rT << " Angström.\n";
#endif
		// Calculate external sphere
		parameters.group_radii[nr] = GroupRadius(group,parameters.group_elements,parameters.group_centers[nr]);
#if DEBUG_LEVEL>1
		cout << "-> Effective radius of bounding sphere able to contain group is: " << parameters.group_radii[nr] << " Angström\n";
#endif
		if(group->Type->LOD){
			cout << "-> Check for user-specified components:\n";
			double* components=get_flex_tupel("components",conf);
			count=0;
			if(components){
				group->Type->LOD->nr_components=0;
				while(components[count]>0){
					count+=(unsigned int)components[count]+1;
					group->Type->LOD->nr_components++;
				}
				if(group->Type->LOD->nr_components==0){
					cout << "No components specified after keyword.\n";
					exit(4);
				}
				if(group->Type->LOD->nr_components>max_nr_components) max_nr_components=group->Type->LOD->nr_components;
				group->Type->LOD->nr_elements=new unsigned int[group->Type->LOD->nr_components];
				group->Type->LOD->component_elements=new unsigned int*[group->Type->LOD->nr_components];
				group->Type->LOD->component_show_dipole=new bool[group->Type->LOD->nr_components];
				group->Type->LOD->component_color=new Vec3[group->Type->LOD->nr_components];
				group->Type->LOD->component_transparency=new double[group->Type->LOD->nr_components];
				group->Type->LOD->component_texture=new string[group->Type->LOD->nr_components];
				group->Type->LOD->component_texture_is_dipole=new bool[group->Type->LOD->nr_components];
				group->Type->LOD->element_in_component=new int[group->nr_elements];
				memset(group->Type->LOD->element_in_component,0xFF,group->nr_elements*sizeof(int)); // puts -1 in each field
				count=0;
				int current=0;
				while(components[count]>0){
					SetParam(group->Type->LOD->component_show_dipole[current],("component_show_dipole."+int2str(current+1)).c_str(),conf,false);
					SetParam(group->Type->LOD->component_color[current].vec,("component_color."+int2str(current+1)).c_str(),conf,Vec3(1.0));
					SetParam(group->Type->LOD->component_transparency[current],("component_transparency."+int2str(current+1)).c_str(),conf,-1.0);
					SetParam(group->Type->LOD->component_texture[current],("component_texture."+int2str(current+1)).c_str(),conf,"");
					SetParam(group->Type->LOD->component_texture_is_dipole[current],("component_texture_is_dipole."+int2str(current+1)).c_str(),conf,group->Type->LOD->component_show_dipole[current]);
					group->Type->LOD->nr_elements[current]=(unsigned int)components[count];
					group->Type->LOD->component_elements[current]=new unsigned int[group->Type->LOD->nr_elements[current]];
					for(i=0; i<group->Type->LOD->nr_elements[current]; i++){
						group->Type->LOD->component_elements[current][i]=(unsigned int)components[count+i+1]-1;
						int curr=(int)components[count+i+1]-1;
						if((curr<0) || (curr>(int)group->nr_elements)){
							cout << "Specified elements need to be between 1 and the total number of elements in the group.\n";
							exit(2);
						}
						group->Type->LOD->element_in_component[(unsigned int)curr]=current;
					}
					current++;
					count+=(unsigned int)components[count]+1;
				}
				delete[] components;
				// If components are specified we can also use them to specify which dipole to use for order calculation
				SetParam(group->Type->LOD->order_dipole_component,"order_dipole_component",conf,0); // default is 0 (whole group)
				for(i=0; i<group->nr_elements; i++){
					Element_Type* element_type=parameters.group_elements[lodgroup->elements[group->Type->LOD->element_in_ellipsoid[lodgroup->levelofdetail-1][i]]].MyType;
					current=group->Type->LOD->element_in_component[i];
					if(current>=0){
						element_type->color=group->Type->LOD->component_color[current];
						if(group->Type->LOD->component_texture[current]!="") element_type->texture=group->Type->LOD->component_texture[current];
						if(group->Type->LOD->component_texture_is_dipole[current]) element_type->texture_is_dipole=true; else element_type->texture_is_dipole=false;
						if(group->Type->LOD->component_transparency[current]>=0.0)
							element_type->transparency+=group->Type->LOD->component_transparency[current];
					}
				}
			} else{
				group->Type->LOD->order_dipole_component=0;
				group->Type->LOD->nr_components=0;
				group->Type->LOD->nr_elements=NULL;
				group->Type->LOD->component_elements=NULL;
				group->Type->LOD->component_show_dipole=NULL;
				group->Type->LOD->component_color=NULL;
				group->Type->LOD->component_transparency=NULL;
				group->Type->LOD->component_texture=NULL;
			}
			cout << "\t-> " << group->Type->LOD->nr_components << " component(s) defined.\n";
			lodgroup->nr_potentials=0;
			lodgroup->potentials=NULL;
			SetGroupRangeField(lodgroup,parameters.group_elements);
			DetermineGroupRings(lodgroup,parameters.group_elements);
			// Calculate center of group
			parameters.group_centers=(Vec3*)realloc(parameters.group_centers,parameters.num_groups*sizeof(Vec3));
			if(!parameters.group_centers){
				cout << "Not enough memory for new group centers.\n";
				exit(2);
			}
			parameters.group_centers[lodgroup->type] = GroupCenter(lodgroup,parameters.group_elements);
#if DEBUG_LEVEL>1
			cout << "-> Center (cartesian average of element centers) of LOD group is: (" << parameters.group_centers[lodgroup->type].vec[0] << ", " << parameters.group_centers[lodgroup->type].vec[1] << ", " << parameters.group_centers[lodgroup->type].vec[2] << ")\n";
#endif
			lodgroup->Type->Volume=-1;
			// Calculate external sphere
			parameters.group_radii=(double*)realloc(parameters.group_radii,parameters.num_groups*sizeof(double));
			if(!parameters.group_radii){
				cout << "Not enough memory for new group radii.\n";
				exit(2);
			}
			parameters.group_radii[lodgroup->type] = GroupRadius(lodgroup,parameters.group_elements,parameters.group_centers[lodgroup->type]);
#if DEBUG_LEVEL>1
			cout << "-> Effective radius of bounding sphere able to contain LOD group is: " << parameters.group_radii[lodgroup->type] << " Angström\n";
#endif
			group->Type->LOD->placement_deltas[0]=parameters.group_centers[lodgroup->type]-parameters.group_centers[nr];
#if DEBUG_LEVEL>1
			cout << "-> Placement difference to original is: (" << group->Type->LOD->placement_deltas[0].vec[0] << ", " << group->Type->LOD->placement_deltas[0].vec[1] << ", " << group->Type->LOD->placement_deltas[0].vec[2] << ")\n";
#endif
			AssignFixedElements(lodgroup);
			SetParam(lodgroup->Type->rand_elements,"rand_elements",conf,(unsigned int)round((double)group->Type->rand_elements/(double)group->Type->nr_movable*lodgroup->Type->nr_movable));
#if DEBUG_LEVEL>1
			cout << "-> Number of elements randomized for LOD group: " << lodgroup->Type->rand_elements << "\n";
#endif
			SetParam(lodgroup->Type->nr_rand_per_cycle,"nr_rand_per_cycle_lod",conf,sqrt(lodgroup->Type->nr_movable));
#if DEBUG_LEVEL>1
			cout << "-> Number of elements randomized on average per cycle for LOD group: " << lodgroup->Type->nr_rand_per_cycle << "\n";
#endif
			if(lodgroup->Type->nr_movable>0){
				lodgroup->Type->inv_sqrt_nrpc = sqrt(1.0/lodgroup->Type->nr_rand_per_cycle);
				if(lodgroup->Type->inv_sqrt_nrpc<EPS) lodgroup->Type->inv_sqrt_nrpc=1.0;
			} else lodgroup->Type->inv_sqrt_nrpc=1.0;
			dipole.V3Zeros();
			color.V3Zeros();
			for(i=0; i<lodgroup->nr_elements; i++){
				color+=parameters.group_elements[lodgroup->elements[i]].MyType->color; // calculate average color while we're here
				dipole+=parameters.group_elements[lodgroup->elements[i]].dipole;
				for(j=0; j<parameters.group_elements[lodgroup->elements[i]].MyType->nr_charges; j++){
					dipole+=parameters.group_elements[lodgroup->elements[i]].center*parameters.group_elements[lodgroup->elements[i]].MyType->q[j];
				}
			}
			color/=group->nr_elements;
			SetParam(lodgroup->Type->group_dipole_color.vec,"group_dipole_color",conf,color); // also needs to be done for super groups
#if DEBUG_LEVEL>1
			cout << "-> LOD group's dipole moment is (" << dipole.V3Str(',') << ") Debye (magnitude: " << dipole.V3Norm() << " D)\n";
#endif
		}
	}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return finished;
} // groups with respective elements are now created - interaction potentials will be created next ...

void MC_Config::LoadPotential(char* conf, int nr, string specifier)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	// First, separate specifier into group name and potential name (which may be empty)
	unsigned int i,j,count;
	bool found;
	Element_Group* group=NULL;
	parameters.potentials[nr]=new Interaction_Potential;
	Interaction_Potential* potential = parameters.potentials[nr];
	bool separate=false;
	bool whitespace=true;
	
	string groupname="";
	string potentialname="";
	for(i=0; i<specifier.length(); i++){
		if(specifier[i]==','){
			if(separate){ // bail out if there are too many commas
				cout << "Please do not use a comma in a group or interaction potential name.\n";
				exit(3);
			}
			separate=true; // separate at the comma
		} else if(separate){
			if((!whitespace) || ((specifier[i]!=' ') && (specifier[i]!='\t'))){
				whitespace=false;
				potentialname+=specifier[i];
			}
		} else groupname+=specifier[i];
	}
	// take care of whitespace after groupname
	j=groupname.length()-1;
	while((groupname[j]==' ') || (groupname[j]=='\t')){ j--; }
	j++;
	groupname.erase(j,groupname.length()-j);
	
	potential->name=new string;
	if(separate){ // interaction potential is defined between elements in group (syntax: [Interaction: <group>, <potential>])
		*potential->name=potentialname;
		cout << "Reading interaction potential: " << *potential->name << "\n";
		
		// find group the interaction potential belongs to
		found=false;
		for(i=0; i<parameters.num_groups; i++){
			group=parameters.groups[i];
			if(compare_strings(group->Type->name.c_str(),groupname.c_str())){
				found=true;
				break;
			}
		} // group points to respective group if found is true (to last group if not)
		if(!found){ // bail out if group archetype specified does not exist
			cout << "Group <" << groupname << "> not found.\n";
			exit(3);
		} else{
#if DEBUG_LEVEL>0
			cout << "Found group: " << group->Type->name << "\n";
#endif
		}
	} else{ // interactions is global between defined element types
		potentialname=groupname;
		groupname="";
		*potential->name=potentialname;
		cout << "Reading interaction potential: " << *(*parameters.potentials[nr]).name << " (global)\n";
	}
	
	// from now on: if groupname is empty interaction potential is global, group-local is not
	string type;
	SetParam(type,"type",conf);
	
	if(groupname==""){ // global interaction between element types
		// <TODO>
		cout << "-> Patience young padawan. Feature you are looking for not yet written it is.\n";
	} else{ // group-local interaction
		double* params = get_flex_tupel("parameters",conf);
		if(params==NULL){
			cout << "Interaction potential needs parameters.\n";
			exit(1);
		}
		if(params[0]>200){
			cout << "Too many parameters specified.\n";
			exit(2);
		}
		// look for recognized interaction potential types (exit if unknown is found)
		if(compare_strings(type.c_str(),"stretch")){
			potential->type=stretch_potential;
			potential->n=2;
			potential->nr_parameters=2;
		} else{
			if(compare_strings(type.c_str(),"bend")){
				potential->type=bend_potential;
				potential->n=3;
				potential->nr_parameters=2;
			} else{
				if(compare_strings(type.c_str(),"dihedral")){
					potential->type=dihedral_potential;
					potential->n=4;
					potential->nr_parameters=(unsigned short)params[0];
					if(potential->nr_parameters%2){
						cout << "Need pairs of parameters for dihedral potential.\n";
						exit(3);
					}
					if(potential->nr_parameters==0) potential->nr_parameters=2; // just to have a meaningful error message a couple lines later
				} else{
					if(compare_strings(type.c_str(),"improper_dihedral")){
						potential->type=dihedral_potential;
						potential->n=4;
						potential->nr_parameters=2;
					} else{
						if(compare_strings(type.c_str(),"spring")){
							potential->type=spring_potential;
							potential->n=2;
							potential->nr_parameters=4;
						} else{
							if(compare_strings(type.c_str(),"chainspring")){
								potential->type=chainspring_potential;
								potential->n=2;
								potential->nr_parameters=8;
							} else{
								cout << "Type <" << type << "> is not a recognized interaction potential type.\n";
								cout << "Current types are: stretch, bend, dihedral, improper_dihedral, and spring.\n";
								exit(3);
							}
						}
					}
				}
			}
		}
		
		// Load parameters
		if((int)params[0]>(int)potential->nr_parameters) cout << "Too many parameters specified. Will only use first " << potential->nr_parameters << " parameters.\n";
		potential->parameters = new double[potential->nr_parameters];
		if((int)params[0]<(int)potential->nr_parameters){
			if((potential->type==chainspring_potential) && ((int)params[0]>=4)){
				cout << "Using default values for missing parameters of the chainspring potential (";
				if((int)params[0]==4) cout << "a=0.5 Ang, l_a=1.5 Ang, m=0, V_0=0";
				if((int)params[0]==5) cout << "l_a=1.5 Ang, m=0, V_0=0";
				if((int)params[0]==6) cout << "m=0, V_0=0";
				if((int)params[0]==7) cout << "V_0=0";
				cout << ").\n";
				potential->parameters[4]=0.5;
				potential->parameters[5]=1.5;
				potential->parameters[6]=0.0;
				potential->parameters[7]=0.0;
			} else{
				cout << "Not enough parameters specified, " << potential->nr_parameters << " parameters are needed.\n";
				exit(3);
			}
		}
		for(j=0; (int)j<(int)params[0]; j++) potential->parameters[j]=params[j+1]; // first element is element count
		// Set result array to NULL, so it can be initialized by the potential (it doesn't have to for all potentials ...)
		
		double* elements = get_flex_tupel("elements",conf);
		if(elements==NULL){ // If that error occurs I'd probably switch on the improbability drive and get out ...
			cout << "The Vorgons succeeded. The mice are unhappy.\n";
			exit(42);
		}
		Element* element;
		i=0;
		count=0;
		Interaction_Potential* current = potential; // current interaction potential being worked on (to allow for creation of more)
		while((int)elements[i]>=0){ // first create additional potentials
			count++;
			if((int)elements[i]!=current->n){
				cout << "Type <" << type << "> requires " << current->n << " elements - element group number " << count << " only has " << (int)elements[i] << ".\n";
				cout << "Syntax: elements = {<element group 1>|<element group 2>|...|<element group n>} ; <element group> = 1,2,...,m\n";
				exit(3);
			}
			if(count>1){ // need to create more interaction potentials dynamically
				parameters.num_potentials++;
				parameters.potentials = (Interaction_Potential**)realloc(parameters.potentials,parameters.num_potentials*sizeof(Interaction_Potential*));
				if(parameters.potentials==NULL){ // ran out of memory
					cout << "Could not find a memory block large enough to hold a total of " << parameters.num_potentials << " group interaction potentials.\n";
					exit(1);
				}
				potential=parameters.potentials[nr];
				parameters.potentials[parameters.num_potentials-1]=new Interaction_Potential;
				current=parameters.potentials[parameters.num_potentials-1];
				current->name=potential->name;
				current->type=potential->type;
				current->n=potential->n;
				current->parameters=potential->parameters;
			}
			// assign partners
			current->partners=new unsigned int[current->n];
#if DEBUG_LEVEL>0
			cout << "Setting up potential between ";
#endif
			for(j=0; j<current->n; j++){
				if(((int)elements[i+j+1]<1) || ((int)elements[i+j+1]>(int)group->nr_elements)){
					cout << "Cannot use group element number " << (int)elements[i+j+1] << ". Group <" << group->Type->name << "> only has elements number 1-" << group->nr_elements << ".\n";
					exit(3);
				}
				current->partners[j]=(unsigned int)elements[i+j+1]-1;
#if DEBUG_LEVEL>0
				if(parameters.group_elements[group->elements[current->partners[j]]].MyType) cout << parameters.group_elements[group->elements[current->partners[j]]].MyType->name; else cout << group->Type->name << ":" << current->partners[j]+1;
				if(j+1<(unsigned int)elements[i]) cout << ", "; else cout << "\n";
#endif
			} // current->partners now contains valid list of group elements belonging to this potential
			// assign potential to the interaction between elements involved
			current->links=new unsigned int[current->n*current->n];
			for(j=0; j<current->n; j++){ // go through elements belonging to this potential ...
				// ... to find the other elements on the other side of (start at 0 so both direction interactions are properly assigned) ...
				for(unsigned int k=0; k<current->n; k++){
					if(k!=j){
						element=&parameters.group_elements[group->elements[current->partners[j]]]; // points to element j
						for(unsigned int l=0; l<element->nr_interactions; l++){ // ... the interaction between these two elements j and k
							if(element->interactions[l].partner==current->partners[k]){ // interaction between j and k found - it's the l-th interaction on element j
#if DEBUG_LEVEL>0
								cout << "-> Found interaction: ";
								if(element->MyType) cout << element->MyType->name; else cout << group->Type->name << ":" << current->partners[j]+1;
								cout << "->";
								if(parameters.group_elements[group->elements[current->partners[k]]].MyType) cout << parameters.group_elements[group->elements[current->partners[k]]].MyType->name; else cout << group->Type->name << ":" << current->partners[k]+1;
								cout << " (" << current->partners[j] << "-" << current->partners[k] << ")";
#endif
								current->links[j*current->n+k]=l;
								element->interactions[l].nr_potentials++;
#if DEBUG_LEVEL>0
								cout << " " << element->interactions[l].nr_potentials;
#endif
								element->interactions[l].potentials=(Interaction_Potential**)realloc(element->interactions[l].potentials,element->interactions[l].nr_potentials*sizeof(Interaction_Potential*));
								if(!element->interactions[l].potentials){
									cout << "\nCould not get enough memory to store interaction potential pointers for interaction.\n";
									exit(3);
								}
								// assign just created new potential pointer the current potential
								element->interactions[l].potentials[element->interactions[l].nr_potentials-1]=current;
#if DEBUG_LEVEL>0
								cout << " - " << type << " potential assigned.\n";
#endif
								break; // we found it, so no need to go through other l's, right?
							}
						}
					}
				}
			}
			// assign potential to group it belongs to
			group->nr_potentials++;
			group->potentials=(Interaction_Potential**)realloc(group->potentials,group->nr_potentials*sizeof(Interaction_Potential*));
			if(group->potentials==NULL){
				cout << "\nCould not get enough memory to store interaction potential pointers in adjoined group.\n";
				exit(3);
			}
			// assign just created new potential pointer the current potential
			group->potentials[group->nr_potentials-1]=current;
#if DEBUG_LEVEL>0
			cout << "<- Potential registered with group.\n";
#endif
			i+=(int)elements[i]+1; // next counter is nr of elements+1 away
		} // count is now the total number of interaction potentials sharing the same type and parameters
		// clean up
		delete[] elements;
	}
#if DEBUG_LEVEL>0
	cout << "<- Created interaction potential " << *potential->name;
	if(groupname=="") cout << " (global).\n"; else cout << " (group: " << group->Type->name << ").\n";
#endif
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

bool MC_Config::GetAllElementProperties(Traj_EP** elements, unsigned int* steps, double* V, double* Xm, double* time, unsigned int &nr_steps, unsigned int start_step)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	nr_steps=parameters.last_step/parameters.grfreq+(parameters.last_step%parameters.grfreq>0)+1-(start_step+(start_step%parameters.grfreq>1)*parameters.grfreq-(start_step%parameters.grfreq))/parameters.grfreq;
#if DEBUG_LEVEL>1
	cout << "-> Reading element coordinates for steps "<< start_step << "-" << parameters.last_step << "\n";
#endif
	bool success=true;
	bool warn=false;
	bool readprev=false;
	bool step_not_found=false;
	string basedon=parameters.trajectoryfile;
	readfile trajfile;
	trajfile.filename=parameters.trajectorydir+basedon;
	trajfile.file.open(trajfile.filename.c_str(),ios::binary);
	unsigned int stored=0;
	if(!elements){
		cout << "Not enough memory to store step data.\n";
		exit(2);
	}
	for(unsigned int i=0; i<parameters.n_oids; i++){
		elements[i]=new Traj_EP[nr_steps];
		if(!elements[i]){
			cout << "Not enough memory to store step data.\n";
			exit(2);
		}
	}
	cout << "\t-> Progress: 0%";
	cout.flush();
	unsigned int percentage=0;
	unsigned int stepnr=0;
	do{
		readprev=false;
		unsigned int position=0;
		unsigned int inc_pos=0;
		unsigned int first_step=0;
		unsigned last_step=0;
		char* first=GetSection(&trajfile,"General",position);
		if(first){
			SetParam(first_step,"first_step",first,0);
			SetParam(last_step,"last_step",first);
			SetParam(basedon,"configuration",first);
		} else{
			cout << "ERROR: File \"" << trajfile.filename << "\" is not a trajectory file.\n";
			success=false;
		}
		delete[] first;
		
		unsigned int step=first_step;
		if(first_step==1) step=0; // done as a safety net
		if(first_step<=start_step){
			stepnr=0;
		} else stepnr=(unsigned int)ceil(double(step-(start_step+(start_step%parameters.grfreq>1)*parameters.grfreq-(start_step%parameters.grfreq)))/double(parameters.grfreq));
		
		char* step_section=NULL;
		char* vel_section=NULL;
		string subname="";
		position=0;
		inc_pos=0;
		while((stepnr<nr_steps) && (success && (step_section=GetSection(&trajfile,"Step",position,inc_pos,subname,false,false)))){
			from_string(step,subname);
			if(parameters.time){
				vel_section=GetSection(&trajfile,"Velocity",position,inc_pos,subname,false,false);
			}
			if(step>(start_step+(start_step%parameters.grfreq>1)*parameters.grfreq-(start_step%parameters.grfreq))){ // Take care of "holes" in the trajectory file
				stepnr=(unsigned int)ceil(double(step-(start_step+(start_step%parameters.grfreq>1)*parameters.grfreq-(start_step%parameters.grfreq)))/double(parameters.grfreq));
			}
			double newV;
			double t=0.0;
			if(step>=(start_step+(start_step%parameters.grfreq>1)*parameters.grfreq-(start_step%parameters.grfreq))){
				double percent=100.0*double(stored+1)/double(nr_steps);
				if(percent>=percentage+5){
					percentage=(unsigned int)floor(percent/5.0)*5;
					if(percentage%20==0) cout << percentage << "%"; else cout << ".";
					cout.flush();
				}
				if(step_section){
					SetParam(t,"simulation_time_in_ps",step_section,0.0,true);
					if(!parameters.transition){
						SetParam(newV,"volume",step_section,parameters.targetV,true);
					} else SetParam(newV,"volume",step_section,parameters.transition_target_volume,true);
					if(parameters.LJwall_calc && parameters.LJwall_fixed){
						double newXm;
						SetParam(newXm,"LJwall_xm",step_section,parameters.LJwall_xm,true);
						double scale=sqrt((newV*parameters.boxlength[0])/(parameters.V*2.0*newXm));
						parameters.LJwall_xm=newXm;
						parameters.boxlength[0]=2.0*newXm;
						parameters.boxlength[1]*=scale; parameters.boxlength[2]*=scale;
						parameters.nndist*=scale;
						parameters.V=parameters.boxlength[0]*parameters.boxlength[1]*parameters.boxlength[2];
					} else update_volume(&parameters,newV);
					SetParam(parameters.kT,"kT",step_section,parameters.kT,true);
					parameters.beta=1.0/parameters.kT;
#ifndef USE_CMWC4096
					if(!parameters.restart_calculations) SetParam(*parameters.idum,"idum",step_section);
#endif
				} else{
					cout << "ERROR: Requested step " << step << " is not stored in specified trajectory file.\n";
					step_not_found=true;
					success=false;
				}
			}
			if(step>=(start_step+(start_step%parameters.grfreq>1)*parameters.grfreq-(start_step%parameters.grfreq))){
				double* step_positions;
				unsigned int element_nr=0;
				while((stepnr<nr_steps) && (success && (step_positions=get_flex_tupel(int2str(element_nr+1).c_str(),step_section)))){
					if((unsigned int)step_positions[0]!=2){
						cout << "ERROR: First field does not consist of two numbers.\n";
						success=false;
					} else{
						if((int)step_positions[3]!=3){
							cout << "ERROR: Position 3-vector expected as second field.\n";
							success=false;
						} else{
							if((int)step_positions[7]!=4){
								cout << "ERROR: Rotation axis and angle 4-vector expected as third field\n";
								success=false;
							} else{
								elements[element_nr][stepnr].element_type=(unsigned int)step_positions[1];
								elements[element_nr][stepnr].group_type=(int)step_positions[2];
								elements[element_nr][stepnr].position.vec[0]=step_positions[4];
								elements[element_nr][stepnr].position.vec[1]=step_positions[5];
								elements[element_nr][stepnr].position.vec[2]=step_positions[6];
								elements[element_nr][stepnr].rotation_vector.vec[0]=step_positions[8];
								elements[element_nr][stepnr].rotation_vector.vec[1]=step_positions[9];
								elements[element_nr][stepnr].rotation_vector.vec[2]=step_positions[10];
								elements[element_nr][stepnr].rotation_vector.vec[3]=step_positions[11];
							}
						}
					}
					element_nr++;
					if(element_nr>parameters.n_oids){
						cout << "ERROR: Trajectory file contains more elements than specified in configuration file.\n";
						success=false;
					}
					delete[] step_positions;
				}
				if(vel_section){
					double* vel_velocities;
					element_nr=0;
					while((stepnr<nr_steps) && (success && (vel_velocities=get_flex_tupel(int2str(element_nr+1).c_str(),vel_section)))){
						if((unsigned int)vel_velocities[0]!=2){
							cout << "ERROR: First field does not consist of two numbers.\n";
							success=false;
						} else{
							if((int)vel_velocities[3]!=1){
								cout << "ERROR: Speed scalar expected as second field.\n";
								success=false;
							} else{
								if((int)vel_velocities[5]!=3){
									cout << "ERROR: Velocity 3-vector expected as third field\n";
									success=false;
								} else{
									elements[element_nr][stepnr].speed=vel_velocities[4];
									elements[element_nr][stepnr].velocity.vec[0]=vel_velocities[6];
									elements[element_nr][stepnr].velocity.vec[1]=vel_velocities[7];
									elements[element_nr][stepnr].velocity.vec[2]=vel_velocities[8];
								}
							}
						}
						element_nr++;
						if(element_nr>parameters.n_oids){
							cout << "ERROR: Trajectory file contains more elements than specified in configuration file.\n";
							success=false;
						}
						delete[] vel_velocities;
					}
				}
				if((element_nr<parameters.n_oids) && (!step_not_found)){
					cout << "ERROR: Trajectory file contains less elements than specified in configuration file.\n";
					success=false;
				}
				
				time[stepnr]=t;
				V[stepnr]=parameters.V;
				if(Xm) Xm[stepnr]=parameters.LJwall_xm;
				steps[stepnr]=step;
				step+=parameters.grfreq;
				stored++;
				if(stepnr>=nr_steps){ // what if there are too many steps in trajectory file (which shouldn't happen)
					cout << "WARNING: Too many steps, only read " << nr_steps << " steps in accordance with configuration file.\n";
					if(success) warn=true;
					success=false;
				}
				stepnr++;
			}
			delete[] step_section;
			if(vel_section) delete[] vel_section;
			subname="";
		}
		if(warn){
			success=true; // was just a warning
			warn=false;
		}
		if(success && ((first_step>start_step) && (first_step>1))){ // need to read previous trajectory file
			trajfile.file.close();
			trajfile.filename=parameters.trajectorydir+basedon;
			GetDirectory(trajfile);
			trajfile.file.open((trajfile.directory+trajfile.filename).c_str(),ios::binary);
			if(trajfile.file.fail()){
				cout << "ERROR: Cannot open connected trajectory file <" << basedon << ">.\n";
				success=false;
			}
			readprev=true;
		}
	} while(readprev && success);
	
	trajfile.file.close();
	if(stepnr>0) nr_steps=stepnr;
	
	
#if DEBUG_LEVEL>1
	if(success) cout << "\n<- Done.\n";
#endif
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return success;
}

bool MC_Config::GetElementProperties(const unsigned step, Traj_EP* elements)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
#if DEBUG_LEVEL>1
	cout << "Reading element coordinates for step " << step << ".\n";
#endif
	bool success=true;
	bool readprev=false;
	bool step_not_found=false;
	string basedon=parameters.trajectoryfile;
	readfile trajfile;
	trajfile.filename=parameters.trajectorydir+basedon;
	trajfile.file.open(trajfile.filename.c_str(),ios::binary);
	unsigned int element_nr=0;
	do{
		unsigned int position=0;
		unsigned int inc_pos=0;
		unsigned int first_step=0;
		unsigned last_step=0;
		char* first=GetSection(&trajfile,"General",position);
		if(first){
			SetParam(first_step,"first_step",first,0);
			SetParam(last_step,"last_step",first);
			SetParam(basedon,"configuration",first);
		} else{
			cout << "ERROR: Requested step " << step << " is not in current file. Trajectory file is expected as parameter <configuration>.\n";
			success=false;
		}
		delete[] first;
		
		if(success && (((first_step<=step) && (last_step>=step)) || ((first_step==1) && (step==0)))){
			position=0;
			inc_pos=0;
			string subname=int2str(step);
			char* step_section=NULL;
			char* vel_section=NULL;
			step_section=GetSection(&trajfile,"Step",position,inc_pos,subname,false,false);
			double newV;
			if(step_section){
				SetParam(newV,"volume",step_section,parameters.targetV,true);
				if(parameters.restart_volume){
					parameters.targetV=newV;
				} else{
					if(parameters.LJwall_calc && parameters.LJwall_fixed){
						double newXm;
						SetParam(newXm,"LJwall_xm",step_section,parameters.LJwall_xm,true);
						double scale=sqrt((newV*parameters.boxlength[0])/(parameters.V*2.0*newXm));
						parameters.LJwall_xm=newXm;
						parameters.boxlength[0]=2.0*newXm;
						parameters.boxlength[1]*=scale; parameters.boxlength[2]*=scale;
						parameters.nndist*=scale;
						parameters.V=parameters.boxlength[0]*parameters.boxlength[1]*parameters.boxlength[2];
					} else update_volume(&parameters,newV);
				}
				SetParam(parameters.kT,"kT",step_section,parameters.kT,true);
				parameters.beta=1.0/parameters.kT;
#ifndef USE_CMWC4096
				if(!parameters.restart_calculations) SetParam(*parameters.idum,"idum",step_section);
#endif
			} else{
				cout << "ERROR: Requested step " << step << " is not stored in specified trajectory file.\n";
				step_not_found=true;
				success=false;
			}
			if(parameters.time){
				vel_section=GetSection(&trajfile,"Velocity",position,inc_pos,subname,false,false);
			}
			double* step_positions;
			while(success && (step_positions=get_flex_tupel(int2str(element_nr+1).c_str(),step_section))){
				if((unsigned int)step_positions[0]!=2){
					cout << "ERROR: First field does not consist of two numbers.\n";
					success=false;
				} else{
					if((int)step_positions[3]!=3){
						cout << "ERROR: Position 3-vector expected as second field.\n";
						success=false;
					} else{
						if((int)step_positions[7]!=4){
							cout << "ERROR: Rotation axis and angle 4-vector expected as third field\n";
							success=false;
						} else{
							elements[element_nr].element_type=(unsigned int)step_positions[1];
							elements[element_nr].group_type=(int)step_positions[2];
							elements[element_nr].position.vec[0]=step_positions[4];
							elements[element_nr].position.vec[1]=step_positions[5];
							elements[element_nr].position.vec[2]=step_positions[6];
							elements[element_nr].rotation_vector.vec[0]=step_positions[8];
							elements[element_nr].rotation_vector.vec[1]=step_positions[9];
							elements[element_nr].rotation_vector.vec[2]=step_positions[10];
							elements[element_nr].rotation_vector.vec[3]=step_positions[11];
						}
					}
				}
				element_nr++;
				if(element_nr>parameters.n_oids){
					cout << "ERROR: Trajectory file contains more elements than specified in configuration file.\n";
					success=false;
				}
				delete[] step_positions;
			}
			if(vel_section){
				double* vel_velocities;
				element_nr=0;
				while(success && (vel_velocities=get_flex_tupel(int2str(element_nr+1).c_str(),vel_section))){
					if((unsigned int)vel_velocities[0]!=2){
						cout << "ERROR: First field does not consist of two numbers.\n";
						success=false;
					} else{
						if((int)vel_velocities[3]!=1){
							cout << "ERROR: Speed scalar expected as second field.\n";
							success=false;
						} else{
							if((int)vel_velocities[5]!=3){
								cout << "ERROR: Velocity 3-vector expected as third field\n";
								success=false;
							} else{
								elements[element_nr].speed=vel_velocities[4];
								elements[element_nr].velocity.vec[0]=vel_velocities[6];
								elements[element_nr].velocity.vec[1]=vel_velocities[7];
								elements[element_nr].velocity.vec[2]=vel_velocities[8];
							}
						}
					}
					element_nr++;
					if(element_nr>parameters.n_oids){
						cout << "ERROR: Trajectory file contains more elements than specified in configuration file.\n";
						success=false;
					}
					delete[] vel_velocities;
				}
			}
			delete[] step_section;
			if(vel_section) delete[] vel_section;
			readprev=false;
		} else{
			if(success){
				if(last_step>=step){ // need to read previous trajectory file
					trajfile.file.close();
					trajfile.filename=parameters.trajectorydir+basedon;
					trajfile.file.open(trajfile.filename.c_str(),ios::binary);
					if(trajfile.file.fail()){
						cout << "ERROR: Cannot open connected trajectory file <" << basedon << ">.\n";
						success=false;
					}
					readprev=true;
				} else{
					cout << "ERROR: Specified step " << step << " is not stored in trajectory file <" << basedon << ">.\n";
					success=false;
				}
			}
		}
	} while(readprev && success);
	
	if((element_nr<parameters.n_oids) && (!step_not_found)){
		cout << "ERROR: Trajectory file contains less elements than specified in configuration file.\n";
		success=false;
	}
	
	trajfile.file.close();
#if DEBUG_LEVEL>1
	if(success) cout << "<- Done.\n";
#endif
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return success;
}

bool MC_Config::GetStatistics(const unsigned int initialstep, const unsigned int uptostep, Vstore* potentials, Vec3* dipoles, Vec3* cosmeans, double* Vs, double* msmoved, unsigned int* accepted, unsigned int* tries)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
#if DEBUG_LEVEL>1
	if(uptostep-initialstep>1){
		if(initialstep>0){
			cout << "Reading statistics between steps " << initialstep << "-" << uptostep << ".\n";
		} else cout << "Reading statistics upto step " << uptostep << ".\n";
	} else cout << "Reading statistics for step " << initialstep << ".\n";
#endif
	bool success=true;
	bool readprev=false;
	bool progress_first=false;
	unsigned int prev_first=uptostep;
	string basedon=parameters.trajectoryfile;
	readfile trajfile;
	trajfile.filename=parameters.trajectorydir+basedon;
	trajfile.file.open(trajfile.filename.c_str(),ios::binary);
	if(uptostep<1){
		cout << "WARNING: Statistics does not exist for initial step.\n-> Statistics not used.\n";
		success=false;
		return success;
	}
	bool volume=true;
	bool stepsizes=true;
	do{
		unsigned int include_pos=0;
		unsigned int position=0;
		unsigned int first_step=0;
		char* first=GetSection(&trajfile,"General",position);
		if(first){
			SetParam(first_step,"first_step",first,0);
#if DEBUG_LEVEL>1
			if(success){
				if(first_step<prev_first) cout << "-> Found steps " << first_step << " to " << prev_first << " in trajectory file <" << basedon << ">.\n"; else cout << "-> Trajectory file <" << basedon << "> does not contain relevant steps.\n";
			}
#endif
			SetParam(basedon,"configuration",first);
			delete[] first;
		} else{
			cout << "ERROR: To read previous calculations a trajectory file is expected as parameter <configuration>.\n";
			success=false;
		}
		
		bool show_progress=true;
		char* stat=NULL;
		string nrstring="";
		int stepnr=first_step-1;
		
		if(initialstep>0){
			if((first_step<initialstep) && (prev_first>initialstep)){
				if(parameters.chkcoord){
					cout << "\t-> Jumping to statistics for step " << initialstep << ".\n";
					nrstring=int2str((initialstep/parameters.grfreq-(initialstep%parameters.grfreq==0))*parameters.grfreq);
					stat=GetSection(&trajfile,"Statistics",position,nrstring);
					if(!stat){
						cout << "ERROR: Could not find statistics for step " << initialstep << ".\n";
						exit(3);
					}
					delete[] stat;
					nrstring="";
				}
				stepnr=initialstep;
			}
		}
		if(!((stepnr<(int)uptostep) && (uptostep-initialstep>1))) show_progress=false;
		cout.flush();
		unsigned int percentage=0;
		while((stat=GetSection(&trajfile,"Statistics",position,include_pos,nrstring,false,false)) && (success && (stepnr<(int)uptostep))){
			if(nrstring==int2str(uptostep)){ // read msmoved and accepted arrays
				double* moved=get_flex_tupel("msmoved",stat);
				if(moved){
					if((int)moved[0]!=2){
						cout << "WARNING: <msmoved> is expected to be a 2-vector.\n-> Use (0.0,0.0) as default.\n";
						msmoved[0]=0.0;
						msmoved[1]=0.0;
					} else{
						msmoved[0]=moved[1];
						msmoved[1]=moved[2];
					}
					delete[] moved;
				} else{
					cout << "WARNING: <msmoved> not found.\n-> Use (0.0,0.0) as default.\n";
					msmoved[0]=0.0;
					msmoved[1]=0.0;
				}
				double* accept=get_flex_tupel("acceptance",stat);
				if(accept){
					if((int)accept[0]!=3){
						cout << "WARNING: <acceptance> expected to be 3-vector.\n-> Use (0,0,0) as default.\n";
						accepted[0]=0;
						accepted[1]=0;
						accepted[2]=0;
					} else{
						accepted[0]=(int)accept[1];
						accepted[1]=(int)accept[2];
						accepted[2]=(int)accept[3];
					}
					delete[] accept;
				} else{
					cout << "WARNING: <acceptance> not found.\n-> Use (0,0,0) as default.\n";
					accepted[0]=0;
					accepted[1]=0;
					accepted[2]=0;
				}
				double* nr_tries=get_flex_tupel("nr_tries",stat);
				if(nr_tries){
					if((int)nr_tries[0]!=3){
						cout << "WARNING: <nr_tries> expected to be 3-vector.\n-> Use (0,0,0) as default.\n";
						tries[0]=0;
						tries[1]=0;
						tries[2]=0;
					} else{
						tries[0]=(unsigned int)nr_tries[1];
						tries[1]=(unsigned int)nr_tries[2];
						tries[2]=(unsigned int)nr_tries[3];
					}
					delete[] nr_tries;
				} else{
					cout << "WARNING: <nr_tries> not found.\n-> Use (0,0,0) as default.\n";
					tries[0]=0;
					tries[1]=0;
					tries[2]=0;
				}
			}
			double* step_stat;
			while((stepnr<(int)uptostep) && (success && (step_stat=get_flex_tupel(int2str(stepnr+1).c_str(),stat)))){
				double percent=100.0*double(stepnr-first_step+1)/double(prev_first-first_step);
				if(initialstep>0){
					if((first_step<initialstep) && (prev_first>initialstep)) percent=100.0*double(stepnr-initialstep+1)/double(prev_first-initialstep);
				}
				if(show_progress && (percent>=percentage+5)){
					percentage=(unsigned int)floor(percent/5.0)*5;
					if(percentage%20==0) cout << percentage << "%"; else cout << ".";
					cout.flush();
				}
				if((unsigned int)step_stat[0]!=NUM_V_STORE){
					cout << "WARNING, STEP " << stepnr+1 << ": Trajectory file records " << (unsigned int)step_stat[0] << " potential energies. However, this version uses " << NUM_V_STORE << ".\n-> Statistics not used.\n";
					success=false;
				} else{
					if((int)step_stat[NUM_V_STORE+1]!=3){
						cout << "WARNING, STEP " << stepnr+1 << ": Dipole 3-vector expected as second field.\n-> Statistics not used.\n";
						success=false;
					} else{
						if((int)step_stat[NUM_V_STORE+5]<3){
							cout << "WARNING, STEP " << stepnr+1 << ": Need at least first three cos^n moments.\n-> Statistics not used.\n";
							success=false;
						} else{
							if((int)step_stat[NUM_V_STORE+9]<1){
								if(volume) cout << "\t-> No volume for each step, using configuration file volume " << parameters.targetV << " Angström³.\n";
								if(stepsizes) cout << "\t   No stepsizes in trajectory either, using configuration file values: " << parameters.maxtrans << " Angström, " << parameters.maxrot << " radians.\n";
								volume=false;
								stepsizes=false;
								
							} else{ // test if stepsize averages exist (was introduced after volume, so can only be there if volume exists
								if((int)step_stat[NUM_V_STORE+11]<1){
									if(stepsizes) cout << "\t-> No stepsizes in trajectory, using configuration file values: " << parameters.maxtrans << " Angström, " << parameters.maxrot/pi << "*pi radians.\n";
									stepsizes=false;
								}
							}
							for(unsigned int j=0; j<NUM_V_STORE; j++) potentials[stepnr-initialstep][j]=step_stat[j+1];
							dipoles[stepnr-initialstep].vec[0]=step_stat[NUM_V_STORE+2];
							dipoles[stepnr-initialstep].vec[1]=step_stat[NUM_V_STORE+3];
							dipoles[stepnr-initialstep].vec[2]=step_stat[NUM_V_STORE+4];
							cosmeans[stepnr-initialstep].vec[0]=step_stat[NUM_V_STORE+6];
							cosmeans[stepnr-initialstep].vec[1]=step_stat[NUM_V_STORE+7];
							cosmeans[stepnr-initialstep].vec[2]=step_stat[NUM_V_STORE+8];
							if(volume) Vs[stepnr-initialstep]=step_stat[NUM_V_STORE+10]; else Vs[stepnr-initialstep]=parameters.targetV;
							if(stepsizes){
							} else{
							}
						}
					}
				}
				if(show_progress && !progress_first){
					cout << "\t-> Progress: 0%";
					progress_first=true;
				}
				stepnr++;
				delete[] step_stat;
			}
			// clean up
			delete[] stat;
			nrstring="";
		}
		if(initialstep<first_step-1){ // need to read previous trajectory file
			trajfile.file.close();
#if DEBUG_LEVEL>1
			prev_first=first_step;
#endif
			trajfile.filename=parameters.trajectorydir+basedon;
			trajfile.file.open(trajfile.filename.c_str(),ios::binary);
			if(trajfile.file.fail()){
				cout << "WARNING: Could not open connected trajectory file <" << basedon << ">.\n-> Statistics not used.\n";
				success=false;
			}
			readprev=true;
		} else readprev=false;
		if(show_progress) cout << "\n";
	} while(readprev && success);
#if DEBUG_LEVEL>1
	cout << "<- Done.\n";
#endif
	trajfile.file.close();
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return success;
}

void MC_Config::GetFromFile(const char* filename)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	unsigned int size, position;
	bool was_trajectory=false;
	bool is_trajectory=false;
	string conf_filename;
	string adj_fn=filename;
	readfile conffile;
	conffile.filename=filename;
	conffile.directory="";
	GetDirectory(conffile);
	parameters.resume_step=0; // default value is start from step 0
	parameters.trajectorynr=0;
	
	do{
		parameters.trajectorynr++;
		if(is_trajectory){
			conffile.filename=conf_filename;
			GetDirectory(conffile);
			conffile.file.open((conffile.directory+conffile.filename).c_str(),ios::binary);
		} else{
			conffile.filename=filename;
			conffile.directory="";
			GetDirectory(conffile);
			conffile.file.open((conffile.directory+conffile.filename).c_str(),ios::binary);
		}
		if(conffile.file.fail()){
			// Try .conf extension (could be that the user didn't specify an extension)
			conffile.file.open((conffile.directory+conffile.filename+".conf").c_str(),ios::binary);
			if(conffile.file.fail()){
				// Try .traj extension (could be that the user didn't specify an extension)
				conffile.file.open((conffile.directory+conffile.filename+".traj").c_str(),ios::binary);
				if(conffile.file.fail()){
					cout << "Could not open file " << conffile.directory+conffile.filename << ".\n";
					exit(1);
				} else{
					conffile.filename+=".traj";
					if(parameters.trajectorynr==1) adj_fn+=".traj";
				}
			} else{
				conffile.filename+=".conf";
				if(parameters.trajectorynr==1) adj_fn+=".conf";
			}
		}
		cout << "-> Opened file: " << conffile.directory+conffile.filename << "\n";
		// Get file size
		conffile.file.seekg(0,ifstream::end);
		size=conffile.file.tellg();
		conffile.file.seekg(0);
		// Sanity checks
		if(size==0){
			cout << "File has no content.\n";
			exit(1);
		}
		
		// Test for trajectory file
		position=0;
		char* first=GetSection(&conffile,"General",position);
		conf_filename="";
		if(first) SetParam(conf_filename,"configuration",first,"");
		was_trajectory=is_trajectory;
		is_trajectory=(conf_filename!="");
		
		if((size>MAXCONFSIZE) && (!is_trajectory)){ // config files should not be bigger than MAXCONFSIZE/(1024*1024) MB (only exception: trajectory files)
			cout << "Configuration file is too big (>" << MAXCONFSIZE/(1024*1024) << " MB).\n";
			exit(1);
		}
		if(is_trajectory){
			conffile.file.close();
			if(!was_trajectory) cout << "\nSpecified file " << conffile.directory+conffile.filename << " is a trajectory file.\n-> Looking for original configuration file.\n";
		} else if(was_trajectory) cout << "<- Found: " << conffile.filename << "\n";
		if(first && (is_trajectory && !was_trajectory)){
			SetParam(parameters.restart_calculations,"restart_calculations",first,false);
			SetParam(parameters.restart_volume,"restart_volume",first,false);
			SetParam(parameters.resume_step,"resume_step",first,-1);
			SetParam(parameters.placement_rngseed,"rngseed",first);
			SetParam(parameters.last_step,"last_step",first);
			if(parameters.resume_step==-1) parameters.resume_step=parameters.last_step;
			if(parameters.resume_step>(int)parameters.last_step){
				cout << "Cannot read step " << parameters.resume_step << ". Trajectory file only contains data upto step " << parameters.last_step << ".\n";
				exit(1);
			}
		}
		if(first) delete[] first;
	} while(is_trajectory);
	
	cout << "\nReading configuration file " << conffile.directory+conffile.filename << " (" << size << " bytes)\n";
	
	parameters.use_trajectory=was_trajectory;
	if(parameters.use_trajectory){
		parameters.trajectoryfile=adj_fn;
		parameters.trajectorydir=GetDirectory(adj_fn);
	}
	
	parameters.configfile=conffile.filename;
	parameters.configdir=conffile.directory;
	
	// First, get simulation parameters ...
	position=0;
	constant_block=GetSection(&conffile,"Simulation Parameters",position);
	if(!constant_block){
		cout << "Section [Simulation Parameters] does not exist in configuration file.\n";
		exit(1);
	}
	LoadSimParams(constant_block);
	
	// ... then get element types
	position=0;
	parameters.num_element_types=GetNrSections(&conffile,"Element",position);
	
	// Allocate memory
	parameters.element_types = (Element_Type**)malloc(parameters.num_element_types*sizeof(Element_Type*));
	if(parameters.element_types==NULL){
		cout << "Cough, cough. Not enough memory to create element types.\n";
		exit(1);
	}
	
	char* elementtype;
	parameters.offctrmu=false;
	for(unsigned int i=0; i<parameters.num_element_types; i++){
		parameters.element_types[i]=new Element_Type;
		elementtype=GetSection(&conffile,"Element",position,parameters.element_types[i]->name); // subname is element name
		// need to check if element has already been found (this section will also be included then)
		bool found=true;
		while(found){
			found=false;
			for(unsigned int j=0; j<i; j++){
				if(compare_strings(parameters.element_types[i]->name.c_str(),parameters.element_types[j]->name.c_str())){
					found=true;
					break;
				}
			}
			if(found){
				delete[] elementtype;
				parameters.element_types[i]->name="";
				elementtype=GetSection(&conffile,"Element",position,parameters.element_types[i]->name); // subname is element name
			}
		}
		cout << "Reading element type: ";
		if(parameters.element_types[i]->name=="") parameters.element_types[i]->name="<"+arbitrary_element_names[i%arbitrary_element_names_nr]+">"; // empty names are not good ...
		cout << parameters.element_types[i]->name << "\n";
		LoadElementTypes(elementtype,i);
		delete[] elementtype; // delete when done with it (a new one will be created by next GetSection ...)
		cout << "<- Done.\n";
	}
	
	// Get element groups
	position=0;
	parameters.num_groups=GetNrSections(&conffile,"Group",position);
	
	// Allocate memory
	parameters.groups = (Element_Group**)malloc(parameters.num_groups*sizeof(Element_Group*));
	if(parameters.groups==NULL){
		cout << "Cough, cough. Not enough memory to create group archetypes.\n";
		exit(1);
	}
	parameters.group_centers = (Vec3*)malloc(parameters.num_groups*sizeof(Vec3));;
	if(parameters.group_centers==NULL){
		cout << "Not enough memory to create group archetype center vectors.\n";
		exit(1);
	}
	parameters.group_radii = (double*)malloc(parameters.num_groups*sizeof(double));
	if(parameters.group_centers==NULL){
		cout << "Not enough memory to create group archetype bounding sphere radii.\n";
		exit(1);
	}
	
	char* group;
	char** rerun_groups=NULL;
	unsigned int nr_reruns=0;
	unsigned int* reruns=NULL;
	bool success;
	for(unsigned int i=0; i<parameters.num_groups; i++){
		parameters.groups[i] = new Element_Group;
		parameters.groups[i]->type=i;
		parameters.groups[i]->Type = new Element_Group_Type;
		group=GetSection(&conffile,"Group",position,parameters.groups[i]->Type->name); // subname is element name
		// need to check if group has already been found (this section will also be included then)
		bool found=true;
		while(found){
			found=false;
			for(unsigned int j=0; j<i; j++){
				if(compare_strings(parameters.groups[i]->Type->name.c_str(),parameters.groups[j]->Type->name.c_str())){
					found=true;
					break;
				}
			}
			if(found){
				delete[] group;
				parameters.groups[i]->Type->name="";
				group=GetSection(&conffile,"Group",position,parameters.groups[i]->Type->name); // subname is element name
			}
		}
		cout << "Reading group: ";
		if(parameters.groups[i]->Type->name=="") parameters.groups[i]->Type->name="<"+arbitrary_group_names[i%arbitrary_group_names_nr]+">"; // empty names are not good ...
		cout << parameters.groups[i]->Type->name << "\n";
		success=LoadGroup(group,i,false);
		if(success){
			delete[] group; // delete text when done with it (a new one will be created by next GetSection ...)
			cout << "<- Created group " << parameters.groups[i]->Type->name << ".\n";
		} else{
			nr_reruns++;
			reruns=(unsigned int*)realloc(reruns,sizeof(unsigned int)*nr_reruns);
			if(!reruns){
				cout << "ERROR: Not enough memory to create rerun table.\n";
				exit(2);
			}
			rerun_groups=(char**)realloc(rerun_groups,sizeof(char*)*nr_reruns);
			if(!rerun_groups){
				cout << "ERROR: Not enough memory to store rerun table.\n";
				exit(2);
			}
			reruns[nr_reruns-1]=i;
			rerun_groups[nr_reruns-1]=group;
			cout << "<- Group " << parameters.groups[i]->Type->name << " will be created later.\n";
		}
	}
	
	// Get interaction potentials
	position=0;
	parameters.num_potentials=GetNrSections(&conffile,"Interaction",position);
	
	// Allocate initial interaction potential memory (more can be appended later)
	parameters.potentials = (Interaction_Potential**)malloc(parameters.num_potentials*sizeof(Interaction_Potential*));
	if(parameters.potentials==NULL){
		cout << "Cough, cough. Not enough memory to create interaction potentials.\n";
		exit(1);
	}
	
	char* potential;
	string specifier;
	unsigned int nr_potentials=parameters.num_potentials; // parameters.num_potentials can be extended by LoadPotential, need to use current state as it reflects the config file
	for(unsigned int i=0; i<nr_potentials; i++){
		specifier="";
		potential=GetSection(&conffile,"Interaction",position,specifier); // subname is <group name, potential name> *or* <potential name> -- treated in LoadPotential
		// need to check if potential has already been found (this section will also be included then)
		bool found=true;
		while(found){
			found=false;
			for(unsigned int j=0; j<i; j++){
				if(compare_strings(specifier.c_str(),parameters.potentials[j]->name->c_str())){
					found=true;
					break;
				}
			}
			if(found){
				delete[] potential;
				specifier="";
				potential=GetSection(&conffile,"Interaction",position,specifier); // subname is <group name, potential name> *or* <potential name> -- treated in LoadPotential
			}
		}
		LoadPotential(potential,i,specifier);
		delete[] potential; // delete when done with it (a new one will be created by next GetSection ...)
	}
	
	// Get Level of Detail sections
	position=0;
	parameters.num_levelofdetail=GetNrSections(&conffile,"Detail",position);
	// Allocate memory
	parameters.lods = (Level_of_Detail**)malloc(parameters.num_levelofdetail*sizeof(Level_of_Detail*));
	if(!parameters.lods){
		cout << "Cough, cough. Not enough memory to create level of detail definitions.\n";
		exit(1);
	}
	
	char* levelofdetail;
	Element_Group* current_group=NULL;
	for(unsigned int i=0; i<parameters.num_levelofdetail; i++){
		string groupname="";
		levelofdetail=GetSection(&conffile,"Detail",position,groupname); // subname is group name
		bool found=false;
		unsigned int j;
		// Find out if group exists
		for(j=0; j<parameters.num_groups; j++){
			current_group=parameters.groups[j];
			if(compare_strings(groupname.c_str(),current_group->Type->name.c_str())){
				found=true;
				break;
			}
		}
		cout << "Reading details for group: " << groupname << "\n";
		if(!found){
			cout << "ERROR: Group does not exist.\n";
			exit(1);
		}
		for(unsigned int k=0; k<nr_reruns; k++){
			if(reruns[k]==j){
				cout << "ERROR: Super groups (groups of groups) are (currently) not allowed to have user-specified LOD groups.\n";
				exit(3);
			}
		}
		parameters.lods[i]=new Level_of_Detail;
		LoadLevelofDetail(levelofdetail,current_group,i);
		delete[] levelofdetail;
		cout << "<- Done.\n";
	}
	
	// Take care of super groups aka reruns
	for(unsigned int i=0; i<nr_reruns; i++){
		cout << "Creating super group: ";
		cout << parameters.groups[reruns[i]]->Type->name << "\n";
		success=LoadGroup(rerun_groups[i],reruns[i],true);
		delete[] rerun_groups[i];
	}
	
	if(reruns) free(reruns);
	if(rerun_groups) free(rerun_groups);
	
	// automatically assign potentials for groups requiring it
	for(unsigned int i=0; i<parameters.num_groups; i++){
		current_group=parameters.groups[i];
		AutoPotentials(current_group);
	}
	
	// Sum up how things and assign residual charges
	parameters.n_oids=0;
	unsigned int nr_moving=0;
	parameters.ions_present=false;
	for(unsigned int i=0; i<parameters.num_element_types; i++){
		parameters.n_oids+=parameters.element_types[i]->number;
		nr_moving+=parameters.element_types[i]->number;
		parameters.element_types[i]->charge=0.0;
		for(unsigned int j=0; j<parameters.element_types[i]->nr_charges; j++) parameters.element_types[i]->charge+=parameters.element_types[i]->q[j];
		if((parameters.element_types[i]->number>0) && (fabs(parameters.element_types[i]->charge)>EPS)) parameters.ions_present=true;
	}
	parameters.n_groups=0;
	parameters.n_group_oids=0;
	for(unsigned int i=0; i<parameters.num_groups; i++){
		parameters.n_oids+=parameters.groups[i]->number*parameters.groups[i]->nr_elements;
		nr_moving+=(unsigned int)(parameters.groups[i]->number*parameters.groups[i]->Type->nr_rand_per_cycle);
		if((unsigned int)(parameters.groups[i]->number*parameters.groups[i]->Type->nr_rand_per_cycle)==0) nr_moving+=parameters.groups[i]->number;
		parameters.n_group_oids+=parameters.groups[i]->number*parameters.groups[i]->nr_elements;
		parameters.n_groups+=parameters.groups[i]->number;
		parameters.groups[i]->Type->charge=0.0;
		parameters.groups[i]->Type->mincharge=0.0;
		parameters.groups[i]->Type->maxcharge=0.0;
		for(unsigned int j=0; j<parameters.groups[i]->nr_elements; j++){
			double eq=0.0;
			for(unsigned int k=0; k<parameters.group_elements[parameters.groups[i]->elements[j]].MyType->nr_charges; k++){
				parameters.groups[i]->Type->charge+=parameters.group_elements[parameters.groups[i]->elements[j]].MyType->q[k];
				eq+=parameters.group_elements[parameters.groups[i]->elements[j]].MyType->q[k];
			}
			if(eq<parameters.groups[i]->Type->mincharge) parameters.groups[i]->Type->mincharge=eq;
			if(eq>parameters.groups[i]->Type->maxcharge) parameters.groups[i]->Type->maxcharge=eq;
		}
		if((parameters.groups[i]->number>0) && (fabs(parameters.groups[i]->Type->charge)>EPS)) parameters.ions_present=true;
	}
	// Find out if the user wants to have a neutral RF sphere
	cout << "-> Are ions present in simulation volume?\n";
	if(!parameters.ions_present){
		cout << "\t-> No ions are present in simulation box.\n";
	} else{
		cout << "\t-> Ions are present.\n";
		if(parameters.dyneps){
//			cout << "\t\t-> Switch group_dipole on for ionic groups to allow proper epsRF calculation:\n";
			cout << "\t\t-> Switch group_dipole off for groups to allow proper epsRF calculation:\n";
			for(unsigned int i=0; i<parameters.num_groups; i++){
/*				if((parameters.groups[i]->number>0) && (fabs(parameters.groups[i]->Type->charge)>EPS)){
					cout << "\t\t\t-> " << parameters.groups[i]->Type->name << "\n";
					parameters.groups[i]->Type->group_dipole=true;
				}*/
				if(parameters.groups[i]->number>0){
					cout << "\t\t\t-> " << parameters.groups[i]->Type->name << "\n";
					parameters.groups[i]->Type->group_dipole=false;
				}
			}
		} else cout << "\t\t-> No changes because dyneps is off.\n";
	}
	cout << "<- Done.\n";
	
	parameters.N=parameters.n_oids-parameters.n_group_oids+parameters.n_groups;
	parameters.N_NpT=parameters.N;
	if(parameters.NpT_lnV_change) parameters.N_NpT++;
	SetSystemVolume();
	// Now is a good time to set these
	if(parameters.NpT) SetParam(parameters.NpT_move_nr,"NpT_move_nr",constant_block,nr_moving);
	if(parameters.transition){
		SetParam(parameters.transition_start_volume,"transition_start_volume",constant_block,parameters.V);
		SetParam(parameters.transition_target_volume,"transition_target_volume",constant_block,parameters.targetV);
		if(parameters.transition_start_volume>parameters.V) parameters.V=parameters.transition_start_volume;
		if(fabs(parameters.transition_start_volume-parameters.targetV)>EPS) parameters.targetV=parameters.transition_start_volume; // make sure transition method start at specified start volume
	}
	SetSystemProperties();
	conffile.file.close();
	delete[] constant_block;
	if(parameters.use_trajectory) parameters.trajfile=filename; // in case trajectory file is used it becomes the configuration file
		else parameters.trajfile=parameters.configdir+parameters.configfile;
	cout << "\nConfiguration file read and parsed.\n";
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

MC_Config global_Config;

MC_Config* GetConfig()
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return &global_Config;
}

