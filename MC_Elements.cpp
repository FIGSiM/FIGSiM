/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * updated by AT, Jan 29, 2011
 * - implemented faster calculation of MC_Elements::Rotate
 *
 * updated by AT, Feb 14, 2011
 * - fixed "funny little" bug leading to non-positive rotation matrix (in Yaw, rot[1][0] was set to Rb[2][0], instead Rb[1][0])
 *
 * updated by AT, Feb 22, 2011
 * - fixed bug in offcenter correction, now use corresponding rot matrix per dipole
 * - merged most external functions modifying elements into class
 */

#include "MC_Elements.h"

// Fit variables from alpha(gamma) plot for variable size LJ
#define alpha_a 0.000231713
#define alpha_b 0.107135
#define alpha_c 2.231724
#define alpha_d 0.017576
#define alpha_e 0.601184
#define alpha_f 2.975710

inline double alpha(double gamma)
{
	return ((alpha_a*gamma+alpha_b)*gamma+alpha_c)/((alpha_d*gamma+alpha_e)*gamma+alpha_f)+0.25;
}

#define MC_ELEMENTS

/// Standard constructor
MC_Elements::MC_Elements(Config_Data* conf)
{
	configuration=conf;
	LJlambda=configuration->LJlambda;
	selfLJlambda=1.0;
	smallest_volume=configuration->V;
	// decide which Lennar-Jones pair function to use
	switch (configuration->vdwtype){
		default:
		case 1:
			LJpair = &MC_Elements::simpleVLJ_notouch;
			break;
		case 2:
			LJpair = &MC_Elements::VLJ_notouch;
			break;
		case 3:
			LJpair = &MC_Elements::VSS_notouch;
			break;
		case 4:
			LJpair = &MC_Elements::simpleVLJ_touch;
			break;
		case 5:
			LJpair = &MC_Elements::VSE_touch;
			break;
		case 6:
			LJpair = &MC_Elements::modulate_touch;
			break;
		case 7:
			LJpair = &MC_Elements::nmVLJ_notouch;
			break;
	}
	max_charges=0;
	
	// Handle signals manually
	exitnow=false;
	for(int i=0; i<32; i++){
		switch(i){
#ifdef SIGKILL
			case SIGKILL: // uncatcheable anyway ...
#endif
#ifdef SIGSTOP
			case SIGSTOP:
#endif
#ifdef SIGILL
			case SIGILL: // do not catch the following
#endif
#ifdef SIGTRAP
			case SIGTRAP:
#endif
#ifdef SIGIOT
			case SIGIOT:
#endif
#ifdef SIGBUS
			case SIGBUS:
#endif
#ifdef SIGFPE
			case SIGFPE:
#endif
#ifdef SIGSEGV
			case SIGSEGV:
#endif
#ifdef SIGPIPE
			case SIGPIPE:
#endif
#ifdef SIGSTKFLT
			case SIGSTKFLT:
#endif
#ifdef SIGPROF
			case SIGPROF:
#endif
				break;
			default:
				signal(i,MC_Elements::exit_not);
		}
	}
#ifdef SIGCONT
	signal(SIGCONT, MC_Elements::continue_execute);
#endif
#ifdef SIGINT
	signal(SIGINT, MC_Elements::exit_gracefully);
#endif
#ifdef SIGHUP
	signal(SIGHUP, MC_Elements::exit_gracefully);
#endif
#ifdef SIGQUIT
	signal(SIGQUIT, MC_Elements::exit_gracefully);
#endif
#ifdef SIGABRT
	signal(SIGABRT, MC_Elements::exit_gracefully);
#endif
#ifdef SIGTERM
	signal(SIGTERM, MC_Elements::exit_gracefully);
#endif
#ifdef SIGPWR
	signal(SIGPWR, MC_Elements::exit_gracefully);
#endif
	
	traj_laststep=0;
	traj_created=false;
	
	elements_allocated=configuration->n_oids; // Start with number of elements as suggested by configuration
	Elements=(Element*)malloc(elements_allocated*sizeof(Element));
	if(!Elements){ // Houston control, we have a problem ...
		cout << "Could not find a memory block large enough to hold " << elements_allocated << " elements.\n";
		exit(1);
	}
	
	Groups=NULL;
	Rmus=NULL; // is created with "new" so shouldn't matter, but is tested at end of calculation hence why
	Rdist2=NULL; // is created with "new" so shouldn't matter, but is tested at end of calculation hence why
	Rdist2LJ=NULL; // is created with "new" so shouldn't matter, but is tested at end of calculation hence why
	
	N=0;
	max_group_nrelements=0;
	nr_elements=0;
	nr_group_elements=0;
	nr_groups=0;
	nr_N_to_kk=0;
	rf_correction=0.0;
	rf_neutralization=0.0;
	Evec=Vec3(0.0);
	
	epsRF=configuration->epsilon;
}

/// Destructor
MC_Elements::~MC_Elements()
{
	// Clean up element memory
	if(elements_allocated>0){
		for(unsigned int i=0; i<nr_elements; i++){
			if(Elements[i].nr_interactions>0){
				for(unsigned int j=0; j<Elements[i].nr_interactions; j++){
					if(Elements[i].interactions[j].location) delete Elements[i].interactions[j].location;
					if(Elements[i].interactions[j].normal) delete Elements[i].interactions[j].normal;
					if(Elements[i].interactions[j].tangent) delete Elements[i].interactions[j].tangent;
				}
				delete[] Elements[i].interactions;
			}
			if(Elements[i].MyType->nr_charges>0) delete[] Elements[i].q_pos;
		}
		free(Elements); // finally free ...
	}
	// Clean up group memory
	if(max_group_nrelements>0) FreeGroupStorage();
	if(nr_groups>0){
		for(unsigned int i=0; i<nr_groups; i++){
			if(Groups[i]->nr_elements>0) delete[] Groups[i]->elements;
			delete Groups[i];
		}
		free(Groups);
	}
}

/// Add a new Element to class and return number of element
unsigned int MC_Elements::AddElement(const unsigned int type)
{
	N++;
	if(nr_elements>=elements_allocated){ // Need more space (> is just a safety measure, should not occur)
		elements_allocated+=16; // reallocate memory in chunks of 16 elements (blocks of 16*160 Bytes = 2560 Bytes)
		Elements=(Element*)realloc(Elements,elements_allocated*sizeof(Element));
		if(!Elements){ // Houston control, we have a problem ...
			cout << "Could not find a memory block large enough to hold " << elements_allocated << " elements.\n";
			exit(1);
		}
	}
	Elements[nr_elements].mytype = type;
	Elements[nr_elements].MyType = configuration->element_types[type];
	Elements[nr_elements].delta=1.0;
	Elements[nr_elements].fixed=false;
	Elements[nr_elements].gamma_quarter=0.0;
	if(Elements[nr_elements].MyType->nr_charges>0){
		if(Elements[nr_elements].MyType->nr_charges>max_charges) max_charges=Elements[nr_elements].MyType->nr_charges;
		Elements[nr_elements].q_pos=new Vec3[Elements[nr_elements].MyType->nr_charges];
		for(unsigned int i=0; i<Elements[nr_elements].MyType->nr_charges; i++) Elements[nr_elements].q_pos[i]=Elements[nr_elements].MyType->q_pos[i];
	}
	Elements[nr_elements].dipole = configuration->element_types[type]->initial_dipole;
	Elements[nr_elements].rot = Mat33();
	Elements[nr_elements].center = Vec3(0.0);
	Elements[nr_elements].ds = Vec3(0.0);
	Elements[nr_elements].current_s = Vec3(0.0);
	Elements[nr_elements].dtheta = Vec3(0.0);
	
	Elements[nr_elements].nr_interactions = 0;
	Elements[nr_elements].interactions = NULL;
	Elements[nr_elements].group = NULL;
	Elements[nr_elements].VE=0.0;
	double V=4.0*pi/3.0*Elements[nr_elements].MyType->saxes.vec[0]*Elements[nr_elements].MyType->saxes.vec[1]*Elements[nr_elements].MyType->saxes.vec[2];
	if(V<smallest_volume) smallest_volume=V;
	
	nr_elements++;
	return nr_elements-1;
}

void MC_Elements::CreateGroupStorage()
{
	rebuild_storage.element_storage=new unsigned int[2*max_group_nrelements];
	rebuild_storage.visited=new unsigned int[max_group_nrelements];
	rebuild_storage.to_link=new unsigned int[max_group_nrelements];
	rebuild_storage.randomized=new unsigned int[max_group_nrelements];
	rebuild_storage.en_route=new bool[max_group_nrelements];
	rebuild_storage.delta_trans=new Vec3[max_group_nrelements];
	rebuild_storage.delta_rot=new Mat33[max_group_nrelements];
	rebuild_storage.theta2=new double[max_group_nrelements];
	group_storage=new group_save*[configuration->num_groups];
	Element* element;
	for(unsigned int i=0; i<configuration->num_groups; i++){
		group_storage[i]=new group_save[configuration->groups[i]->nr_elements];
		for(unsigned int j=0; j<configuration->groups[i]->nr_elements; j++){
			element=&configuration->group_elements[configuration->groups[i]->elements[j]];
			if(element->nr_interactions>0){
				group_storage[i][j].interaction_location=new Vec3[element->nr_interactions];
				group_storage[i][j].interaction_normal=new Vec3[element->nr_interactions];
				group_storage[i][j].interaction_tangent=new Vec3[element->nr_interactions];
			}
			if(element->MyType->nr_charges>0) group_storage[i][j].q_pos=new Vec3[element->MyType->nr_charges];
		}
	}
}

void MC_Elements::FreeGroupStorage()
{
	Element* element;
	for(unsigned int i=0; i<configuration->num_groups; i++){
		for(unsigned int j=0; j<configuration->groups[i]->nr_elements; j++){
			element=&configuration->group_elements[configuration->groups[i]->elements[j]];
			if(element->nr_interactions>0){
				delete[] group_storage[i][j].interaction_location;
				delete[] group_storage[i][j].interaction_normal;
				delete[] group_storage[i][j].interaction_tangent;
			}
			if(element->MyType->nr_charges>0) delete[] group_storage[i][j].q_pos;
		}
		delete[] group_storage[i];
	}
	delete[] group_storage;
	delete[] rebuild_storage.element_storage;
	delete[] rebuild_storage.visited;
	delete[] rebuild_storage.to_link;
	delete[] rebuild_storage.randomized;
	delete[] rebuild_storage.en_route;
	delete[] rebuild_storage.delta_trans;
	delete[] rebuild_storage.delta_rot;
	delete[] rebuild_storage.theta2;
}

/// Add a new element group to class
Element_Group* MC_Elements::AddGroup(const unsigned int type)
{
	N++;
	nr_groups++;
	Element_Group* archetype = configuration->groups[type];
	// reallocate memory for rebuild storage if more needed
	if(archetype->nr_elements>max_group_nrelements) max_group_nrelements=archetype->nr_elements;
	// Create new (blank) group element
	Groups=(Element_Group**)realloc(Groups,nr_groups*sizeof(Element_Group*));
	if(!Groups){
		cout << "Not enough memory to create group elements.\n";
		exit(3);
	}
	Element_Group* group = new Element_Group;
	Groups[nr_groups-1] = group;
	group->number=nr_groups-1; // running counter of groups created
	group->type=type;
	group->levelofdetail=archetype->levelofdetail;
	group->Type=archetype->Type; // pointer ...
	group->nr_elements=archetype->nr_elements;
	group->inv_nr_elements=1.0/archetype->nr_elements;
	group->elements = new unsigned int[group->nr_elements];
	group->nr_potentials=archetype->nr_potentials;
	group->potentials=archetype->potentials; // pointer ...
	group->fixed=false;
	// Add group elements
	Element* current;
	Element* element;
	group->center.vec[0]=0.0; group->center.vec[1]=0.0; group->center.vec[2]=0.0;
	if(archetype->Type->Volume<smallest_volume) smallest_volume=archetype->Type->Volume;
	double save_sv=smallest_volume;
	for(unsigned int i=0; i<group->nr_elements; i++){
		element=&configuration->group_elements[archetype->elements[i]];
		group->elements[i]=AddElement(element->MyType->thistype);
		N--;
		nr_group_elements++;
		// set element properties
		current=&Elements[group->elements[i]];
		RotateElement(current,element->rot);
		current->group=group; // this element belongs to this group
		current->center=element->center-configuration->group_centers[type]; // adjust parameters of created group so bounding sphere center is (0,0,0)
		group->center+=current->center; // recalculate geometric center on the fly
		current->nr_interactions=element->nr_interactions;
		if(element->nr_interactions>0) current->interactions = new Interaction[current->nr_interactions];
		for(unsigned int j=0; j<current->MyType->nr_charges; j++) current->q_pos[j]=current->rot*current->MyType->q_pos[j];
		for(unsigned int j=0; j<current->nr_interactions; j++){
			current->interactions[j].fixed=element->interactions[j].fixed;
			current->interactions[j].allow_bond_stretch=element->interactions[j].allow_bond_stretch;
			current->interactions[j].allow_bond_bend=element->interactions[j].allow_bond_bend;
			current->interactions[j].bond_length=element->interactions[j].bond_length;
			current->interactions[j].potential_map[0]=0.0; // since my memory left me too,
			current->interactions[j].potential_map[1]=0.0; // potential_map is for mapping
			current->interactions[j].potential_map[2]=0.0; // arbitrary potentials to link restrictions (if wanted)
			current->interactions[j].nr_potentials=element->interactions[j].nr_potentials;
			current->interactions[j].potentials=element->interactions[j].potentials; // pointer to array of pointers to potentials will be filled later
			current->interactions[j].partner=element->interactions[j].partner;
			current->interactions[j].back_link=element->interactions[j].back_link;
			
			current->interactions[j].initial_location=element->interactions[j].initial_location; // pointer
			current->interactions[j].initial_normal=element->interactions[j].initial_normal; // pointer
			current->interactions[j].initial_tangent=element->interactions[j].initial_tangent; // pointer
			
			if(element->interactions[j].location){
				current->interactions[j].location=new Vec3;
				*current->interactions[j].location=*element->interactions[j].location;
			} else current->interactions[j].location=NULL;
			
			if(element->interactions[j].normal){
				current->interactions[j].normal=new Vec3;
				*current->interactions[j].normal=*element->interactions[j].normal;
			} else current->interactions[j].normal=NULL;
			
			if(element->interactions[j].tangent){
				current->interactions[j].tangent=new Vec3;
				*current->interactions[j].tangent=*element->interactions[j].tangent;
			} else current->interactions[j].tangent=NULL;
		}
	}
	nr_N_to_kk=nr_group_elements-nr_groups;
	smallest_volume=save_sv;
	group->center*=group->inv_nr_elements;
	return group;
}

/*!
 * Calculation of Lennard-Jones interactions. Notouch variants can only handle spheres, Touch variants handle anisotropy
 * using the method of Perram and Wertheim (J. Comp. Phys, 1985). Several functions are available, providing tradeoffs
 * between generality and speed
 * Variable defitions:
 * 	kk - the index of the currently molecule
 * 	i - the index of the molecule it is interacting with
 * 	distance - the norm of the distance vector between the molecules
 * Return value is the LJ energy for one pair of ellipsoids
 *
 * Soft sphere - only calculate nuclear repulsion term of LJ potential
 */
double MC_Elements::VSS_notouch(const unsigned int kk, const unsigned int i)
{
	unsigned int ikk = Elements[i].mytype*configuration->num_element_types+Elements[kk].mytype;
	//Calculate energy
	return configuration->pre_eps[ikk]*configuration->r*qqpwr(configuration->pre_sigma2[ikk]/Rdist2LJ[i],configuration->LJexp[0]);
}

/// Simple Lennard-Jones potential (nuclear repulsion power is twice that of dispersion, e.g. LJ 12-6 potential)
double MC_Elements::simpleVLJ_notouch(const unsigned int kk, const unsigned int  i)
{
	unsigned int ikk = Elements[i].mytype*configuration->num_element_types+Elements[kk].mytype;
	double disp = qqpwr(configuration->pre_sigma2[ikk]/Rdist2LJ[i],configuration->LJexp[0]>>1);
#if DEBUG_LEVEL>3
		cout << configuration->pre_eps[ikk] << ", " << configuration->r << ", " << sqrt(configuration->pre_sigma2[ikk]) << ", " << sqrt(Rdist2LJ[i]) << ", " << LJlambda << " = " << configuration->pre_eps[ikk]*(configuration->r*disp*(disp-LJlambda)) << "\n";
#endif
	//Calculate energy
	return configuration->pre_eps[ikk]*(configuration->r*disp*(disp-LJlambda));
}

/// n-m Lennard-Jones potential
double MC_Elements::nmVLJ_notouch(const unsigned int kk, const unsigned int  i)
{
	unsigned int ikk = Elements[i].mytype*configuration->num_element_types+Elements[kk].mytype;
	double disp = sqrt(configuration->pre_sigma2[ikk]/Rdist2LJ[i]);
	//Calculate energy
	return configuration->pre_eps[ikk]*configuration->r*(qqpwr(disp,configuration->LJexp[1])-LJlambda*qqpwr(disp,configuration->LJexp[0]));
}


/// General LJ potential including optional Bruce correction (attenuating Gaussian function in bottom of well)
double MC_Elements::VLJ_notouch(const unsigned int kk, const unsigned int i)
{
	unsigned int ikk = Elements[i].mytype*configuration->num_element_types+Elements[kk].mytype;
	
	//Calculate energy
	double rmrmin = sqrt(Rdist2LJ[i]) - configuration->Solvent[2];
	double sx = configuration->pre_sigma2[ikk]/Rdist2LJ[i];
	
	//LJ + Bruce correction around r=rmin. Tanh functions have been replaced with a Gaussian for speed
	return configuration->pre_eps[ikk]*configuration->r*(qqpwr(sx,configuration->LJexp[1]>>1)-LJlambda*qqpwr(sx,configuration->LJexp[0]>>1)) + configuration->Solvent[1]*exp(configuration->Solvent[0]*rmrmin*rmrmin);
}

/// Simplified VLJ, but using Touch algorithm to handle non-spherical ellipsoids. This ONLY works for even powers!
/// -- extended to be able to use square root LJ epsilon texture if present for either party
double MC_Elements::simpleVLJ_touch(const unsigned int kk, const unsigned int i)
{
	double Rx; //reduced distance (r0/r)^2
	unsigned int ikk = Elements[i].mytype*configuration->num_element_types+Elements[kk].mytype;
	double epsilon = configuration->pre_eps[ikk];
	
	//Determine whether touch should be run or the average radius should be used
	if((configuration->pre_touch[ikk]) && ((!configuration->touchtrunc) || (Rdist2LJ[i]<configuration->tcut2[ikk]))){
		Rx = touch(kk, i); // Square of reduced distance corrected based on closest contact
		//Calculate energy
		Element_Type* type_i=Elements[i].MyType;
		Element_Type* type_kk=Elements[kk].MyType;
		Vec3 rmu_i_frame=Elements[i].rot.TransMulVec(Rmus[i]*(-1.0)); // Rmu is vector from element kk to element i
		Vec3 rmu_kk_frame=Elements[kk].rot.TransMulVec(Rmus[i]);
		if(configuration->LJ_adjust_width){
			double wi=type_i->avg_width;
			double wkk=type_kk->avg_width;
			if(wi<EPS) wi=EllipsoidRmin(Rmus[i],type_i->saxes,Elements[i].rot);
			if(wkk<EPS) wkk=EllipsoidRmin(Rmus[i],type_kk->saxes,Elements[kk].rot);
			double s=wi+wkk;
			double rab=sqrt(Rdist2LJ[i]);
			double r=rab*(1.0-sqrt(Rx))+s;
			if(r<0.0) r=0.0;
			Rx=s/r;
			Rx*=Rx;
		}
		double phi, cost;
		unsigned int t,p;
		if(type_i->eps_texture){
			cost=VecDist2Phi(rmu_i_frame,sqrt(Rdist2LJ[i]),phi);
			t=(unsigned int)((1.0-cost)*theta_res_half-0.5); // looks dangerous, but (unsigned int)(-0.5) is 0
			p=(phi_res+(int)(phi*phi_res_per_tau))%phi_res;
			epsilon*=type_i->eps_texture[t*phi_res+p];
		} else if(configuration->LJ_interaction_area) epsilon*=IA(type_i,rmu_i_frame);
		if(type_kk->eps_texture){
			cost=VecDist2Phi(rmu_kk_frame,sqrt(Rdist2LJ[i]),phi);
			t=(unsigned int)((1.0-cost)*theta_res_half-0.5);
			p=(phi_res+(int)(phi*phi_res_per_tau))%phi_res;
			epsilon*=type_kk->eps_texture[t*phi_res+p];
		} else if(configuration->LJ_interaction_area) epsilon*=IA(type_kk,rmu_kk_frame);
	} else{
		Rx = configuration->pre_sigma2[ikk]/Rdist2LJ[i]; // (r0/r)^2
	}
	//Calculate parameters
	double disp = qqpwr(Rx,configuration->LJexp[0]>>1); // Since Rx is already squared, the small term must be halved (bitshift right by one is integer/2 - AT)
	return epsilon*(configuration->r*disp*(disp-LJlambda));
}

/// Modulated VLJ (adjustable sigma) using Touch algorithm to handle non-spherical ellipsoids.
/// -- extended to be able to use square root LJ epsilon texture if present for either party
double MC_Elements::modulate_touch(const unsigned int kk, const unsigned int i)
{
	double Rx; //reduced distance (r0/r)^2
	Element* element_i=&Elements[i];
	Element* element_kk=&Elements[kk];
	Element_Type* type_i=element_i->MyType;
	Element_Type* type_kk=element_kk->MyType;
	
	unsigned int ikk = element_i->mytype*configuration->num_element_types+element_kk->mytype;
	double epsilon = configuration->pre_eps[ikk];
	
	//Determine whether touch should be run or the average radius should be used
	if((configuration->pre_touch[ikk]) && ((!configuration->touchtrunc) || (Rdist2LJ[i]<configuration->tcut2[ikk]))){
		Rx = touch(kk, i); // Square of reduced distance corrected based on closest contact
		//Calculate energy
		Element_Type* type_i=Elements[i].MyType;
		Element_Type* type_kk=Elements[kk].MyType;
		Vec3 rmu_i_frame=Elements[i].rot.TransMulVec(Rmus[i]*(-1.0)); // Rmu is vector from element kk to element i
		Vec3 rmu_kk_frame=Elements[kk].rot.TransMulVec(Rmus[i]);
		if(configuration->LJ_adjust_width){
			double wi=type_i->avg_width;
			double wkk=type_kk->avg_width;
			if(wi<EPS) wi=EllipsoidRmin(Rmus[i],type_i->saxes,Elements[i].rot);
			if(wkk<EPS) wkk=EllipsoidRmin(Rmus[i],type_kk->saxes,Elements[kk].rot);
			double s=wi+wkk;
			double r=sqrt(Rdist2LJ[i])-(sqrt(Rx*Rdist2LJ[i])-s); // r - (sigma - width)
			if(r<0.0) r=0.0;
			Rx=s/r;
			Rx*=Rx;
		}
		double phi, cost;
		unsigned int t,p;
		if(type_i->eps_texture){
			cost=VecDist2Phi(rmu_i_frame,sqrt(Rdist2LJ[i]),phi);
			t=(unsigned int)((1.0-cost)*theta_res_half-0.5); // looks dangerous, but (unsigned int)(-0.5) is 0
			p=(phi_res+(int)(phi*phi_res_per_tau))%phi_res;
			epsilon*=type_i->eps_texture[t*phi_res+p];
		} else if(configuration->LJ_interaction_area) epsilon*=IA(type_i,rmu_i_frame);
		if(type_kk->eps_texture){
			cost=VecDist2Phi(rmu_kk_frame,sqrt(Rdist2LJ[i]),phi);
			t=(unsigned int)((1.0-cost)*theta_res_half-0.5);
			p=(phi_res+(int)(phi*phi_res_per_tau))%phi_res;
			epsilon*=type_kk->eps_texture[t*phi_res+p];
		} else if(configuration->LJ_interaction_area) epsilon*=IA(type_kk,rmu_kk_frame);
	} else{
		double sigma=element_i->delta*average(type_i->saxes.vec,3)+element_kk->delta*average(type_kk->saxes.vec,3);
		Rx = sigma*sigma/Rdist2LJ[i]; // (r0/r)^2
	}
	//Calculate parameters
	double gamma=element_i->gamma_quarter+element_kk->gamma_quarter;
	epsilon*=alpha(gamma);
	double disp = qqpwr(Rx,configuration->LJexp[0]>>1); // Since Rx is already squared, the small term must be halved (bitshift right by one is integer/2 - AT)
	return epsilon*(configuration->r*disp*(disp*pow(Rx,gamma)-LJlambda));
}

/// Soft-sphere potential, but using Touch algorithm to handle non-spherical ellipsoids
double MC_Elements::VSE_touch(const unsigned int kk, const unsigned int i)
{
	double Rx; //reduced distance (r0/r)^2
	unsigned int ikk = Elements[i].mytype*configuration->num_element_types+Elements[kk].mytype;
	double epsilon = configuration->pre_eps[ikk];
	
	//Determine whether touch should be run or the average radius should be used
	if((configuration->pre_touch[ikk]) && ((!configuration->touchtrunc) || (Rdist2LJ[i]<configuration->tcut2[ikk]))){
		Rx = touch(kk, i); // Square of reduced distance corrected based on closest contact
		//Calculate energy
		Element_Type* type_i=Elements[i].MyType;
		Element_Type* type_kk=Elements[kk].MyType;
		Vec3 rmu_i_frame=Elements[i].rot.TransMulVec(Rmus[i]*(-1.0)); // Rmu is vector from element kk to element i
		Vec3 rmu_kk_frame=Elements[kk].rot.TransMulVec(Rmus[i]);
		if(configuration->LJ_adjust_width){
			double wi=type_i->avg_width;
			double wkk=type_kk->avg_width;
			if(wi<EPS) wi=EllipsoidRmin(Rmus[i],type_i->saxes,Elements[i].rot);
			if(wkk<EPS) wkk=EllipsoidRmin(Rmus[i],type_kk->saxes,Elements[kk].rot);
			double s=wi+wkk;
			double rab=sqrt(Rdist2LJ[i]);
			double r=rab*(1.0-sqrt(Rx))+s;
			if(r<0.0) r=0.0;
			Rx=s/r;
			Rx*=Rx;
		}
		double phi, cost;
		unsigned int t,p;
		if(type_i->eps_texture){
			cost=VecDist2Phi(rmu_i_frame,sqrt(Rdist2LJ[i]),phi);
			t=(unsigned int)((1.0-cost)*theta_res_half-0.5); // looks dangerous, but (unsigned int)(-0.5) is 0
			p=(phi_res+(int)(phi*phi_res_per_tau))%phi_res;
			epsilon*=type_i->eps_texture[t*phi_res+p];
		} else if(configuration->LJ_interaction_area) epsilon*=IA(type_i,rmu_i_frame);
		if(type_kk->eps_texture){
			cost=VecDist2Phi(rmu_kk_frame,sqrt(Rdist2LJ[i]),phi);
			t=(unsigned int)((1.0-cost)*theta_res_half-0.5);
			p=(phi_res+(int)(phi*phi_res_per_tau))%phi_res;
			epsilon*=type_kk->eps_texture[t*phi_res+p];
		} else if(configuration->LJ_interaction_area) epsilon*=IA(type_kk,rmu_kk_frame);
	} else{
		Rx = configuration->pre_sigma2[ikk]/Rdist2LJ[i]; // (r0/r)^2
	}
	//Calculate energy
	return configuration->pre_eps[ikk]*configuration->r*qqpwr(Rx,configuration->LJexp[0]); // since Rx is already squared, this gives the larger (repulsive term)
}

void MC_Elements::Translate(Element_Group* group, Vec3 delta) ///< Translate group by delta
{
	// easy :-)
	for(unsigned int i=0; i<group->nr_elements; i++) Elements[group->elements[i]].center+=delta;
	group->center+=delta;
}

double MC_Elements::Rotate(Element_Group* group, Vec3 &center) ///< Rotate group around center by (r)oll, (y)aw, (p)itch angles
{
	//Random number for picking an axis and for the move itself
#ifdef USE_CMWC4096
	unsigned int choice=CMWC4096()%3;
	double theta = configuration->maxrot*(2.0*ranQ()-1.0)*group->Type->inv_sqrt_nrpc;
#else
	unsigned int choice=(unsigned int)ran2_int(*configuration->idum)%3;
	double theta = configuration->maxrot*(2.0*ran2(*configuration->idum)-1.0)*group->Type->inv_sqrt_nrpc;
#endif
	double co = cos(theta);
	double si = sin(theta);
	
	Mat33 rot;
	switch(choice){
		case 0: /* yaw: x-rot
			 * (1   0   0)
			 * (0  co  si)
			 * (0 -si  co)
			 */
			rot.mat[0][0] = 1;	rot.mat[0][1] = 0;	rot.mat[0][2] = 0;
			rot.mat[1][0] = 0;	rot.mat[1][1] = co;	rot.mat[1][2] = si;
			rot.mat[2][0] = 0;	rot.mat[2][1] = -si;	rot.mat[2][2] = co;
			break;
		case 1: /* pitch: y-rot
			 * ( co  0  si)
			 * ( 0   1   0)
			 * (-si  0  co)
			 */
			rot.mat[0][0] = co;	rot.mat[0][1] = 0;	rot.mat[0][2] = si;
			rot.mat[1][0] = 0;	rot.mat[1][1] = 1;	rot.mat[1][2] = 0;
			rot.mat[2][0] = -si;	rot.mat[2][1] = 0;	rot.mat[2][2] = co;
			break;
		case 2: /* roll: z-rot
			 * ( co  si  0)
			 * (-si  co  0)
			 * ( 0   0   1)
			 */
			rot.mat[0][0] = co;	rot.mat[0][1] = si;	rot.mat[0][2] = 0;
			rot.mat[1][0] = -si;	rot.mat[1][1] = co;	rot.mat[1][2] = 0;
			rot.mat[2][0] = 0;	rot.mat[2][1] = 0;	rot.mat[2][2] = 1;
			break;
	}
	// Apply to group elements
	Rotate(group,center,rot);
	theta*=group->nr_elements;
	
	return theta*theta;
}

void MC_Elements::Rotate(Element_Group* group, Vec3 &center, double r, double y, double p) ///< Rotate group around center by (r)oll, (y)aw, (p)itch angles
{
	// Get rotation matrix
	double co,si;
	Mat33 roll,yaw,pitch;
	/* roll:
	 * ( co  si  0)
	 * (-si  co  0)
	 * ( 0   0   1)
	 */
	co = cos(r); si = sin(r);
	roll.mat[0][0] = co;	roll.mat[0][1] = si;	roll.mat[0][2] = 0;
	roll.mat[1][0] = -si;	roll.mat[1][1] = co;	roll.mat[1][2] = 0;
	roll.mat[2][0] = 0;	roll.mat[2][1] = 0;	roll.mat[2][2] = 1;
	/* yaw:
	 * (1   0   0)
	 * (0  co  si)
	 * (0 -si  co)
	 */
	co = cos(y); si = sin(y);
	yaw.mat[0][0] = 1;	yaw.mat[0][1] = 0;	yaw.mat[0][2] = 0;
	yaw.mat[1][0] = 0;	yaw.mat[1][1] = co;	yaw.mat[1][2] = si;
	yaw.mat[2][0] = 0;	yaw.mat[2][1] = -si;	yaw.mat[2][2] = co;
	/* pitch:
	 * ( co  0  si)
	 * ( 0   1   0)
	 * (-si  0  co)
	 */
	co = cos(p); si = sin(p);
	pitch.mat[0][0] = co;	pitch.mat[0][1] = 0;	pitch.mat[0][2] = si;
	pitch.mat[1][0] = 0;	pitch.mat[1][1] = 1;	pitch.mat[1][2] = 0;
	pitch.mat[2][0] = -si;	pitch.mat[2][1] = 0;	pitch.mat[2][2] = co;

	Mat33 rot=roll*yaw*pitch;

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

void MC_Elements::Translate(const unsigned int idx, Vec3 delta)
{
	Elements[idx].center+=delta;
	apply_PBCs(Elements[idx].center);
}


/*!
 * Regular translate method for cubic cell with PBCs in an arbitrary number of directions.
 * If a move in a non-PBC'd direction would cause the particle to leave the box, the component
 * of the move along the offending axis is ignored, but the components along other axes
 * are retained. Return value is the squared distance of the move.
 */
double MC_Elements::Translate(const unsigned int idx)
{
	double trans, transtemp;
	double moved = 0.0;
	
	for (int i = 0; i < 3; i++){
		//Get a random number
#ifdef USE_CMWC4096
		trans = (configuration->maxtrans)*(2.0*ranQ()-1.0);
#else
		trans = (configuration->maxtrans)*(2.0*ran2(*configuration->idum)-1.0);
#endif
		
		//Make a move appropriate to the PBCs along that axis
		if(configuration->PBCs[i]){
			Elements[idx].center.vec[i] += trans;
			Elements[idx].ds.vec[i] += trans;
			if(fabs(Elements[idx].center.vec[i])>configuration->boxlength[i]*0.5) Elements[idx].center.vec[i] -= copysign(configuration->boxlength[i],Elements[idx].center.vec[i]);
			moved += trans*trans;
		} else{
			transtemp = Elements[idx].center.vec[i] + trans;
			if (abs(2*transtemp) < configuration->boxlength[i]){
				Elements[idx].center.vec[i] = transtemp;
				Elements[idx].ds.vec[i] += trans;
				moved += trans*trans;
			}
		}
	}
	return moved;
}

/*!
 * Translate method for isolated spherical cell (no PBCs)
 * If a move would cause a particle to leave the sphere, the old coords are retained.
 * Return value is the squared distance of the move
 */
double MC_Elements::Translate_Sphere(const unsigned int idx)
{
	double randmov;
	double moved = 0.0;
	Vec3 transtemp(Elements[idx].center);

	//Try a random move
	for (int i = 0; i<3; i++){
#ifdef USE_CMWC4096
		randmov = (configuration->maxtrans)*(2.0*ranQ()-1.0);
#else
		randmov = (configuration->maxtrans)*(2.0*ran2(*configuration->idum)-1.0);
#endif
		transtemp.vec[i] += randmov;
		moved += randmov*randmov;
	}

	double rtemp = transtemp.V3Norm();

	//Accept the move only if the particle is still in the sphere
	if(rtemp < configuration->spherecylr){
		Elements[idx].center = transtemp;
	}
	else moved = 0.0;

	return moved;
}

/*!
 * Translate method for cylindrical cell (PBCs currently only available on z-axis)
 * If a move would cause a particle to leave the cylinder, the component of the
 * move along the offending axis (r or z) is discarded. Moves along other axes
 * are retained. Return value is the squared distance of the move.
 */
double MC_Elements::Translate_Cylinder(const unsigned int idx)
{

	double trans1, trans2, trans3, chkr, pbcadj, transtemp;
	double moved = 0.0;
	
	//Get all of the random numbers
#ifdef USE_CMWC4096
	trans1 = (configuration->maxtrans)*(2.0*ranQ()-1.0);
	trans2 = (configuration->maxtrans)*(2.0*ranQ()-1.0);
	trans3 = (configuration->maxtrans)*(2.0*ranQ()-1.0);
#else
	trans1 = (configuration->maxtrans)*(2.0*ran2(*configuration->idum)-1.0);
	trans2 = (configuration->maxtrans)*(2.0*ran2(*configuration->idum)-1.0);
	trans3 = (configuration->maxtrans)*(2.0*ran2(*configuration->idum)-1.0);
#endif
	//Check to see if the r componentof the move is ok
	chkr = sqrt(qpwr(trans1+Elements[idx].center.vec[0], 2) + qpwr(trans2+Elements[idx].center.vec[1],2));
	if (chkr < configuration->spherecylr)
	{
		Elements[idx].center.vec[0] += trans1;
		Elements[idx].center.vec[1] += trans2;
		moved = trans1*trans1+trans2*trans2;
	}

	//Make a move appropriate to the PBCs along the z-axis
	if (configuration->PBCs[2] == 1)
	{
		Elements[idx].center.vec[2] += trans3;
		pbcadj = Elements[idx].center.vec[2]*configuration->inv_boxlength[2];
		Elements[idx].center.vec[2] -= configuration->boxlength[2]*specialround(pbcadj);
		moved += trans3*trans3;
	}
	else
	{
		transtemp = Elements[idx].center.vec[2] + trans3;
		if (abs(2*transtemp) < configuration->boxlength[2])
		{
			Elements[idx].center.vec[2] = transtemp;
			moved += trans3*trans3;
		}
	}
	return moved;
}

void MC_Elements::Rotate(const unsigned int idx, double r, double y, double p) ///< Rotate element around center by (r)oll, (y)aw, (p)itch angles
{
	// Get rotation matrix
	double co,si;
	Mat33 roll,yaw,pitch;
	/* roll:
	 * ( co  si  0)
	 * (-si  co  0)
	 * ( 0   0   1)
	 */
	co = cos(r); si = sin(r);
	roll.mat[0][0] = co;	roll.mat[0][1] = si;	roll.mat[0][2] = 0;
	roll.mat[1][0] = -si;	roll.mat[1][1] = co;	roll.mat[1][2] = 0;
	roll.mat[2][0] = 0;	roll.mat[2][1] = 0;	roll.mat[2][2] = 1;
	/* yaw:
	 * (1   0   0)
	 * (0  co  si)
	 * (0 -si  co)
	 */
	co = cos(y); si = sin(y);
	yaw.mat[0][0] = 1;	yaw.mat[0][1] = 0;	yaw.mat[0][2] = 0;
	yaw.mat[1][0] = 0;	yaw.mat[1][1] = co;	yaw.mat[1][2] = si;
	yaw.mat[2][0] = 0;	yaw.mat[2][1] = -si;	yaw.mat[2][2] = co;
	/* pitch:
	 * ( co  0  si)
	 * ( 0   1   0)
	 * (-si  0  co)
	 */
	co = cos(p); si = sin(p);
	pitch.mat[0][0] = co;	pitch.mat[0][1] = 0;	pitch.mat[0][2] = si;
	pitch.mat[1][0] = 0;	pitch.mat[1][1] = 1;	pitch.mat[1][2] = 0;
	pitch.mat[2][0] = -si;	pitch.mat[2][1] = 0;	pitch.mat[2][2] = co;

	Mat33 rot=roll*yaw*pitch;

	// Apply to group elements
	Vec3 delta;
	Element* element=&Elements[idx];
	// Now take care of element properties
	element->rot=rot*element->rot; // takes care of things like mupos, charge position ...
	// new dipole vector
	element->dipole=element->rot*element->MyType->initial_dipole;
	// charges
	for(unsigned int j=0; j<element->MyType->nr_charges; j++) element->q_pos[j]=element->rot*element->MyType->q_pos[j];
	// link related properties (just in case ...)
	for(unsigned int j=0; j<element->nr_interactions; j++){
		if(element->interactions[j].location) *(element->interactions[j].location)=element->rot*(*(element->interactions[j].initial_location));
		if(element->interactions[j].normal) *(element->interactions[j].normal)=element->rot*(*(element->interactions[j].initial_normal));
		if(element->interactions[j].tangent) *(element->interactions[j].tangent)=element->rot*(*(element->interactions[j].initial_tangent));
	}
}

/*!
 * Rotation method for all simulation types
 * Rotates around a single axis, picking at random, limited by configuration->maxrot.
 * Return value is the squared angular displacement of the move.
 */
double MC_Elements::Rotate(const unsigned int idx)
{
	Mat33 Rb;
	Mat33* rot;
	Element* element=&Elements[idx];
	rot = &element->rot; // to keep code at least a bit readable, rot is now the rotation matrix of our Element
	//Random number for picking an axis
#ifdef USE_CMWC4096
	unsigned int choice = CMWC4096()%3; // it so happens that the range is divisible by 3 ...
	double theta = configuration->maxrot*(2.0*ranQ()-1.0);
#else
	unsigned int choice = (unsigned int)ran2_int(*configuration->idum)%3;
	double theta = configuration->maxrot*(2.0*ran2(*configuration->idum)-1.0);
#endif
	//Random number for the move itself
	double co = cos(theta);
	double si = sin(theta);
#ifdef LAB_FRAME_ROTATION
	switch(choice){
		case 0: /* yaw, lab frame:
			 * (1   0   0)   (a00 a01 a02)   (     a00            a01            a02     )
			 * (0  co  si) * (a10 a11 a12) = (co*a10+si*a20  co*a11+si*a21  co*a12+si*a22)
			 * (0 -si  co)   (a20 a21 a22)   (co*a20-si*a10  co*a21-si*a11  co*a22-si*a12)
			 */
			Rb.mat[1][0] = co*rot->mat[1][0]+si*rot->mat[2][0]; Rb.mat[1][1] = co*rot->mat[1][1]+si*rot->mat[2][1]; Rb.mat[1][2] = co*rot->mat[1][2]+si*rot->mat[2][2];
			Rb.mat[2][0] = co*rot->mat[2][0]-si*rot->mat[1][0]; Rb.mat[2][1] = co*rot->mat[2][1]-si*rot->mat[1][1]; Rb.mat[2][2] = co*rot->mat[2][2]-si*rot->mat[1][2];
			rot->mat[1][0] = Rb.mat[1][0]; rot->mat[1][1] = Rb.mat[1][1] ; rot->mat[1][2] = Rb.mat[1][2];
			rot->mat[2][0] = Rb.mat[2][0]; rot->mat[2][1] = Rb.mat[2][1] ; rot->mat[2][2] = Rb.mat[2][2];
			break;
		case 1: /* pitch
			 * ( co  0  si)   (a00 a01 a02)   (co*a00+si*a20  co*a01+si*a21  co*a02+si*a22)
			 * ( 0   1   0) * (a10 a11 a12) = (     a10            a11            a12     )
			 * (-si  0  co)   (a20 a21 a22)   (co*a20-si*a00  co*a21-si*a01  co*a22-si*a02)
			 */
			Rb.mat[0][0] = co*rot->mat[0][0]+si*rot->mat[2][0]; Rb.mat[0][1] = co*rot->mat[0][1]+si*rot->mat[2][1]; Rb.mat[0][2] = co*rot->mat[0][2]+si*rot->mat[2][2];
			Rb.mat[2][0] = co*rot->mat[2][0]-si*rot->mat[0][0]; Rb.mat[2][1] = co*rot->mat[2][1]-si*rot->mat[0][1]; Rb.mat[2][2] = co*rot->mat[2][2]-si*rot->mat[0][2];
			rot->mat[0][0] = Rb.mat[0][0]; rot->mat[0][1] = Rb.mat[0][1] ; rot->mat[0][2] = Rb.mat[0][2];
			rot->mat[2][0] = Rb.mat[2][0]; rot->mat[2][1] = Rb.mat[2][1] ; rot->mat[2][2] = Rb.mat[2][2];
			break;
		case 2: /* roll, lab frame:
			 * ( co  si  0)   (a00 a01 a02)   (co*a00+si*a10  co*a01+si*a11  co*a02+si*a12)
			 * (-si  co  0) * (a10 a11 a12) = (co*a10-si*a00  co*a11-si*a01  co*a12-si*a02)
			 * ( 0   0   1)   (a20 a21 a22)   (     a20            a21            a22     )
			 */
			Rb.mat[0][0] = co*rot->mat[0][0]+si*rot->mat[1][0]; Rb.mat[0][1] = co*rot->mat[0][1]+si*rot->mat[1][1]; Rb.mat[0][2] = co*rot->mat[0][2]+si*rot->mat[1][2];
			Rb.mat[1][0] = co*rot->mat[1][0]-si*rot->mat[0][0]; Rb.mat[1][1] = co*rot->mat[1][1]-si*rot->mat[0][1]; Rb.mat[1][2] = co*rot->mat[1][2]-si*rot->mat[0][2];
			rot->mat[0][0] = Rb.mat[0][0]; rot->mat[0][1] = Rb.mat[0][1] ; rot->mat[0][2] = Rb.mat[0][2];
			rot->mat[1][0] = Rb.mat[1][0]; rot->mat[1][1] = Rb.mat[1][1] ; rot->mat[1][2] = Rb.mat[1][2];
			break;
	}
#else
	switch(choice){
		case 0: /* yaw, element frame:
			 * (a00 a01 a02)    (1   0   0)   (     a00       co*a01-si*a02  si*a01+co*a02)
			 * (a10 a11 a12)  * (0  co  si) = (     a10       co*a11-si*a12  si*a11+co*a12)
			 * (a20 a21 a22)    (0 -si  co)   (     a20       co*a21-si*a22  si*a21+co*a22)
			 */
			Rb.mat[0][1] = co*rot->mat[0][1]-si*rot->mat[0][2]; Rb.mat[0][2] = si*rot->mat[0][1]+co*rot->mat[0][2];
			Rb.mat[1][1] = co*rot->mat[1][1]-si*rot->mat[1][2]; Rb.mat[1][2] = si*rot->mat[1][1]+co*rot->mat[1][2];
			Rb.mat[2][1] = co*rot->mat[2][1]-si*rot->mat[2][2]; Rb.mat[2][2] = si*rot->mat[2][1]+co*rot->mat[2][2];
			rot->mat[0][1] = Rb.mat[0][1]; rot->mat[0][2] = Rb.mat[0][2];
			rot->mat[1][1] = Rb.mat[1][1]; rot->mat[1][2] = Rb.mat[1][2];
			rot->mat[2][1] = Rb.mat[2][1]; rot->mat[2][2] = Rb.mat[2][2];
			break;
		case 1: /* pitch
			 * (a00 a01 a02)   ( co  0  si)   (co*a00-si*a02    a01    si*a00+co*a02)
			 * (a10 a11 a12) * ( 0   1   0) = (co*a10-si*a12    a11    si*a10+co*a12)
			 * (a20 a21 a22)   (-si  0  co)   (co*a20-si*a22    a21    si*a20+co*a22)
			 */
			Rb.mat[0][0] = co*rot->mat[0][0]-si*rot->mat[0][2]; Rb.mat[0][2] = si*rot->mat[0][0]+co*rot->mat[0][2];
			Rb.mat[1][0] = co*rot->mat[1][0]-si*rot->mat[1][2]; Rb.mat[1][2] = si*rot->mat[1][0]+co*rot->mat[1][2];
			Rb.mat[2][0] = co*rot->mat[2][0]-si*rot->mat[2][2]; Rb.mat[2][2] = si*rot->mat[2][0]+co*rot->mat[2][2];
			rot->mat[0][0] = Rb.mat[0][0]; rot->mat[0][2] = Rb.mat[0][2];
			rot->mat[1][0] = Rb.mat[1][0]; rot->mat[1][2] = Rb.mat[1][2];
			rot->mat[2][0] = Rb.mat[2][0]; rot->mat[2][2] = Rb.mat[2][2];
			break;
		case 2: /* roll, element frame:
			 * (a00 a01 a02)  ( co  si  0)   (co*a00-si*a01  si*a00+co*a01  a02)
			 * (a10 a11 a12) *(-si  co  0) = (co*a10-si*a11  si*a10+co*a11  a12)
			 * (a20 a21 a22)  ( 0   0   1)   (co*a20-si*a21  si*a20+co*a21  a22)
			 */
			Rb.mat[0][0] = co*rot->mat[0][0]-si*rot->mat[0][1]; Rb.mat[0][1] = si*rot->mat[0][0]+co*rot->mat[0][1];
			Rb.mat[1][0] = co*rot->mat[1][0]-si*rot->mat[1][1]; Rb.mat[1][1] = si*rot->mat[1][0]+co*rot->mat[1][1];
			Rb.mat[2][0] = co*rot->mat[2][0]-si*rot->mat[2][1]; Rb.mat[2][1] = si*rot->mat[2][0]+co*rot->mat[2][1];
			rot->mat[0][0] = Rb.mat[0][0]; rot->mat[0][1] = Rb.mat[0][1];
			rot->mat[1][0] = Rb.mat[1][0]; rot->mat[1][1] = Rb.mat[1][1];
			rot->mat[2][0] = Rb.mat[2][0]; rot->mat[2][1] = Rb.mat[2][1];
			break;
	}
#endif
	element->dipole=element->rot*element->MyType->initial_dipole;
	// charges
	for(unsigned int j=0; j<element->MyType->nr_charges; j++) element->q_pos[j]=element->rot*element->MyType->q_pos[j];
	
	return theta*theta;
}

double MC_Elements::calc_charge_interaction()
{
	//Zero energy before calculation
	double Ves = 0.0;
	double selfVes = 0.0; // need to have a separate self ES energy (otherwise 1/n^2 is applied twice)
	double V;
	Vec3 rmu;
	for (unsigned int kk = 0; kk < nr_groups; kk++){
		selfInteraction(Groups[kk],Groups[kk]->VES,NULL,false); // don't calculate VLJ, but calculate charge interaction
		selfVes+=Groups[kk]->VES+0.5*Groups[kk]->Type->charge*Groups[kk]->Type->charge*rf_neutralization*configuration->in2;
	}
	bool calc_dipole;
	Vec3 position;
	Element_Group* kk_group;
	Vec3 RF, RFmirror; // reaction field
	double RFp, RFpmirror; // reaction field potential (interaction starts from charge)
	Vec3 Emu, Emirror;
	double phiq;
	bool updateRF=(fabs(rf_correction)>EPS);
	unsigned int nr_charges;
	double kk_charge=0.0;
	unsigned int start_n=0;
	unsigned int end_n=0;
	unsigned int start;
	unsigned int i_kk, end_kk;
	unsigned int i_g, neighbor, end_nr;
	Element_Group* n_group;
	unsigned int index;
	double* energy;
	for(unsigned int jj = 0; jj<N-1; jj++){
		start=jj*configuration->max_neighbors;
		if(jj<nr_groups){
			kk_group=Groups[jj];
			i_kk=kk_group->elements[0]; // first group element index
			end_kk=i_kk+kk_group->nr_elements;
			start=kk_group->number*configuration->max_neighbors;
			kk_charge=kk_group->Type->charge*rf_neutralization;
			start_n=kk_group->nr_neighbors_above;
			end_n=kk_group->nr_neighbors;
		} else{
			i_kk=jj+nr_N_to_kk;
			end_kk=i_kk+1;
			kk_group=NULL;
		}
		for(unsigned int kk = i_kk; kk < end_kk; kk++){ // Loop over 2-body interactions (interaction partners start at kk+1, hence need only go over nr_elements-1)
			//Calculate the distance array
			get_distances_above_index(kk); // fills vector array Rmus with distances from element kk to all other elements
			calc_dipole=Elements[kk].MyType->hasmu;
			if(calc_dipole){
				if(configuration->offctrmu){
					position=Elements[kk].rot*Elements[kk].MyType->mu_pos; // correction in case of offcenter dipoles
				} else position.V3Zeros();
			}
			nr_charges=Elements[kk].MyType->nr_charges;
			if(!kk_group){
				start=(nr_groups+kk-nr_group_elements)*configuration->max_neighbors;
				kk_charge=Elements[kk].MyType->charge*rf_neutralization;
				start_n=Elements[kk].nr_neighbors_above;
				end_n=Elements[kk].nr_neighbors;
			}
			// Loop over all 2-body neighbor interactions
			for(unsigned int j=start_n; j<end_n; j++){
				index=start+j;
				neighbor=neighbors[index];
				energy=&pair_energies[energy_indices[index]+1]; // go straight to ES energies
				if(kk==i_kk) *energy=0.0;
				if(neighbor<nr_groups){
					n_group=Groups[neighbor];
					i_g=n_group->elements[0]; // first group element index
					end_nr=i_g+n_group->nr_elements;
				} else{
					i_g=neighbor+nr_N_to_kk;
					end_nr=i_g+1;
				}
				if((Rdist2[i_g]<configuration->escut2) || (!configuration->cut)){
					V=0.0;
					for(unsigned int i=i_g; i<end_nr; i++){
						if(calc_dipole){
							rmu=Rmus[i]-position; // rmu is vector from element kk to element i
							Emu=E_mu(i,rmu,RF,updateRF);
							if(LJwallMirror[i]){ // consider image charges from dielectric walls here
								rmu.vec[0]=Rmirror[i]+position.vec[0];
								Emirror=E_mu(i_g,rmu,RFmirror,updateRF);
								Emirror.vec[0]*=-1.0;
								RFmirror.vec[0]*=-1.0;
								Emu+=Emirror*configuration->LJwall_mirror_factor;
								RF+=RFmirror*configuration->LJwall_mirror_factor;
							}
							V+=(Emu-RF*rf_correction)*Elements[kk].dipole;
						}
						if(nr_charges){
							if(kk==i_kk) V -= kk_charge*Elements[i].MyType->charge;
							for(unsigned int m=0; m<nr_charges; m++){
								rmu=Rmus[i]-Elements[kk].q_pos[m];
								phiq=phi_q(i,rmu,RFp,updateRF);
								if(LJwallMirror[i]){ // consider image charges from dielectric walls here
									if(kk==i_kk) V -= kk_charge*configuration->LJwall_mirror_factor*Elements[i].MyType->q[m];
									rmu.vec[0]=Rmirror[i]+Elements[kk].q_pos[m].vec[0];
									phiq+=phi_q(i,rmu,RFpmirror,updateRF)*configuration->LJwall_mirror_factor;
									RFp+=RFpmirror*configuration->LJwall_mirror_factor;
								}
								V+=(phiq+RFp*rf_correction)*Elements[kk].MyType->q[m];
							}
						}
					}
					*energy += V;
					Ves += V;
				}
			} // done calculating electrostatics and reaction field as well as Lennard-Jones
		}
	}
	return Ves*configuration->in2+selfVes;
}

/*!
 * Function for changing E and T during the run, which evaluates every cycle.
 * Changesteps is an array of the steps where changes take place, trackchange
 * is a counter that moves forwards one every time a change is made. At minimum,
 * this function will set the values before the simulation starts, so all of the
 * arrays must have at least one element. This function also controls simulated
 * annealing at the beginning of the simulation if it is used.
 * Inputs:
 * 	step - current simulation step
 * Return value is a flag to indicate whether or not dipole-field interactions need to be recalculated
 */
bool MC_Elements::setEandT(unsigned int step)
{
	// VmuE must be recalculated if E changes
	bool recalcuE = false;
	double e_scale = MV_per_m_to_perg_per_Debye;
	if((configuration->realfield) && (step>=configuration->randsteps || configuration->rand_Efield)){
		if(configuration->dyneps) e_scale *= (3.0*epsRF)/(2.0*epsRF+configuration->n2);
			else e_scale *= (2.0*configuration->epsilon + epsRF)/(2.0*configuration->epsilon + configuration->n2);
		Evec = configuration->Efield*e_scale;
		recalcuE = true;
	}
	// Simulated annealing for randomization phase
	// This is currently a simple model with the initial temperature dependent on the number of randsteps - LEJ 10/23/2009
	// Updated to use two phases, one in the initial randomization, and one right when electrostatics turn on - LEJ 02/09/2010
	
	// Critical steps
	unsigned int cr12 = configuration->randsteps/2;
	unsigned int cr2 = configuration->randsteps*2;
	// Phase 1 - randomization with no electrostatics. Lasts for first half of randsteps.
	if((configuration->anneal) && (step<=cr12)){
		// Temperature increment in K.
		double tinc = 0.20;
		// Enable annealing, reduce temperature, or disable annealing, depending on step.
		if(step==0){
			cout << "Using simulated annealing for pre-minimization.\n";
			cout << "Initial annealing temp is " << configuration->T + (double)(configuration->randsteps/10) << " K\nTemperature increment is " << tinc << " K\n\n";
		} else{
			if(step==1){ // Enable annealing and expand move sizes
				configuration->kT += kB*(double)(configuration->randsteps/10);
				configuration->beta = 1.0/configuration->kT;
				configuration->maxtrans *= 1.6;
				configuration->maxrot *= 1.6;
				adjust_cut2plus(configuration);
			} else{
				if(step<cr12){ // Reduce the temperature by tinc K per step during the annealing period
					configuration->kT -= tinc*kB;
					configuration->beta = 1.0/configuration->kT;
				} else{
					if(step==cr12){ // Reset move sizes once annealing period is over
						configuration->kT = kB*configuration->T;
						configuration->beta = 1.0/configuration->kT;
						configuration->maxtrans *= 0.625;
						configuration->maxrot *= 0.625;
						adjust_cut2plus(configuration);
						cout << "\nSystem has cooled to final temp of " << configuration->T << " K\n";
					}
				}
			}
		}
	} else if((configuration->anneal) && (step<=cr2)){
		// Phase 2 - randomization with electrostatics. Lasts for randsteps after randomization phase finishes.
		// Functionally identical to the first simulated annealing phase, but cooling slower for longer
		// For this function, the simulation must not start averaging until at least step 2*randsteps
		// Temperature increment in K.
		double tinc = 0.10;
		// Enable annealing, reduce temperature, or disable annealing, depending on step.
		if(step==configuration->randsteps){
			cout << "\nUsing simulated annealing on for next " << cr2 << " steps.\n";
			cout << "Initial annealing temp is " << configuration->T + (double)(configuration->randsteps/10) << " K\nTemperature increment is " << tinc << " K\n";
		} else{
			if(step == (configuration->randsteps + 1)){
				configuration->kT += kB*(double)(configuration->randsteps/10);
				configuration->beta = 1.0/configuration->kT;
				configuration->maxtrans *= 1.6;
				configuration->maxrot *= 1.6;
				adjust_cut2plus(configuration);
			} else{
				if((step>configuration->randsteps) && (step<cr2)){
					configuration->kT -= tinc*kB;
					configuration->beta = 1.0/configuration->kT;
				} else{
					if(step==cr2){
						configuration->kT = kB*configuration->T;
						configuration->beta = 1.0/configuration->kT;
						configuration->maxtrans *= 0.625;
						configuration->maxrot *= 0.625;
						adjust_cut2plus(configuration);
						cout << "\nSystem has cooled to final temp of " << configuration->T << " K\n";
					}
				}
			}
		}
	}
	// DEBUG - arrays have been converted to single-element variables! BEWARE! -LEJ
	// Make change if on a step flagged for a change
	if((step==configuration->changesteps) || ((int)step==configuration->resume_step)){ // aka step==0
		// Update kT
		double oldkT = configuration->kT;
		configuration->kT = (kB*configuration->T);
		configuration->beta = 1.0/configuration->kT;
		if ((fabs(oldkT-configuration->kT)>EPS) && (step!=0)) cout << "T is now " << configuration->T << " K\n";
		
		Vec3 oldE = Evec;
		if(!configuration->realfield) Evec=configuration->Efield*MV_per_m_to_perg_per_Debye;
		if(((oldE-Evec).V3Norm()>EPS) && (step!=0)) recalcuE=true;
		// Increment change counter
		configuration->trackchange++;
	}
	return recalcuE;
}

bool MC_Elements::selfInteraction(Element_Group* group, double &VES, double* VLJ, bool skip_charges)
{
	VES=0.0;
	double VRF=0.0;
	double cLJl=LJlambda;
	LJlambda*=selfLJlambda;
	if(VLJ) *VLJ=0.0;
	
	unsigned int bond_range=group->Type->bond_range;
	bool calculation;
	
	__uint32_t* range_field=group->Type->range_field;
	
	unsigned int element_nr;
	
	if(bond_range>0){
		Element* i_element;
		Element* kk_element;
		
		unsigned int g_i, g_kk;
		
		bool updateRF=((fabs(rf_correction)>EPS) && configuration->selfrf);
		Vec3 RF(0.0);
		double RFp=0.0;
		Vec3 rmu;
		for(unsigned int kk=0; kk<group->nr_elements-1; kk++){
			g_kk=group->elements[kk];
			kk_element=&Elements[g_kk];
			
			bool calc_dipole=kk_element->MyType->hasmu;
			Vec3 position(0.0);
			if(calc_dipole){
				if(configuration->offctrmu) position=kk_element->rot*kk_element->MyType->mu_pos; // correction in case of offcenter dipoles
			}
			unsigned int nr_charges=kk_element->MyType->nr_charges;
			
			Vec3 Emu;
			double phiq=0.0;
			double factor;
			bool calc_VLJ=(VLJ);
			
			unsigned int m=0;
			do{
				if(calc_dipole){
					RF.vec[0]=0.0; RF.vec[1]=0.0; RF.vec[2]=0.0;
					Emu.vec[0]=0.0; Emu.vec[1]=0.0; Emu.vec[2]=0.0;
				} else{
					if(nr_charges>0) position=kk_element->q_pos[m];
					RFp=0.0;
					phiq=0.0;
				}
				for(unsigned int i=kk+1; i<group->nr_elements; i++){
					// addressing is i*(i-1)/2+kk for kk<i
					element_nr=((i*(i-1))>>1)+kk;
					__uint32_t storage=range_field[element_nr>>(5-SPECIAL_DISTANCES_BIT_EXPONENT)];
					unsigned int shift=((element_nr&31)<<SPECIAL_DISTANCES_BIT_EXPONENT); // getting SPECIAL_DISTANCES bits (nr&31==nr%32) -> shift 1 left by nr%32
					// magic happens here: it's a smoke screen (essentially) shifting 1111 (for four bits) to respective location to only get those bits, then shift them down for readout
					__uint32_t special=(storage&(SPECIAL_DISTANCES_BIT_SCREEN<<shift))>>shift;
					calculation=(special>0) || (bond_range==1);
#if DEBUG_LEVEL>4
					if(special>1) cout << "Factor for interaction between " << kk+1 << " and " << i+1 << " (" << group->Type->range_factors[special-2].distance << " bonds away): " << group->Type->range_factors[special-2].factor << "\n";
#endif
					factor=1.0;
					if(special>1) factor=group->Type->range_factors[special-2].factor;
					if(calculation){
						g_i=group->elements[i];
						i_element=&Elements[g_i];
						//Calculate distances
						Rmus[g_i]=i_element->center-kk_element->center; // rmu is vector from element kk to element i
						rmu=Rmus[g_i];
						Rdist2LJ[g_i]=rmu*rmu;
						if(calc_VLJ){
							// VdW interactions (pick whether you want to use the soft-sphere, simple, extended, or Touch LJ function)
							if((Rdist2LJ[g_i]<configuration->ljcut2) || (!configuration->cut_internal)) *VLJ += factor*(this->*LJpair)(g_kk,g_i);
						}
						if(!skip_charges){
							rmu-=position;
							double cut_dist2=rmu.vec[0]*rmu.vec[0]+rmu.vec[1]*rmu.vec[1]+rmu.vec[2]*rmu.vec[2];
							//Bomb out if a zero distance is detected to avoid dividing by zero
							if(cut_dist2 < EPS*EPS) return false;
							if((cut_dist2<configuration->escut2) || (!configuration->cut_internal)){
								if(calc_dipole){
									Emu+=self_E_mu(g_i,rmu,RF,updateRF)*factor;
								} else if(nr_charges>0) phiq+=self_phi_q(g_i,rmu,RFp,updateRF)*factor;
							}
						}
					}
				}
				if(!skip_charges){
					if(calc_dipole){ // around a dipole
						VES += kk_element->dipole*Emu;
						if(updateRF) VRF -= RF*kk_element->dipole;
					} else{ // around a charge
						if(nr_charges>0){
							VES += kk_element->MyType->q[m]*phiq;
							if(updateRF) VRF += RFp*kk_element->MyType->q[m];
						}
						m++;
					}
				}
				// done with dipole and VLJ (both done in first round)
				calc_VLJ=false;
				calc_dipole=false;
			} while(((calc_dipole) || (m<nr_charges)) && (!skip_charges));
			// done calculating electrostatics, VLJ interactions, and reaction field
		}
	}
	if(!skip_charges){
		VES += VRF*rf_correction-group->Type->charge*group->Type->charge*rf_neutralization;
		VES*=configuration->in2;
	}
	LJlambda=cLJl;
	return true;
}

/*!
 * OneStep computes a single step of all of the interactions between one particle (kk from main loop) and every other particle (i) in the system.
 * Inputs:
 * 	kk - the index of the active ellipsoid (ellipsoid being moved)
 * 	Ves - passthrough for mu-mu energy (actually an output)
 * 	VLJ - passthrough for Lennard-Jones energy (actually an output)
 * 	skip_charges - flag used to control whether charge/dipole interactions are calculated
 */
bool MC_Elements::OneStep(const unsigned int kk, double &Ves, double &VLJ, bool skip_charges)
{
	//Calculate the distance array
	bool success=get_distances(kk); // fills vector array Rmus with distances from element kk to all other elements
	if(!success) return success;
	
	//Zero energies before calculation
	Ves = 0.0;
	VLJ = 0.0;
	
	unsigned int start=(nr_groups+kk-nr_group_elements)*configuration->max_neighbors;
	Element_Group* n_group;
	unsigned int i_g, neighbor, end_nr;
	double energy;
	// Loop over all 2-body neighbor interactions
	for(unsigned int index=start; index<start+Elements[kk].nr_neighbors; ++index){
		neighbor=neighbors[index];
		if(neighbor<nr_groups){
			n_group=Groups[neighbor];
			i_g=n_group->elements[0]; // first group element index
			end_nr=i_g+n_group->nr_elements;
		} else{
			i_g=neighbor+nr_N_to_kk;
			end_nr=i_g+1;
		}
		// VdW interactions (pick whether you want to use the soft-sphere, simple, extended, or Touch LJ function)
		energy=0.0;
		for( ; i_g<end_nr; ++i_g) if((Rdist2LJ[i_g]<configuration->ljcut2) || (!configuration->cut)) energy += (this->*LJpair)(kk,i_g);
		VLJ+=energy;
		pair_energies[energy_indices[index]]=energy;
	}
	if(!skip_charges){
		bool calc_dipole=Elements[kk].MyType->hasmu;
		Vec3 position(0.0);
		if(calc_dipole){
			if(configuration->offctrmu) position=Elements[kk].rot*Elements[kk].MyType->mu_pos; // correction in case of offcenter dipoles
		}
		unsigned int nr_charges=Elements[kk].MyType->nr_charges;
		double kk_charge=Elements[kk].MyType->charge*rf_neutralization;
		bool updateRF=(fabs(rf_correction)>EPS);
		Vec3 rmu;
		Vec3 RF, RFmirror; // reaction field
		Vec3 Emu, Emirror;
		double RFp, RFpmirror; // reaction field potential
		double phiq;
		
		for(unsigned int index=start; index<start+Elements[kk].nr_neighbors; ++index){
			neighbor=neighbors[index];
			if(neighbor<nr_groups){
				n_group=Groups[neighbor];
				i_g=n_group->elements[0]; // first group element index
				end_nr=i_g+n_group->nr_elements;
			} else{
				i_g=neighbor+nr_N_to_kk;
				end_nr=i_g+1;
			}
			energy=0.0;
			if((Rdist2[i_g]<configuration->escut2) || (!configuration->cut)){
				for( ; i_g<end_nr; ++i_g){
					if(calc_dipole){
						rmu=Rmus[i_g]-position; // rmu is vector from element kk to element i
						Emu=E_mu(i_g,rmu,RF,updateRF);
						if(LJwallMirror[i_g]){ // consider image charges from dielectric walls here
							rmu.vec[0]=Rmirror[i_g]+position.vec[0];
							Emirror=E_mu(i_g,rmu,RFmirror,updateRF);
							Emirror.vec[0]*=-1.0;
							RFmirror.vec[0]*=-1.0;
							Emu+=Emirror*configuration->LJwall_mirror_factor;
							RF+=RFmirror*configuration->LJwall_mirror_factor;
						}
						energy+=(Emu-RF*rf_correction)*Elements[kk].dipole;
					}
					if(nr_charges){
						energy -= Elements[i_g].MyType->charge*kk_charge;
						for(unsigned int m=0; m<nr_charges; m++){
							rmu=Rmus[i_g]-Elements[kk].q_pos[m];
							phiq=phi_q(i_g,rmu,RFp,updateRF);
							if(LJwallMirror[i_g]){ // consider image charges from dielectric walls here
								energy -= kk_charge*configuration->LJwall_mirror_factor*Elements[i_g].MyType->q[m];
								rmu.vec[0]=Rmirror[i_g]+Elements[kk].q_pos[m].vec[0];
								phiq+=phi_q(i_g,rmu,RFpmirror,updateRF)*configuration->LJwall_mirror_factor;
								RFp+=RFpmirror*configuration->LJwall_mirror_factor;
							}
							energy+=(phiq+RFp*rf_correction)*Elements[kk].MyType->q[m];
						}
					}
				}
			}
			Ves+=energy;
			pair_energies[energy_indices[index]+1]=energy;
		}
	}
	// done calculating electrostatics, VLJ interactions, and reaction field
	
	Ves*=configuration->in2;
	return true;
}

/// OneStep (in case kk belongs to a group)
bool MC_Elements::OneStepGroup(const unsigned int kk, double &Ves, double &VLJ, bool skip_charges)
{
	//Calculate the distance array
	bool success=get_distances(kk); // fills vector array Rmus with distances from element kk to all other elements
	if(!success) return false;
	
	//Zero energies before calculation
	Ves = 0.0;
	VLJ = 0.0;
	
	Element_Group* kk_group=Elements[kk].group; // element kk has to be part of a group -- otherwise OneStepGroup would not have been called
	Element_Group* n_group;
	
	bool start_kk=(kk==kk_group->elements[0]);
	unsigned int i_g, neighbor, end_nr;
	double energy;
	unsigned int start=kk_group->number*configuration->max_neighbors;
	// Loop over all 2-body neighbor interactions
	for(unsigned int index=start; index<start+kk_group->nr_neighbors; ++index){
		neighbor=neighbors[index];
		if(start_kk) pair_energies[energy_indices[index]]=0.0;
		if(neighbor<nr_groups){
			n_group=Groups[neighbor];
			i_g=n_group->elements[0]; // first group element index
			end_nr=i_g+n_group->nr_elements;
		} else{
			i_g=neighbor+nr_N_to_kk;
			end_nr=i_g+1;
		}
		//Calculate distances
		// VdW interactions (pick whether you want to use the soft-sphere, simple, extended, or Touch LJ function)
		energy=0.0;
		for( ; i_g<end_nr; ++i_g) if((Rdist2LJ[i_g]<configuration->ljcut2) || (!configuration->cut)) energy += (this->*LJpair)(kk,i_g);
		VLJ += energy;
		pair_energies[energy_indices[index]] += energy;
	}
	if(!skip_charges){
		bool calc_dipole=Elements[kk].MyType->hasmu;
		Vec3 position(0.0);
		if(calc_dipole){
			if(configuration->offctrmu) position=Elements[kk].rot*Elements[kk].MyType->mu_pos; // correction in case of offcenter dipoles
		}
		unsigned int nr_charges=Elements[kk].MyType->nr_charges;
		bool updateRF=(fabs(rf_correction)>EPS);
		double kk_charge=kk_group->Type->charge;
		kk_charge*=rf_neutralization;
		Vec3 rmu;
		Vec3 RF, RFmirror; // reaction field
		Vec3 Emu, Emirror;
		double RFp, RFpmirror; // reaction field potential
		double phiq;
		
		for(unsigned int index=start; index<start+kk_group->nr_neighbors; ++index){
			neighbor=neighbors[index];
			if(start_kk) pair_energies[energy_indices[index]+1]=0.0;
			if(neighbor<nr_groups){
				n_group=Groups[neighbor];
				i_g=n_group->elements[0]; // first group element index
				end_nr=i_g+n_group->nr_elements;
			} else{
				i_g=neighbor+nr_N_to_kk;
				end_nr=i_g+1;
			}
			if((Rdist2[i_g]<configuration->escut2) || (!configuration->cut)){
				energy=0.0;
				for( ; i_g<end_nr; ++i_g){
					if(calc_dipole){
						rmu=Rmus[i_g]-position; // rmu is vector from element kk to element i
						Emu=E_mu(i_g,rmu,RF,updateRF);
						if(LJwallMirror[i_g]){ // consider image charges from dielectric walls here
							rmu.vec[0]=Rmirror[i_g]+position.vec[0];
							Emirror=E_mu(i_g,rmu,RFmirror,updateRF);
							Emirror.vec[0]*=-1.0;
							RFmirror.vec[0]*=-1.0;
							Emu+=Emirror*configuration->LJwall_mirror_factor;
							RF+=RFmirror*configuration->LJwall_mirror_factor;
						}
						energy+=(Emu-RF*rf_correction)*Elements[kk].dipole;
					}
					if(nr_charges){
						if(start_kk) energy -= kk_charge*Elements[i_g].MyType->charge;
						for(unsigned int m=0; m<nr_charges; m++){
							rmu=Rmus[i_g]-Elements[kk].q_pos[m];
							phiq=phi_q(i_g,rmu,RFp,updateRF);
							if(LJwallMirror[i_g]){ // consider image charges from dielectric walls here
								if(start_kk) energy -= kk_charge*configuration->LJwall_mirror_factor*Elements[i_g].MyType->q[m];
								rmu.vec[0]=Rmirror[i_g]+Elements[kk].q_pos[m].vec[0];
								phiq+=phi_q(i_g,rmu,RFpmirror,updateRF)*configuration->LJwall_mirror_factor;
								RFp+=RFpmirror*configuration->LJwall_mirror_factor;
							}
							energy+=(phiq+RFp*rf_correction)*Elements[kk].MyType->q[m];
						}
					}
				}
				Ves += energy;
				pair_energies[energy_indices[index]+1] += energy;
			}
		}
	}
	// done calculating electrostatics, VLJ interactions, and reaction field
	//Scale dipole-dipole interactions by background dielectric constant (1 unless implicit polarizability is used)
	Ves *= configuration->in2;
	return true;
}

inline double MC_Elements::NpT_step(double* Vstep, double* save_groupVG, unsigned int step, unsigned int eqsteps, bool &adjust_internal)
{
	NpT_tries++;
	double tmovtemp=0.0;
	double NewVstep[NUM_V_STORE];
	double initialV=configuration->V;
#if DEBUG_LEVEL>2
	cout << "NpT step " << step+1 << ": " << configuration->V;
#endif
	// randomly change volume
	double volscale;
#ifdef USE_CMWC4096
	if(configuration->NpT_lnV_change){
		volscale=exp(log(initialV)+(1.0-2.0*ranQ())*configuration->maxdV)/initialV;
	} else volscale=(initialV+(1.0-2.0*ranQ())*configuration->maxdV)/initialV;
#else
	if(configuration->NpT_lnV_change){
		volscale=exp(log(initialV)+(1.0-2.0*ran2(*configuration->idum))*configuration->maxdV)/initialV;
	} else volscale=(initialV+(1.0-2.0*ran2(*configuration->idum))*configuration->maxdV)/initialV;
#endif
	double rfc_save=configuration->rfcut;
	double esc2_save=configuration->escut2;
	double OldV = 0.0;
	if((step>=configuration->randsteps) && (configuration->cut && configuration->NpT_adjust_EScut)){ // in case electrostatics are on, adjust cutoff and recalculate energies
#if DEBUG_LEVEL>2
		cout << "adjusting ES cutoff: " << configuration->rfcut << " -> ";
#endif
		// escut2, the square of the electrostatics cutoff, is the only thing used during runtime (for distances, r*r without the additional sqrt saves time)
		configuration->rfcut=configuration->escut*cbrt(volscale*configuration->V/configuration->targetV);
		configuration->escut2=qqpwr(configuration->rfcut,2);
		adjust_cut2plus(configuration);
		find_neighbors();
#if DEBUG_LEVEL>2
		cout << configuration->rfcut << "\n(";
		for(unsigned int jj=0; jj<NUM_V_STORE; jj++) cout << NewVstep[jj] << "\t";
#endif
		recalcall(NewVstep,false); // dipole/charges are calculated
#if DEBUG_LEVEL>2
		for(unsigned int jj=0; jj<NUM_V_STORE; jj++) cout << NewVstep[jj] << "\t";
		cout << ")\n";
#endif
		for(unsigned int jj=0; jj<NUM_V_STORE; jj++) OldV += NewVstep[jj];
	} else{
		// old energy is in Vstep
		for(unsigned int jj=0; jj<NUM_V_STORE; jj++) OldV += Vstep[jj];
	}
	double scale=1.0;
	if(configuration->LJwall_calc && configuration->LJwall_fixed){
		scale=sqrt(volscale); // how much each dimension needs to change
		switch(configuration->latticetype){
			case 1:
			case 2: // Rectangular box
				configuration->boxlength[0]=2.0*configuration->LJwall_xm;
				configuration->inv_boxlength[0]=1.0/configuration->boxlength[0];
				configuration->boxlength[1]*=scale;
				configuration->inv_boxlength[1]=1.0/configuration->boxlength[1];
				configuration->boxlength[2]*=scale;
				configuration->inv_boxlength[2]=1.0/configuration->boxlength[2];
				break;
			default: // Default condition
				cout << "Cannot change volume for LJ wall configuration. Lattice type " << configuration->latticetype << " not recognized. Exiting.\n";
				exit(1);
		}
		configuration->nndist*=scale;
		configuration->V*=volscale;
	} else update_volume(configuration,configuration->V*volscale);
	if(configuration->LJwall_calc && !configuration->LJwall_fixed) configuration->LJwall_xm=0.5*configuration->boxlength[0];
	double dV = configuration->V-initialV;
	Vec3 delta_trans;
	// scale center position of independent element centers according to change of volume
	if(configuration->NpT_relative_coords){ // always true for LJ wall
		if(!configuration->LJwall_fixed) scale=cbrt(volscale); // in case of the fixed LJ wall use the scale factor defined above (sqrt(volscale))
		for(unsigned int i=nr_group_elements; i<nr_elements; i++){
			delta_trans=Elements[i].center*(-1.0);
			if(!Elements[i].fixed){
				if(!configuration->LJwall_fixed) Elements[i].center.vec[0]*=scale;
				Elements[i].center.vec[1]*=scale; Elements[i].center.vec[2]*=scale;
			}
			delta_trans+=Elements[i].center;
//			Elements[i].ds+=delta_trans;
			tmovtemp += delta_trans*delta_trans;
		}
		// scale center position of groups (move whole group)
		for(unsigned int i=0; i<nr_groups; i++){
			delta_trans=Groups[i]->center*(-1.0);
			if(!Groups[i]->fixed){
				if(!configuration->LJwall_fixed) Groups[i]->center.vec[0]*=scale;
				Groups[i]->center.vec[1]*=scale; Groups[i]->center.vec[2]*=scale;
			}
			delta_trans+=Groups[i]->center;
			tmovtemp += (delta_trans*delta_trans)*Groups[i]->nr_elements;
			for(unsigned int j=0; j<Groups[i]->nr_elements; j++) Elements[Groups[i]->elements[j]].center+=delta_trans;
		}
		if(!configuration->NpT_adjust_EScut) find_neighbors();
	} else find_neighbors();
	if(configuration->LJwall_calc){
		for(unsigned int i=0; i<nr_groups; i++) save_groupVG[i]=Groups[i]->VG;
	}
	// calculate energy change
	bool success=recalcall(NewVstep,(configuration->randsteps>step)); // dipole/charges are calculated
	double NewV = 0.0;
	for(unsigned int jj=0; jj<NUM_V_STORE; jj++) NewV += NewVstep[jj];
	// Metropolis Conditional: dE = dU + pdV - (N+1)kT dlnV
	double deltaV = NewV - OldV + configuration->pext*dV; // for NpT ensemble add pdV (pext is in pErg/Angström^3)
#ifdef USE_CMWC4096
	double rand = ranQ();
#else
	double rand = ran2(*configuration->idum);
#endif
#if DEBUG_LEVEL>3
	cout << " (try volume: " << configuration->V << " (dV=" << dV << ", OldV=" << OldV << ", NewV=" << NewV << ", pext*dV=" << configuration->pext*dV << ", (N+1)ln(V'/V)=" << (N+1)*log(configuration->V/initialV) << ", -dE/kT=" << -deltaV*configuration->beta+(N+1)*log(configuration->V/initialV) << ", rand=" << rand << ")";
#endif
	// Metropolis conditional
	if(success && (exp(-deltaV*configuration->beta+configuration->N_NpT*log(volscale)) > rand)){ // for whatever reason log is actually ln in C++
#if DEBUG_LEVEL>3
		cout << " -> used)";
#endif
		// acceptance tracking [vdw preeq eq]
		if(step>=eqsteps){
			accept[2]++;
			NpT_accept[2]++;
		} else{
			if(step>=configuration->randsteps){
				accept[1]++;
				NpT_accept[1]++;
			}else{
				accept[0]++;
				NpT_accept[0]++;
			}
		}
		if(!configuration->NpT_relative_coords){
			// scale center position of independent element centers according to change of volume
			for(unsigned int i=nr_group_elements; i<nr_elements; i++){
//				delta_trans=Elements[i].center*(-1.0);
				apply_PBCs(Elements[i].center);
//				Elements[i].ds+=delta_trans;
			}
			// scale center position of groups (move whole group)
			for(unsigned int i=0; i<nr_groups; i++){
				delta_trans=Groups[i]->center*(-1.0);
				apply_PBCs(Groups[i]->center);
				delta_trans+=Groups[i]->center;
				for(unsigned int j=0; j<Groups[i]->nr_elements; j++) Elements[Groups[i]->elements[j]].center+=delta_trans;
			}
		}
		// Take care of internal energy adjustments if necessary
		if(adjust_internal){
			selfLJlambda=calc_lambda();
			if(fabs(selfLJlambda-1.0)<EPS){
				adjust_internal=false;
				selfLJlambda=1.0;
			}
			recalcall(NewVstep,(configuration->randsteps>step));
			if(!adjust_internal) cout << "\t\t\tInternal LJ dispersion adjustment turned off at " << configuration->V << " Angström\n";
		}
		// Record instantaneous value of total energy
		for(unsigned int jj=0; jj<NUM_V_STORE; jj++) Vstep[jj] = NewVstep[jj];
	} else{ // Revert things back if move is rejected
#if DEBUG_LEVEL>3
		cout << " -> not used)";
#endif
		configuration->rfcut=rfc_save;
		configuration->escut2=esc2_save;
		adjust_cut2plus(configuration);
		double invscale=1.0/scale;
		if(configuration->LJwall_calc && configuration->LJwall_fixed){
			switch(configuration->latticetype){
				case 1:
				case 2: // Rectangular box
					configuration->boxlength[0]=2.0*configuration->LJwall_xm;
					configuration->inv_boxlength[0]=1.0/configuration->boxlength[0];
					configuration->boxlength[1]*=invscale;
					configuration->inv_boxlength[1]=1.0/configuration->boxlength[1];
					configuration->boxlength[2]*=invscale;
					configuration->inv_boxlength[2]=1.0/configuration->boxlength[2];
					break;
				default: // Default condition
					cout << "Cannot change volume for LJ wall configuration. Lattice type " << configuration->latticetype << " not recognized. Exiting.\n";
					exit(1);
			}
			configuration->nndist*=invscale;
			configuration->V=initialV;
		} else update_volume(configuration,initialV);
		if(configuration->LJwall_calc && !configuration->LJwall_fixed) configuration->LJwall_xm=0.5*configuration->boxlength[0];
		tmovtemp=0.0;
		if(configuration->NpT_relative_coords){
			// scale center position of independent element centers according to change of volume
			for(unsigned int i = nr_group_elements; i < nr_elements; i++){
				if(!Elements[i].fixed){
					if(!configuration->LJwall_fixed) Elements[i].center.vec[0]*=invscale;
					Elements[i].center.vec[1]*=invscale; Elements[i].center.vec[2]*=invscale;
				}
			}
			// scale center position of groups (move whole group)
			for(unsigned int i=0; i<nr_groups; i++){
				delta_trans=Groups[i]->center*(-1.0);
				if(!Groups[i]->fixed){
					if(!configuration->LJwall_fixed) Groups[i]->center.vec[0]*=invscale;
					Groups[i]->center.vec[1]*=invscale; Groups[i]->center.vec[2]*=invscale;
				}
				delta_trans+=Groups[i]->center;
				for(unsigned int j=0; j<Groups[i]->nr_elements; j++) Elements[Groups[i]->elements[j]].center+=delta_trans;
			}
		}
		if(configuration->LJwall_calc){
			for(unsigned int i=0; i<nr_groups; i++) Groups[i]->VG=save_groupVG[i];
		}
		if(configuration->NpT_relative_coords && !configuration->NpT_adjust_EScut) find_neighbors();
		recalcall(Vstep,(configuration->randsteps>step));
	}
#if DEBUG_LEVEL>3
	cout << "\n";
#endif
	return tmovtemp;
}

/// Calculate potential energy for one group
double MC_Elements::GroupPotentialEnergy(Element_Group* group)
{
#if DEBUG_LEVEL>3
	cout << "Potential calculation\n";
#endif
	double Vg=0.0;
	Interaction_Potential* potential;
	Vec3 a,b,delta;
	double theta;
	Interaction* link;
	unsigned int i, k, l, m;
	double r, cosine, sine;
	for(unsigned int j=0; j<group->nr_potentials; j++){
		potential=group->potentials[j];
		switch(potential->type){
			case stretch_potential:
						// get element number
						i=group->elements[potential->partners[0]];
						k=group->elements[potential->partners[1]];
						// distance vector between elements
						delta=(Elements[i].center-Elements[k].center);
						// correct for individual connection starting points
						link=&Elements[i].interactions[potential->links[0*2+1]]; // connection start on A of A->B
						if(link->location) delta+=*(link->location);
						link=&Elements[k].interactions[potential->links[1*2+0]]; // connection start on B of B->A
						if(link->location) delta-=*(link->location); // minus because delta=A-B
						// calculate potential: V=k*(r-r0)^2
						Vg+=potential->parameters[0]*qqpwr(delta.V3Norm()-potential->parameters[1],2);
#if DEBUG_LEVEL>3
						cout << "stretch: " << Vg << ", " << delta.V3Norm() << "\n";
#endif
						break;
			case bend_potential:
						// get element numbers
						i=group->elements[potential->partners[0]]; // element A
						k=group->elements[potential->partners[1]]; // center element
						l=group->elements[potential->partners[2]]; // element B
						// distances between center and other two elements
						a=(Elements[i].center-Elements[k].center); // vector center->A
						b=(Elements[l].center-Elements[k].center); // vector center->B
#if DEBUG_LEVEL>3
						cout << "bend: (<" << a.vec[0] << ", " << a.vec[1] << ", " << a.vec[2] << ">, <" << b.vec[0] << ", " << b.vec[1] << ", " << b.vec[2] << ">), ";
#endif
						// correct for individual connection starting points
						link=&Elements[k].interactions[potential->links[1*3+0]]; // connection start on center of center->A
						if(link->location) a-=*(link->location); // minus because a=A-center
						link=&Elements[k].interactions[potential->links[1*3+2]]; // connection start on center of center->B
						if(link->location) a-=*(link->location); // minus because b=B-center
						link=&Elements[i].interactions[potential->links[0*3+1]]; // connection start on A of A->center
						if(link->location) a+=*(link->location);
						link=&Elements[l].interactions[potential->links[2*3+1]]; // connection start on B of B->center
						if(link->location) b+=*(link->location);
						// calculate theta (this calculation saves one sqrt ...)
#if DEBUG_LEVEL>3
						cout << "<" << a.vec[0] << ", " << a.vec[1] << ", " << a.vec[2] << ">, <" << b.vec[0] << ", " << b.vec[1] << ", " << b.vec[2] << ">, ";
#endif
						cosine=(a*b)/sqrt((a*a)*(b*b));
						if(cosine>1.0){ // the Cauchy-Schwartz inequality can be broken by computers ...
							theta=0.0;
						} else{
							if(cosine<-1.0){
								theta=pi;
							} else theta=acos(cosine);
						}
						// calculate potential: V=k*(theta-theta0)^2
						Vg+=potential->parameters[0]*qqpwr(theta-potential->parameters[1],2);
#if DEBUG_LEVEL>3
						cout << Vg << ", " << cosine << ", " << theta << "\n";
#endif
						break;
			case dihedral_potential:
			case improper_dihedral_potential:
						// get element numbers (elements should be sorted in row)
						i=group->elements[potential->partners[0]]; // element A_1
						k=group->elements[potential->partners[1]]; // element A_2
						l=group->elements[potential->partners[2]]; // element B_1
						m=group->elements[potential->partners[3]]; // element B_2
						// distances between center and other two elements
						a=(Elements[k].center-Elements[i].center); // vector A_1->A_2
						b=(Elements[m].center-Elements[l].center); // vector B_1->B_2
						delta=(Elements[l].center-Elements[k].center); // vector A_2->B_1 (rotation axis)
						// correct for individual connection starting points
						link=&Elements[i].interactions[potential->links[0*4+1]]; // connection start on A_1 of A_1->A_2
						if(link->location) a-=*(link->location); // minus because a=A_2-A_1
						link=&Elements[k].interactions[potential->links[1*4+0]]; // connection start on A_2 of A_2->A_1
						if(link->location) a+=*(link->location);
						link=&Elements[l].interactions[potential->links[2*4+3]]; // connection start on B_1 of B_1->B_2
						if(link->location) b-=*(link->location); // minus because b=B_2-B_1
						link=&Elements[m].interactions[potential->links[3*4+2]]; // connection start on B_2 of B_2->B_1
						if(link->location) b+=*(link->location);
						link=&Elements[k].interactions[potential->links[1*4+2]]; // connection start on A_2 of A_2->B_1
						if(link->location) delta-=*(link->location); // minus because delta=B_1-A_2
						link=&Elements[l].interactions[potential->links[2*4+1]]; // connection start on B_1 of B_1->A_2
						if(link->location) delta+=*(link->location);
						
						b.V3Cross(delta); // b is normal on plane spanned by B_1->B_2 and A_2->B_1 vectors
						sine=(a*b)*delta.V3Norm(); // sin(theta) = [(B_1->B_2) x (A_2->B_1)]*(A_1->A_2) * |A_2->B_1| / (|(A_1->A_2) x (A_2->B_1)|*|(B_1->B_2) x (A_2->B_1)|)
						a.V3Cross(delta); // a is normal on plane spanned by A_1->A_2 and A_2->B_1 vectors
						cosine=a*b; // cos(theta) = [(A_1->A_2) x (A_2->B_1)]*[(B_1->B_2) x (A_2->B_1)] / (|(A_1->A_2) x (A_2->B_1)|*|(B_1->B_2) x (A_2->B_1)|)
						// denominator is omitted for sine and cosine:
						// theta=atan(sine/cosine) <- it cancels, right?
						if(abs(cosine)>1E-12) theta=atan(sine/cosine);
							else theta=0.0;
						switch(potential->type){
							case dihedral_potential:
							// calculate potential: V=V_1/2*(1+cos(theta+1*f_1)) + V_2/2*(1+cos(2*theta+2*f_2)) + V_3/2*(1+cos(3*theta+3*f_3))
										for(short n=0; n<potential->nr_parameters>>2; n++)
											Vg+=0.5*potential->parameters[n<<1]*(1.0+cos(n*(theta+potential->parameters[(n<<1)+1])));
#if DEBUG_LEVEL>3
										cout << "dihedral: " << Vg << ", " << sine << ", " << cosine << ", " << theta << "\n";
#endif
										break;
							case improper_dihedral_potential:
										Vg+=0.5*potential->parameters[0]*qqpwr(theta-potential->parameters[1],2);
#if DEBUG_LEVEL>3
										cout << "improper dihedral: " << Vg << ", " << sine << ", " << cosine << ", " << theta << "\n";
#endif
										break;
						}
						break;
			case chainspring_potential:
			case spring_potential:
						// get element number
						i=group->elements[potential->partners[0]];
						k=group->elements[potential->partners[1]];
						// distance vector between elements ultimately points from partner 2 to 1
						delta=(Elements[i].center-Elements[k].center);
						// correct for individual connection starting points
						link=&Elements[i].interactions[potential->links[0*2+1]]; // connection start on A of A->B
						if(link->location) delta+=*(link->location);
						link=&Elements[k].interactions[potential->links[1*2+0]]; // connection start on B of B->A
						if(link->location) delta-=*(link->location); // minus because delta=A-B
						// calculate stretch part of the potential: V=k*(r-r0)^2
						r=delta.V3Norm();
						if(potential->type==spring_potential){ // Hook's law: k*(l-l_0)^2
							Vg+=potential->parameters[0]*qqpwr(r-potential->parameters[1],2);
						} else{ // potential well fit: exp((l_a-l)/a) + exp((l-l_b)/b) + m*(l-l_a) + V_0
							theta=(potential->parameters[5]-r); // l_a-l
							Vg+=exp((r-potential->parameters[1])/potential->parameters[0])+exp(theta/potential->parameters[4])-potential->parameters[6]*theta+potential->parameters[7];
						}
#if DEBUG_LEVEL>3
						cout << "bondlength: " << Vg << ", " << delta.V3Norm() << "\n";
#endif
						// calculate bend part from partner 1
						a=*Elements[i].interactions[potential->links[0*2+1]].normal;
						cosine=-(a*delta)/sqrt((a*a)*(delta*delta));
						if(cosine>1.0){ // the Cauchy-Schwartz inequality can be broken by computers ...
							theta=0.0;
						} else{
							if(cosine<-1.0){
								theta=pi;
							} else theta=acos(cosine);
						}
						// calculate potential: V=k*(theta-theta0)^2
						Vg+=potential->parameters[2]*qqpwr(theta,2);
#if DEBUG_LEVEL>3
						cout << "bend A: " << Vg << ", " << cosine << ", " << theta << "\n";
#endif
						// calculate bend part from partner 2
						a=*Elements[k].interactions[potential->links[1*2+0]].normal;
						cosine=(a*delta)/sqrt((a*a)*(delta*delta));
						if(cosine>1.0){ // the Cauchy-Schwartz inequality can be broken by computers ...
							theta=0.0;
						} else{
							if(cosine<-1.0){
								theta=pi;
							} else theta=acos(cosine);
						}
						// calculate potential: V=k*(theta-theta0)^2
						Vg+=potential->parameters[3]*qqpwr(theta,2);
#if DEBUG_LEVEL>3
						cout << "bend B: " << Vg << ", " << cosine << ", " << theta << "\n";
#endif
						break;
		}
	}
	Vg*=configuration->energy_scale[3];
	return Vg;
}

/// Initial potential energy calculation - Provides a reference energy to update as the simulation progresses. Inefficient and N^2, but only called once.
bool MC_Elements::recalcall(double* totalVs, bool skip_charges)
{
	bool success=true;
	double selfVES=0.0;
	double selfVLJ=0.0;
	
	for(unsigned j=0; j<NUM_V_STORE; j++) totalVs[j] = 0.0;
	// loop over elements which are not in groups
	for(unsigned int k = nr_group_elements; k < nr_elements; k++){
		// MuE energy
		if(!skip_charges || configuration->rand_Efield) totalVs[0] += VmuE(k);
		if(configuration->LJwall_calc) Elements[k].VE=Vwall(k,skip_charges,success);
		totalVs[3] += Elements[k].VE;
	}
	if(!success) return false;
	// loop over groups
	Element_Group* group;
	for(unsigned int k = 0; k < nr_groups; k++){
		group=Groups[k];
		success&=selfInteraction(group,group->VES,&(group->VLJ),skip_charges);
		if(!skip_charges) selfVES+=group->VES+0.5*group->Type->charge*group->Type->charge*rf_neutralization*configuration->in2;
		selfVLJ+=group->VLJ;
		// Group potential
		group->VG=GroupPotentialEnergy(group);
		for(unsigned int i=0; i<group->nr_elements; i++){
			// MuE energy
			if(!skip_charges || configuration->rand_Efield) totalVs[0] += VmuE(group->elements[i]);
			if(configuration->LJwall_calc) group->VG+=Vwall(group->elements[i],skip_charges,success);
		}
		totalVs[3] += group->VG;
	}
	if(!success) return false;
	double V;
	Vec3 rmu, position;
	double Ves=0.0;
	double VLJ=0.0;
	bool updateRF=(fabs(rf_correction)>EPS);
	Vec3 RF, RFmirror; // reaction field
	double RFp, RFpmirror; // reaction field potential (interaction starts from charge)
	Vec3 Emu, Emirror; // electrostatics vector potential
	double phiq; // electrostatics scalar potential
	bool calc_dipole;
	Element_Group* kk_group;
	unsigned int nr_charges;
	double kk_charge=0.0;
	unsigned int start_n=0;
	unsigned int end_n=0;
	unsigned int start;
	unsigned int i_kk, end_kk;
	unsigned int i_g, neighbor, end_nr;
	unsigned int index;
	double* energy;
	Element_Group* n_group;
	for(unsigned int jj = 0; jj<N-1; jj++){
		if(jj<nr_groups){
			kk_group=Groups[jj];
			i_kk=kk_group->elements[0]; // first group element index
			end_kk=i_kk+kk_group->nr_elements;
			start=kk_group->number*configuration->max_neighbors;
			kk_charge=kk_group->Type->charge*rf_neutralization;
			start_n=start+kk_group->nr_neighbors_above;
			end_n=start+kk_group->nr_neighbors;
		} else{
			i_kk=jj+nr_N_to_kk;
			end_kk=i_kk+1;
			start=jj*configuration->max_neighbors;
			kk_charge=Elements[i_kk].MyType->charge*rf_neutralization;
			start_n=start+Elements[i_kk].nr_neighbors_above;
			end_n=start+Elements[i_kk].nr_neighbors;
		}
		for(unsigned int kk = i_kk; kk < end_kk; kk++){ // Loop over jj's element(s)
			//Calculate the distance array
			get_distances_above_index(kk); // fills vector array Rmus with distances from element kk to all other elements
			calc_dipole=Elements[kk].MyType->hasmu;
			if(calc_dipole){
				if(configuration->offctrmu){
					position=Elements[kk].rot*Elements[kk].MyType->mu_pos; // correction in case of offcenter dipoles
				} else position.V3Zeros();
			}
			nr_charges=Elements[kk].MyType->nr_charges;
			// Loop over all 2-body neighbor interactions
			for(index=start_n; index<end_n; ++index){
				energy=&pair_energies[energy_indices[index]]; // go straight to ES energies
				if(kk==i_kk){
					*energy=0.0;
					*(energy+1)=0.0;
				}
				neighbor=neighbors[index];
				if(neighbor<nr_groups){
					n_group=Groups[neighbor];
					i_g=n_group->elements[0]; // first group element index
					end_nr=i_g+n_group->nr_elements;
				} else{
					i_g=neighbor+nr_N_to_kk;
					end_nr=i_g+1;
				}
				//Calculate distances
				// VdW interactions (pick whether you want to use the soft-sphere, simple, extended, or Touch LJ function)
				V=0.0;
				for(unsigned int i=i_g; i<end_nr; i++) if((Rdist2LJ[i]<configuration->ljcut2) || (!configuration->cut)) V += (this->*LJpair)(kk,i);
				*energy += V;
				VLJ += V;
				if(!skip_charges){
					if((Rdist2[i_g]<configuration->escut2) || (!configuration->cut)){
						V=0.0;
						for(unsigned int i=i_g; i<end_nr; i++){
							if(calc_dipole){
								rmu=Rmus[i]-position; // rmu is vector from element kk to element i
								Emu=E_mu(i,rmu,RF,updateRF);
								if(LJwallMirror[i]){ // consider image charges from dielectric walls here
									rmu.vec[0]=Rmirror[i]+position.vec[0];
									Emirror=E_mu(i,rmu,RFmirror,updateRF);
									Emirror.vec[0]*=-1.0;
									RFmirror.vec[0]*=-1.0;
									Emu+=Emirror*configuration->LJwall_mirror_factor;
									RF+=RFmirror*configuration->LJwall_mirror_factor;
								}
								V+=(Emu-RF*rf_correction)*Elements[kk].dipole;
							}
							if(nr_charges){
								if(kk==i_kk) V -= kk_charge*Elements[i].MyType->charge;
								for(unsigned int m=0; m<nr_charges; m++){
									rmu=Rmus[i]-Elements[kk].q_pos[m];
									phiq=phi_q(i,rmu,RFp,updateRF);
									if(LJwallMirror[i]){ // consider image charges from dielectric walls here
										if(kk==i_kk) V -= kk_charge*configuration->LJwall_mirror_factor*Elements[i].MyType->q[m];
										rmu.vec[0]=Rmirror[i]+Elements[kk].q_pos[m].vec[0];
										phiq+=phi_q(i,rmu,RFpmirror,updateRF)*configuration->LJwall_mirror_factor;
										RFp+=RFpmirror*configuration->LJwall_mirror_factor;
									}
									V+=(phiq+RFp*rf_correction)*Elements[kk].MyType->q[m];
								}
							}
						}
						*(energy+1) += V;
						Ves += V;
					}
				}
			} // done calculating electrostatics and reaction field as well as Lennard-Jones
		}
	}
	totalVs[1] += Ves*configuration->in2+selfVES;
	totalVs[2] += VLJ+selfVLJ;
	
	return true;
}

void MC_Elements::find_bound_neighbors(unsigned int kk, unsigned int level)
{
	// update pair energies and test for bound things
	unsigned int start=kk*configuration->max_neighbors;
	unsigned int curr_nr_neighbors;
	if(kk<nr_groups){
		curr_nr_neighbors=Groups[kk]->nr_neighbors;
	} else{
		curr_nr_neighbors=Elements[kk-nr_groups+nr_group_elements].nr_neighbors;
	}
	unsigned int e_idx, bd;
	double ilj, ies;
	if(level==0){
		memset(bound_neighbors,0,configuration->N*sizeof(bool));
		nr_bound_neighbors=0;
	}
	// Loop over all neighbors
	for(unsigned int j=0; j<curr_nr_neighbors; j++){
		e_idx=energy_indices[start+j];
		ilj=pair_energies[e_idx];
		ies=pair_energies[e_idx+1];
		if(ilj+ies*configuration->in2<-10.0*configuration->kT){
			bd=neighbors[start+j];
			if(bound_neighbors[bd]==false){
				bound_neighbors[bd] = true;
				nr_bound_neighbors++;
				find_bound_neighbors(bd,level+1);
			}
		}
	}
}

void MC_Elements::evolve_system(bool benchmark)
{
	// Now that everything is set up, create temporary storage space for groups
	if(max_group_nrelements>0) CreateGroupStorage();
	// Initialize clock and output
	cout << "\nSetup complete. Starting simulation...\n";
	cout.flush();
	
	double tstart = (double)clock();
	configuration->Rand_Time=-1.0;
	configuration->storetime=false;
	// Allocate memory
	cosmeans = new Vec3[configuration->steps]; // Stores order parameters
	Ms = new Vec3[configuration->steps]; // Stores total dipole moment of system
	Vs = new double[configuration->steps]; // Stores system volume
	totalV = new Vstore[configuration->steps]; // Stores potential energy components for systems
	
	msmoved[0] = 0.0; msmoved[1] = 0.0;
	accept[0] = 0; accept[1] = 0; accept[2] = 0;
	NpT_accept[0] = 0; NpT_accept[1] = 0; NpT_accept[2] = 0;
	unsigned int tries[3];
	tries[0]=0; tries[1]=0; tries[2]=0;
	nr_tries = 0;
	nr_tries_initial=0;
	nr_tries_equi=0;
	NpT_tries = 0;
	NpT_tries_initial=0;
	NpT_tries_equi=0;
	
	dt=0.0;
	if(configuration->time) dt=hbar_perg_ps/(2*configuration->kT); // set first time-step according to Heisenberg uncertainty principle
	simulation_time=0.0;
	
	Vec3 overallM;
	Vec3 Tmpmu;
	double Vstep[NUM_V_STORE];
	double oldmuE = 0.0, muE = 0.0;
	double oldES = 0.0 , VES = 0.0;
	double oldVLJ  = 0.0, VLJ = 0.0;
	double tmovtemp = 0.0;
	double rmovtemp = 0.0;
	double tmovstep, rmovstep;
	double OldV, NewV, deltaV, totVstep, nMstep;
	double moving_oids, rotating_oids;
	double iNkT; // 1/NkT
	unsigned int ebf_size=(nr_elements>>5)+(nr_elements%32>0); // right shift by five is division by 32
	__uint32_t* element_bitfield=new __uint32_t[ebf_size];
	for(unsigned int i=0; i<ebf_size; i++) element_bitfield[i]=0;
	unsigned int move_array_size=2*nr_elements; // that should most likely be enough ...
	unsigned int* move_elements=(unsigned int*)malloc(move_array_size*sizeof(unsigned int)); // ... but, just in case, keep the array dynamic
	unsigned int max_nr_moves_per_one=0;
	unsigned int max_nr_neighbor_moves=NEIGHBOR_UPDATE_FREQ;
	nr_individual_elements=nr_elements-nr_group_elements;
	double* save_groupVG=NULL;
	if(configuration->LJwall_calc) save_groupVG=new double[nr_groups];
	// determine how many neighbors we'll expect maximally
	configuration->max_neighbors=configuration->N-1; // worst case if there's no cutoff
	if(configuration->cut){
		configuration->rfcut=configuration->escut;
		adjust_cut2plus(configuration);
		double cutr=sqrt(configuration->cut2plus);
		if(configuration->NpT && !configuration->NpT_adjust_EScut){
			configuration->max_neighbors=(unsigned int)round(qqpwr(cutr,3)/smallest_volume);
			if(configuration->max_neighbors>configuration->N-1) configuration->max_neighbors=configuration->N-1;
		} else{
			configuration->max_neighbors=(unsigned int)round(4.0*pi/3.0*qqpwr(cutr,3)/configuration->targetV*configuration->N);
			if(configuration->max_neighbors>configuration->N-1) configuration->max_neighbors=configuration->N-1;
		}
		if(configuration->max_neighbors<8) configuration->max_neighbors=8;
		if(configuration->max_neighbors>configuration->N-1) configuration->max_neighbors=configuration->N-1;
		cout << "-> Maximum number of neighbors: " << configuration->max_neighbors << "\n";
	}
	max_nr_neighbors=(configuration->N*configuration->max_neighbors);
	neighbor_switch=false;
	
	neighbor_store=new unsigned int[max_nr_neighbors<<1];
	pair_energy_store=new double[max_nr_neighbors<<1];
	energy_index_store=new unsigned int[max_nr_neighbors<<1];
	
	neighbors=&neighbor_store[(unsigned int)neighbor_switch*max_nr_neighbors];
	pair_energies=&pair_energy_store[(unsigned int)neighbor_switch*max_nr_neighbors];
	energy_indices=&energy_index_store[(unsigned int)neighbor_switch*max_nr_neighbors];
	
	bound_neighbors = new bool[configuration->N];
	memset(bound_neighbors,0,configuration->N*sizeof(bool));
	nr_bound_neighbors = 0;
	Vec3* center_store = new Vec3[configuration->N];
	double* VE_store = new double[configuration->N];
	
	initialize_neighbors();
	
	Vec3 v(0.0);
	Vec3 ds(0.0);
	double m=0.0;
	double max_step=0.0;
	double qterms=0.0;
	double mds2=0.0;
	double initial_mds2=0.0;
	double avg_speed=0.0;
	if(dt>EPS){
		for(unsigned int i=0; i<nr_elements; i++){ // set initial speeds according to kT
			Vec3 ran_vec(0.0);
			if(Elements[i].MyType->mass*amu_to_perg_ps2_per_Ang2>EPS){
				double boltzmann_speed=sqrt(configuration->kT/(Elements[i].MyType->mass*amu_to_perg_ps2_per_Ang2)); // speed is in Angström/ps
				for(unsigned int j=0; j<3; j++) ran_vec.vec[j]=ran_n();
				Elements[i].ds=ran_vec*boltzmann_speed*dt;
			}
			double mds2_one=Elements[i].MyType->mass*amu_to_perg_ps2_per_Ang2*(Elements[i].ds*Elements[i].ds);
			Elements[i].Ekin=0.5*mds2_one/(dt*dt);
			Elements[i].max_step=0.0;
			mds2+=mds2_one;
			initial_mds2+=Elements[i].MyType->mass*amu_to_perg_ps2_per_Ang2*configuration->maxtrans*configuration->maxtrans*(ran_vec*ran_vec);
		}
		if(mds2>EPS){
			double timescale=sqrt(initial_mds2/mds2);
			cout << "\nInitial timestep is " << dt*timescale << " ps\n";
			if(timescale>1.0){
				cout << "WARNING:\n\tInitial timestep is quite large (from Heisenberg's uncertainty one would expect " << dt << " ps).\n";
				cout << "\tIf the system is quite dense and does not hold temperature well, consider adjusting the average step size.\n";
				cout << "\tFor example, from currently " << configuration->maxtrans << " Angström to " << configuration->maxtrans/timescale << " Angström.\n";
			}
			dt*=timescale;
			for(unsigned int i=0; i<nr_elements; i++){
				Elements[i].ds*=timescale;
				avg_speed+=Elements[i].ds.V3Norm()/dt; // Angström/ps
			}
			avg_speed/=nr_elements;
		}
		cout << "-> Average speed is " << avg_speed*100.0 << " m/s\n"; // 1 Angström/ps = 100 m/s
	}
	Mat33 delta_rot;
	Vec3 delta_trans;
	
	// Define flags and constants;
	bool skip_charges = true;
	bool recalc_muE = false;
	unsigned int step = configuration->steps;
	unsigned int start_step = 0;
	unsigned int eqsteps = configuration->steps-configuration->laststep; // Pre-equilibrium cycles before simulation starts recording data
	
	unsigned int nr_entities=nr_elements-nr_group_elements+nr_groups;
	double* avg_trans=NULL;
	double* avg_rot=NULL;
	double halfbox=average(configuration->boxlength,3)/2.0;
	if(configuration->stepsize_average>0){
		avg_trans=new double[nr_entities];
		avg_rot=new double[nr_entities];
		trans_store=new double[nr_entities*configuration->stepsize_average];
		rot_store=new double[nr_entities*configuration->stepsize_average];
		for(unsigned int i=0; i<nr_entities*configuration->stepsize_average; i++){
			trans_store[i]=configuration->maxtrans*configuration->maxtrans;
			rot_store[i]=configuration->maxrot*configuration->maxrot/3.0;
		}
		running_average=new unsigned int[nr_entities];
		for(unsigned int i=0; i<nr_entities; i++){
			avg_trans[i]=configuration->maxtrans*configuration->maxtrans*configuration->stepsize_average;
			avg_rot[i]=configuration->maxrot*configuration->maxrot*configuration->stepsize_average/3.0;
			running_average[i]=0;
		}
	}
	
	// Create Rmus
	Rmus = new Vec3[nr_elements]; // set up distance array
	Rdist2 = new double[nr_elements];
	Rdist2LJ = new double[nr_elements];
	Rmirror = new double[nr_elements];
	LJwallMirror = new bool[nr_elements];
	// Set initial dielectric conditions and calculate initial reaction field
	update_rf_correction();
	// Calculate overall dipole moment
	get_means(cosmeans[configuration->resume_step],overallM);
	
	unsigned int kk,temp;
	// Determine number of oids that can move or rotate
	moving_oids = 0.0;
	rotating_oids = 0.0;
	for(unsigned int i=0; i<configuration->num_element_types; i++){
		if((!configuration->element_types[i]->rot_notrans) && (!configuration->element_types[i]->still)) moving_oids += configuration->element_types[i]->number;
		if(!configuration->element_types[i]->still) rotating_oids += configuration->element_types[i]->number;
	}
	bool adjust_internal=configuration->NpT && configuration->adjust_internal;
	if(adjust_internal){
		selfLJlambda=calc_lambda();
		if(fabs(selfLJlambda-1.0)<EPS){
			adjust_internal=false;
			selfLJlambda=1.0;
		}
	}
	if(configuration->use_trajectory){
		Traj_EP* elements = new Traj_EP[nr_elements];
		if(GetConfig()->GetElementProperties(configuration->resume_step,elements)){
			for(unsigned int i=0; i<nr_elements; i++){
				Element* element=&Elements[i];
				if(elements[i].element_type!=element->mytype){
					cout << "Element types differ between configuration file (" << element->mytype << ") and trajectory file (" << elements[i].element_type << ").\n";
					exit(3);
				}
				if(element->group){
					if((elements[i].group_type==-1) || (elements[i].group_type!=(int)element->group->type)){
						cout << "Group types differ between configuration file and trajectory file.\n";
						exit(3);
					}
				} else{
					if(elements[i].group_type!=-1){
						cout << "Group types differ between configuration file and trajectory file.\n";
						exit(3);
					}
				}
				element->center=elements[i].position;
				element->rot=AxisAngle2Rot(elements[i].rotation_vector);
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
			// calculate group centers
			Element_Group* group;
			for(unsigned int k = 0; k < nr_groups; k++){
				group=Groups[k];
				group->center.vec[0]=0.0; group->center.vec[1]=0.0; group->center.vec[2]=0.0;
				for(unsigned int i=0; i<group->nr_elements; i++) group->center+=Elements[group->elements[i]].center; // calculate geometric group center
				group->center*=group->inv_nr_elements;
			}
#ifndef USE_CMWC4096
			// restore random number generator state
			__int32_t idum=configuration->rngseed;
			if(idum>0){
				cout << "WARNING: Cannot restore random number generator state. Initial seed is positive (needs to be negative).\n";
			} else{
				while(idum!=*configuration->idum) ran2(idum);
			}
#endif
		} else exit(1);
		Evec=configuration->Efield*MV_per_m_to_perg_per_Debye;
		if(!configuration->restart_calculations){
			cout << "\nResuming previous calculation from step " << configuration->resume_step << ".\n";
			traj_laststep=configuration->resume_step;
			if(!GetConfig()->GetStatistics(configuration->resume_step,totalV,Ms,cosmeans,Vs,msmoved,accept,tries)){
				cout << "ERROR: Could not restore previous calculations.\n-> Use <restart_calculations=1> in trajectory file to restart calculations from scratch.\n";
				exit(1);
			}
			update_volume(configuration,Vs[configuration->resume_step-1]);
			if(adjust_internal){
				selfLJlambda=calc_lambda(configuration->resume_step);
				if(fabs(selfLJlambda-1.0)<EPS){
					adjust_internal=false;
					selfLJlambda=1.0;
				}
			}
			if(configuration->transition){
				LJlambda=calc_lambda(configuration->resume_step);
				if(fabs(LJlambda-1.0)<EPS){
					if(configuration->LJwall_calc){
						if(configuration->LJwall_fixed) configuration->LJwall_xm=configuration->LJwall_x; else configuration->LJwall_xm=0.5*configuration->boxlength[0];
					}
					configuration->targetV=configuration->transition_target_volume;
					configuration->transition=false;
					LJlambda=1.0;
				}
				LJlambda*=configuration->transition_LJ_delta;
				LJlambda+=configuration->transition_LJ_start;
			}
			nr_tries=tries[0]+tries[1]+tries[2];
			if(tries[0]>0){
				nr_tries_initial=tries[0];
				if(tries[1]>0) nr_tries_equi=tries[0]+tries[1];
			}
			for(unsigned int j=0; j<NUM_V_STORE; j++) Vstep[j]=totalV[configuration->resume_step-1][j];
			skip_charges=!(configuration->resume_step>(int)configuration->randsteps);
			if(!skip_charges){
				if(configuration->NpT && configuration->cut){
					// escut2, the square of the electrostatics cutoff, is the only thing used during runtime (for distances, r*r without the additional sqrt saves time)
					configuration->rfcut=configuration->escut*cbrt(configuration->V/configuration->targetV);
					configuration->escut2=qqpwr(configuration->rfcut,2);
				}
#if DEBUG_LEVEL>1
				cout << "Electrostatics cutoff is at " << configuration->rfcut << " Angström.\n";
#endif
				if(configuration->dyneps){
					get_means(cosmeans[configuration->resume_step-1],overallM);
					update_eps(overallM,configuration->resume_step); // use the overall dipole from the previous step to update epsilon
					update_rf_correction();
#if DEBUG_LEVEL>3
					cout << "epsRF: " << epsRF << "\n";
#endif
				}
			}
#if DEBUG_LEVEL>1
			find_neighbors();
			cout << "ES in trajectory: " << Vstep[1] << " <-> ES recalculated: " << calc_charge_interaction() << "\n";
			cout << "LJ in trajectory: " << Vstep[2] << "\n";
			recalcall(Vstep,skip_charges);
			cout << "ES recalcall(): " << Vstep[1] << "\n";
			cout << "LJ recalcall(): " << Vstep[2] << "\n";
//			exit(42);
#endif
		}
		if(configuration->restart_volume){
			cout << "-> Increasing volume of resume step system state to match initial bounding-sphere limited volume.\n";
			double factor=cbrt(configuration->V/configuration->targetV);
			// loop over elements which are not in groups
			for(unsigned int k = nr_group_elements; k < nr_elements; k++) Elements[k].center*=factor;
			// loop over groups
			Element_Group* group;
			for(unsigned int k = 0; k < nr_groups; k++){
				group=Groups[k];
				delta_trans=group->center*(factor-1.0); // new position - old position points to new position
				group->center*=factor; // move group center but leave group elements undisturbed
				for(unsigned int i=0; i<group->nr_elements; i++) Elements[group->elements[i]].center+=delta_trans;
			}
		}
	}
	if(configuration->steps-configuration->resume_step<=0){
		cout << "\nNothing to do here, already calculated " << configuration->steps << " steps. (Is steps >= resume_step ?)\n\n";
		exit(0);
	}
	changevolume=false;
	if((fabs(configuration->V-configuration->targetV)>EPS) && (!configuration->NpT || configuration->transition)){
		changevolume=true;
		V_diff=configuration->V-configuration->targetV;
		double sum_invn=1.0;
		for(unsigned int i=2; i<=configuration->randsteps>>1; i++) sum_invn+=1.0/sqrt((double)i);
		V_diff/=sum_invn;
	}
	if(changevolume) cout << "\nADJUSTING VOLUME FROM " << configuration->V << " TO " << configuration->targetV << " ANGSTRÖM³\n\n";
	double LJwall_start_trans=0.5*configuration->boxlength[0];
	if(configuration->transition){
		LJlambda=calc_lambda(configuration->resume_step);
		if(fabs(LJlambda-1.0)<EPS){
			configuration->transition=false;
			LJlambda=1.0;
		}
		LJlambda*=configuration->transition_LJ_delta;
		LJlambda+=configuration->transition_LJ_start;
		if(configuration->LJwall_calc) configuration->LJwall_xm=0.5*configuration->boxlength[0];
		cout << "\nTRANSITIONING FROM " << configuration->transition_start_volume << " TO " << configuration->transition_target_volume << " ANGSTRÖM³\n\n";
	}
	if(configuration->NpT){
		cout << "\nUSING NPT ENSEMBLE CALCULATIONS (" << configuration->pext/atm_in_pErg_per_Ang3 << " ATM, +/- " << configuration->maxdV;
		if(configuration->NpT_lnV_change) cout << " LN V"; else cout << " ANGSTRÖM";
		cout << " CHANGE EVERY " << configuration->NpT_move_nr << " MOVES) TO DETERMINE SYSTEM VOLUME\n\n";
	}
	cout.flush();
	
	
	if(configuration->chkcoord) store_trajectory(configuration->resume_step);
	
	if(adjust_internal) cout << "\nADJUSTING INTERNAL ENERGIES UNTIL SYSTEM REACHES TARGET VOLUME OF " << configuration->targetV << " Angström\n\n";
	
	if(configuration->collective_moves) cout << "\nUSING COLLECTIVE MOVES.\n\n";
	
	cout.setf(ios::fixed,ios::floatfield);
	cout.precision(6);
	bool success=true;
	
	if(benchmark){
		cout << "\n***************************************************************\n";
		cout << "*                                                             *\n";
		cout << "* BENCHMARKING POTENTIAL ENERGY CALCULATION FOR ENTIRE SYSTEM *\n";
		cout << "*                                                             *\n";
		cout << "***************************************************************\n\n";
		for(step=configuration->resume_step; step<configuration->steps; step++){
			iNkT = 1.0/((double)N*configuration->kT);
			// Calculate first step energy if this is the first step
			if((int)step == configuration->resume_step){
				success=recalcall(Vstep,false); // dipole/charges are calculated
				if(!success){
					cout << "Initial position overlap.\n";
					exit(3);
				}
				totVstep=0.0;
				cout << "Initial energy components (NkT):\nVmuE\t\tVES\t\tVLJ\t\tVG\t\tTotal\n";
				for(unsigned int jj=0; jj<NUM_V_STORE; jj++){
					totVstep += Vstep[jj];
					cout << Vstep[jj]*iNkT << "\t";
				}
				cout << totVstep*iNkT << "\n";
#if DEBUG_LEVEL>2
				cout << "Initial energy components (kcal/mol of boxes):\nVmuE\t\tVES\t\tVLJ\t\tVG\t\tTotal\n";
				for(unsigned int jj=0; jj<NUM_V_STORE; jj++){
					cout << Vstep[jj]*perg_to_kcal_per_mol << "\t";
				}
				cout << totVstep*perg_to_kcal_per_mol << "\n";
				cout.flush();
#endif
				// Check to see if randomization is used (not needed if already beyond)
				if(configuration->randsteps>step){
					skip_charges = true;
					cout << "Disabling electrostatics and randomizing system...\n\n";
					// Discard everything except the VdW term, so energies don't get
					// unphysically propagated during the randomization stage
					Vstep[0] = 0.0;
					Vstep[1] = 0.0;
					for(unsigned int jj = 0; jj < ((N*configuration->max_neighbors)>>1); jj++){
						pair_energies[(jj<<1)+1]=0.0;
					}
					for(unsigned int jj = 0; jj < nr_groups; jj++) Groups[jj]->VES=0.0;
				} else{
					skip_charges = false;
					cout << "Randomization phase skipped.\n\n";
				}
				// Set up columns for remainder of log
				cout << "STEP\t\tU (NkT)\n";
				cout.flush();
			}
			unsigned int splus1 = step+1;
			if(step == configuration->randsteps){ // Recalculate total energy once mu-mu interactions are turned on and display energy/timing update
				skip_charges = false;
				configuration->Rand_Time=((double)clock()-tstart)/CLOCKS_PER_SEC;
				cout << "\nInitial randomization complete. Randomization took " << configuration->Rand_Time << " seconds.\n";
				cout << "\nEnabling electrostatics and recalculating total energy.\n";
				Vstep[0] = 0.0;
				// loop over elements which are not in groups
				for(unsigned int k = nr_group_elements; k < nr_elements; k++) Vstep[0] += VmuE(k);
				// loop over groups
				Element_Group* group;
				for(unsigned int k = 0; k < nr_groups; k++){
					group=Groups[k];
					for(unsigned int i=0; i<group->nr_elements; i++) Vstep[0] += VmuE(group->elements[i]);
				}
				Vstep[1] = calc_charge_interaction(); // updates charge/dipole interaction
				cout << "Energy components (NkT):\nVmuE\t\tVES\t\tVLJ\t\tVG\t\tTotal\n";
				totVstep = 0.0;
				for(unsigned int jj=0; jj<NUM_V_STORE; jj++){
					totVstep += Vstep[jj];
					cout << Vstep[jj]*iNkT << "\t";
				}
				cout << totVstep*iNkT << "\n";
				cout << "\nSTEP\t\tU (NkT)\n";
				cout.flush();
			}
			success=recalcall(Vstep, skip_charges); // dipole/charges are not calculated if in randomization
			if(!success){
				cout << "WARNING: Overlap detected.\n";
			}
			// Calculate cos^n, M, V
			get_means(cosmeans[step],overallM);
			Ms[step]=overallM;
			Vs[step]=configuration->V;
			for(unsigned int jj=0; jj<NUM_V_STORE; jj++) totalV[step][jj] = Vstep[jj];
			// check for unexpected program termination
			if(exitnow){
				configuration->Runtime=((double)clock()-tstart)/CLOCKS_PER_SEC;
				configuration->storetime=true;
				store_trajectory(splus1);
				cout << "\n";
				delete this;
				exit(42);
			}
			// Output progress indicator to user, including total energy trace.
			if (splus1 % 100 == 0){
				totVstep = 0.0;
				for(unsigned int jj=0; jj<NUM_V_STORE; jj++) totVstep += Vstep[jj];
				totVstep *= iNkT;
				cout << "\n" << splus1  << "\t\t" << totVstep << "\n";
			}
			// set time if this is the last step
			if(step+1==configuration->steps){
				configuration->storetime=true;
				configuration->Runtime=((double)clock()-tstart)/CLOCKS_PER_SEC;
			}
			// Save intermediate coordinates. Note that this can require large amounts of disk space
			if((configuration->chkcoord) && (splus1 % configuration->grfreq==0)) store_trajectory(splus1);
			// Premature termination code for debugging purposes. Do NOT use for normal operation.
#if DEBUG_EXIT_STEP>=0
			if(step+1>=DEBUG_EXIT_STEP){
				cout << "Terminating run at step " << step+1 << "\n";
				configuration->storetime=true;
				configuration->Runtime=((double)clock()-tstart)/CLOCKS_PER_SEC;
				store_trajectory(splus1);
				exit(12);
			}
#endif
		}
		start_step=configuration->steps; // We are benchmarking and don't want to run the production loop below
	} else start_step=configuration->resume_step;
	
	// This is it! Start the main (steps) loop
	for(step=start_step; step<configuration->steps; step++){
		// Check if we need to add new elements/groups
		if((configuration->add_slowly) && (step%configuration->add_slowly==0)){
			// Place new element/group
			unsigned int added_other_groups=0;
			bool added_something=false;
			Vec3 add_vec;
			add_vec.vec[2]=2.0*ranQ()-1.0; // cos(theta)
			add_vec.vec[0]=sqrt(1.0-add_vec.vec[2]*add_vec.vec[2]); // sin(theta)
			add_vec.vec[1]=add_vec.vec[0];
			double phi=2.0*pi*ranQ();
			add_vec.vec[0]*=cos(phi);
			add_vec.vec[1]*=sin(phi);
			double dist=cbrt(configuration->V)/2.0;
			if(configuration->NpT){
				if(configuration->cut) dist+=sqrt(configuration->ljcut2); else dist*=2.0;
			}
			add_vec*=dist;
			if(nr_elements==nr_group_elements){// only add a group if the last thing added was a group
				for(unsigned int i=0; i<configuration->num_groups; i++){ // check to see if a group should be added
					// if last group type corresponds to current type or is the next one (>=)
					if(i>=Elements[nr_elements-1].group->type){
						if(nr_groups-added_other_groups<configuration->groups[i]->number){
							update_volume(configuration,8.0*qqpwr(dist,3));
							AddGroup(i);
							Translate(Groups[nr_groups-1],add_vec);
							double theta=2.0*pi*ranQ();
							double phi=2.0*pi*ranQ();
							double gamma=2.0*pi*ranQ();
							Rotate(Groups[nr_groups-1],Groups[nr_groups-1]->center,theta,phi,gamma);
							cout << "\tAdded new group to simulation.\n";
							added_something=true;
							break;
						}
					}
					added_other_groups+=configuration->groups[i]->number;
				}
			}
			unsigned int added_other_elements=0;
			if(!added_something){
				for(unsigned int i=0; i<configuration->num_element_types; i++){
					if(i>=Elements[nr_elements-1].mytype){
						if(nr_elements-nr_group_elements-added_other_elements<configuration->element_types[i]->number){
							update_volume(configuration,8.0*qqpwr(dist,3));
							AddElement(i);
							Translate(nr_elements-1,add_vec);
							double theta=2.0*pi*ranQ();
							double phi=2.0*pi*ranQ();
							double gamma=2.0*pi*ranQ();
							Rotate(nr_elements-1,theta,phi,gamma);
							cout << "\tAdded new elements to simulation.\n";
							added_something=true;
							break;
						}
					}
					added_other_elements+=configuration->element_types[i]->number;
				}
			}
			if(!added_something){
				cout << "\n-> Finished adding elements/groups to simulation.\n";
				configuration->add_slowly=0;
			}
			
			// grow arrays
			if((nr_elements>>5)+(nr_elements%32>0)>ebf_size){
				ebf_size=(nr_elements>>5)+(nr_elements%32>0); // right shift by five is division by 32
				delete[] element_bitfield;
				element_bitfield=new __uint32_t[ebf_size];
				for(unsigned int i=0; i<ebf_size; i++) element_bitfield[i]=0;
			}
			move_array_size=2*nr_elements; // that should most likely be enough ...
			move_elements=(unsigned int*)realloc(move_elements,move_array_size*sizeof(unsigned int)); // ... but, just in case, keep the array dynamic
			if(!move_elements){ // big problem ...
				cout << "ERROR: Can not allocate enough memory to hold element numbers for moving.\n";
				exit(7);
			}
			nr_individual_elements=nr_elements-nr_group_elements;
			if(configuration->LJwall_calc){
				delete[] save_groupVG;
				save_groupVG=new double[nr_groups];
			}
			nr_entities=nr_elements-nr_group_elements+nr_groups;
			if(configuration->stepsize_average>0){
				delete[] avg_trans;
				delete[] avg_rot;
				delete[] trans_store;
				delete[] rot_store;
				delete[] running_average;
				avg_trans=new double[nr_entities];
				avg_rot=new double[nr_entities];
				trans_store=new double[nr_entities*configuration->stepsize_average];
				rot_store=new double[nr_entities*configuration->stepsize_average];
				for(unsigned int i=0; i<nr_entities*configuration->stepsize_average; i++){
					trans_store[i]=configuration->maxtrans*configuration->maxtrans;
					rot_store[i]=configuration->maxrot*configuration->maxrot/3.0;
				}
				running_average=new unsigned int[nr_entities];
				for(unsigned int i=0; i<nr_entities; i++){
					avg_trans[i]=configuration->maxtrans*configuration->maxtrans*configuration->stepsize_average;
					avg_rot[i]=configuration->maxrot*configuration->maxrot*configuration->stepsize_average/3.0;
					running_average[i]=0;
				}
			}
			delete[] Rmus;
			delete[] Rdist2;
			delete[] Rdist2LJ;
			delete[] Rmirror;
			delete[] LJwallMirror;
			// Create Rmus
			Rmus = new Vec3[nr_elements]; // set up distance array
			Rdist2 = new double[nr_elements];
			Rdist2LJ = new double[nr_elements];
			Rmirror = new double[nr_elements];
			LJwallMirror = new bool[nr_elements];
			// Set dielectric conditions and calculate reaction field
			update_rf_correction();
			// Calculate overall dipole moment
			get_means(cosmeans[configuration->resume_step],overallM);
			// Determine number of oids that can move or rotate
			moving_oids = 0.0;
			rotating_oids = 0.0;
			for(unsigned int i=0; i<configuration->num_element_types; i++){
				if((!configuration->element_types[i]->rot_notrans) && (!configuration->element_types[i]->still)) moving_oids += configuration->element_types[i]->number;
				if(!configuration->element_types[i]->still) rotating_oids += configuration->element_types[i]->number;
			}
			// Recalculate everything
			success=recalcall(Vstep,false); // dipole/charges are calculated
			if(!success){
				cout << "Position overlap.\n";
				exit(3);
			}
		}
		// Reset displacement counter to zero
		tmovstep = 0.0;
		rmovstep = 0.0;
		// Update dielectric if using dynamic reaction field, and calculate the reaction field coefficient
		// Note that this requires recalculating mu-mu interactions, causing a performance hit
		if((configuration->dyneps) && (!skip_charges)){
			update_eps(overallM,step); // use the overall dipole from the previous step to update epsilon
			update_rf_correction();
			Vstep[1]=calc_charge_interaction();
		}
		// check if E and T need to be updated, and update uE energy if E is updated
		recalc_muE = setEandT(step);
		if((recalc_muE) && (!skip_charges || configuration->rand_Efield)){
			Vstep[0] = 0.0;
			// loop over elements which are not in groups
			for(unsigned int k = nr_group_elements; k < nr_elements; k++) Vstep[0] += VmuE(k);
			// loop over groups
			Element_Group* group;
			for(unsigned int k = 0; k < nr_groups; k++){
				group=Groups[k];
				for(unsigned int i=0; i<group->nr_elements; i++) Vstep[0] += VmuE(group->elements[i]);
			}
			recalc_muE = false;
		}
		iNkT = 1.0/((double)N*configuration->kT);
		// Calculate first step energy if this is the first step
		if((int)step == configuration->resume_step){
			success=recalcall(Vstep,false); // dipole/charges are calculated
			if(!success){
				cout << "Initial position overlap.\n";
				exit(3);
			}
			totVstep=0.0;
			cout << "Initial energy components (NkT):\nVmuE\t\tVES\t\tVLJ\t\tVG\t\tTotal\n";
			for(unsigned int jj=0; jj<NUM_V_STORE; jj++){
				totVstep += Vstep[jj];
				cout << Vstep[jj]*iNkT << "\t";
			}
			cout << totVstep*iNkT << "\n";
#if DEBUG_LEVEL>2
			cout << "Initial energy components (kcal/mol of boxes):\nVmuE\t\tVES\t\tVLJ\t\tVG\t\tTotal\n";
			for(unsigned int jj=0; jj<NUM_V_STORE; jj++){
				cout << Vstep[jj]*perg_to_kcal_per_mol << "\t";
			}
			cout << totVstep*perg_to_kcal_per_mol << "\n";
			cout.flush();
#endif
			// Check to see if randomization is used (not needed if already beyond)
			if(configuration->randsteps>step){
				skip_charges = true;
				cout << "Disabling electrostatics and randomizing system...\n\n";
				// Discard everything except the VdW term, so energies don't get
				// unphysically propagated during the randomization stage
				Vstep[0] = 0.0;
				Vstep[1] = 0.0;
				for(unsigned int jj = 0; jj < ((N*configuration->max_neighbors)>>1); jj++){
					pair_energies[(jj<<1)+1]=0.0;
				}
				for(unsigned int jj = 0; jj < nr_groups; jj++) Groups[jj]->VES=0.0;
			} else{
				skip_charges = false;
				cout << "Randomization phase skipped.\n\n";
			}
			// Set up columns for remainder of log
			cout << "STEP\t\tU (NkT)\t\t\tNORMALIZED ORDER\n";
			cout.flush();
		}
		unsigned int splus1 = step+1;
		if(changevolume){
			ShrinkVolume(splus1-configuration->restart_volume*(configuration->resume_step-1));
			success=recalcall(Vstep, skip_charges); // dipole/charges are not calculated if in randomization
			if(!success){
				cout << "Unrecoverable overlap while shrinking volume.\n";
				exit(3);
			}
		}
		if(step == configuration->randsteps){ // Recalculate total energy once mu-mu interactions are turned on and display energy/timing update
			skip_charges = false;
			nr_tries_initial=nr_tries;
			NpT_tries_initial=NpT_tries;
			configuration->Rand_Time=((double)clock()-tstart)/CLOCKS_PER_SEC;
			cout << "\nInitial randomization complete. Randomization took " << configuration->Rand_Time << " seconds.\n";
			cout << "Initial acceptance rate (%%): " << (double)accept[0]/(double)(nr_tries)*100.0 << "\n";
			cout << "Initial RMS move sizes: " << sqrt(msmoved[0]/configuration->randsteps) << " (A)\t" << sqrt(msmoved[1]/configuration->randsteps) << " (rad)\n";
			cout << "\nEnabling electrostatics and recalculating total energy.\n";
			Vstep[0] = 0.0;
			// loop over elements which are not in groups
			for(unsigned int k = nr_group_elements; k < nr_elements; k++) Vstep[0] += VmuE(k);
			// loop over groups
			Element_Group* group;
			for(unsigned int k = 0; k < nr_groups; k++){
				group=Groups[k];
				for(unsigned int i=0; i<group->nr_elements; i++) Vstep[0] += VmuE(group->elements[i]);
			}
			Vstep[1] = calc_charge_interaction(); // updates charge/dipole interaction
			cout << "Energy components (NkT):\nVmuE\t\tVES\t\tVLJ\t\tVG\t\tTotal\n";
			totVstep = 0.0;
			for(unsigned int jj=0; jj<NUM_V_STORE; jj++){
				totVstep += Vstep[jj];
				cout << Vstep[jj]*iNkT << "\t";
			}
			cout << totVstep*iNkT << "\n";
#if DEBUG_LEVEL>2
			cout << "Energy components (kcal/mol of boxes):\nVmuE\t\tVES\t\tVLJ\t\tVG\t\tTotal\n";
			totVstep = 0.0;
			for(unsigned int jj=0; jj<NUM_V_STORE; jj++){
				totVstep += Vstep[jj];
				cout << Vstep[jj]*perg_to_kcal_per_mol << "\t";
			}
			cout << totVstep*perg_to_kcal_per_mol << "\n\n";
#else
			cout << "\n";
#endif
			cout << "STEP\t\tU (NkT)\t\t\tNORMALIZED ORDER\n";
			cout.flush();
		}
		if(step==eqsteps){
			nr_tries_equi=nr_tries;
			NpT_tries_equi=NpT_tries;
		}
#ifdef USE_CMWC4096
		int nr_elements_movable=0;
		double nr_group_movable=0.0;
		for(unsigned int ii=0; ii<nr_elements; ii++){
			Elements[ii].nr_trials=0;
			Elements[ii].nr_moves=0;
			Elements[ii].current_s=Vec3(0.0);
			// set bitfield to zero
			if((ii&31)==0) element_bitfield[ii>>5]=0; // for b being power of 2: (a%b = a&(b-1)) => ii%32 = ii&31
			bool movable=false;
			if(Elements[ii].group){
				Element_Group* group=Elements[ii].group;
				if(group->fixed){
					movable=false;
				} else{
					if(ii-group->elements[0]==0) nr_group_movable+=group->Type->nr_rand_per_cycle; // at the first element in group
					if((group->Type->nr_movable==0) && (ii-group->elements[0]==0)){
						nr_group_movable++; // ensure that groups with fixed bonds can still be moved around
						movable=true;
					} else{
						if(!group->Type->fixed_elements[ii-group->elements[0]]) movable=true;
					}
				}
			} else{
				if(Elements[ii].fixed){
					movable=false;
				} else{
					if(!Elements[ii].MyType->still){
						nr_elements_movable++;
						movable=true;
					}
				}
			}
			// set bitfield to 1 for elements that can move and to 0 for ones that can't
			element_bitfield[ii>>5]|=movable<<(ii&31);
		} // the bitfield is now set with ones for elements that can move (the still ones are set to zero)
		nr_elements_movable+=(unsigned int)nr_group_movable;
		unsigned int nr_elements_moving=0;
		bool setmore=true;
		while(setmore){
			kk = (unsigned int)((rotating_oids+nr_group_elements)*ranQ()); // rotating_oids: nr of oids at least able to rotate (not still)
			bool accept=false;
			if(kk>=nr_group_elements){ // did not choose a group element
				temp=nr_group_elements; // group elements are created first
				for(unsigned int iii=0; iii<configuration->num_element_types; iii++){
					temp+=configuration->element_types[iii]->number;
					if(configuration->element_types[iii]->still){ // all of that type are still
						kk += configuration->element_types[iii]->number; // this works b/c rotating_oids <= nr of oids in total
					} else{
						if (temp>=kk) break;
					}
				} // kk is now a random element's index which is allowed to rotate (and maybe even to move)
				if(kk<nr_elements){
					if((Elements[kk].nr_moves<1) || (dt<=EPS)) accept=true; // in non-time mode we use normal MC random selection (in time mode only move once)
					if(Elements[kk].fixed) accept=false;
				}
			} else{ // group element chosen
				Element_Group* group=Elements[kk].group;
				if(!group->Type->fixed_elements[kk-group->elements[0]]) accept=true;
				if(group->Type->nr_movable==0) accept=true;
				if(group->fixed) accept=false;
			}
			if(accept){
				Elements[kk].nr_moves++;
				if(Elements[kk].nr_moves>max_nr_moves_per_one) max_nr_moves_per_one=Elements[kk].nr_moves;
				if(nr_elements_moving>=move_array_size){ // auto-expand if needed (should only occur in time mode (every once in a while)
					move_array_size+=nr_elements;
					move_elements=(unsigned int*)realloc(move_elements,move_array_size*sizeof(unsigned int));
					if(!move_elements){ // big problem ...
						cout << "ERROR: Can not allocate enough memory to hold element numbers for moving.\n";
						exit(7);
					}
				}
				move_elements[nr_elements_moving]=kk;
				element_bitfield[kk>>5]&=~(1<<(kk&31)); // unset bit for element kk (if all bits are zero then every element that can move has been chosen)
				nr_elements_moving++;
			}
			if((int)nr_elements_moving>=nr_elements_movable){ // conditions for finishing the loop:
				if(dt>EPS){ // in time mode run until every element is getting randomized at least once
					__uint32_t ui=0;
					for(unsigned int ii=0; ii<ebf_size; ii++) ui|=element_bitfield[ii];
					setmore=(ui>0);
				} else setmore=false; // in non-time mode end after nr_elements elements are randomized
			}
		}
#else
		unsigned int nr_elements_moving=nr_elements; // old behavior for ran2 RNG
#endif
		unsigned int nr_moves=nr_elements_moving; // need this as array boundary to wrap around in case moves get rejected
		// update neighbors
		if(max_nr_moves_per_one>max_nr_neighbor_moves){
			max_nr_neighbor_moves=max_nr_moves_per_one;
			configuration->max_neighbor_trans=configuration->maxtrans*sqrt(max_nr_neighbor_moves)*NEIGHBOR_UPDATE_FREQ+configuration->group_max_centerdist;
			adjust_cut2plus(configuration);
			find_neighbors();
#if DEBUG_LEVEL>2
			cout << "\nMaximum number of single element moves increased to " << max_nr_moves_per_one << ".\n";
#endif
		} else if(step%NEIGHBOR_UPDATE_FREQ==0) find_neighbors();
		// Inner (oids) loop -> do n_oids moves every MC cycle, picking oids randomly
		for(unsigned int ii=0; ii<nr_moves; ii++){
			nr_tries++;
#ifdef USE_CMWC4096
			kk=move_elements[ii%nr_elements_moving];
			Elements[kk].nr_trials++;
#else
			kk = (unsigned int)((rotating_oids+nr_group_elements)*ran2(*configuration->idum));
			// now we need to remap that onto the possible element 1d-space ...
			if(kk>nr_group_elements){ // did not choose a group element
				temp=nr_group_elements; // group elements are created first
				for(unsigned int iii=0; iii<configuration->num_element_types; iii++){
					temp+=configuration->element_types[iii]->number;
					if(configuration->element_types[iii]->still==1){ // all of that type are still
						kk += configuration->element_types[iii]->number;
					} else{
						if (temp>=kk) break;
					}
				} // kk is now a random element's index which is allowed to rotate (and maybe even to move)
			}
#endif
			// update pair energies and test for bound things
			unsigned int start, e_start, curr_nr_neighbors;
			Vec3 start_center;
			e_start=(unsigned int)(!neighbor_switch)*max_nr_neighbors;
			if(Elements[kk].group){
				start=Elements[kk].group->number*configuration->max_neighbors;
				start_center=Elements[kk].group->center;
				curr_nr_neighbors=Elements[kk].group->nr_neighbors;
			} else{
				start=(nr_groups+kk-nr_group_elements)*configuration->max_neighbors;
				start_center=Elements[kk].center;
				curr_nr_neighbors=Elements[kk].nr_neighbors;
			}
			unsigned int e_idx;
			double ilj,ies;
			oldVLJ=0.0;
			oldES=0.0;
			// Loop over all neighbors
			for(unsigned int j=0; j<curr_nr_neighbors; j++){
				e_idx=energy_indices[start+j];
				ilj=pair_energies[e_idx];
				ies=pair_energies[e_idx+1];
				oldVLJ += ilj;
				oldES += ies;
				e_idx += e_start;
				pair_energy_store[e_idx++]=ilj;
				pair_energy_store[e_idx]=ies;
			}
			oldES *= configuration->in2;
			bool collective_move=false;
			if(configuration->collective_moves){
				find_bound_neighbors(start/configuration->max_neighbors);
#ifdef USE_CMWC4096
				collective_move=((nr_bound_neighbors>1) && (nr_bound_neighbors*ranQ()<1.0));
#else
				collective_move=((nr_bound_neighbors>1) && (nr_bound_neighbors*ran2(*configuration->idum)<1.0));
#endif
			}
			if(collective_move){
/*				cout << "Energy components (NkT):\nVmuE\t\tVES\t\tVLJ\t\tVG\t\tTotal\n";
				totVstep = 0.0;
				for(unsigned int jj=0; jj<NUM_V_STORE; jj++){
					totVstep += Vstep[jj];
					cout << Vstep[jj]*iNkT << "\t";
				}
				cout << totVstep*iNkT << "\n";*/
#if DEBUG_LEVEL>3
				cout << "Step " << step+1 << ": Collective move of " << nr_bound_neighbors << " entities\n";
#endif
				unsigned int nr_bound_elements=0;
				Vec3 bound_center(0.0);
				unsigned int count=0;
				for(unsigned int j=0; j<configuration->N; j++){
					if(bound_neighbors[j]){
						if(j<nr_groups){
							VE_store[j]=Groups[j]->VG;
							center_store[j]=Groups[j]->center;
							undo_PBCs(start_center,Groups[j]->center);
							bound_center+=Groups[j]->center;
							nr_bound_elements+=Groups[j]->nr_elements;
						} else{
							VE_store[j]=Elements[j-nr_groups+nr_group_elements].VE;
							center_store[j]=Elements[j-nr_groups+nr_group_elements].center;
							undo_PBCs(start_center,Elements[j-nr_groups+nr_group_elements].center);
							bound_center+=Elements[j-nr_groups+nr_group_elements].center;
							nr_bound_elements++;
						}
						count++;
						if(count==nr_bound_neighbors) break;
					}
				}
				bound_center/=nr_bound_neighbors;
				double inv_sqrt_bound_N=sqrt(1.0/((double)nr_bound_elements));
#ifdef USE_CMWC4096
				delta_trans.vec[0]=(configuration->maxtrans)*(2.0*ranQ()-1.0)*inv_sqrt_bound_N;
				delta_trans.vec[1]=(configuration->maxtrans)*(2.0*ranQ()-1.0)*inv_sqrt_bound_N;
				delta_trans.vec[2]=(configuration->maxtrans)*(2.0*ranQ()-1.0)*inv_sqrt_bound_N;
#else
				delta_trans.vec[0]=(configuration->maxtrans)*(2.0*ran2(*configuration->idum)-1.0)*inv_sqrt_bound_N;
				delta_trans.vec[1]=(configuration->maxtrans)*(2.0*ran2(*configuration->idum)-1.0)*inv_sqrt_bound_N;
				delta_trans.vec[2]=(configuration->maxtrans)*(2.0*ran2(*configuration->idum)-1.0)*inv_sqrt_bound_N;
#endif
				rmovtemp=0.0;
				tmovtemp=0.0;
				//Random number for picking an axis and for the move itself
#ifdef USE_CMWC4096
				unsigned int choice=CMWC4096()%3;
				double theta = configuration->maxrot*(2.0*ranQ()-1.0)*inv_sqrt_bound_N;
#else
				unsigned int choice=(unsigned int)ran2_int(*configuration->idum)%3;
				double theta = configuration->maxrot*(2.0*ran2(*configuration->idum)-1.0)*inv_sqrt_bound_N;
#endif
				double co = cos(theta);
				double si = sin(theta);
				
				Mat33 rot;
				switch(choice){
					case 0: /* yaw: x-rot
						 * (1   0   0)
						 * (0  co  si)
						 * (0 -si  co)
						 */
						rot.mat[0][0] = 1;	rot.mat[0][1] = 0;	rot.mat[0][2] = 0;
						rot.mat[1][0] = 0;	rot.mat[1][1] = co;	rot.mat[1][2] = si;
						rot.mat[2][0] = 0;	rot.mat[2][1] = -si;	rot.mat[2][2] = co;
						break;
					case 1: /* pitch: y-rot
						 * ( co  0  si)
						 * ( 0   1   0)
						 * (-si  0  co)
						 */
						rot.mat[0][0] = co;	rot.mat[0][1] = 0;	rot.mat[0][2] = si;
						rot.mat[1][0] = 0;	rot.mat[1][1] = 1;	rot.mat[1][2] = 0;
						rot.mat[2][0] = -si;	rot.mat[2][1] = 0;	rot.mat[2][2] = co;
						break;
					case 2: /* roll: z-rot
						 * ( co  si  0)
						 * (-si  co  0)
						 * ( 0   0   1)
						 */
						rot.mat[0][0] = co;	rot.mat[0][1] = si;	rot.mat[0][2] = 0;
						rot.mat[1][0] = -si;	rot.mat[1][1] = co;	rot.mat[1][2] = 0;
						rot.mat[2][0] = 0;	rot.mat[2][1] = 0;	rot.mat[2][2] = 1;
						break;
				}
				double omuE=0.0;
				double oVE=0.0;
				double oVLJ=0.0;
				double oES=0.0;
				count=0;
				Vec3 dtr;
				for(unsigned int j=0; j<configuration->N; j++){
					if(bound_neighbors[j]){
						unsigned int cstart=j*configuration->max_neighbors;
						if(j<nr_groups){ // group
							if(configuration->LJwall_calc) oVE+=Groups[j]->VG;
							curr_nr_neighbors=Groups[j]->nr_neighbors;
							for(unsigned int k=0; k<Groups[j]->nr_elements; k++){
								// Calculate current muE energy and total the energies
								if((step >= configuration->randsteps) || configuration->rand_Efield) omuE += VmuE(Groups[j]->elements[k]);
							}
							Rotate(Groups[j],bound_center,rot);
							Translate(Groups[j],delta_trans);
							// apply boundary conditions
							dtr=Groups[j]->center; // initial
							apply_PBCs(Groups[j]->center);
							dtr-=Groups[j]->center; // - final
							// Move whole group according to boundary conditions (can't use Translate because it would change group center)
							for(unsigned int k=0; k<Groups[j]->nr_elements; k++) Elements[Groups[j]->elements[k]].center-=dtr;
						} else{
							if(configuration->LJwall_calc) oVE+=Elements[j-nr_groups+nr_group_elements].VE;
							curr_nr_neighbors=Elements[j-nr_groups+nr_group_elements].nr_neighbors;
							// Calculate current muE energy and total the energies
							if((step >= configuration->randsteps) || configuration->rand_Efield) omuE += VmuE(j-nr_groups+nr_group_elements);
							Rotate(&Elements[j-nr_groups+nr_group_elements],bound_center,rot);
							Elements[j-nr_groups+nr_group_elements].center+=delta_trans;
							apply_PBCs(Elements[j-nr_groups+nr_group_elements].center);
						}
						// Loop over all neighbors
						for(unsigned int l=0; l<curr_nr_neighbors; l++){
							e_idx=energy_indices[cstart+l];
							ilj=pair_energies[e_idx];
							ies=pair_energies[e_idx+1];
							pair_energy_store[e_start+e_idx]=ilj;
							pair_energy_store[e_start+e_idx+1]=ies;
							oVLJ += ilj;
							oES += ies;
						}
						count++;
						if(count==nr_bound_neighbors) break;
					}
				}
				oES *= configuration->in2;
				// Metropolis Conditional:
				bool success=true;
				double nES, nVLJ;
				double muE=0.0;
				double VES=0.0;
				double VLJ=0.0;
				double VE=0.0;
//				find_neighbors();
				count=0;
				for(unsigned int j=0; j<configuration->N; j++){
					if(bound_neighbors[j]){
						if(j<nr_groups){ // group
							if(configuration->LJwall_calc) Groups[j]->VG = GroupPotentialEnergy(Groups[j]);
							for(unsigned int k=0; k<Groups[j]->nr_elements; k++){
								// Calculate current muE energy and total the energies
								if((step >= configuration->randsteps) || configuration->rand_Efield) muE += VmuE(Groups[j]->elements[k]);
								// Calculate mu-mu and LJ energy
								success&=OneStepGroup(Groups[j]->elements[k], nES, nVLJ, skip_charges);
								if(configuration->LJwall_calc) Groups[j]->VG+=Vwall(Groups[j]->elements[k],skip_charges,success);
								VES+=nES;
								VLJ+=nVLJ;
							}
							if(configuration->LJwall_calc) VE+=Groups[j]->VG;
						} else{
							// Calculate current muE energy and total the energies
							if((step >= configuration->randsteps) || configuration->rand_Efield) muE += VmuE(j-nr_groups+nr_group_elements);
							// Calculate mu-mu and LJ energy
							success&=OneStep(j-nr_groups+nr_group_elements, nES, nVLJ, skip_charges);
							if(configuration->LJwall_calc){
								Elements[j-nr_groups+nr_group_elements].VE=Vwall(j-nr_groups+nr_group_elements,skip_charges,success);
								VE+=Elements[j-nr_groups+nr_group_elements].VE;
							}
							VES+=nES;
							VLJ+=nVLJ;
						}
						count++;
						if(count==nr_bound_neighbors) break;
					}
				}
				double deltaV = muE+VES+VLJ+VE - (omuE+oES+oVLJ+oVE);
//				cout << omuE << ", " << oES << ", " << oVLJ << " -> " << muE << ", " << VES << ", " << VLJ << " = " << deltaV << "\n";
#ifdef USE_CMWC4096
				double rand = ranQ();
#else
				double rand = ran2(*configuration->idum);
#endif
				if(success && (exp(-deltaV*configuration->beta) > rand)){
					// Record squared displacement
					tmovstep += tmovtemp;
					rmovstep += rmovtemp;
#if DEBUG_LEVEL>3
					cout << "<- Accepted\n";
#endif
					// acceptance tracking [vdw preeq eq]
					if (step>=eqsteps) accept[2]++;
						else if (step>=configuration->randsteps) accept[1]++; else accept[0]++;
					// Record instantaneous value of total energy
					Vstep[0] += (muE - omuE);
					Vstep[1] += (VES - oES);
					Vstep[2] += (VLJ - oVLJ);
					Vstep[3] += (VE - oVE);
				} else{ // Revert things back if move is rejected
#if DEBUG_LEVEL>3
					cout << "<- Not accepted\n";
#endif
					count=0;
					// reverse changes
					rot=rot.M3Transpose();
					for(unsigned int j=0; j<configuration->N; j++){
						if(bound_neighbors[j]){
							unsigned int cstart=j*configuration->max_neighbors;
							if(j<nr_groups){ // group
								curr_nr_neighbors=Groups[j]->nr_neighbors;
								Translate(Groups[j],center_store[j]-Groups[j]->center);
								Rotate(Groups[j],rot);
								Groups[j]->VG=VE_store[j];
							} else{
								curr_nr_neighbors=Elements[j-nr_groups+nr_group_elements].nr_neighbors;
								Elements[j-nr_groups+nr_group_elements].center=center_store[j];
								Elements[j-nr_groups+nr_group_elements].VE=VE_store[j];
								RotateElement(&Elements[j-nr_groups+nr_group_elements],rot);
							}
							// Loop over all neighbors
							for(unsigned int l=0; l<curr_nr_neighbors; l++){
								e_idx=energy_indices[cstart+l];
								pair_energies[e_idx]=pair_energy_store[e_start+e_idx];
								pair_energies[e_idx+1]=pair_energy_store[e_start+e_idx+1];
							}
							count++;
							if(count==nr_bound_neighbors) break;
						}
					}
//					find_neighbors();
				}
			} else{
				if(!Elements[kk].group){ // individual element
					// create copy of chosen element's center vector and rotation matrix
					if(configuration->stepsize_average>0){
						// <trans^2> = int_-trans^+trans P(e) e^2 de = 2/(2*trans)*int_0^trans e^2 de = 1/3*trans^2
						// <trans^2> = <trans_x^2+trans_y^2+trans_z^2> = 3/3*trans^2 => trans = sqrt(<trans^2>)
						configuration->maxtrans=1.05*sqrt(avg_trans[kk-nr_group_elements]/configuration->stepsize_average);
						// <rot^2> = 1/3*rot^2 => rot = sqrt(3*<rot^2>)
						configuration->maxrot=1.05*sqrt(3.0*avg_rot[kk-nr_group_elements]/configuration->stepsize_average);
						// enforce upper limits
						if(configuration->maxtrans>halfbox) configuration->maxtrans=halfbox;
						if(configuration->maxrot>pi) configuration->maxrot=pi;
//						if(kk==0) cout << configuration->maxtrans << " : " << configuration->maxrot/pi << " pi (" << avg_trans[ikk]/configuration->stepsize_average << ":" << avg_rot[ikk]/configuration->stepsize_average << ")\n";
					}
					Vec3 old_dipole = Elements[kk].dipole;
					Vec3 old_center = Elements[kk].center;
					Vec3 old_s = Elements[kk].current_s;
					Mat33 old_rot = Elements[kk].rot;
					double old_gamma_quarter=Elements[kk].gamma_quarter;
					double oldVE = Elements[kk].VE;
					// Calculate current muE energy and total the energies
					if((step >= configuration->randsteps) || configuration->rand_Efield) oldmuE = VmuE(kk);
					double oldEkin=0.0;
					if(configuration->time){
						// calculate old kinetic energy first (with new time step dt hence the recalculation
						m=Elements[kk].MyType->mass*amu_to_perg_ps2_per_Ang2;
						v=Elements[kk].ds/dt;
						oldEkin = 0.5*m*(v*v);
						// set move size and apply trial move
						ds=Elements[kk].ds/Elements[kk].nr_moves;
						Elements[kk].center+=ds;
						Elements[kk].current_s+=ds;
						tmovtemp=ds*ds;
						// how much can we move within Heisenberg's uncertainty principle
						max_step=hbar_perg_ps/(2.0*m*average(Elements[kk].MyType->saxes.vec,3))*dt;
						// increase max_step by upto the current displacement |ds| the longer we try
						if(Elements[kk].nr_trials>Elements[kk].nr_moves){
							double sigma=atan((Elements[kk].nr_trials-Elements[kk].nr_moves)/10.0)*2.0/pi;
							max_step+=sqrt(tmovtemp)*sigma;
						}
						for(unsigned int iii=0; iii<3; iii++){
							double trans=ran_n()*max_step;
							Elements[kk].center.vec[iii]+=trans;
							ds.vec[iii]+=trans;
							Elements[kk].current_s.vec[iii]+=trans;
							tmovtemp+=trans*trans;
						}
						apply_PBCs(Elements[kk].center);
					} else{
						max_step=0.0;
						if(!Elements[kk].MyType->rot_notrans){
							if (configuration->latticetype == 1) tmovtemp = Translate(kk);
								else if (configuration->latticetype == 3) tmovtemp = Translate_Sphere(kk);
									else if (configuration->latticetype == 4) tmovtemp = Translate_Cylinder(kk);
										else tmovtemp = Translate(kk);
						}
						rmovtemp = Rotate(kk);
					}
					apply_PBCs(Elements[kk].center);
					if(configuration->vdwtype==6) Elements[kk].gamma_quarter=ranQ()*4.0-1.0;
					// Recalculate energy
					success=OneStep(kk, VES, VLJ, skip_charges);
					if(configuration->LJwall_calc) Elements[kk].VE=Vwall(kk,skip_charges,success);
					OldV = oldmuE + oldES + oldVLJ + oldEkin + oldVE;
//					cout << step << ": " << oldES << ", " << VES << " ; " << oldVLJ << ", " << VLJ << " (" << kk << ")\n";
					double Ekin=0.0;
					if(configuration->time){
						v=Elements[kk].current_s/(dt/Elements[kk].nr_moves);
						Ekin=0.5*m*(v*v);
						// in order to not violate detailed balance we need to take care of the change in stepsizes
						qterms = 0.5*m*(Elements[kk].max_step*Elements[kk].max_step - max_step*max_step)/(dt*dt);
					}
					// Calculate new muE energy and total energies
					if((step >= configuration->randsteps) || configuration->rand_Efield) muE = VmuE(kk);
					NewV = muE + VES + VLJ + Ekin + Elements[kk].VE;
					// Metropolis Conditional:
					deltaV = NewV - OldV + qterms;
#ifdef USE_CMWC4096
					double rand = ranQ();
#else
					double rand = ran2(*configuration->idum);
#endif
					if(success && (exp(-deltaV*configuration->beta) > rand)){
						// Record squared displacement
						tmovstep += tmovtemp;
						rmovstep += rmovtemp;
						Elements[kk].max_step=max_step;
						if(configuration->stepsize_average>0){
							unsigned int base=running_average[kk-nr_group_elements]*nr_entities;
							avg_trans[kk-nr_group_elements]+=tmovtemp-trans_store[base+kk-nr_group_elements];
							avg_rot[kk-nr_group_elements]+=rmovtemp-rot_store[base+kk-nr_group_elements];
							trans_store[base+kk-nr_group_elements]=tmovtemp;
							rot_store[base+kk-nr_group_elements]=rmovtemp;
							running_average[kk-nr_group_elements]++;
							running_average[kk-nr_group_elements]%=configuration->stepsize_average;
						}
						// acceptance tracking [vdw preeq eq]
						if(step>=eqsteps) accept[2]++;
							else if (step>=configuration->randsteps) accept[1]++; else accept[0]++;
						// Record instantaneous value of total energy
						Vstep[0] += (muE-oldmuE);
						Vstep[1] += (VES - oldES);
						Vstep[2] += (VLJ - oldVLJ);
						Vstep[3] += (Elements[kk].VE - oldVE);
					} else{ // Revert things if move is rejected
						Elements[kk].dipole = old_dipole;
						Elements[kk].center = old_center;
						Elements[kk].current_s = old_s;
						Elements[kk].rot = old_rot;
						Elements[kk].VE = oldVE;
						Elements[kk].gamma_quarter = old_gamma_quarter;
						for(unsigned int j=0; j<Elements[kk].MyType->nr_charges; j++) Elements[kk].q_pos[j]=Elements[kk].rot*Elements[kk].MyType->q_pos[j];
						if(configuration->time){
							move_elements[nr_moves%nr_elements_moving]=kk; // add this element to the rerun list
							nr_moves++;
						}
						// Loop over all of our neighbors
						for(unsigned int j=start; j<start+Elements[kk].nr_neighbors; j++){
							e_idx=energy_indices[j];
							pair_energies[e_idx]=pair_energy_store[e_start+e_idx];
							pair_energies[e_idx+1]=pair_energy_store[e_start+e_idx+1];
						}
						tmovtemp=0.0;
						rmovtemp=0.0;
					}
				} else{ // group
					/* Special attention is necessary if we found ourselves a molecule (elements connected to each other):
					 * basic idea:
					 * - allow whole molecule to rotate and translate slowly (a.k.a. adiabatically)
					 * - allow n elements of molecule to "adapt" randomly with Metropolis sampling for whole molecule afterwards
					 */
					Element* element;
					Element_Group* group=Elements[kk].group;
					if(configuration->stepsize_average>0){
						// <trans^2> = int_-trans/sqrt(Nr)^+trans/(sqrt(Nr) P(e) e^2 de = 2/(2*trans/sqrt(Nr))*int_0^trans/sqrt(Nr) e^2 de = 1/3*trans^2/Nr
						// <trans^2> = Nr*<trans_x^2+trans_y^2+trans_z^2> = Nr*3/3*trans^2/Nr => trans = sqrt(<trans^2>)
						configuration->maxtrans=1.05*sqrt(avg_trans[group->number]/configuration->stepsize_average);
						// for group rotation: <rot^2> = 1/3*rot^2 => rot = sqrt(3*<rot^2>)
						// for group randomization:
						//     <rot^2> = <sum rot^2*Nr> = Nr*int_-rot/Nrand^rot/Nrand P(e) sum e_i^2 de = Nr^2 2/(2*rot/Nrand) int_0^rot/Nrand e^2 de = Nr^2 1/(3*rot/Nrand) (rot/Nrand)^3
						//             = 1/3*Nr^2*(rot/Nrand)^2 = 1/3*rot^2*(Nr/Nrand)^2 (<- 1/group_rand_frac^2)
						//  => rot = sqrt(3*<rot^2>*group_rand_frac^2)
						configuration->maxrot=1.05*sqrt(3.0*avg_rot[group->number]/configuration->stepsize_average);
						// enforce upper limits
						if(configuration->maxtrans>halfbox) configuration->maxtrans=halfbox;
						if(configuration->maxrot>pi) configuration->maxrot=pi;
//						if(group->number==0) cout << configuration->maxtrans << " : " << configuration->maxrot << " (" << avg_trans[group->number] << ":" << avg_rot[group->number] << ")\n";
					}
					unsigned int grouptype=group->type;
					Vec3 old_groupcenter = group->center;
					// Store old energies and calculate current position dependent portions
					if(adjust_internal) selfInteraction(group,group->VES,&(group->VLJ),skip_charges); // recalculate group-internal potential energies if they depend on volume-based scaling
					double selfVES=group->VES;
					double selfVLJ=group->VLJ;
					double selfVG=group->VG;
					
					oldmuE = 0.0;
					
//					double oES, oVLJ;
					for(unsigned int i=0; i<group->nr_elements; i++){
						element=&Elements[group->elements[i]];
						// First things first: save group state one element at a time
						group_storage[grouptype][i].center=element->center;
						group_storage[grouptype][i].dipole=element->dipole;
						group_storage[grouptype][i].rot=element->rot;
						group_storage[grouptype][i].gamma_quarter=element->gamma_quarter;
						// Calculate and save initial potentials (this works here, because potentials are updated later done the road ...)
						for(unsigned int j=0; j<element->nr_interactions; j++){
							if(element->interactions[j].location) group_storage[grouptype][i].interaction_location[j]=*element->interactions[j].location;
							if(element->interactions[j].normal) group_storage[grouptype][i].interaction_normal[j]=*element->interactions[j].normal;
							if(element->interactions[j].tangent) group_storage[grouptype][i].interaction_tangent[j]=*element->interactions[j].tangent;
						}
						// Save charge locations
						for(unsigned int j=0; j<element->MyType->nr_charges; j++) group_storage[grouptype][i].q_pos[j] = element->q_pos[j];
						// Calculate current muE energy and total the energies
						if((step >= configuration->randsteps) || configuration->rand_Efield) oldmuE += VmuE(group->elements[i]);
					}
					oldES += selfVES;
					oldVLJ += selfVLJ;
					OldV = selfVG + oldmuE + oldES + oldVLJ;
					// decide if we want to do rigid body rotation/translation *xor* internal changes
					bool rigidbody_nointernal=true;
#ifdef USE_CMWC4096
					if(group->Type->rand_elements>0) rigidbody_nointernal=(ranQ()>0.5);
#else
					if(group->Type->rand_elements>0) rigidbody_nointernal=(ran2(*configuration->idum)>0.5);
#endif
					rigidbody_nointernal|=configuration->randsteps_nointernal && (step<configuration->randsteps);
					rigidbody_nointernal&=!group->Type->still;
					Vec3 center=group->center;
					double rf2=1.0;
					if(rigidbody_nointernal){
						rmovtemp=0.0;
#if DEBUG_LEVEL>3
						cout << "Rotate ";
#endif
						rmovtemp=Rotate(group,center);
						/* random translation for whole group scaled by 1/sqrt(# of group elements)
						 * if some of them can't move that's OK (because I say so)
						 * - we don't want to loop over elements to find out how many can move and sqrt again
						 * - the ones that can't move technically also screw up rotation, however, since
						 *   displacements are small energy penalty should nicely take care of that problem
						 * - rotating element whose center is not the rotation center also translates it
						 * - the rest of the group still needs to rot/trans to maybe find lower energy configuration
						 * - that's why ignoring the "immovables" is not the dumbest thing to do ;-)
						 */
						if(!group->Type->rot_notrans){
#if DEBUG_LEVEL>3
							cout << "Translate ";
#endif
#ifdef USE_CMWC4096
							delta_trans.vec[0]=(configuration->maxtrans)*(2.0*ranQ()-1.0)*group->Type->inv_sqrt_nrpc;
							delta_trans.vec[1]=(configuration->maxtrans)*(2.0*ranQ()-1.0)*group->Type->inv_sqrt_nrpc;
							delta_trans.vec[2]=(configuration->maxtrans)*(2.0*ranQ()-1.0)*group->Type->inv_sqrt_nrpc;
#else
							delta_trans.vec[0]=(configuration->maxtrans)*(2.0*ran2(*configuration->idum)-1.0)*group->Type->inv_sqrt_nrpc;
							delta_trans.vec[1]=(configuration->maxtrans)*(2.0*ran2(*configuration->idum)-1.0)*group->Type->inv_sqrt_nrpc;
							delta_trans.vec[2]=(configuration->maxtrans)*(2.0*ran2(*configuration->idum)-1.0)*group->Type->inv_sqrt_nrpc;
#endif
#if DEBUG_LEVEL>3
							cout << "(" << delta_trans.V3Str(',') << ") ";
#endif
							Translate(group,delta_trans);
						}
						rf2/=group->Type->inv_sqrt_nrpc;
						if(configuration->vdwtype==6){
							for(unsigned int i=0; i<group->nr_elements; i++) Elements[group->elements[i]].gamma_quarter=ranQ()*4.0-1.0;
						}
					} else{
#if DEBUG_LEVEL>3
						cout << "Shuffle ";
#endif
						rmovtemp=GroupShuffle(group,kk);
						rf2=group->Type->rand_elements;
					}
					rf2*=group->inv_nr_elements;
					// Apply boundary conditions at latest possible time (before potential calculation)
					group->center.vec[0]=0.0; group->center.vec[1]=0.0; group->center.vec[2]=0.0;
					tmovtemp=0.0;
					for(unsigned int i=0; i<group->nr_elements; i++){
						// calculate translation displacement from new position
						element=&Elements[group->elements[i]];
						Vec3 mt=group_storage[grouptype][i].center-element->center;
						tmovtemp+=mt*mt;
						 // calculate geometric group center while we're here
						group->center+=element->center;
					}
					group->center*=group->inv_nr_elements;
					if(group->Type->rot_notrans || group->Type->still){ // in case we don't want to translate the group center (just rotate) from internal motion ...
						delta_trans=group->center-center; // we need to adjust the new center back to the old one with delta_trans
						group->center=center;
					} else{
						delta_trans=group->center; // initial
						apply_PBCs(group->center);
						delta_trans-=group->center; // - final
					}
					// Move whole group according to boundary conditions (can't use Translate because it would change group center)
					for(unsigned int i=0; i<group->nr_elements; i++) Elements[group->elements[i]].center-=delta_trans;
					
					// Calculate new potential energies (internal ones only need to be recalculated if group has been shuffled ...)
					if(rigidbody_nointernal && !configuration->LJwall_calc){
						success=true;
					} else{
						group->VG = GroupPotentialEnergy(group);
						success=selfInteraction(group,group->VES,&(group->VLJ),skip_charges);
					}
					muE = 0.0;
					VES = group->VES;
					VLJ = group->VLJ;
					
					double nES, nVLJ;
					for(unsigned int i=0; i<group->nr_elements; i++){
						// Calculate current muE energy and total the energies
						if((step >= configuration->randsteps) || configuration->rand_Efield) muE += VmuE(group->elements[i]);
						// Calculate mu-mu and LJ energy
						success&=OneStepGroup(group->elements[i], nES, nVLJ, skip_charges);
						VES += nES;
						VLJ += nVLJ;
						if(configuration->LJwall_calc) group->VG+=Vwall(group->elements[i],skip_charges,success);
					}
					NewV = group->VG + muE + VES + VLJ;
					
					// Metropolis Conditional:
					deltaV = NewV - OldV;
#ifdef USE_CMWC4096
					double rand = ranQ();
#else
					double rand = ran2(*configuration->idum);
#endif
					if(success && (exp(-deltaV*configuration->beta) > rand)){
						// Record squared displacement
						tmovstep += tmovtemp;
						rmovstep += rmovtemp;
#if DEBUG_LEVEL>3
						cout << "\nAccepted move: " << kk+1 << ", " << VLJ << ", " << tmovtemp << ":" << rmovtemp << "\n";
#endif
						if(configuration->stepsize_average>0){
							unsigned int base=running_average[group->number]*nr_entities;
							avg_trans[group->number]+=tmovtemp*rf2-trans_store[base+group->number];
							avg_rot[group->number]+=rmovtemp*rf2-rot_store[base+group->number];
							trans_store[base+group->number]=tmovtemp*rf2;
							rot_store[base+group->number]=rmovtemp*rf2;
							running_average[group->number]++;
							running_average[group->number]%=configuration->stepsize_average;
						}
						// acceptance tracking [vdw preeq eq]
						if (step>=eqsteps) accept[2]++;
							else if (step>=configuration->randsteps) accept[1]++; else accept[0]++;
						// Record instantaneous value of total energy
						Vstep[0] += (muE - oldmuE);
						Vstep[1] += (VES - oldES);
						Vstep[2] += (VLJ - oldVLJ);
						Vstep[3] += (group->VG - selfVG);
					} else{ // Revert things back if move is rejected
#if DEBUG_LEVEL>3
						cout << "\nMove not accepted";
#endif
						group->center=old_groupcenter;
						group->VES=selfVES;
						group->VLJ=selfVLJ;
						group->VG=selfVG;
						for(unsigned int i=0; i<group->nr_elements; i++){
							element=&Elements[group->elements[i]];
							element->center=group_storage[grouptype][i].center;
							element->dipole=group_storage[grouptype][i].dipole;
							element->rot=group_storage[grouptype][i].rot;
							element->gamma_quarter=group_storage[grouptype][i].gamma_quarter;
							for(unsigned int j=0; j<element->nr_interactions; j++){
								if(element->interactions[j].location) *element->interactions[j].location=group_storage[grouptype][i].interaction_location[j];
								if(element->interactions[j].normal) *element->interactions[j].normal=group_storage[grouptype][i].interaction_normal[j];
								if(element->interactions[j].tangent) *element->interactions[j].tangent=group_storage[grouptype][i].interaction_tangent[j];
							}
							for(unsigned int j=0; j<element->MyType->nr_charges; j++) element->q_pos[j]=group_storage[grouptype][i].q_pos[j];
						}
						// Loop over all of our neighbors
						for(unsigned int j=start; j<start+group->nr_neighbors; ++j){
							e_idx=energy_indices[j];
							pair_energies[e_idx]=pair_energy_store[e_start+e_idx];
							pair_energies[e_idx+1]=pair_energy_store[e_start+e_idx+1];
						}
						tmovtemp=0.0;
						rmovtemp=0.0;
					}
				}
			}
		}// end oids loop.
		if(configuration->transition){
			LJlambda=calc_lambda(step);
			if(fabs(LJlambda-1.0)<EPS){
				configuration->transition=false;
				LJlambda=1.0;
			}
			double newvolume=configuration->transition_start_volume-LJlambda*(configuration->transition_start_volume-configuration->transition_target_volume);
			LJlambda*=configuration->transition_LJ_delta;
			LJlambda+=configuration->transition_LJ_start;
			if((fabs(newvolume-configuration->V)>EPS) && (step>=configuration->transition_start)){
				double factor;
				double xfactor;
				if(configuration->LJwall_calc && configuration->LJwall_fixed){
					xfactor=1.0/configuration->LJwall_xm;
					configuration->LJwall_xm=configuration->LJwall_x*LJlambda+LJwall_start_trans*(1.0-LJlambda);
					xfactor*=configuration->LJwall_xm;
					factor=sqrt((newvolume*configuration->boxlength[0])/(configuration->V*2.0*configuration->LJwall_xm));
				} else{
					factor=cbrt(newvolume/configuration->V);
					xfactor=factor;
					configuration->LJwall_xm=0.5*configuration->boxlength[0]*xfactor;
				}
				// make changes "official"
				switch(configuration->latticetype){
					case 1:
					case 2: // Rectangular box
						if(configuration->LJwall_calc && configuration->LJwall_fixed) configuration->boxlength[0]=2.0*configuration->LJwall_xm; else configuration->boxlength[0]*=factor;
						configuration->inv_boxlength[0]=1.0/configuration->boxlength[0];
						configuration->boxlength[1]*=factor;
						configuration->inv_boxlength[1]=1.0/configuration->boxlength[1];
						configuration->boxlength[2]*=factor;
						configuration->inv_boxlength[2]=1.0/configuration->boxlength[2];
						break;
					default: // Default condition
						cout << "Cannot change volume. Lattice type " << configuration->latticetype << " not recognized. Exiting.\n";
						exit(1);
				}
				// loop over elements which are not in groups
				for(unsigned int k = nr_group_elements; k < nr_elements; k++){
					Elements[k].center.vec[0]*=xfactor;
					Elements[k].center.vec[1]*=factor; Elements[k].center.vec[2]*=factor;
				}
				// loop over groups
				Element_Group* group;
				for(unsigned int k = 0; k < nr_groups; k++){
					group=Groups[k];
					delta_trans=group->center;
					delta_trans.vec[0]*=xfactor;
					delta_trans.vec[1]*=factor; delta_trans.vec[2]*=factor;
					delta_trans-=group->center; // new position - old position points to new position
					group->center+=delta_trans; // move group center but leave group elements undisturbed
					Vec3 safety_move(0.0);
					for(unsigned int i=0; i<group->nr_elements; i++){
						Elements[group->elements[i]].center+=delta_trans;
						for(unsigned int j=0; j<3; j++){
							if(!configuration->PBCs[j]){
								double db=0.5*configuration->boxlength[j]-fabs(Elements[group->elements[i]].center.vec[j]);
								if(db<0.0){
									if(fabs(db)>safety_move.vec[j]) safety_move.vec[j]=fabs(db);
								}
							}
						}
					}
					if(safety_move*safety_move>0.0){
						for(unsigned int i=0; i<group->nr_elements; i++){
							for(unsigned int j=0; j<3; j++){
								if(Elements[group->elements[i]].center.vec[j]<0.0) Elements[group->elements[i]].center.vec[j]+=safety_move.vec[j]+1.0; else Elements[group->elements[i]].center.vec[j]-=safety_move.vec[j]-1.0;
							}
						}
					}
				}
				configuration->nndist*=factor;
				configuration->V=newvolume;
				// escut2, the square of the electrostatics cutoff, is the only thing used during runtime (for distances, r*r without the additional sqrt saves time)
				configuration->rfcut=configuration->escut*cbrt(newvolume/configuration->transition_target_volume);
				configuration->escut2=qqpwr(configuration->rfcut,2);
				configuration->in2=qqpwr(configuration->mucorr,2)*configuration->energy_scale[1]/(configuration->n2+(1.0-LJlambda)*configuration->transition_delta_n2); // only adjust inverse n^2 (leave eps_RF alone) to scale all electrostatics
				success=recalcall(Vstep, skip_charges); // dipole/charges are not calculated if in randomization
				if(!success){
					cout << "Unrecoverable overlap while changing volume.\n";
					exit(3);
				}
			}
			if(!configuration->transition){
				if(configuration->LJwall_calc && configuration->LJwall_fixed) configuration->LJwall_xm=configuration->LJwall_x;
				configuration->targetV=configuration->transition_target_volume;
				cout << "\t\t\tTransition completed at " << configuration->V << " Angström\n";
			}
		} else{ // NpT move ?
			bool NpT_move=false;
			if((!configuration->transition && configuration->NpT) && (!configuration->randsteps_novolume || (configuration->randsteps_novolume && (step>=configuration->randsteps)))){
#ifdef USE_CMWC4096
				NpT_move=(ranQ()*configuration->NpT_move_nr<nr_moves);
#else
				NpT_move=(ran2(*configuration->idum)*configuration->NpT_move_nr<nr_moves);
#endif
				if(NpT_move){
					tmovstep+=NpT_step(Vstep, save_groupVG, step, eqsteps, adjust_internal);
				}
			}
		}
		// Output section. Saves data for processing and dumps messages to standard out.
		// Calculate cos^n, M, V
		get_means(cosmeans[step],overallM);
		Ms[step]=overallM;
		Vs[step]=configuration->V;
		for(unsigned int jj=0; jj<NUM_V_STORE; jj++) totalV[step][jj] = Vstep[jj];
		if(configuration->time){
			mds2=0.0;
			for(unsigned int i=0; i<nr_elements; i++){
				mds2+=Elements[i].MyType->mass*amu_to_perg_ps2_per_Ang2*(Elements[i].current_s*Elements[i].current_s);
				Elements[i].ds=Elements[i].current_s;
				Elements[i].current_s=Vec3(0.0);
			}
			if(mds2>EPS){
				dt=sqrt(mds2/(2.0*configuration->Ekin));
#if DEBUG_LEVEL>3
				cout << "dt = " << dt << "\n";
#endif
				simulation_time+=dt;
			}
		}
		// Calculate mean squared displacement for step
		msmoved[0] += tmovstep/(moving_oids+nr_group_elements);
		msmoved[1] += rmovstep/(rotating_oids+nr_group_elements);
		// check for unexpected program termination
		if(exitnow){
			configuration->Runtime=((double)clock()-tstart)/CLOCKS_PER_SEC;
			configuration->storetime=true;
			store_trajectory(splus1);
			cout << "\n";
			delete this;
			exit(42);
		}
		// Output progress indicator to user, including total energy trace.
		if (splus1 % 100 == 0){
			nMstep=0.0;
			if(fabs(configuration->max_dipole)>EPS) nMstep = overallM.V3Norm()/(configuration->max_dipole*configuration->mucorr);
			totVstep = 0.0;
			for(unsigned int jj=0; jj<NUM_V_STORE; jj++) totVstep += Vstep[jj];
			totVstep *= iNkT;
			cout << "\n" << splus1  << "\t\t" << totVstep << "\t\t\t" << nMstep << "\n";
#if DEBUG_LEVEL>2
			cout << "(" << Vstep[0]*iNkT << "\t" << Vstep[1]*iNkT << "\t" << Vstep[2]*iNkT << "\t" << Vstep[3]*iNkT << "\t" << overallM.V3Norm() << ")\n";
#endif
#if DEBUG_LEVEL>2
			double Ves=Vstep[1];
			recalcall(Vstep,skip_charges);
			if(fabs(Ves-Vstep[1])>EPS) cout << "(" << Vstep[0]*iNkT << "\t" << Vstep[1]*iNkT << "\t" << Vstep[2]*iNkT << "\t" << Vstep[3]*iNkT << "\t" << overallM.V3Norm() << ")\n";
			for(unsigned int jj=0; jj<NUM_V_STORE; jj++) Vstep[jj]=totalV[step][jj];
			cout << "\t\t(" << calc_charge_interaction()*iNkT << ")\n";
#endif
			if (configuration->anneal == 1 && step < configuration->randsteps*2){
				cout << "\t\t\tT = " << configuration->kT/kB << " K\n";
			}
			if (configuration->dyneps == 1 && step > configuration->randsteps){
				cout << "\t\t\tEpsilon = " << epsRF << "\n";
			}
			if (configuration->realfield && (step>configuration->randsteps) && (Evec.V3Norm()>EPS)){
				cout << "\t\t\t|E_0| = " << Evec.V3Norm()/MV_per_m_to_perg_per_Debye << " V/µm\n";
			}
			if(changevolume || configuration->NpT || configuration->transition) cout << "\t\t\tV = " << configuration->V << " Angström³\n";
#if DEBUG_LEVEL>1
			if(configuration->NpT && ((step>=configuration->randsteps) && (configuration->cut && configuration->NpT_adjust_EScut))) cout << "\t\t\tElectrostatics cutoff is at " << configuration->rfcut << " Angström\n";
#endif
			if(configuration->transition) cout << "\t\t\tSolvent transition lambda = " << LJlambda << "\n";
			
			if(configuration->NpT && adjust_internal) cout << "\t\t\tInternal LJ dispersion adjustment lambda = " << selfLJlambda*LJlambda << "\n";
			
			cout.flush();
		}
		// set time if this is the last step
		if(step+1==configuration->steps){
			configuration->storetime=true;
			configuration->Runtime=((double)clock()-tstart)/CLOCKS_PER_SEC;
		}
		// Save intermediate coordinates. Note that this can require large amounts of disk space
		if((configuration->chkcoord) && (splus1 % configuration->grfreq==0)) store_trajectory(splus1);
		
		// adjust time step (dt) and trial move sizes (ds) - needs to happen after data output b/c output uses ds and dt
		if(dt>EPS){
			double avg_v=0.0;
			for(unsigned int i=0; i<nr_elements; i++) avg_v+=Elements[i].ds.V3Norm()/dt; // speed in Angström/ps
			avg_v/=nr_elements;
			double A=0.0;
			for(unsigned int i=0; i<nr_elements; i++){
				double vs=Elements[i].ds.V3Norm()/dt-avg_v; // speed in Angström/ps
				A+=exp(-0.5*Elements[i].MyType->mass*amu_to_perg_ps2_per_Ang2*vs*vs/configuration->kT);
			}
			double timescale=sqrt(initial_mds2/mds2);
			double dv=(avg_speed-avg_v)/A*nr_elements;
			for(unsigned int i=0; i<nr_elements; i++){
				double v_mag=Elements[i].ds.V3Norm()/dt; // speed in Angström/ps
				double vs=v_mag-avg_v;
				double new_v=v_mag+exp(-0.5*Elements[i].MyType->mass*amu_to_perg_ps2_per_Ang2*vs*vs/configuration->kT)*dv;
				Elements[i].ds*=new_v/v_mag;
			}
			dt*=timescale;
		}
		
		// Rotate system so total dipole moment is in direction of external field
		if((configuration->rotate_Efield_steps>0) && (configuration->rotate_Efield_steps==splus1)){
			cout << "\n\tRotating external electric field to match overall dipole moment.\n";
			cout << "\t(" << (Evec/MV_per_m_to_perg_per_Debye).V3Str(',') << ")";
			Mat33 rot=RotAtoB(Evec,overallM);
			Evec=rot*Evec;
			cout << " -> (" << (Evec/MV_per_m_to_perg_per_Debye).V3Str(',') << ") V/µm\n";
			cout.flush();
		}
		
		// Premature termination code for debugging purposes. Do NOT use for normal operation.
#if DEBUG_EXIT_STEP>=0
		if(step+1>=DEBUG_EXIT_STEP){
			cout << "Terminating run at step " << step+1 << "\n";
			configuration->storetime=true;
			configuration->Runtime=((double)clock()-tstart)/CLOCKS_PER_SEC;
			store_trajectory(splus1);
			exit(12);
		}
#endif
	}// end steps loop
	
	// Get final positions and orientations of ellipsoids, then output data to simulation log
	if(traj_laststep!=configuration->steps) store_trajectory(configuration->steps); // do not output twice
	
	cout << "\nFinished calculating " << configuration->steps-configuration->resume_step << " steps with " << configuration->n_oids << " elements after " << configuration->Runtime << " seconds (" << configuration->Runtime/((configuration->steps-configuration->resume_step)*configuration->n_oids)*1E6 << " µs/step-N).\n";
	unsigned int which_acceptance=0+(nr_tries_initial>0)+(nr_tries_equi>0);
	if(!benchmark) cout << "Final acceptance rate (%%): " << (double)accept[which_acceptance]/(double)(nr_tries-(nr_tries_equi==0)*nr_tries_initial-nr_tries_equi)*100.0 << "\n\n";
	
	//clean up
	free(move_elements);
	if(configuration->LJwall_calc) delete[] save_groupVG;
	delete[] element_bitfield;
	delete[] neighbor_store;
	delete[] energy_index_store;
	delete[] cosmeans;
	delete[] Ms;
	delete[] Vs;
	delete[] totalV;
	delete[] pair_energy_store;
	delete[] center_store;
	delete[] VE_store;
	if(Rmus) delete[] Rmus; // for the "resume-at-the-end" case with user-specified trajectory file
	if(Rdist2) delete[] Rdist2; // for the "resume-at-the-end" case with user-specified trajectory file
	if(Rdist2LJ) delete[] Rdist2LJ; // for the "resume-at-the-end" case with user-specified trajectory file
	if(Rmirror) delete[] Rmirror; // for the "resume-at-the-end" case with user-specified trajectory file
	if(LJwallMirror) delete[] LJwallMirror; // for the "resume-at-the-end" case with user-specified trajectory file
}

void MC_Elements::store_trajectory(unsigned int step)
{
	string filename = configuration->fileout + "_" + int2str(configuration->runid);
	string temp_string=filename+".traj";
	if(configuration->trajectorynr>1){ // don't overwrite trajectory file user specifies
		filename+="_"+int2str(configuration->trajectorynr);
	}
	filename+=".traj";

	if((int)step>configuration->resume_step){
		if(step==configuration->steps) cout << "\n\t\t\tRecording final trajectory\n"; else cout << "\t\t\tRecording trajectory\n";
	} else cout << "\t\t\tRecording initial coordinates\n\n";

	fstream trajfile;
	if(traj_created) trajfile.open(filename.c_str(),ios::in|ios::out); else trajfile.open(filename.c_str(),ios::out|ios::trunc);
	if (trajfile.fail()){
		cout << "Could not open trajectory file.\n";
		exit(1);
	}
	string output="";
	unsigned int end=0;
	if(!traj_created){ // write header
		output+="[General]\n";
		output+="configuration = "+configuration->trajfile+"\n";
		output+="N = "+int2str(configuration->n_oids)+"\n";
		output+="restart_calculations = 0 // if 1, potential calculations (as well as statistics) are started from scratch (default: 0)\n";
		output+="restart_volume = 0 // if 1, volume upon restart will be increased to initial bounding-sphere limited volume\n";
		output+="rngseed = "+int2str(configuration->rngseed)+" // initial random number generator seed\n";
		output+="first_step = ";
		if(configuration->trajectorynr>1) output+=int2str(traj_laststep+1); else output+=int2str(configuration->resume_step+1);
		output+=" // step this trajectory is recorded from\n";
		output+="// Please choose the step from which to resume simulation (if set to -1 last step is used)\n";
		output+="// WARNING: steps after selected resume step will be overwritten\n";
		output+="resume_step = -1\n\n";
		output+="// last step included:\n";
		output+="last_step = ";
		step_pos=output.length();
#ifdef IN_WINDOWS
		step_pos+=10; // 10 lines, 10 additional carriage returns in Windows ...
#endif
		output+=int2str_fixed(step,9)+"\n\n";
		output+="// storage formats:\n";
		output+="// [Step: <step number>]\n";
		if(configuration->vdwtype==6) output+="// element nr = [ element type, group type (-1 if not in group) | position 3-vector | rotation axis, angle | gamma ]\n"; else output+="// element nr = [ element type, group type (-1 if not in group) | position 3-vector | rotation axis, angle ]\n";
		if(dt>EPS){
			output+="// [Velocity: <step number>]\n";
			output+="// element nr = [ element type, group type (-1 if not in group) | speed in m/s | velocity vector m/s ]\n";
		}
		output+="// [Statistics: <upto step number>]\n";
		output+="// step nr = { individual potential energies (comma separated, in picoerg) | total dipole moment 3-vector (in Debye) | cos, cos^2, cos^3 | system volume in Angström³ }\n\n";
	} else{ // update header
		trajfile.seekg(0,ios::end);
		end=trajfile.tellg();
		trajfile.seekp(step_pos); // go to where last step is stored
		trajfile.write(int2str_fixed(step,9).c_str(),9); // overwrite with current step
		trajfile.seekp(end); // move to end of file
	}
	traj_created=true;
	output+="[Step: "+int2str(step)+"]";
	if(step==0) output+=" // Initial coordinates";
	output+="\n";
	if(simulation_time>EPS) output+="simulation_time_in_ps = "+double2str(simulation_time)+"\n";
#ifndef USE_CMWC4096
	output+="idum = "+int2str(*configuration->idum)+"\n";
#endif
	output+="rancount = "+to_string(ran_count())+"\n";
	if(fabs(configuration->V-configuration->targetV)>EPS) output+="volume = "+double2str(configuration->V)+"\n";
	if(configuration->LJwall_calc && configuration->LJwall_fixed) output+="LJwall_xm = "+double2str(configuration->LJwall_xm)+"\n";
	if(fabs(kB*configuration->T-configuration->kT)>EPS) output+="kT = "+double2str(configuration->kT)+"\n";
	int grouptype;
	for(unsigned int i=0; i<nr_elements; i++){
		grouptype=-1;
		if(Elements[i].group) grouptype=Elements[i].group->type;
		Vec4 AxisAngle=Rot2AxisAngle(Elements[i].rot);
		output+=int2str(i+1)+"=["+int2str(Elements[i].mytype)+","+int2str(grouptype)+"|"+double2str(Elements[i].center.vec[0])+","+double2str(Elements[i].center.vec[1])+","+double2str(Elements[i].center.vec[2])+"|"+double2str(AxisAngle.vec[0])+","+double2str(AxisAngle.vec[1])+","+double2str(AxisAngle.vec[2])+","+double2str(AxisAngle.vec[3]);
		if(configuration->vdwtype==6) output+="|"+double2str(Elements[i].gamma_quarter*4);
		output+=+"]\n";
	}
	if(dt>EPS){
		output+="\n[Velocity: "+int2str(step)+"]";
		if(step==0) output+=" // Initial velocity";
		output+="\n";
		if(simulation_time>EPS) output+="simulation_time_in_ps = "+double2str(simulation_time)+"\n";
		output+="dt = "+double2str(dt)+" // ps\n";
		double avg_v=0.0;
		for(unsigned int i=0; i<nr_elements; i++) avg_v+=Elements[i].ds.V3Norm()/dt*100.0;
		avg_v/=nr_elements;
		output+="average_speed = "+double2str(avg_v)+" // m/s\n";
		if(fabs(kB*configuration->T-configuration->kT)>EPS) output+="kT = "+double2str(configuration->kT)+"\n";
		for(unsigned int i=0; i<nr_elements; i++){
			grouptype=-1;
			Vec3 vel=Elements[i].ds/dt*100.0; // 100 Angström/ps = m/s
			double v=vel.V3Norm(); // speed in m/s
			if(Elements[i].group) grouptype=Elements[i].group->type;
			output+=int2str(i+1)+"=("+int2str(Elements[i].mytype)+","+int2str(grouptype)+"|"+double2str(v)+"|"+double2str(vel.vec[0])+","+double2str(vel.vec[1])+","+double2str(vel.vec[2])+")\n";
		}
	}
	output+="\n";
	if((int)step>configuration->resume_step){
		output+="[Statistics: "+int2str(step)+"]\n";
		if(configuration->storetime){
			output+="randtime = "+double2str(configuration->Rand_Time)+"\n";
			output+="runtime = "+double2str(configuration->Runtime)+"\n";
		}
		output+="msmoved = ("+double2str(msmoved[0])+", "+double2str(msmoved[1])+")\n";
		output+="acceptance = ("+int2str(accept[0])+", "+int2str(accept[1])+", "+int2str(accept[2])+")\n";
		output+="nr_tries = (";
		if(nr_tries_initial>0) output+=int2str(nr_tries_initial); else output+=int2str(nr_tries);
		output+=", ";
		if(nr_tries_equi>0) output+=int2str(nr_tries_equi-nr_tries_initial)+", "+int2str(nr_tries-nr_tries_equi);
			else if(nr_tries_initial>0) output+=int2str(nr_tries-nr_tries_initial)+", 0"; else output+="0, 0";
		output+=")\n";
		if(configuration->NpT){
			output+="NpT_acceptance = ("+int2str(NpT_accept[0])+", "+int2str(NpT_accept[1])+", "+int2str(NpT_accept[2])+")\n";
			output+="NpT_tries = (";
			if(nr_tries_initial>0) output+=int2str(NpT_tries_initial); else output+=int2str(NpT_tries);
			output+=", ";
			if(nr_tries_equi>0) output+=int2str(NpT_tries_equi-NpT_tries_initial)+", "+int2str(NpT_tries-NpT_tries_equi);
				else if(nr_tries_initial>0) output+=int2str(NpT_tries-NpT_tries_initial)+", 0"; else output+="0, 0";
			output+=")\n";
		}
		for(unsigned int i=traj_laststep; i<step; i++){
			output+=int2str(i+1)+"={";
			for(unsigned int j=0; j<NUM_V_STORE; j++){
				output+=double2str(totalV[i][j]);
				if(j==NUM_V_STORE-1) output+="|"; else output+=",";
			}
			output+=double2str(Ms[i].vec[0])+","+double2str(Ms[i].vec[1])+","+double2str(Ms[i].vec[2])+"|";
			output+=double2str(cosmeans[i].vec[0])+","+double2str(cosmeans[i].vec[1])+","+double2str(cosmeans[i].vec[2])+"|";
			output+=double2str(Vs[i])+"}\n";
		}
		output+="\n";
		traj_laststep=step;
	}
	trajfile.write(output.c_str(),output.length());
	// Project size of trajectory after *next* step and create a new one in case it exceeds the size limit
	trajfile.seekg(0,ios::end);
	if(2*trajfile.tellg()-end>TRAJECTORY_SIZE){
		configuration->trajectorynr++;
		configuration->trajfile=filename; // link to old trajectory file (directory stays the same)
		traj_created=false;
	}
	trajfile.close();
	cout.flush();
}

