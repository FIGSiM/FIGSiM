/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * updated by AT on Jan 24, 2011
 * - can use multiple molecule types now
 */

#include "Latticeoids3d.h"

/*! This is a constructor that populates the BeOids object with molecules, with the BeOids.configuration set by BeOids.configuration->latticetype.
 * Currently, can only build simple cubic lattices or read in a coordinate file.
 * Inputs:
 * 	BeOids - the entire ensemble of molecules (passthrough)
 */
LatticeOids3d::LatticeOids3d(MC_Elements &BeOids)
{
	long rngsave=BeOids.configuration->rngseed;
	if(BeOids.configuration->use_trajectory){
		BeOids.configuration->rngseed=BeOids.configuration->placement_rngseed;
#ifdef USE_CMWC4096
		init_CMWC4096(BeOids.configuration->rngseed);
		zigset();
#else
		configuration->idum = &BeOids.configuration->rngseed; // idum changes each time ran2() is called
#endif
	}
	created_groups=NULL;
	nr_groups=0;
	
	if(BeOids.configuration->add_slowly>0){ // Just add the first thing in the center and don't allow it to translate
		if(BeOids.configuration->n_groups>0){
			for(unsigned int i=0; i<BeOids.configuration->num_groups; i++){
				if(BeOids.configuration->groups[i]->number>0){
					BeOids.AddGroup(i)->fixed=true;
					break;
				}
			}
		} else{
			for(unsigned int i=0; i<BeOids.configuration->num_element_types; i++)
				if(BeOids.configuration->element_types[i]->number>0){
					BeOids[BeOids.AddElement(i)]->fixed=true;
					break;
				}
		}
	} else{
		if(BeOids.configuration->n_groups>0){
			// Create groups first, then regular elements
			if(BeOids.configuration->loader==1){
				SeqLoadGroups(BeOids);
				SeqLoad(BeOids);
			} else if(BeOids.configuration->loader==2){
					RandLoadGroups(BeOids);
					RandLoad(BeOids);
			} else{
				cout << "Only loader 1 and 2 are allowed for groups at the moment.\n";
				exit(1);
			}
			// Arrange groups and elements so they don't overlap (also take care of randomly rotating elements there)
			GroupPacking(BeOids);
		} else{ // business as usual
			// Create array of molecules -- BeOids. Can be done with either a sequential or a random boxloader.
			if (BeOids.configuration->loader == 1) SeqLoad(BeOids);
			else if (BeOids.configuration->loader == 2) RandLoad(BeOids);
				else{
					cout << "Invalid boxloader specified.\n";
					exit(1);
				}
			SimpleCubic(BeOids);
			
			//Randomize orientations
			for (unsigned int i=0;i<BeOids.configuration->n_oids;i++){
#ifdef USE_CMWC4096
				double theta=2.0*pi*ranQ();
				double phi=2.0*pi*ranQ();
				double gamma=2.0*pi*ranQ();
#else
				double theta=2.0*pi*ran2(*BeOids.configuration->idum);
				double phi=2.0*pi*ran2(*BeOids.configuration->idum);
				double gamma=2.0*pi*ran2(*BeOids.configuration->idum);
#endif
				BeOids.Rotate(i,theta,phi,gamma);
			}
			
			//Randomize positions. Using this routine or latticetype == 2 is not essential, as randomizesteps handles
			//most of the work for breaking up the lattice. However, randomizesteps incorporates vdw interactions, so
			//if you want to start with a less energetically favored spatial distribution of molecules, use latticetype 2, and
			//if you want the vdw interactions to dictate the starting BeOids.configuration, use randomizesteps. Do not modify
			//this routine for other latticetypes without making sure you are calling the correct Translate function.
			if (BeOids.configuration->latticetype==2){
				for(unsigned int i=0; i<BeOids.configuration->num_element_types; i++){
					if(BeOids.configuration->element_types[i]->rot_notrans==1){
						cout << "Translating particles that shouldn't be able to move is not recommended.\n";
						exit(1);
					}
				}
				
				for (unsigned int k=0;k<10;k++)
					for (unsigned int i=0;i<BeOids.configuration->n_oids;i++)
						for (unsigned int j=0;j<3;j++) BeOids.Translate(i);
			}
		}
	}
	if(BeOids.configuration->use_trajectory) BeOids.configuration->rngseed=rngsave;
}

/*! Destructor for LatticeOids3D. While it looks tempting to replace this loop with a delete[], it will
 *  result in a segfault if you run a batch job. - LEJ
 *-- looks much better now, right :-) (AT, Feb 18, 2011)
 */
LatticeOids3d::~LatticeOids3d()
{
	if(nr_groups>0) free(created_groups); // it's a class which takes care of it's own deconstruction
}

/// Pack groups as close as possible
void LatticeOids3d::GroupPacking(MC_Elements &BeOids)
{
	unsigned int nr_elements=BeOids.configuration->n_oids-BeOids.configuration->n_group_oids;
#if DEBUG_LEVEL>0
	cout << "\nPlacing " << nr_groups << " groups and " << nr_elements << " elements in " << BeOids.configuration->fcc[0] << "x" << BeOids.configuration->fcc[1] << "x" << BeOids.configuration->fcc[2] << " = " << BeOids.configuration->nr_fcc_cells << " fcc unit cells\n";
	cout << "-> Each unit cell holds up to 4 groups and " << BeOids.configuration->elements_per_fcc << " elements.\n";
#endif
	// start with groups
	unsigned int fcc_nr, fcc_range; // running number of all fcc spaces, available numbers till next element's space
	unsigned int fcc_x, fcc_y, fcc_z;
	double fraction=0.0;
	if(nr_groups>0) fraction=(double)4*BeOids.configuration->nr_fcc_cells/nr_groups; // fraction of fcc spots each group can occupy (fraction>=1)
	Vec3 position;
	for(unsigned int i=0; i<nr_groups; i++){
		// Rotate first while each groups center is still the archetypes one
#ifdef USE_CMWC4096
		double theta=2.0*pi*ranQ();
		double phi=2.0*pi*ranQ();
		double gamma=2.0*pi*ranQ();
#else
		double theta=2.0*pi*ran2(*BeOids.configuration->idum);
		double phi=2.0*pi*ran2(*BeOids.configuration->idum);
		double gamma=2.0*pi*ran2(*BeOids.configuration->idum);
#endif
		position=Vec3(0.0);
		BeOids.Rotate(created_groups[i],position,theta,phi,gamma); // group was placed so minimum radius center is at (0,0,0)
		// Now decide on which fcc lattice space to put group
		fcc_nr=(unsigned int)round(fraction*i);
		fcc_range=(unsigned int)round(fraction*(i+1))-fcc_nr;
#if DEBUG_LEVEL>2
		cout << "fcc cell: " << fcc_nr/4 << "-" << (fcc_nr+fcc_range-1)/4;
#endif
#ifdef USE_CMWC4096
		if(fcc_range>1) fcc_nr+=(unsigned int)round(ranQ()*(fcc_range-1));
#else
		if(fcc_range>1) fcc_nr+=(unsigned int)round(ran2(*BeOids.configuration->idum)*(fcc_range-1));
#endif
#if DEBUG_LEVEL>2
		cout << " - choose: " << fcc_nr/4 << " (spot: " << fcc_nr-(fcc_nr/4)*4 << ")";
#endif
		// calculate address of fcc unit cell in x,y,z
		fcc_x=(fcc_nr/4)%BeOids.configuration->fcc[0];
		fcc_y=((fcc_nr/4)/BeOids.configuration->fcc[0])%BeOids.configuration->fcc[1];
		fcc_z=(fcc_nr/4)/(BeOids.configuration->fcc[0]*BeOids.configuration->fcc[1]);
#if DEBUG_LEVEL>2
		cout << ", address: " << fcc_x << ", " << fcc_y << ", " << fcc_z;
#endif
		// middle of volume is 0,0,0
		position.vec[0]=(double)(fcc_x-BeOids.configuration->fcc[0]/2.0)*(BeOids.configuration->scale[0]*BeOids.configuration->a_group);
		position.vec[1]=(double)(fcc_y-BeOids.configuration->fcc[1]/2.0)*(BeOids.configuration->scale[1]*BeOids.configuration->a_group);
		position.vec[2]=(double)(fcc_z-BeOids.configuration->fcc[2]/2.0)*(BeOids.configuration->scale[2]*BeOids.configuration->a_group);
		switch((unsigned int)(fcc_nr-(fcc_nr/4)*4)){
		// case 0: position currently points to it (front bottom left)
			case 1: // middle of (bottom) x-y face
				position.vec[0]+=BeOids.configuration->scale[0]*BeOids.configuration->a_group/2.0;
				position.vec[1]+=BeOids.configuration->scale[1]*BeOids.configuration->a_group/2.0;
				break;
			case 2: // middle of (front) x-z face
				position.vec[0]+=BeOids.configuration->scale[0]*BeOids.configuration->a_group/2.0;
				position.vec[2]+=BeOids.configuration->scale[2]*BeOids.configuration->a_group/2.0;
				break;
			case 3: // middle of (left) y-z face
				position.vec[1]+=BeOids.configuration->scale[1]*BeOids.configuration->a_group/2.0;
				position.vec[2]+=BeOids.configuration->scale[2]*BeOids.configuration->a_group/2.0;
				break;
		}
#if DEBUG_LEVEL>2
		cout << ", coords: " << position.vec[0] << ", " << position.vec[1] << ", " << position.vec[2] << "\n";
#endif
		BeOids.apply_PBCs(position);
		BeOids.Translate(created_groups[i],position); // this works because center is at location 0,0,0 initially
	}
#if DEBUG_LEVEL>0
	cout << "<- Finished placing and randomly rotating groups.\n";
#endif
	// place elements
	if(nr_elements>0) fraction=(double)(BeOids.configuration->elements_per_fcc*BeOids.configuration->nr_fcc_cells/nr_elements); // fraction of spots availabe per element (fraction>=1)
	unsigned int sc_spot, sc_x, sc_y, sc_z;
	for(unsigned int i=0; i<nr_elements; i++){ // groups are created first
		// Rotate
#ifdef USE_CMWC4096
		double theta=2.0*pi*ranQ();
		double phi=2.0*pi*ranQ();
		double gamma=2.0*pi*ranQ();
#else
		double theta=2.0*pi*ran2(*BeOids.configuration->idum);
		double phi=2.0*pi*ran2(*BeOids.configuration->idum);
		double gamma=2.0*pi*ran2(*BeOids.configuration->idum);
#endif
		BeOids.Rotate(i+BeOids.configuration->n_group_oids,theta,phi,gamma);
		// Now decide on which fcc lattice space to put element
		fcc_nr=(unsigned int)round(fraction*i);
		fcc_range=(unsigned int)round(fraction*(i+1))-fcc_nr;
#if DEBUG_LEVEL>2
		cout << "nr: " << i+BeOids.configuration->n_group_oids << ", fcc cell: " << fcc_nr/BeOids.configuration->elements_per_fcc << "-" << (fcc_nr+fcc_range-1)/BeOids.configuration->elements_per_fcc;
#endif
#ifdef USE_CMWC4096
		if(fcc_range>1) fcc_nr+=(unsigned int)round(ranQ()*(fcc_range-1));
#else
		if(fcc_range>1) fcc_nr+=(unsigned int)round(ran2(*BeOids.configuration->idum)*(fcc_range-1));
#endif
#if DEBUG_LEVEL>2
		cout << " - choose: " << fcc_nr/BeOids.configuration->elements_per_fcc << " (fcc spot: " << (fcc_nr-(fcc_nr/BeOids.configuration->elements_per_fcc)*BeOids.configuration->elements_per_fcc)/(BeOids.configuration->elements_per_fcc/3) << ", sc spot: " << (fcc_nr-(fcc_nr/BeOids.configuration->elements_per_fcc)*BeOids.configuration->elements_per_fcc)%(BeOids.configuration->elements_per_fcc/3) << ")";
#endif
		// calculate address of fcc unit cell in x,y,z
		fcc_x=(fcc_nr/BeOids.configuration->elements_per_fcc)%BeOids.configuration->fcc[0];
		fcc_y=((fcc_nr/BeOids.configuration->elements_per_fcc)/BeOids.configuration->fcc[0])%BeOids.configuration->fcc[1];
		fcc_z=(fcc_nr/BeOids.configuration->elements_per_fcc)/(BeOids.configuration->fcc[0]*BeOids.configuration->fcc[1]);
		
		sc_spot=(fcc_nr-(fcc_nr/BeOids.configuration->elements_per_fcc)*BeOids.configuration->elements_per_fcc)%(BeOids.configuration->elements_per_fcc/3);
		sc_x=sc_spot%BeOids.configuration->sc[0];
		sc_y=(sc_spot/BeOids.configuration->sc[0])%BeOids.configuration->sc[1];
		sc_z=sc_spot/(BeOids.configuration->sc[0]*BeOids.configuration->sc[1]);
#if DEBUG_LEVEL>2
		cout << ", address: " << fcc_x << ", " << fcc_y << ", " << fcc_z << ", sc #" << (fcc_nr-(fcc_nr/BeOids.configuration->elements_per_fcc)*BeOids.configuration->elements_per_fcc)/(BeOids.configuration->elements_per_fcc/3) << ", sc address: " << sc_x << ", " << sc_y << ", " << sc_z;
#endif
		// middle of volume is 0,0,0
		position.vec[0]=(double)(fcc_x-BeOids.configuration->fcc[0]/2.0)*(BeOids.configuration->scale[0]*BeOids.configuration->a_group);
		position.vec[1]=(double)(fcc_y-BeOids.configuration->fcc[1]/2.0)*(BeOids.configuration->scale[1]*BeOids.configuration->a_group);
		position.vec[2]=(double)(fcc_z-BeOids.configuration->fcc[2]/2.0)*(BeOids.configuration->scale[2]*BeOids.configuration->a_group);
		switch((unsigned int)(fcc_nr-(fcc_nr/BeOids.configuration->elements_per_fcc)*BeOids.configuration->elements_per_fcc)/(BeOids.configuration->elements_per_fcc/3)){ // choose fcc spot
			case 0: // bottom front middle
				position.vec[0]+=BeOids.configuration->rmax_g;
				position.vec[1]-=BeOids.configuration->c_group[1];
				position.vec[2]-=BeOids.configuration->c_group[2];
				break;
			case 1: // left front middle
				position.vec[0]-=BeOids.configuration->c_group[0];
				position.vec[1]-=BeOids.configuration->c_group[1];
				position.vec[2]+=BeOids.configuration->rmax_g;
				break;
			case 2: // left bottom middle
				position.vec[0]-=BeOids.configuration->c_group[0];
				position.vec[1]+=BeOids.configuration->rmax_g;
				position.vec[2]-=BeOids.configuration->c_group[2];
				break;
		}
		position.vec[0]+=(double)(sc_x+1)*BeOids.configuration->rmax_e; // need to shift sc lattice by (r,r,r) in order for spheres to fully be inside box
		position.vec[1]+=(double)(sc_y+1)*BeOids.configuration->rmax_e;
		position.vec[2]+=(double)(sc_z+1)*BeOids.configuration->rmax_e;
		
		// now assign spot in simple cubic lattice inside group fcc lattice
#if DEBUG_LEVEL>2
		cout << ", coords: " << position.vec[0] << ", " << position.vec[1] << ", " << position.vec[2] << "\n";
#endif
		BeOids.apply_PBCs(position);
		BeOids.Translate(i+BeOids.configuration->n_group_oids,position); // groups elements are created first
	}
#if DEBUG_LEVEL>0
	cout << "<- Finished placing and randomly rotating elements.\n";
	cout << "<- Done.\n";
#endif
}

/// Simple Cubic Lattice generator. Assigns coordinates to a pre-constructed ensemble of molecules to form a simple cubic lattice.
void LatticeOids3d::SimpleCubic(MC_Elements &BeOids)
{
	unsigned int nr_elements=BeOids.configuration->n_oids-BeOids.configuration->n_group_oids;
#if DEBUG_LEVEL>0
	cout << "\nPlacing " << nr_elements << " elements in " << BeOids.configuration->sc[0] << "x" << BeOids.configuration->sc[1] << "x" << BeOids.configuration->sc[2] << " = " << BeOids.configuration->nr_sc_cells << " simple cubic unit cells\n";
#endif
	// place elements
	Vec3 position;
	unsigned int sc_nr, sc_range;
	double fraction=(double)(BeOids.configuration->nr_sc_cells/nr_elements); // fraction of spots availabe per element (fraction>=1)
	unsigned int sc_x, sc_y, sc_z;
	for(unsigned int i=0; i<nr_elements; i++){
		// Now decide on which sc lattice space to put element
		sc_nr=(unsigned int)round(fraction*i);
		sc_range=(unsigned int)round(fraction*(i+1))-sc_nr;
#if DEBUG_LEVEL>2
		cout << "nr: " << i << ", sc cell: " << sc_nr << "-" << (sc_nr+sc_range-1);
#endif
#ifdef USE_CMWC4096
		if(sc_range>1) sc_nr+=(unsigned int)round(ranQ()*(sc_range-1));
#else
		if(sc_range>1) sc_nr+=(unsigned int)round(ran2(*BeOids.configuration->idum)*(sc_range-1));
#endif
#if DEBUG_LEVEL>2
		cout << " - choose: " << sc_nr;
#endif
		// to stay consistent with previous builds filling direction is z->y->x
		sc_x=sc_nr/(BeOids.configuration->sc[1]*BeOids.configuration->sc[2]);
		sc_y=(sc_nr/BeOids.configuration->sc[2])%BeOids.configuration->sc[1];
		sc_z=sc_nr%BeOids.configuration->sc[2];
#if DEBUG_LEVEL>2
		cout << ", address: " << sc_x << ", " << sc_y << ", " << sc_z;
#endif
		position.vec[0]=((double)sc_x-BeOids.configuration->sc[0]/(2.0)+0.5)*(BeOids.configuration->scale[0]*BeOids.configuration->a_element);
		position.vec[1]=((double)sc_y-BeOids.configuration->sc[1]/(2.0)+0.5)*(BeOids.configuration->scale[1]*BeOids.configuration->a_element);
		position.vec[2]=((double)sc_z-BeOids.configuration->sc[2]/(2.0)+0.5)*(BeOids.configuration->scale[2]*BeOids.configuration->a_element);
		// now assign spot in simple cubic lattice
#if DEBUG_LEVEL>2
		cout << ", coords: " << position.vec[0] << ", " << position.vec[1] << ", " << position.vec[2] << "\n";
#endif
		BeOids.Translate(i,position);
	}
#if DEBUG_LEVEL>0
	cout << "<- Done.\n";
#endif
}

/*!
Sequential loader for groups. Should only be used for homogeneous systems or systems where you *really* want
to start with two distinct layers, as it loads all A oids, then all B oids. Called if BeOids.configuration->loader == 1
*/
void LatticeOids3d::SeqLoadGroups(MC_Elements &BeOids)
{
	for(unsigned int i=0; i<BeOids.configuration->num_groups; i++){
		for(unsigned int j=0; j<BeOids.configuration->groups[i]->number; j++){
			nr_groups++;
			created_groups=(Element_Group**)realloc(created_groups,sizeof(Element_Group*)*nr_groups);
			if(created_groups==NULL){
				cout << "Tiny problem here, not enough memory.\n";
				exit(3);
			}
			created_groups[nr_groups-1]=BeOids.AddGroup(i);
		}
	}
}

/*!
Random loader for groups
Used with loader 2
*/
void LatticeOids3d::RandLoadGroups(MC_Elements &BeOids)
{
	unsigned int created[BeOids.configuration->num_groups];
	double fraction[BeOids.configuration->num_groups];
	created[0] = 0;
	fraction[0] = 0.0;
	for(unsigned int i=1; i<BeOids.configuration->num_groups; i++){
		created[i]=0; // we haven't created any groups yet
		fraction[i]=fraction[i-1]+(double)BeOids.configuration->groups[i-1]->number/(double)BeOids.configuration->n_groups;
	}
	double r;
	unsigned int i,j;
	i=0;
	while(i<BeOids.configuration->n_groups){
#ifdef USE_CMWC4096
		r=ranQ();
#else
		r=ran2(*BeOids.configuration->idum);
#endif
		// which type of group did we randomly (and proportionally distributed) choose?
		j=BeOids.configuration->num_groups;
		while(j>0){
			j--;
			if (r>fraction[j]) break;
		}
		if(created[j]<BeOids.configuration->groups[j]->number){
			nr_groups++;
			// Need to be able to easily access all groups so keep track of them ...
			created_groups=(Element_Group**)realloc(created_groups,sizeof(Element_Group*)*nr_groups);
			if(created_groups==NULL){
				cout << "Tiny problem here, not enough memory.\n";
				exit(3);
			}
			created_groups[nr_groups-1]=BeOids.AddGroup(j);
			created[j]++;
			i++;
		}
	}
}

/*!
Sequential loader. Should only be used for homogeneous systems or systems where you *really* want
to start with two distinct layers, as it loads all A oids, then all B oids. Called if BeOids.configuration->loader == 1
*/
void LatticeOids3d::SeqLoad(MC_Elements &BeOids)
{
	for(unsigned int i=0; i<BeOids.configuration->num_element_types; i++)
		for(unsigned int j=0; j<BeOids.configuration->element_types[i]->number; j++) BeOids.AddElement(i);
}

/*!
Random loader for mixtures
Used with loader 2
*/
void LatticeOids3d::RandLoad(MC_Elements &BeOids)
{
	unsigned int nr_elements=BeOids.configuration->n_oids-BeOids.configuration->n_group_oids;
	unsigned int created[BeOids.configuration->num_element_types];
	double fraction[BeOids.configuration->num_element_types];
	created[0] = 0;
	fraction[0] = 0.0;
	for(unsigned int i=1; i<BeOids.configuration->num_element_types; i++){
	    created[i]=0; // we haven't created any molecules yet
	    fraction[i]=fraction[i-1]+(double)BeOids.configuration->element_types[i-1]->number/(double)BeOids.configuration->n_oids;
	}
	double r;
	unsigned int i,j;
	i=0;
	while(i<nr_elements){
#ifdef USE_CMWC4096
		r=ranQ();
#else
		r=ran2(*BeOids.configuration->idum);
#endif
		// which type of molecule did we randomly (and proportionally distributed) choose?
		j=BeOids.configuration->num_element_types;
		while(j>0){
			j--;
			if (r>fraction[j]) break;
		}
		if(created[j]<BeOids.configuration->element_types[j]->number){
			BeOids.AddElement(j);
			created[j]++;
			i++;
		}
	}
}

