/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * Robinson group Monte Carlo simulation code Trajectory file to MOL2 converter
 * Currently requires the entire simulation code for configuration trajectory
 * parsing as well as the simulation configuration (.conf) file.
 * Command line utility with two arguments, (1) the trajectory file, and (2) the step (optional, uses last step if omitted)
 * Usage example: traj2mol2 test0.traj 1000
 * Created by Andreas Tillack, June 2012
 */

#include "MC_Config.h"
#include "MC_Elements.h"
#include "ScalarMat.h"
#include "VecMat.h"
#include "setup.h"

string clean_lastdot(string &input);
string clean_numbers(string &input, string &num);
string clean_numbers(string &input){ string num; return clean_numbers(input,num); };

string clean_lastdot(string &input)
{
	unsigned int dot=input.length();
	while((dot>0) && (input[dot]!='.')) dot--;
	if(input[dot]!='.') dot=input.length();
	string out="";
	for(unsigned int i=0; i<dot; i++) out+=input[i];
	return out;
}

string clean_numbers(string &input, string &number)
{
	unsigned int il=input.length();
	string out="";
	char number_key[]="1234567890.";
	unsigned int nr=strcspn(input.c_str(),number_key);
	for(unsigned int i=0; i<nr; i++) out+=input[i];
	
	number="";
	if(il!=nr){
		if(input[nr]=='.') nr++;
		for(;nr<il; nr++) number+=input[nr];
	}
	return out;
}

int main(int argc, char* argv[])
{
	cout << "\nRobinson Group Trajectory File Converter\n";
	cout << "Compiled " << __DATE__ << " (Build " << BUILD << ", " << VERSION << ")\n";
	
	MC_Config config;
	Config_Data* configuration = &config.parameters;
	
	string conffile="";
	// Check for command line parameters
	if(argc>1){ // yes, there are some -- only paramter accepted is configuration filename
		conffile=argv[1]; // easy enough
	} else{
		cout << "Syntax: " << argv[0] << " <trajectory file> <optional: step nr>\n\n";
		exit(1);
	}
	
	// load configuration from file
#ifndef USE_CMWC4096
	configuration->idum = new __int32_t; // so we don't get segfaults ...
#endif
	config.GetFromFile(conffile.c_str());
	phys_configuration(configuration);
	
	unsigned int step=configuration->last_step;
	if(argc>2){ // second parameter: requested step
		int id;
		bool fail=from_string(id,argv[2]);
		if(!fail) if(id<0) fail=true;
		if(!fail){
			step=(unsigned int)id;
		} else{
			cout << "Second command line parameter (step) needs to be an integer greater or equal to zero.\n";
			exit(1);
		}
	}
	
	// get step data
	Traj_EP* elements = NULL;
	if(configuration->n_oids>0) elements = new Traj_EP[configuration->n_oids];
	if(elements){
		cout << "-> Reading element positions (step " << step << ")\n";
		if(config.GetElementProperties(step,elements)){
			// Check for inconsistencies
			for(unsigned int i=0; i<configuration->n_oids; i++){
				if(i<configuration->n_group_oids){
					if((elements[i].group_type<0) || (elements[i].group_type>=(int)configuration->num_groups)){
						cout << "ERROR: Trajectory file content is inconsistent with configuration file (specified groups do not exist).\n";
						exit(1);
					}
				} else{
					if(elements[i].element_type>=configuration->num_element_types){
						cout << "ERROR: Trajectory file content is inconsistent with configuration file (specified elements do not exist).\n";
						exit(1);
					}
				}
			}
			cout << "-> Generating MOL2 file:\n";
			
			string atoms="";
			string bonds="";
			unsigned int atom_counter=0;
			unsigned int bond_counter=0;
			unsigned int current_group_start=0;
			unsigned int full_model_delta=0;
			unsigned int group_count=1;
			string name;
			for(unsigned int i=0; i<configuration->n_oids; i++){
				Element_Type* element_type=configuration->element_types[elements[i].element_type];
				name=element_type->name;
				if(element_type->archetype>=0) name=configuration->element_types[element_type->archetype]->name;
				if(elements[i].group_type<0){ // element
					bool zerocharge=false;
					if((element_type->nr_charges!=1 || element_type->hasmu) || ((element_type->saxes.vec[0]!=element_type->saxes.vec[1]) || (element_type->saxes.vec[0]!=element_type->saxes.vec[2]))){
						if(element_type->nr_charges!=1 || element_type->hasmu){
							cout << "WARNING: Sphere (element type " << elements[i].element_type << ") contains more than one charge or a dipole, mol2 allows only one charge - outputting value of zero for it.\n";
							zerocharge=true;
						} else{
							cout << "Mol2 files cannot contain user-specified ellipsoids. Aborting.\n";
							exit(2);
						}
					}
					atom_counter++;
					atoms+="\t"+int2str(atom_counter)+"\t"+clean_numbers(name)+"\t"+elements[i].position.V3Str('\t')+"\t"+name+"\t"+"1"+"\t"+"LIG1"+"\t";
					if(zerocharge) atoms+="0.0\n"; else atoms+=double2str(element_type->q[0]/e_in_esu)+"\n";
				} else{ // group
					int* distances=NULL;
					Element_Group* group=configuration->groups[elements[i].group_type];
					if(group->levelofdetail==0){ // fully-atomistic model
						if(current_group_start==i){ // at beginning of group
							int* element_storage = new int[group->nr_elements<<2];
							distances = new int[group->nr_elements+1];
							unsigned int steps=0;
							unsigned int rings=0;
							memset(distances,0xFF,group->nr_elements*sizeof(int)); // puts -1 in each field
							follow_links2(element_storage,distances,group,configuration->group_elements,0,group->nr_elements-1,steps,NULL,rings);
							delete[] element_storage;
						}
						Element* element=&configuration->group_elements[group->elements[i-current_group_start]];
						atom_counter++;
						atoms+="\t"+int2str(atom_counter)+"\t"+clean_numbers(name)+"\t"+elements[i].position.V3Str('\t')+"\t"+name+"\t"+int2str(group_count)+"\t"+group->Type->name+"\t"+double2str(element_type->q[0]/e_in_esu)+"\n";
						// add bonds from this particular molecule (if there are any)
						for(unsigned int j=0; j<element->nr_interactions; j++){
							if(element->interactions[j].partner>i-current_group_start){ // only record links we haven't recorded yet
								string bo;
								if(element->interactions[j].bond_order>0) bo=double2str(element->interactions[j].bond_order); else bo=bondorder_types[(unsigned int)(-element->interactions[j].bond_order)];
								bond_counter++;
								bonds+=int2str(bond_counter)+"\t"+int2str(atom_counter)+"\t"+int2str(current_group_start+full_model_delta+element->interactions[j].partner+1)+"\t"+bo+"\n";
							}
						}
					} else{ // LOD model -> we want underlying fully-atomistic model
						Element_Group* full_group=group->Type->LOD->groups[0];
						for(unsigned int j=0; j<full_group->nr_elements; j++){
							Element* element=&configuration->group_elements[full_group->elements[j]];
							if((element->MyType->nr_charges>1 || element->MyType->hasmu) || ((element->MyType->saxes.vec[0]!=element->MyType->saxes.vec[1]) || (element->MyType->saxes.vec[0]!=element->MyType->saxes.vec[2]))){
								cout << "Underlying fully-atomistic structure cannot contain user-specified ellipsoids. Aborting.\n";
								exit(3);
							}
							unsigned int ellipsoid_nr=group->Type->LOD->element_in_ellipsoid[group->levelofdetail-1][j]; // fully-atomistic level (0) has no ellipsoids
							Element* ellipsoid=&configuration->group_elements[group->elements[ellipsoid_nr]];
							Vec3 el_pos=elements[ellipsoid_nr+current_group_start].position;
							Mat33 el_rot=AxisAngle2Rot(elements[ellipsoid_nr+current_group_start].rotation_vector);
							Vec3 full_pos=el_pos+el_rot*(ellipsoid->rot.M3Transpose()*(element->center-ellipsoid->center));
							name=element->MyType->name;
							if(element->MyType->archetype>=0) name=configuration->element_types[element->MyType->archetype]->name;
							atom_counter++;
							double qel=0.0;
							if(element->MyType->nr_charges>0) qel=element->MyType->q[0]/e_in_esu;
							atoms+="\t"+int2str(atom_counter)+"\t"+clean_numbers(name)+"\t"+full_pos.V3Str('\t')+"\t"+name+"\t"+int2str(group_count)+"\t"+group->Type->name+"\t"+double2str(qel)+"\n";
							// add bonds from this particular molecule (if there are any)
							for(unsigned int k=0; k<element->nr_interactions; k++){
								if(element->interactions[k].partner>j){ // only record links we haven't recorded yet
									string bo;
									if(element->interactions[k].bond_order>0) bo=double2str(element->interactions[k].bond_order); else bo=bondorder_types[(unsigned int)(-element->interactions[k].bond_order)];
									bond_counter++;
									bonds+=int2str(bond_counter)+"\t"+int2str(atom_counter)+"\t"+int2str(current_group_start+full_model_delta+element->interactions[k].partner+1)+"\t"+bo+"\n";
								}
							}
						}
						full_model_delta+=full_group->nr_elements-group->nr_elements;
						i+=group->nr_elements-1; // need to jump over remaining elements in LOD group since we went over the underlying fully-atomistic model
					}
					// check if we reached the end of the current group
					if(i-current_group_start+1==group->nr_elements){
						current_group_start=i+1;
						if(i+1<configuration->n_oids){
							if(elements[i].group_type==elements[i+1].group_type) group_count++; else group_count=1;
						}
					}
				}
			}
			// Create file
			ofstream mol2out;
			string mol2name = clean_lastdot(configuration->trajectoryfile)+"-"+int2str(step)+".mol2";
			mol2out.open(mol2name.c_str());
			if(mol2out.fail()){
				cout << "Unable to open output file.\n";
				exit(1);
			}
			string header="# System state of "+configuration->fileout+" simulation at step "+int2str(step)+"\n";
			// Get ending time
			struct tm * currenttime;
			time_t rawtime;
			time (&rawtime);
			currenttime = localtime (&rawtime);
			string tstr=asctime(currenttime);
			unsigned int pbcs=(unsigned int)configuration->PBCs[0]+(unsigned int)configuration->PBCs[1]+(unsigned int)configuration->PBCs[2];
			header+="# System volume: "+double2str(configuration->V)+" Angström³";
			if(pbcs>0){
				header+=" (periodic boundary conditions in ";
				unsigned int c=0;
				for(unsigned int i=0; i<3; i++){
					if(configuration->PBCs[i]){
						c++;
						switch(i){
							case 0: header+="x";
								break;
							case 1: header+="y";
								break;
							case 2: header+="z";
								break;
						}
						if(c<pbcs) header+=",";
					}
				}
				header+=" direction)";
			}
			header+="\n";
			header+="# Created "+tstr;
			mol2out << header << "\n";
			mol2out << "@<TRIPOS>MOLECULE\n" << clean_lastdot(conffile) << "\n" << " " << atom_counter << " " << bond_counter << " 0 0 0\nSMALL\nUSER_CHARGES\n\n";
			mol2out << "@<TRIPOS>ATOM\n" << atoms;
			mol2out << "@<TRIPOS>BOND\n" << bonds << "\n";
			mol2out.close();
			cout << "<- MOL2 file succesfully created.\n";
		} else exit(1);
		delete[] elements;
	} else{
		cout << "No elements specified in configuration file => No output.\n";
		exit(42);
	}
	return 0;
}

