/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * Robinson group Monte Carlo simulation code Trajectory file to X3D converter
 * Currently requires the entire simulation code for configuration trajectory
 * parsing as well as the simulation configuration (.conf) file.
 * Command line utility with two arguments, (1) the trajectory file, and (2) the step (optional, uses last step if omitted)
 * Usage example: traj2x3d test0.traj 1000
 * Created by Andreas Tillack and Lewis Johnson, May 2011
 */

#include "traj2dat.h"
#include "x3d_functions.h"

#define calc_distance 0
#define group_beginning 0
#define group_end 31

string name_vizfile(string &trajfile);
void output_element(ofstream &vizout, Traj_EP &an_element, Element_Type* type, const char* DEFname, const char* extra, bool animate, Traj_EP** movements, unsigned int current_nr, unsigned int group_start, unsigned int steps, Config_Data* configuration, double* V);

///Auto-generate x3d filename from trajectory filename
string name_vizfile(string &trajfile)
{
	//Get length of streng
	unsigned int trajf_len = trajfile.length();
	
	//Copy filename without extension
	trajfile.begin();
	string vizfile = "";
	for (unsigned int i = 0; i < (trajf_len - 5); i++) {
		vizfile += trajfile[i];
	}
	
	//Assemble new filename
	vizfile += ".x3d";
	return vizfile;
}

///Output element to x3d file
void output_element(ofstream &vizout, Traj_EP &an_element, Element_Type* type, const char* DEFname, const char* extra, bool animate, Traj_EP** movements, unsigned int current_nr, unsigned int group_start, unsigned int steps, Config_Data* configuration, double* V, double* Xm)
{
	//Convert rotation from unit vector to angle and axis
	string groupname="";
	if(an_element.group_type>=0){
		if(configuration->groups[an_element.group_type]->Type->visual_PBCs) apply_PBCs(an_element.position,configuration);
		groupname=" ("+configuration->groups[an_element.group_type]->Type->name+")";
	}
	vizout << extra << "\t\t<Transform>\n";
	vizout << extra <<  "\t\t\t<ProtoInstance name='" << type->name+groupname << "' DEF='" << DEFname << "'>\n";
	vizout << extra <<  "\t\t\t\t<fieldValue name='position' value='" << an_element.position.V3Str(' ') << "'/>\n";
	vizout << extra <<  "\t\t\t\t<fieldValue name='rotation' value='" << an_element.rotation_vector.V4Str(' ') << "'/>\n";
	double transparency=type->transparency;
	if(an_element.group_type>=0){
		if(configuration->groups[an_element.group_type]->Type->transparency>=0.0) transparency=configuration->groups[an_element.group_type]->Type->transparency;
		if(configuration->groups[an_element.group_type]->Type->LOD && (configuration->groups[an_element.group_type]->levelofdetail>0)){
			double inside_transparency=configuration->groups[an_element.group_type]->Type->LOD->inside_transparency[configuration->groups[an_element.group_type]->levelofdetail-1];
			if(inside_transparency>EPS) vizout << extra <<  "\t\t\t\t<fieldValue name='inside_transparency' value='" << inside_transparency << "'/>\n";
		}
		if(configuration->groups[an_element.group_type]->Type->LOD && (configuration->groups[an_element.group_type]->levelofdetail==0)){
			if((configuration->groups[an_element.group_type]->Type->LOD->nr_components>0) && (configuration->groups[an_element.group_type]->Type->LOD->element_in_component[current_nr-group_start]>=0)){
				if(configuration->groups[an_element.group_type]->Type->LOD->component_transparency[configuration->groups[an_element.group_type]->Type->LOD->element_in_component[current_nr-group_start]]>=0.0){
					transparency=configuration->groups[an_element.group_type]->Type->LOD->component_transparency[configuration->groups[an_element.group_type]->Type->LOD->element_in_component[current_nr-group_start]];
				}
			}
		}
	}
	vizout << extra <<  "\t\t\t\t<fieldValue name='transparency' value='" << transparency << "'/>\n";
	if(animate && movements){
		string key="";
		string positions="";
		string rotations="";
		double ismo=1.0/(steps-1);
		Vec3 pos, oldpos, groupcenter, oldgroupcenter;
		Vec4 rot;
		bool visual_PBCs=true;
		if(an_element.group_type>=0){
			if(!configuration->groups[an_element.group_type]->Type->visual_PBCs) visual_PBCs=false;
		}
		for(unsigned int i=0; i<steps; i++){
			if(configuration->LJwall_calc && configuration->LJwall_fixed){
				double scale=sqrt((V[i]*configuration->boxlength[0])/(configuration->V*2.0*Xm[i]));
				configuration->LJwall_xm=Xm[i];
				configuration->boxlength[0]=2.0*Xm[i];
				configuration->boxlength[1]*=scale; configuration->boxlength[2]*=scale;
				configuration->nndist*=scale;
				configuration->V=configuration->boxlength[0]*configuration->boxlength[1]*configuration->boxlength[2];
			} else update_volume(configuration,V[i]);
			pos=movements[current_nr][i].position;
			if(an_element.group_type>=0){
				if(configuration->groups[an_element.group_type]->Type->visual_PBCs) apply_PBCs(pos,configuration);
			}
			rot=movements[current_nr][i].rotation_vector;
#if calc_distance>0
			if((current_nr==group_end) && (i>=0.8*steps)){
				cout << i*configuration->grfreq << "\t" << (pos-movements[group_beginning][i].position).V3Norm() << "\n";
			}
#endif
			if(!visual_PBCs){ // only available for groups
				groupcenter=Vec3(0.0);
				for(unsigned int j=0; j<configuration->groups[an_element.group_type]->nr_elements; j++){
					groupcenter+=movements[j+group_start][i].position; // only works because group elements are created sequentially
				}
				groupcenter/=configuration->groups[an_element.group_type]->nr_elements;
			}
			if(i>0){
				bool jump;
				Vec3 dpos=pos-oldpos;
				Vec3 dposPBCs=Vec3(0.0);
				if(configuration->smooth_animation){
					dposPBCs=dpos;
					apply_PBCs(dposPBCs,configuration);
					if(!visual_PBCs){ // only available for groups
						Vec3 dgpos=groupcenter-oldgroupcenter;
						Vec3 dgposPBCs=dgpos;
						apply_PBCs(dgposPBCs,configuration);
						oldpos=pos-dgpos+dgposPBCs;
						dposPBCs=Vec3(0.0);
						jump=(dgposPBCs!=dgpos);
					} else jump=(dposPBCs!=dpos);
				} else jump=true;
				if(jump){ // Element center or group center "jumped" across a boundary condition (which we don't want to smoothly animate)
					// easy fix: move smoothly over boundary and then jump over in this step
					key+=" "+double2str(i*ismo);
					positions+=" "+(oldpos+dposPBCs).V3Str(' ');
					if(configuration->smooth_animation) rotations+=" "+rot.V4Str(' ');
				}
				key+=" ";
				positions+=" ";
				rotations+=" ";
			}
			key+=double2str(i*ismo);
			positions+=pos.V3Str(' ');
			rotations+=rot.V4Str(' ');
			if(!configuration->smooth_animation) rotations+=" "+rot.V4Str(' ');
			oldpos=pos;
			if(!visual_PBCs) oldgroupcenter=groupcenter;
		}
		vizout << extra <<  "\t\t\t\t<fieldValue name='key' value='" << key << "'/>\n";
		vizout << extra <<  "\t\t\t\t<fieldValue name='position_values' value='" << positions << "'/>\n";
		vizout << extra <<  "\t\t\t\t<fieldValue name='rotation_values' value='" << rotations << "'/>\n";
	}
	vizout << extra <<  "\t\t\t</ProtoInstance>\n";
	vizout << extra <<  "\t\t</Transform>\n";
}

void output_animated_dipole(cfparams_t &params, ofstream &vizout, unsigned int group_type, Traj_EP** movements, unsigned int current_nr, unsigned int group_start, unsigned int steps, Config_Data* configuration, double* V, double* Xm)
{
	string key="";
	string positions="";
	string rotations="";
	string lengths="";
	double ismo=1.0/(steps-1);
	Vec3 oldpos, groupcenter, oldgroupcenter;
	Vec3* pos=new Vec3; 
	Vec3* dipole=new Vec3; 
	Vec3* cp_maxdipole=new Vec3; 
	Mat33* cp_rot=new Mat33;
	Vec4 rot;
	bool visual_PBCs=true;
	if(group_type>=0){
		if(!configuration->groups[group_type]->Type->visual_PBCs) visual_PBCs=false;
	}
	Vec3 y_dir(0.0,1.0,0.0);
	double length_ref=1.0;
	for(unsigned int i=0; i<steps; i++){
		unsigned int nr=1;
		GetPosDipole(params,0,i,pos,dipole,cp_maxdipole,cp_rot,nr,current_nr,true);
		if(configuration->LJwall_calc && configuration->LJwall_fixed){
			double scale=sqrt((V[i]*configuration->boxlength[0])/(configuration->V*2.0*Xm[i]));
			configuration->LJwall_xm=Xm[i];
			configuration->boxlength[0]=2.0*Xm[i];
			configuration->boxlength[1]*=scale; configuration->boxlength[2]*=scale;
			configuration->nndist*=scale;
			configuration->V=configuration->boxlength[0]*configuration->boxlength[1]*configuration->boxlength[2];
		} else update_volume(configuration,V[i]);
		if(configuration->groups[group_type]->Type->visual_PBCs) apply_PBCs(*pos,configuration);
		rot=Rot2AxisAngle(RotAtoB(y_dir,*dipole));
		if(i==0){
			length_ref=dipole->V3Norm()/configuration->dipoleQ;
			vizout << "\t\t\t\t<fieldValue name='length' value='" << length_ref << "'/>\n";
			vizout << "\t\t\t\t<fieldValue name='tip_vector' value='0 " << 0.5*length_ref << " 0'/>\n";
			vizout << "\t\t\t\t<fieldValue name='position' value='" << pos->V3Str(' ') << "'/>\n";
			vizout << "\t\t\t\t<fieldValue name='rotation_from_y' value='" << rot.V4Str(' ') << "'/>\n";
		}
		if(!visual_PBCs){ // only available for groups
			groupcenter=Vec3(0.0);
			for(unsigned int j=0; j<configuration->groups[group_type]->nr_elements; j++){
				groupcenter+=movements[j+group_start][i].position; // only works because group elements are created sequentially
			}
			groupcenter/=configuration->groups[group_type]->nr_elements;
		}
		if(i>0){
			bool jump;
			Vec3 dpos=*pos-oldpos;
			Vec3 dposPBCs=Vec3(0.0);
			if(configuration->smooth_animation){
				dposPBCs=dpos;
				apply_PBCs(dposPBCs,configuration);
				if(!visual_PBCs){ // only available for groups
					Vec3 dgpos=groupcenter-oldgroupcenter;
					Vec3 dgposPBCs=dgpos;
					apply_PBCs(dgposPBCs,configuration);
					oldpos=*pos-dgpos+dgposPBCs;
					dposPBCs=Vec3(0.0);
					jump=(dgposPBCs!=dgpos);
				} else jump=(dposPBCs!=dpos);
			} else jump=true;
			if(jump){ // Element center or group center "jumped" across a boundary condition (which we don't want to smoothly animate)
				// easy fix: move smoothly over boundary and then jump over in this step
				key+=" "+double2str(i*ismo);
				positions+=" "+(oldpos+dposPBCs).V3Str(' ');
				lengths+=" 1 "+double2str((dipole->V3Norm()/configuration->dipoleQ)/length_ref)+" 1";
				if(configuration->smooth_animation) rotations+=" "+rot.V4Str(' ');
			}
			key+=" ";
			positions+=" ";
			rotations+=" ";
			lengths+=" ";
		}
		key+=double2str(i*ismo);
		positions+=pos->V3Str(' ');
		rotations+=rot.V4Str(' ');
		lengths+="1 "+double2str((dipole->V3Norm()/configuration->dipoleQ)/length_ref)+" 1";
		if(!configuration->smooth_animation) rotations+=" "+rot.V4Str(' ');
		oldpos=*pos;
		if(!visual_PBCs) oldgroupcenter=groupcenter;
	}
	vizout << "\t\t\t\t<fieldValue name='key' value='" << key << "'/>\n";
	vizout << "\t\t\t\t<fieldValue name='position_values' value='" << positions << "'/>\n";
	vizout << "\t\t\t\t<fieldValue name='rotation_values' value='" << rotations << "'/>\n";
	vizout << "\t\t\t\t<fieldValue name='length_values' value='" << lengths << "'/>\n";
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
	
	bool animate=false;
	unsigned int step=configuration->last_step;
	if(argc>2){ // second parameter: requested step
		if(compare_strings(argv[2],"all")){
			animate=true;
			step=0;
		} else{
			int id;
			bool fail=from_string(id,argv[2]);
			if(!fail) if(id<0) fail=true;
			if(!fail){
				step=(unsigned int)id;
			} else{
				cout << "Second command line parameter (step to be plotted) needs to be an integer greater or equal to zero.\n";
				exit(1);
			}
		}
		if(argc>3){
			if(compare_strings(argv[3],"all")){
				animate=true;
			} else{
				int id;
				bool fail=from_string(id,argv[3]);
				if(!fail) if(id<0) fail=true;
				if(!fail){
					if(animate){ // third parameter describes start step
						step=(unsigned int)id;
					} else{ // third parameter describes end step
						configuration->last_step=(unsigned int)id;
						animate=true;
					}
				} else{
					cout << "Third command line parameter (last step to be animated) needs to be an integer greater or equal to initial step (" << step << ").\n";
					exit(1);
				}
			}
		}
	}
	
	// get step data
	Traj_EP* elements = NULL;
	unsigned int* individual_steps=NULL;
	double* V=NULL;
	double* Xm=NULL;
	double* time=NULL;
	Traj_EP** movements=NULL;
	if(configuration->n_oids>0) elements = new Traj_EP[configuration->n_oids];
	if(elements){
		cout << "-> Reading";
		if(animate) cout << " initial";
		cout << " element positions (step " << step << ")\n";
		if(config.GetElementProperties(step,elements)){
			unsigned int steps=(configuration->last_step-step)/configuration->grfreq+(configuration->last_step%configuration->grfreq>0)+1;
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
			cout << "-> Generating X3D file:\n";
			// Create file
			ofstream vizout;
			string vizname = name_vizfile(configuration->trajectoryfile);
			vizout.open(vizname.c_str());
			if (vizout.fail()) {
				cout << "Unable to open output file.\n";
				exit(1);
			}
			double averageV=configuration->V;
			if(animate){
				cout << "\t-> Reading in data for animation (may take a while) ...\n";
				individual_steps=new unsigned int[steps];
				V=new double[steps];
				if(configuration->LJwall_calc) Xm=new double[steps];
				time=new double[steps];
				movements=new Traj_EP*[configuration->n_oids];
				if(!config.GetAllElementProperties(movements,individual_steps,V,Xm,time,steps,step)) exit(1);
				averageV=0.0;
				double wsum=0.0;
				for(unsigned int i=0; i<steps; i++){
					double w=1.0;
					if(i+step<configuration->randsteps) w=0.1;
					averageV+=w*V[i];
					wsum+=w;
				}
				averageV/=wsum;
			}
			double position_scale=cbrt(averageV/configuration->V);
			cfparams_t* cf_params=new cfparams_t;
			cf_params->startstep = step;
			cf_params->endstep = configuration->last_step;
			cf_params->start_frame_in_data = step/configuration->grfreq;
			cf_params->configuration=configuration;
			cf_params->individual_output = true;
			if(animate){
				cf_params->trajectory=movements;
				cf_params->V=V;
				cf_params->Xm=Xm;
				cf_params->time=time;
			} else{
				cf_params->trajectory=new Traj_EP*[configuration->n_oids];
				for(unsigned int i=0; i<configuration->n_oids; i++){
					cf_params->trajectory[i]=new Traj_EP;
					cf_params->trajectory[i][0]=elements[i];
				}
			}
			cf_params->startframe=cf_params->startstep/configuration->grfreq;
			cf_params->endframe = configuration->last_step/configuration->grfreq+(configuration->last_step%configuration->grfreq>0)+1;
			cf_params->n_frames = cf_params->endframe-cf_params->startframe;
			cf_params->items_involved=1;
			cf_params->items=new Corr_Item;
			
			cout << "\t-> Defining element prototypes ...\n";
			// Configure viewpoint -- want to have left-handed coordinate system (x to the right, y to the back, z up)
			viewpoint_t viewpoint;
			viewpoint.Center = Vec3(0,0,0);
			double lookdownangle=-pi/7.854;
			viewpoint.Orientation = Vec4(1,0,0,pi*0.5+lookdownangle);
			viewpoint.Position = Vec3(0.0,position_scale*0.75*configuration->boxlength[1]/lookdownangle,position_scale*0.75*configuration->boxlength[2]);
			viewpoint.Description = "Default viewpoint";
			
			// Output file header info
			output_header(vizout, vizname, viewpoint,animate,steps);
			
			// Load and define element types
			vizout << "\t\t<!-- Define element types -->\n";
			for(unsigned int i = 0; i < configuration->num_element_types; i++){
				define_element(vizout,configuration->element_types[i],NULL,animate,steps,configuration); // define elements independent if they are used or not
				for(unsigned int j=0; j<configuration->n_group_oids; j++){
					if(elements[j].element_type==i){
						define_element(vizout,configuration->element_types[i],configuration->groups[elements[j].group_type],animate,steps,configuration,true);
						break;
					}
				}
			}
			vizout << "\t\t<!-- done -->\n";
			
			cout << "\t-> Defining independent elements ...\n";
			vizout << "\t\t<!-- Independent elements -->\n";
			unsigned int* element_numbers = new unsigned int[configuration->num_element_types];
			for(unsigned int i=0; i<configuration->num_element_types; i++) element_numbers[i]=0;
			for(unsigned int i = configuration->n_group_oids; i < configuration->n_oids; i++){
				if(animate) output_element(vizout, elements[i], configuration->element_types[elements[i].element_type],(configuration->element_types[elements[i].element_type]->name+"_"+int2str(element_numbers[elements[i].element_type]+1)).c_str(),"",animate,movements,i,0,steps,configuration,V,Xm);
					else output_element(vizout, elements[i], configuration->element_types[elements[i].element_type],(configuration->element_types[elements[i].element_type]->name+"_"+int2str(element_numbers[elements[i].element_type]+1)).c_str(),"",animate,NULL,i,0,steps,configuration,V,Xm);
				element_numbers[elements[i].element_type]++;
			}
			if(element_numbers) delete[] element_numbers;
			vizout << "\t\t<!-- ... -->\n";
			cout << "\t-> Defining groups ...\n";
			vizout << "\t\t<!-- Groups -->\n";
			unsigned int* group_numbers = new unsigned int[configuration->num_groups];
			if(!group_numbers){
				cout << "Not enough memory to store group numbers.\n";
				exit(2);
			}
			for(unsigned int i=0; i<configuration->num_groups; i++) group_numbers[i]=0;
			unsigned int group_start=0;
			bool has_components=false;
			bool component_dipole_show=false;
			for(unsigned int i=0; i<configuration->n_group_oids; i++){
				Element_Group* group=configuration->groups[elements[i].group_type];
				if(group->Type->LOD){
					if(group->Type->LOD->nr_components>0){
						has_components=true;
						for(unsigned int j=0; j<group->Type->LOD->nr_components; j++){
							if(group->Type->LOD->component_show_dipole[j]){
								component_dipole_show=true;
								break;
							}
						}
					}
				} else has_components=false;
				if(group_numbers[elements[i].group_type]%group->nr_elements==0){ // this only works because group elements are created sequentially
					group_start=i;
#if calc_distance>0
					if(animate) cout << group->Type->name << " initial end-to-end distance (between " << configuration->group_elements[group->elements[group_beginning]].MyType->name << " and " << configuration->group_elements[group->elements[group_end]].MyType->name << "): " << (movements[group_end][0].position-movements[group_beginning][0].position).V3Norm() << " Angstrom\n";
#endif
					if(has_components){ // output component dipoles first if they are to be drawn
						Vec3* cp_pos=new Vec3;
						Vec3* cp_dipole=new Vec3;
						Vec3* cp_maxdipole=new Vec3;
						Mat33* cp_rot=new Mat33;
						Vec3 y_dir(0.0,1.0,0.0);
						for(unsigned int j=0; j<group->Type->LOD->nr_components; j++){
							if(group->Type->LOD->component_show_dipole[j]){
								cf_params->items->is_element=false;
								cf_params->items->is_group=true;
								cf_params->items->group_idx=elements[i].group_type;
								cf_params->items->is_LOD=false;
								cf_params->items->LOD_nr=-1;
								cf_params->items->is_component=true;
								cf_params->items->element_idx=j;
								cf_params->items->type_idx=j;
								unsigned int nr=1;
								if(animate){
									vizout << "\t\t<Transform>\n";
									vizout << "\t\t\t<ProtoInstance name='Vector3D' DEF='" << configuration->groups[elements[i].group_type]->Type->name << ".Component:" << j << " dipole'>\n";
									vizout << "\t\t\t\t<fieldValue name='color' value='" << group->Type->LOD->component_color[j].V3Str(' ') << "'/>\n";
									output_animated_dipole(*cf_params,vizout,elements[i].group_type,movements,i,group_start,steps,configuration,V,Xm);
									vizout << "\t\t\t</ProtoInstance>\n";
									vizout << "\t\t</Transform>\n";
								} else{
									GetPosDipole(*cf_params,0,0,cp_pos,cp_dipole,cp_maxdipole,cp_rot,nr,i,true);
									vizout << "\t\t<Transform>\n";
									vizout << "\t\t\t<ProtoInstance name='Vector3D' DEF='" << configuration->groups[elements[i].group_type]->Type->name << ".Component:" << j << " dipole'>\n";
									vizout << "\t\t\t\t<fieldValue name='length' value='" << cp_dipole->V3Norm()/configuration->dipoleQ << "'/>\n";
									vizout << "\t\t\t\t<fieldValue name='tip_vector' value='0 " << 0.5*cp_dipole->V3Norm()/configuration->dipoleQ << " 0'/>\n";
									vizout << "\t\t\t\t<fieldValue name='color' value='" << group->Type->LOD->component_color[j].V3Str(' ') << "'/>\n";
									vizout << "\t\t\t\t<fieldValue name='position' value='" << cp_pos->V3Str(' ') << "'/>\n";
									vizout << "\t\t\t\t<fieldValue name='rotation_from_y' value='" << Rot2AxisAngle(RotAtoB(y_dir,*cp_dipole)).V4Str(' ') << "'/>\n";
									vizout << "\t\t\t</ProtoInstance>\n";
									vizout << "\t\t</Transform>\n";
								}
							}
						}
						delete cp_pos;
						delete cp_dipole;
						delete cp_maxdipole;
					}
					if(group->Type->show_just_dipoles && !component_dipole_show){ // fallback to whole group's dipole in case user didn't specify which dipole to show
						Vec3* cp_pos=new Vec3;
						Vec3* cp_dipole=new Vec3;
						Vec3* cp_maxdipole=new Vec3;
						Mat33* cp_rot=new Mat33;
						Vec3 y_dir(0.0,1.0,0.0);
						cf_params->items->is_element=false;
						cf_params->items->is_group=true;
						cf_params->items->group_idx=elements[i].group_type;
						cf_params->items->is_LOD=false;
						cf_params->items->LOD_nr=-1;
						cf_params->items->is_component=false;
						cf_params->items->element_idx=-1;
						cf_params->items->type_idx=-1;
						unsigned int nr=1;
						if(animate){
							vizout << "\t\t<Transform>\n";
							vizout << "\t\t\t<ProtoInstance name='Vector3D' DEF='" << configuration->groups[elements[i].group_type]->Type->name << " dipole'>\n";
							vizout << "\t\t\t\t<fieldValue name='color' value='" << group->Type->group_dipole_color.V3Str(' ') << "'/>\n";
							output_animated_dipole(*cf_params,vizout,elements[i].group_type,movements,i,group_start,steps,configuration,V,Xm);
							vizout << "\t\t\t</ProtoInstance>\n";
							vizout << "\t\t</Transform>\n";
						} else{
							GetPosDipole(*cf_params,0,0,cp_pos,cp_dipole,cp_maxdipole,cp_rot,nr,i,true);
							vizout << "\t\t<Transform>\n";
							vizout << "\t\t\t<ProtoInstance name='Vector3D' DEF='" << configuration->groups[elements[i].group_type]->Type->name << " dipole'>\n";
							vizout << "\t\t\t\t<fieldValue name='color' value='" << group->Type->group_dipole_color.V3Str(' ') << "'/>\n";
							vizout << "\t\t\t\t<fieldValue name='length' value='" << cp_dipole->V3Norm()/configuration->dipoleQ << "'/>\n";
							vizout << "\t\t\t\t<fieldValue name='tip_vector' value='0 " << 0.5*cp_dipole->V3Norm()/configuration->dipoleQ << " 0'/>\n";
							vizout << "\t\t\t\t<fieldValue name='position' value='" << cp_pos->V3Str(' ') << "'/>\n";
							vizout << "\t\t\t\t<fieldValue name='rotation_from_y' value='" << Rot2AxisAngle(RotAtoB(y_dir,*cp_dipole)).V4Str(' ') << "'/>\n";
							vizout << "\t\t\t</ProtoInstance>\n";
							vizout << "\t\t</Transform>\n";
						}
						delete cp_pos;
						delete cp_dipole;
						delete cp_maxdipole;
					}
					if(!group->Type->show_just_dipoles) vizout << "\t\t<Group DEF='" << group->Type->name << "_" << group_numbers[elements[i].group_type]/group->nr_elements+1 << "'>\n";
				}
				if(!group->Type->show_just_dipoles){
					if(animate) output_element(vizout, elements[i], configuration->element_types[elements[i].element_type],(configuration->groups[elements[i].group_type]->Type->name+"_"+int2str(group_numbers[elements[i].group_type]/configuration->groups[elements[i].group_type]->nr_elements+1)+":"+int2str(group_numbers[elements[i].group_type]%group->nr_elements+1)).c_str(),"\t",animate,movements,i,group_start,steps,configuration,V,Xm);
						else output_element(vizout, elements[i], configuration->element_types[elements[i].element_type],(configuration->groups[elements[i].group_type]->Type->name+"_"+int2str(group_numbers[elements[i].group_type]/configuration->groups[elements[i].group_type]->nr_elements+1)+":"+int2str(group_numbers[elements[i].group_type]%group->nr_elements+1)).c_str(),"\t",animate,NULL,i,group_start,steps,configuration,V,Xm);
				}
				group_numbers[elements[i].group_type]++;
				if(group_numbers[elements[i].group_type]%configuration->groups[elements[i].group_type]->nr_elements==0){
					if(!group->Type->show_just_dipoles) vizout << "\t\t</Group>\n";
				}
			}
			if(group_numbers) delete[] group_numbers;
			vizout << "\t\t<!-- done -->\n";
			if(configuration->LJwall_calc){
				vizout << "\t\t<!-- Draw LJ wall-->\n";
				if(!animate){
					double pos=configuration->boxlength[0]/2.0;
					if(configuration->LJwall_fixed) pos=configuration->LJwall_xm;
					vizout << "\t\t<Transform translation='" << pos << " 0.0 0.0'>\n";
					vizout << "\t\t\t<Shape>\n";
					vizout << "\t\t\t\t<Box size='0.1 " << configuration->boxlength[1] << " " << configuration->boxlength[2] << "'/>\n";
					vizout << "\t\t\t\t<Appearance>\n";
					vizout << "\t\t\t\t\t<Material diffuseColor='0.5 0.5 1.0' emissiveColor='0.5 0.5 1.0' transparency='0.4'/>\n";
					vizout << "\t\t\t\t</Appearance>\n";
					vizout << "\t\t\t</Shape>\n";
					vizout << "\t\t</Transform>\n";
					vizout << "\t\t<Transform translation='" << -pos << " 0.0 0.0'>\n";
					vizout << "\t\t\t<Shape>\n";
					vizout << "\t\t\t\t<Box size='0.1 " << configuration->boxlength[1] << " " << configuration->boxlength[2] << "'/>\n";
					vizout << "\t\t\t\t<Appearance>\n";
					vizout << "\t\t\t\t\t<Material diffuseColor='0.5 0.5 1.0' emissiveColor='0.5 0.5 1.0' transparency='0.4'/>\n";
					vizout << "\t\t\t\t</Appearance>\n";
					vizout << "\t\t\t</Shape>\n";
					vizout << "\t\t</Transform>\n";
				} else{
					string key="";
					string posValues="";
					string nposValues="";
					string BoxScales="";
					for(unsigned int i=0; i<steps; i++){
						if(configuration->LJwall_calc && configuration->LJwall_fixed){
							double scale=sqrt((V[i]*configuration->boxlength[0])/(configuration->V*2.0*Xm[i]));
							configuration->LJwall_xm=Xm[i];
							configuration->boxlength[0]=2.0*Xm[i];
							configuration->boxlength[1]*=scale; configuration->boxlength[2]*=scale;
							configuration->nndist*=scale;
							configuration->V=configuration->boxlength[0]*configuration->boxlength[1]*configuration->boxlength[2];
						} else{
							update_volume(configuration,V[i]);
							configuration->LJwall_xm=configuration->boxlength[0]/2.0;
						}
						if(i>0){
							key+=" ";
							posValues+=" ";
							nposValues+=" ";
							BoxScales+=" ";
						}
						key+=double2str((double)i/(steps-1.0));
						posValues+=double2str(configuration->LJwall_xm)+" 0.0 0.0";
						nposValues+=double2str(-configuration->LJwall_xm)+" 0.0 0.0";
						BoxScales+="0.0 "+double2str(configuration->boxlength[1])+" "+double2str(configuration->boxlength[2]);
					}
					vizout << "\t\t<TimeSensor DEF='SimGrStep' cycleInterval='" << steps << "' loop='true'/>\n";
					vizout << "\t\t<PositionInterpolator DEF='positions' key='" << key << "' keyValue='" << posValues << "'/>\n";
					vizout << "\t\t<PositionInterpolator DEF='npositions' key='" << key << "' keyValue='" << nposValues << "'/>\n";
					vizout << "\t\t<PositionInterpolator DEF='BoxScale' key='" << key << "' keyValue='" << BoxScales << "'/>\n";
					vizout << "\t\t<Transform DEF='pLJwall'>\n";
					vizout << "\t\t\t<Shape>\n";
					vizout << "\t\t\t\t<Box size='0.1 1.0 1.0'/>\n";
					vizout << "\t\t\t\t<Appearance>\n";
					vizout << "\t\t\t\t\t<Material diffuseColor='0.5 0.5 1.0' emissiveColor='0.5 0.5 1.0' transparency='0.4'/>\n";
					vizout << "\t\t\t\t</Appearance>\n";
					vizout << "\t\t\t</Shape>\n";
					vizout << "\t\t</Transform>\n";
					vizout << "\t\t<Transform DEF='nLJwall'>\n";
					vizout << "\t\t\t<Shape>\n";
					vizout << "\t\t\t\t<Box size='0.1 1.0 1.0'/>\n";
					vizout << "\t\t\t\t<Appearance>\n";
					vizout << "\t\t\t\t\t<Material diffuseColor='0.5 0.5 1.0' emissiveColor='0.5 0.5 1.0' transparency='0.4'/>\n";
					vizout << "\t\t\t\t</Appearance>\n";
					vizout << "\t\t\t</Shape>\n";
					vizout << "\t\t</Transform>\n";
					vizout << "\t\t<ROUTE fromNode='SimGrStep' fromField='fraction_changed' toNode='positions' toField='set_fraction'/>\n";
					vizout << "\t\t<ROUTE fromNode='SimGrStep' fromField='fraction_changed' toNode='npositions' toField='set_fraction'/>\n";
					vizout << "\t\t<ROUTE fromNode='SimGrStep' fromField='fraction_changed' toNode='BoxScale' toField='set_fraction'/>\n";
					vizout << "\t\t<ROUTE fromNode='positions' fromField='value_changed' toNode='pLJwall' toField='translation'/>\n";
					vizout << "\t\t<ROUTE fromNode='npositions' fromField='value_changed' toNode='nLJwall' toField='translation'/>\n";
					vizout << "\t\t<ROUTE fromNode='BoxScale' fromField='value_changed' toNode='pLJwall' toField='scale'/>\n";
					vizout << "\t\t<ROUTE fromNode='BoxScale' fromField='value_changed' toNode='nLJwall' toField='scale'/>\n";
				}
				vizout << "\t\t<!-- done -->\n";
			}
			cout << "\t-> Closing statements ...\n";
			
			//Configure footer - replace this with zvis read and defaults
			axis_t axis;
			axis.DrawAxes = 1; //draw axis markers if true
			axis.ColorX = Vec3(1,0,0); //color of x axis marker
			axis.ColorY = Vec3(0.1,0.5,1); //color of y axis marker
			axis.ColorZ = Vec3(0,1,0); //you've probably already guessed what this is (green? -- AT)
			axis.length = Vec3(1.0);
			
			//Write footer
			output_footer(vizout,axis);
			
			vizout.close();
			cout << "<- X3D file succesfully created.\n";
			
			if(animate){
				for(unsigned int j=0; j<configuration->n_oids; j++) delete[] movements[j];
				delete[] movements;
				delete[] individual_steps;
				delete[] V;
				delete[] time;
			}
		} else exit(1);
		delete[] elements;
	} else{
		cout << "No elements specified in configuration file => No output.\n";
		exit(42);
	}
	return 0;
}

