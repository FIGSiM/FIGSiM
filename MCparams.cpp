/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * Robinson group Monte Carlo simulation parameter extraction tool
 *
 * Command line utility with three arguments, (1) the configuration file, (2) the section name, and (3) the parameter name
 * Usage example: traj2stat acetonitrile.conf "Simulation Parameters" "randsteps"
 * Created by Andreas Tillack, June 2012
 */

#include "MC_Config.h"

int main (int argc, char * const argv[]) {
	MC_Config config;
	
	string filename="";
	string section="";
	char* param;
	bool get_filelist=false;
	bool want_trajectory=false;
	bool want_number=false;
	
	// Check to make sure there are enough command line parameters
	if(argc<4){
		cout << "Syntax: " << argv[0] << " <configuration file> <section name> <parameter name>\n";
		exit(1);
	} else{
		// Configuration file name
		filename=argv[1];
		section=argv[2];
		param=argv[3];
		if(compare_strings(param,"files_used")){
			get_filelist=true;
		}
		if(argc>4){
			if(compare_strings(argv[4],"want_traj")) want_trajectory=true;
			if(compare_strings(argv[4],"want_number")) want_number=true;
		}
	}
	
	unsigned int size, position;
	bool is_trajectory=false;
	string conf_filename;
	readfile conffile;
	string files_used="";
	conffile.filename=filename;
	conffile.directory="";
	GetDirectory(conffile);
	
	do{
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
			cout << "Could not open file " << conffile.directory+conffile.filename << ".\n";
			exit(1);
		}
		// Get file size
		conffile.file.seekg(0,ifstream::end);
		size=conffile.file.tellg();
		conffile.file.seekg(0);
		// Sanity checks
		if(size==0){
			cout << "File has no content.\n";
			exit(1);
		}
		
		if(!want_trajectory){
			// Test for trajectory file
			position=0;
			unsigned int include_pos=0;
			string subname="";
			char* first=config.GetSection(&conffile,"General",position,include_pos,subname,false,false);
			conf_filename="";
			if(first){
				config.SetParam(conf_filename,"configuration",first,"");
				delete[] first;
			}
			is_trajectory=(conf_filename!="");
			
			if((size>MAXCONFSIZE) && (!is_trajectory)){ // config files should not be bigger than MAXCONFSIZE/(1024*1024) MB (only exception: trajectory files)
				cout << "Configuration file is too big (>" << MAXCONFSIZE/(1024*1024) << " MB).\n";
				exit(1);
			}
			if(is_trajectory){
				conffile.file.close();
			}
		}
		if(files_used!="") files_used+="\n";
		if(conffile.directory=="") files_used+="./";
		files_used+=conffile.directory+"\t"+conffile.filename;
	} while(is_trajectory && !want_trajectory);
	
	// First, get simulation parameters ...
	position=0;
	
	string subname;
	unsigned int include_position=0;
	
	char* constant_block=config.GetSection(&conffile,section,position,include_position,subname,true,!want_trajectory,&files_used);
	if(get_filelist){
		position=0;
		unsigned int num_groups=config.GetNrSections(&conffile,"Group",position);
		char* group;
		for(unsigned int i=0; i<num_groups; i++){
			subname="";
			group=config.GetSection(&conffile,"Group",position,subname); // subname is element name
			if(group){
				readfile mol2;
				config.SetParam(mol2.filename,"structure_mol2",group,"");
				if(mol2.filename!=""){
					GetDirectory(mol2);
					if(files_used!="") files_used+="\n";
					mol2.directory=conffile.directory+mol2.directory;
					if(mol2.directory=="") files_used+="./";
					files_used+=mol2.directory+"\t"+mol2.filename;
				}
				delete[] group;
			}
		}
		cout << files_used << "\n";
	} else{
		if(!constant_block){
			cout << "Section [" << section << "] does not exist in configuration file.\n";
			exit(1);
		}
		if(want_number){
			double answer;
			config.SetParam(answer,param,constant_block);
			cout << answer << "\n";
		} else{
			string answer;
			config.SetParam(answer,param,constant_block);
			cout << answer << "\n";
		}
	}
	if(constant_block) delete[] constant_block;
}
