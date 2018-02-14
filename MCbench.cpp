/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * Robinson group C++ Monte Carlo simulation code.
 * Created by Dr. Robin Barnes and Lewis Johnson, September 2008
 * Based on the group's older Matlab codebase.
 * Current build 01/16/2011
 * While currently limited in capability compared to the Matlab code, this code is much faster (up to 10000x) and can be
 * compiled under Win32, OSX, or Linux without modification. It currently only requires the standard C++ library, and has
 * been tested with gcc 2.96, 3.43, 4.01, 4.12, 4.2, LLVM-gcc 4.2, LLVM-clang 1.6, icc 9, icc 10.1, and icc 11.1.
 * Settings for a simulation are currently declared in the main file.
 * Linear algebra methods are declared in VecMat.h, stats methods are declared in stdafx.h
 * Parameters that do not need to be changed by the user, and most of the error-handling code are in Setup.cpp
 * Output is generated as text files. Some are easier for Matlab to read than others.
 * Simulation parameters are punched at beginning of run for ease of monitoring a simulation.
 * Premature exit codes are (1) - initial parameters error, (2) - singular matrix in Touch, (3) Fully overlapped particle
 * (invalid distance), (4) - coordinate file loading error, and (12) - User initiated exit for debugging purposes.
 * Note that error type (3) is often the result of a corrupted build; try cleaning and rebuilding before digging into the code.
 *
 * updated Jan 24, 2011 by Andreas Tillack
 * - added more flexible molecule type system
 *
 * update on Feb 9, 2011 by Andreas Tillack
 * - added configuration file reading
 */

#include "setup.h"
#include "Latticeoids3d.h"

EXTERN_OBJLOAD(Source_tar_gz)

int main(int argc, char* argv[])
{
	cout << "\nRobinson Group Classical Monte Carlo Molecular Simulation Benchmark Tool\n";
	cout << "Compiled " << __DATE__ << " (Build " << BUILD << ", " << VERSION << ")\n";
	MC_Config* global_config = GetConfig(); // Access global simulation memory space
	Config_Data* configuration = &global_config->parameters;
	
	string conffile=""; // argv[0] is program name
	// Check for command line parameters
	if(argc>1){ // yes, there are some -- only paramter accepted is configuration filename
		conffile=argv[1]; // easy enough
	} else{
		cout << "Syntax: " << argv[0] << " <configuration file> <optional: run id>\n\n";
		exit(1);
	}
	 
	// load configuration from file
	global_config->GetFromFile(conffile.c_str());
	
	if(argc>2){ // second parameter: runid
		int id;
		bool fail=from_string(id,argv[2]);
		if(!fail) if(id<0) fail=true;
		if(!fail){
			configuration->runid=(unsigned int)id;
		} else{
			cout << "Second command line parameter (run id) needs to be an integer greater or equal to zero.\n";
			exit(1);
		}
	}
	
	// Throw out the random warning ...
	if((configuration->dogr) && (configuration->latticetype>3)){
		cout << "g(r) is currently untested for latticetype > 3\n";
		exit(1);
	}
	// Run setup and validation functions for simulation parameters. If everything passes,
	// then we'll start loading molecules.
	cout << "\nSetting up the simulation...\n";
	phys_configuration(configuration);
	setup_vdw(configuration);
	setup_electrostatics(configuration);
	check_box(configuration);
	final_validation(configuration); // One important check in there: n_oids = volume spanned by lin_dim (otherwise things will break badly, like really badly ...)
	cout << "Setup complete. Creating simulation ...\n";
	
	unsigned int SRC_SIZE=EXTERN_OBJLENGTH(Source_tar_gz);
	// Output code if wanted (default: true)
	if((configuration->output_code) && (SRC_SIZE>0)){
		string codefn = configuration->fileout + "_" + int2str(configuration->runid);
		string temp_string=codefn+".traj";
		if(configuration->trajectorynr>1){ // don't overwrite trajectory file user specifies
			codefn+="_"+int2str(configuration->trajectorynr);
		}
		codefn+="_code.tar.gz";
		cout << "\n-> Outputting simulation source code to\n\t" << codefn << "\n";
		FILE* codefile = fopen(codefn.c_str(),"w");
		const unsigned char* p=EXTERN_OBJDATA(Source_tar_gz);
		for(unsigned int i=0; i<SRC_SIZE; i++) fputc(*p++,codefile);
		fclose(codefile);
		cout << "<- Done.\n\n";
	}
	
	// Done with loading parameters. Prepare to run!
	MC_Elements* BeOids = new MC_Elements(configuration);
	
	// Get current time and initialize memory for RNG seed
	cout << "Initializing RNG seed...\n";
#ifdef USE_CMWC4096
	__uint32_t tempseed;
#else
	__int32_t tempseed;
#endif
	time_t rawtime;
	struct tm * batchstarttime;
	time (&rawtime);
	batchstarttime = localtime (&rawtime);
	
	// Open batch summary file
	cout << "Opening batch summary file...\n";
	open_summary(configuration);
	cout << "Done. Batch started on " << asctime(batchstarttime) << "\n";
	
	// Execute
	configuration->trackchange = false; // reset E and T to initial.
	if (!configuration->test){
	// get RNG seed from time. This won't repeat until 2038 :D
#ifdef USE_CMWC4096
		tempseed = time(NULL);
#else
		tempseed = -1*time(NULL);
#endif
		tempseed /= configuration->runid+1; // make RNGs unique even if two are started at the same time
		tempseed += configuration->runid; // in case we do lapse there can still be some variation
	} else{ // Locked seed for test mode
		cout << "DIAGNOSTIC MODE ON\n";
#ifdef USE_CMWC4096
		tempseed = 8;
#else
		tempseed = -7;
#endif
	}
	configuration->rngseed = tempseed; // record the seed
#ifdef USE_CMWC4096
	init_CMWC4096(tempseed);
	zigset();
#else
	configuration->idum = &tempseed; // idum changes each time ran2() is called
#endif
	
	// Run simulation
	LatticeOids3d* lattice = new LatticeOids3d(*BeOids); // This constructor loads molecules into BeOids.
	if(configuration->add_slowly<0) check_oids(*BeOids);
	BeOids->benchmark_system();
	delete lattice;
	delete BeOids;
	
	// Get ending time
	struct tm * batchendtime;
	time (&rawtime);
	batchendtime = localtime (&rawtime);
	
	// Clean up
	string summfn = configuration->fileout + "_summary_" + int2str(configuration->runid) + ".dat";
	FILE* summ = fopen(summfn.c_str(), "a");
	fprintf(summ, "\n#END\n");
	fclose(summ);
	
	// Normal termination
	cout << "Batch finished on " << asctime(batchendtime) << "\nDone.\n";
	return 0;
}

