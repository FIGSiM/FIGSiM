/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * Robinson group Monte Carlo simulation code statistics analysis tool
 *
 * Command line utility with three arguments, (1) the trajectory file, and, optionally, (2) the starting step, (3) the ending step
 * Usage example: traj2stat test0.traj 1000
 * Created by Lewis Johnson and Andreas Tillack, November 2011
 * - minor bugfixes since original inception by LEJ:
 *   - initial step handling works correctly now
 *   - ending step is now properly taken from third command line parameter (was second)
 *   - standard deviation for epsilon and epsilon_RF are now sqrt of variance (was variance)
 *   - cleaned up most NaNs by making sure the variances are >=EPS (otherwise they are set to 0.0)
 */

#define BLOCK_AVERAGES 0 // if set to anything but 0, block averaged values are used for derived quantities (e.g. total dipole moment^2) instead of per-step values (*should not be used*)
#define AUTOCORR_CUTOFF 0.95 // fraction of frames after which to stop autocorrelation calculation
#define CALC_ENTROPY_dG 1 // if set to anything greater than 0, calculate entropy and delta G using total system energy and ideal gas solution (experimental and most likely wrong, default should be 0)

#include "MC_Config.h"
#include "ScalarMat.h"
#include "VecMat.h"
#include "setup.h"

///Auto-generate statistics filename from trajectory filename
//Used for consistency with traj2x3d, this may be unnecessary depending on how configuration->fileout
//is currently used, and one convention should probably be used in the future
string name_statfile(string &trajfile)
{
	//Get length of streng
	unsigned int trajf_len = trajfile.length();
	
	//Copy filename without extension
	trajfile.begin();
	string coorfile = "";
	for (unsigned int i = 0; i < (trajf_len - 5); i++) {
		coorfile += trajfile[i];
	}
	
	//Assemble new filename
	coorfile += ".stat";
	return coorfile;
}


///Storage struct for stats
struct stats_t {
	//Filename and length parameters
	string stats_fn;
	unsigned int start_step;
	unsigned int end_step;
	unsigned int nr_steps;
	unsigned int block_size;
	unsigned int running_average;
	unsigned int num_blocks;
	string potentials[NUM_V_STORE];
	
	//Thermodynamic statistics
	Vstore MeanU_NkT; //Dimensionless energy components
	Vstore MeanU_SI; //Energy components in KJ/mol
	Vstore StdU_NkT; //Standard error of dimensionless energy components
	Vstore StdU_SI; //Standard error of SI energy components
	double TotalU_NkT[2]; //Total energy and its standard error (dimensionless)
	double TotalU_SI[2]; //Total energy and its standard error (SI)
#if CALC_ENTROPY_dG>0
	double S[2]; // Entropy
	double JarzynskiG[2]; // Total Jarzynski energy average
#endif
	double VarTotalU_NkT[2]; //Variance of energy and its standard error (dimensionless)
	double VarTotalU_SI[2];
	double Heat_Capacity[2]; //Cv or Cp, depending on ensemble and its stdev
	double Mass;
	double Volume[2]; //Volume of the box in A^3
	double Number_Density[2]; //Number density (concentration) in particles/A^3
	double Molarity[2]; //Concentration in mol/L
	bool NpT;
	
	//Electrostatic statistics
	Mat33 CovarM; //Covariance matrix of the total dipole moment
	Vec3 MeanM; //Total dipole moment in D
	Vec3 StdM; //Stderr of dipole moment in D
	double NormM[2]; //Dipole magnitude and its standard error
	double M2[2]; //Squared dipole moment and its standard error (used for calculating dielectric constant from fluctuations)
	double ME[2]; //Dipole moment in field direction and its variance
	double VarM[2]; //Scalar variance of dipole moment and its standard error (used for calculating dielectric constant from fluctuations)
	double Epsilon[2]; //Calculated dielectric constant from simulation and its stderr
	double Epsilon_s[2]; //Reaction field dielectric constant and its stderr
	double MaxwellField[2]; //Externally-measured average field in material and its stderr (V/um)
	double CavityField[2]; //Local field E0 inside the reaction field cavity and its stderr (V/um)
	
	//Order statistics
	Vec3 MeanCos; //<cos^n(theta)> wrt z-axis
	Vec3 StdCos; //Standard deviation thereof
	double Var_Skew[2]; //Variance and skewness of <cos(theta)>
	Vec3 P_n_theta; //<P_n(theta)> wrt z-axis
	Vec3 Std_Pn; //Standard deviation thereof
	
	//Acceptance statistics
	double RMS_Move[2]; //RMS translational and rotational displacements in A
	Vec3 Acceptance;
};


///Print statistics to file
void print_statistics(stats_t Statistics, Config_Data* configuration)
{
	ofstream statfn;
	statfn.open(Statistics.stats_fn.c_str());
	if(statfn.fail()) {
		cout << "ERROR: Unable to open output file.\n";
		exit(1);
	}
	statfn.setf(ios::fixed,ios::floatfield);
	statfn.precision(6);
	cout << "\t-> Outputting statistics data to " << Statistics.stats_fn << "\n";
	
	// Write out tab-delimited data with descriptor in first column
	statfn << "# Statistics from step " << Statistics.start_step << " to " << Statistics.end_step << " using a block size of " << Statistics.block_size << "\n\n";
	statfn << "[Thermodynamics]\n\n";
	statfn << "Average Volume:\t" << Statistics.Volume[0] << "\t+/-\t" << Statistics.Volume[1] << "\tAngström^3\n";
	statfn.setf(ios::scientific,ios::floatfield);
	statfn << "Number Density:\t" << Statistics.Number_Density[0] << "\t+/-\t" << Statistics.Number_Density[1] << "\tcm^(-3)\n";
	statfn.setf(ios::fixed,ios::floatfield);
	statfn << "Average Molarity:\t" << Statistics.Molarity[0] << "\t+/-\t" << Statistics.Molarity[1] << "\tM\n";
	statfn << "Average Reduced Dipole Density:\t" << 4.0*pi*Statistics.Number_Density[0]*1E-24*qqpwr(configuration->max_dipole/configuration->N,2)/(9.0*configuration->kT) << "\t+/-\t" << 4.0*pi*Statistics.Number_Density[1]*1E-24*qqpwr(configuration->max_dipole/configuration->N,2)/(9.0*configuration->kT) << "\n"; // 10^20/cm^3 = 10^20/(10^8 Ang)^3 = 10^-4 / Ang^3
	// Calculate average density if system mass exists (1 Angström^3 = 1E-24 cm^3
	if(Statistics.Mass>EPS) statfn << "Average Density:\t" << Statistics.Mass/(NA*1E-24*Statistics.Volume[0]) << "\t+/-\t" << Statistics.Mass/(NA*1E-24*Statistics.Volume[0]*Statistics.Volume[0])*Statistics.Volume[1] << "\tg/cm^3\n";
	statfn << "\n";
	statfn << "# Potential energies in NkT units\n";
	for (unsigned int i = 0; i < NUM_V_STORE; i++) statfn << "<" << Statistics.potentials[i] << ">\t" << Statistics.MeanU_NkT[i] << "\t+/-\t" << Statistics.StdU_NkT[i] << "\tNkT\n";
	statfn << "# Potential energies in SI units\n";
	for (unsigned int i = 0; i < NUM_V_STORE; i++) statfn << "<" << Statistics.potentials[i] << ">\t" << Statistics.MeanU_SI[i] << "\t+/-\t" << Statistics.StdU_SI[i] << "\tkJ/mol\n";
	if(Statistics.NpT){
		statfn << "\n# Total enthalpy in NkT units\n";
		statfn << "<H>\t" << Statistics.TotalU_NkT[0] << "\t+/-\t" << Statistics.TotalU_NkT[1] << "\tNkT\n";
		statfn << "Var(H)\t" << Statistics.VarTotalU_NkT[0] << "\t+/-\t" << Statistics.VarTotalU_NkT[1] << "\t(NkT)^2\n";
		statfn << "# Total enthalpy in SI units\n";
		statfn << "<H>\t" << Statistics.TotalU_SI[0] << "\t+/-\t" << Statistics.TotalU_SI[1] << "\tkJ/mol\n";
		statfn << "Var(H)\t" << Statistics.VarTotalU_SI[0] << "\t+/-\t" << Statistics.VarTotalU_SI[1] << "\t(kJ/mol)^2\n\n";
		statfn << "# Heat capacity\n";
		statfn << "<Cp>";
	} else{
		statfn << "\n# Total energy in NkT units\n";
		statfn << "<U>\t" << Statistics.TotalU_NkT[0] << "\t+/-\t" << Statistics.TotalU_NkT[1] << "\tNkT\n";
		statfn << "Var(U)\t" << Statistics.VarTotalU_NkT[0] << "\t+/-\t" << Statistics.VarTotalU_NkT[1] << "\t(NkT)^2\n";
		statfn << "# Total energy in SI units\n";
		statfn << "<U>\t" << Statistics.TotalU_SI[0] << "\t+/-\t" << Statistics.TotalU_SI[1] << "\tkJ/mol\n";
		statfn << "Var(U)\t" << Statistics.VarTotalU_SI[0] << "\t+/-\t" << Statistics.VarTotalU_SI[1] << "\t(kJ/mol)^2\n\n";
		statfn << "# Heat capacity\n";
		statfn << "<Cv>";
	}
	statfn << "\t" << Statistics.Heat_Capacity[0] << "\t+/-\t" << Statistics.Heat_Capacity[1] << "\tJ/mol-K\n\n\n";
#if CALC_ENTROPY_dG>0
	statfn << "Entropy\t" << Statistics.S[0] << "\tJ/mol-K\n\n\n";
	statfn << "# Jarzynski equality related values\n";
	statfn << "dG\t" << Statistics.JarzynskiG[0] << "\t+/-\t" << Statistics.JarzynskiG[1] << "\tkJ/mol\n";
	statfn << "RTln(V/Angström)\t" << configuration->kT*perg_to_kJ_per_mol*log(Statistics.Volume[0]) << "\t+/-\t" << configuration->kT*perg_to_kJ_per_mol*Statistics.Volume[1]/Statistics.Volume[0] << "\tkJ/mol\n\n\n";
#endif
	statfn << "[Electrostatics]\n\n";
	statfn << "<M_x>\t" << Statistics.MeanM.vec[0] << "\t+/-\t" << Statistics.StdM.vec[0] << "\tD\n";
	statfn << "<M_y>\t" << Statistics.MeanM.vec[1] << "\t+/-\t" << Statistics.StdM.vec[1] << "\tD\n";
	statfn << "<M_z>\t" << Statistics.MeanM.vec[2] << "\t+/-\t" << Statistics.StdM.vec[2] << "\tD\n\n";
	statfn << "<M^2>\t" << Statistics.M2[0] << "\t+/-\t" << Statistics.M2[1] << "\tD^2\n";
	statfn << "<|M|>\t" << Statistics.NormM[0] << "\t+/-\t" << Statistics.NormM[1] << "\tD\n";
	statfn << "var(M)\t" << Statistics.VarM[0] << "\t+/-\t" << Statistics.VarM[1] << "\tD^2\n\n";
	statfn << "covar(M) (D^2):\n";
	statfn << Statistics.CovarM.M3Str() << "\n\n";
	if(Statistics.MaxwellField[0]>EPS){
		statfn << "<M_E>\t" << Statistics.ME[0] << "\tD\n";
		statfn << "<M_E^2> - <M_E>^2\t" << Statistics.ME[1] << "\tD^2\n\n";
	}
	statfn << "<epsilon>\t" << Statistics.Epsilon[0] << "\t+/-\t" << Statistics.Epsilon[1] << "\n";
	statfn << "<epsilon_RF>\t" << Statistics.Epsilon_s[0] << "\t+/-\t" << Statistics.Epsilon_s[1] << "\n\n";
	statfn << "Maxwell Field <E>\t" << Statistics.MaxwellField[0] << "\t+/-\t" << Statistics.MaxwellField[1] << "\tV/µm\n";
	statfn << "Cavity Field <E_0>\t" << Statistics.CavityField[0] << "\t+/-\t" << Statistics.CavityField[1] << "\tV/µm\n\n\n";
	statfn << "[Order]\n\n";
	statfn << "<cos(theta)>\t" << Statistics.MeanCos.vec[0] << "\t+/-\t" << Statistics.StdCos.vec[0] << "\n";
	statfn << "<cos^2(theta)>\t" << Statistics.MeanCos.vec[1] << "\t+/-\t" << Statistics.StdCos.vec[1] << "\n";
	statfn << "<cos^3(theta)>\t" << Statistics.MeanCos.vec[2] << "\t+/-\t" << Statistics.StdCos.vec[2] << "\n\n";
	statfn << "<P_1(theta)>\t" << Statistics.P_n_theta.vec[0] << "\t+/-\t" << Statistics.Std_Pn.vec[0] << "\n";
	statfn << "<P_2(theta)>\t" << Statistics.P_n_theta.vec[1] << "\t+/-\t" << Statistics.Std_Pn.vec[1] << "\n";
	statfn << "<P_3(theta)>\t" << Statistics.P_n_theta.vec[2] << "\t+/-\t" << Statistics.Std_Pn.vec[2] << "\n\n";
	statfn << "var(cos(theta)):\t" << Statistics.Var_Skew[0] << "\n";
	statfn << "skew(cos(theta)):\t" << Statistics.Var_Skew[1] << "\n\n\n";
	statfn << "[Acceptance]\n\n";
	statfn << "Acceptance (%):\t" << Statistics.Acceptance.V3Str() << "\n";
	statfn << "RMS Translation (Angström):\t" << Statistics.RMS_Move[0] << "\n";
	statfn << "RMS Rotation (rad):\t" << Statistics.RMS_Move[1] << "\n\n";
	statfn.close();
}

///Calculate average dielectric constant and electric fields
// - depends on Volume (+stderr), StdM, M2 (+stderr), VarM (+stderr), MeanM
// - returns epsilon, epsilon_RF, MaxwellField, and CavityField (+stderrs)
void average_dielectric(stats_t& Statistics, Config_Data* configuration, Vec3 edirection)
{
	// Allocate memory
	double epsilon, Var_epsilon, q, qs, Ecorr, deps_dV, deps_dME, dq_dV, dq_dME, deps_dq, deps_ds, dE_deps, dE_ds;
	double eps_s = configuration->epsilon;
	double Var_eps_s = 0.0;
	double n2 = configuration->n2;
	double n4 = n2*n2;
	double n_blocks = (double)(Statistics.num_blocks);
	
	// Unpack standard errors needed for error propagation and convert back to variances
	double Var_V = Statistics.Volume[1]*Statistics.Volume[1]*n_blocks;
	double Var_ME = qqpwr(Statistics.StdM*edirection,2)*n_blocks;
	double Var_M2 = Statistics.M2[1]*Statistics.M2[1]*n_blocks;
	double Var_Var_M = Statistics.VarM[1]*Statistics.VarM[1]*n_blocks;
	
	// Dipole density scaling factor
	double y_scale = 4*pi/(9*Statistics.Volume[0]*configuration->kT); // Conveniently also the derivative of the RHS of the Kirkwood equation wrt M2 *and* V
	double dy_dV = -4*pi/(9*Statistics.Volume[0]*Statistics.Volume[0]*configuration->kT);
	
	double inv_n = 1.0/(double)(Statistics.num_blocks);
	
	// Calculate scalar dielectric constant based on external field
	if(!configuration->noEfield){ // Constant external potential (E)
		if(configuration->realfield){
			epsilon = 4.0*pi*(Statistics.MeanM*edirection)/(configuration->Efield.V3Norm()*MV_per_m_to_perg_per_Debye*Statistics.Volume[0]) + configuration->n2; // [D/(perg/D*Ang^3)] = D^2/Ang^3 / perg = perg/perg = 1
			deps_dV = -4.0*pi*(Statistics.MeanM*edirection)/(configuration->Efield.V3Norm()*MV_per_m_to_perg_per_Debye*Statistics.Volume[0]*Statistics.Volume[0]);
			deps_dME = 4.0*pi/(configuration->Efield.V3Norm()*MV_per_m_to_perg_per_Debye*Statistics.Volume[0]);
			Var_epsilon = (deps_dV*deps_dV)*Var_V + (deps_dME*deps_dME)*Var_ME;
			if(configuration->dyneps){
				eps_s = epsilon;
				Var_eps_s = Var_epsilon;
			}
		} else{ // Constant effective/internal field (E0)
			q = (4.0*pi*Statistics.MeanM.vec[2])/(3.0*configuration->Efield.V3Norm()*MV_per_m_to_perg_per_Debye*Statistics.Volume[0]);
			dq_dV = (-4.0*pi*(Statistics.MeanM*edirection))/(3.0*configuration->Efield.V3Norm()*MV_per_m_to_perg_per_Debye*Statistics.Volume[0]*Statistics.Volume[0]);
			dq_dME = (4.0*pi)/(3.0*configuration->Efield.V3Norm()*MV_per_m_to_perg_per_Debye*Statistics.Volume[0]);
			double Var_q = (dq_dME*dq_dME)*Var_ME + (dq_dV)*(dq_dV)*Var_V;
			if(!configuration->dyneps){
				deps_dq = (6.0*eps_s)/(n2+2.0*eps_s-3.0*q) + 3.0*(6.0*q*eps_s+2.0*n2*eps_s+n4)/qqpwr((n2+2.0*eps_s-3.0*q),2);
				epsilon = (6.0*q*eps_s+2.0*n2*eps_s+n4)/(n2+2.0*eps_s-3.0*q);
				Var_epsilon = (deps_dq*deps_dq)*Var_q; //Epsilon_s is constant so not included in propagation
			} else{
				epsilon = 0.25*n2+2.25*q+0.75*sqrt(n4+2.0*n2*q+9.0*q*q);
				deps_dq = 2.25 + (0.75*n2+6.75*q)/sqrt(n4+2.0*n2*q+9.0*q*q);
				Var_epsilon = (deps_dq*deps_dq)*Var_q;
				eps_s = epsilon;
				Var_eps_s = Var_epsilon;
			}
		}
	} else{ // Calculate dielectric constant from variance of total dipole moment (Kirkwood method)
		q = y_scale*(Statistics.VarM[0]); // Calculate product of dipole density and Var(M), typically called y
		double Var_q = (dy_dV*dy_dV)*Var_V + (y_scale*y_scale)*Var_Var_M; // Very convenient for error propagation...
		if(configuration->dyneps){ // Calculate average epsilon_s (reaction field dielectric constant) from mean square M
			qs = y_scale*Statistics.M2[0];
//			qs = y_scale*(Statistics.M2[0]-Statistics.MeanM*Statistics.MeanM);
			double Var_qs = (dy_dV*dy_dV)*Var_V + (y_scale*y_scale)*Var_M2;
			eps_s = 0.25*n2+2.25*qs+0.75*sqrt(n4+2.0*n2*qs+9.0*qs*qs);
			deps_dq = 2.25 + (0.75*n2+6.75*qs)/sqrt(n4+2.0*n2*qs+9.0*qs*qs);
			Var_eps_s = (deps_dq*deps_dq)*Var_qs; // Variance of reaction field dielectric constant
		}
		epsilon = 0.25*n2+2.25*q+0.75*sqrt(n4+2.0*n2*q+9.0*q*q);
		deps_dq = (6.0*eps_s)/(n2+2.0*eps_s-3.0*q) + 3.0*(6.0*q*eps_s+2.0*n2*eps_s+n4)/qqpwr((n2+2.0*eps_s-3.0*q),2);
		deps_ds = (6.0*q+2.0*n2)/(n2+2.0*eps_s-3.0*q)+ 2.0*(6.0*q*eps_s+2.0*n2*eps_s+n4)/qqpwr((n2+2.0*eps_s-3.0*q),2);
		Var_epsilon = (deps_dq*deps_dq)*Var_q + (deps_ds*deps_ds)*Var_q; // Epsilon_s not constant here
	}
	
	Statistics.Epsilon[0] = epsilon;
	Statistics.Epsilon[1] = sqrt(Var_epsilon*inv_n);
	Statistics.Epsilon_s[0] = eps_s;
	Statistics.Epsilon_s[1] = sqrt(Var_eps_s*inv_n);
	
	// Calculate cavity field and externally measured (Maxwell) field
	if(configuration->realfield){
		Ecorr = (2.0*eps_s+epsilon)/(2.0*eps_s + configuration->n2);
		dE_deps = 1.0/(2.0*eps_s+1.0);
		dE_ds = 2.0/(2.0*eps_s+1.0) - 2.0*(2.0*eps_s+epsilon)/qqpwr((2.0*eps_s+1.0),2);
		Statistics.MaxwellField[0] = configuration->Efield.V3Norm();
		Statistics.MaxwellField[1] = 0.0;
		Statistics.CavityField[0] = configuration->Efield.V3Norm()*Ecorr;
		Statistics.CavityField[1] = sqrt(((dE_deps*dE_deps)*Var_epsilon + (dE_ds*dE_ds)*Var_eps_s)/n_blocks)*configuration->Efield.V3Norm(); // directly as standard error
	} else{
		Ecorr = (2.0*eps_s + configuration->n2)/(2.0*eps_s + epsilon);
		dE_deps = -1.0*(2.0*eps_s + 1.0)/qqpwr((2.0*eps_s+epsilon), 2);
		dE_ds = 2.0/(2.0*eps_s+epsilon) + (4.0*eps_s+2.0)/qqpwr((2.0*eps_s+epsilon),2);
		Statistics.CavityField[0] = configuration->Efield.V3Norm();
		Statistics.CavityField[1] = 0.0;
		Statistics.MaxwellField[0] =configuration->Efield.V3Norm()*Ecorr;
		Statistics.MaxwellField[1] = sqrt(((dE_deps*dE_deps)*Var_epsilon + (dE_ds*dE_ds)*Var_eps_s)/n_blocks)*configuration->Efield.V3Norm(); // directly as standard error
	}
}


///Calculate average properties and other statistics
void calc_means(Config_Data* configuration, stats_t& Statistics, Vstore* energy_components, Vec3* total_dipole, Vec3* cos_n, double* V, double* msmoved, unsigned int* accepted, unsigned int* tries, bool calc_dipole_autocorr, Vec3 edirection)
{	
	//Temporary accumulators
	double UCompTemp, TotalUTemp, StdUTemp, VarTotalUTemp, ScalarM2Temp;
	Vec3 CosnTemp, MMTemp, StdMTemp;
	
	
	// Allocate memory for block averages
	double* BlockTotalU = new double[Statistics.num_blocks];
	double* BlockTotalU2 = new double[Statistics.num_blocks];
	double* BlockTotalU4 = new double[Statistics.num_blocks];
	double* BlockV = new double[Statistics.num_blocks];
	double* BlockME = new double[Statistics.num_blocks];
	Vec3* BlockM = new Vec3[Statistics.num_blocks];
	Mat33* BlockM2 = new Mat33[Statistics.num_blocks];
	Vec3* BlockCosn = new Vec3[Statistics.num_blocks];
	Vstore* BlockU = new Vstore[Statistics.num_blocks];
	for(unsigned int i = 0; i < Statistics.num_blocks; i++){
		BlockTotalU[i] = 0.0;
		BlockTotalU2[i] = 0.0;
		BlockTotalU4[i] = 0.0;
		BlockV[i] = 0.0;
		BlockME[i] = 0.0;
		BlockM[i].V3Zeros();
		BlockM2[i].M3Zeros();
		BlockCosn[i].V3Zeros();
		for(unsigned int j = 0; j < NUM_V_STORE; j++) BlockU[i][j] = 0.0;
	}
	

	//Block average loop - calculate block sums and convert to averages once each block
	cout << "\tSeparating data into blocks of " << Statistics.block_size << " steps...\n";
	//Denominators for averages (1/n_blocks, 1/NkT, etc.)
	double inv_bs = 1.0/(double)(Statistics.block_size);
	double NkT = (double)(configuration->N)*configuration->kT;
	double inv_NkT = 1.0/NkT;
	double kJ_mol = perg_to_kJ_per_mol/((double)configuration->N); //Precalculated conversion factor from pico-ergs to kJ/mol
	double vartemp;
	
#if CALC_ENTROPY_dG>0
	double inv_nrsteps=1.0/Statistics.nr_steps;
	bool firstJG=true;
	double U0=0.0;
	long double JG=0.0;
	long double JG2=0.0;
	long double Usum=0.0;
	long double Z=0.0;
#endif
	
	for(unsigned int step = 0; step < Statistics.nr_steps; step++){
		//Get block index
		unsigned int i = step/Statistics.block_size;
		
		//Energy components
		TotalUTemp = 0.0;
		for(unsigned int j = 0; j < NUM_V_STORE; j++){
			UCompTemp = energy_components[step][j];
			BlockU[i][j] += UCompTemp;
			TotalUTemp += UCompTemp;
		}
		// calculate enthalpy when in NpT
		if(Statistics.NpT) TotalUTemp+=configuration->pext*V[step];
#if CALC_ENTROPY_dG>0
		if(firstJG){
			U0=TotalUTemp;
			firstJG=false;
		}
		long double ej=expl(-(TotalUTemp-U0)/configuration->kT);
		Usum += ej*TotalUTemp;
		Z += ej;
		// Jarzynski energy average
		JG += inv_nrsteps*ej;
		JG2 += inv_nrsteps*ej*ej;
#endif
		
		//Total Energy
		BlockTotalU[i] += TotalUTemp;
#if BLOCK_AVERAGES==0
		BlockTotalU2[i] += (TotalUTemp*TotalUTemp);
		BlockTotalU4[i] += qqpwr(TotalUTemp,4);
#endif
		
		//Volume
		BlockV[i] += V[step];
		
		//Total Dipole (M)
		BlockM[i] += total_dipole[step];
		BlockME[i] += total_dipole[step]*edirection;
#if BLOCK_AVERAGES==0
		BlockM2[i] += total_dipole[step].V3TensProd(total_dipole[step]);
#endif
		
		//cos^n
		BlockCosn[i] += cos_n[step];
		
		if((step+1)%Statistics.block_size==0){
			for(unsigned int j = 0; j < NUM_V_STORE; j++) BlockU[i][j] *= inv_bs;
			BlockTotalU[i] *= inv_bs;
#if BLOCK_AVERAGES==0
			BlockTotalU2[i] *= inv_bs;
			BlockTotalU4[i] *= inv_bs;
#else
			BlockTotalU2[i] = qqpwr(BlockTotalU[i],2);
			BlockTotalU4[i] = qqpwr(BlockTotalU[i],4);
#endif
			BlockV[i] *= inv_bs;
			BlockME[i] *= inv_bs;
			BlockM[i] *= inv_bs;
#if BLOCK_AVERAGES==0
			BlockM2[i] *= inv_bs;
#else
			BlockM2[i] = BlockM[i].V3TensProd(BlockM[i]);
#endif
			BlockCosn[i] *= inv_bs;
		}
	}
#if CALC_ENTROPY_dG>0
	Usum/=Z;
	Z*=inv_nrsteps;
	Statistics.S[0]=kB*perg_to_kJ_per_mol*1000*(Usum/configuration->kT+logl(Z)-U0/configuration->kT);
	// Jarzynski delta G goes here
	Statistics.JarzynskiG[0]=-configuration->kT*(logl(JG)-U0/configuration->kT)*kJ_mol;
	Statistics.JarzynskiG[1]=fabs(configuration->kT*sqrt(JG2-JG*JG)/JG)*kJ_mol;
#endif
	//Accumulators for final averages - first element is sum, second is sum of squares
	Vstore SumU[2];
	for(unsigned int i=0; i<NUM_V_STORE; i++){
		SumU[0][i]=0.0;
		SumU[1][i]=0.0;
	}
	double SumTotalU[2] = {0.0,0.0};
	double SumTotalU2[2] = {0.0,0.0};
	double SumV[2] = {0.0,0.0};
	double SumME[2] = {0.0,0.0};
	double SumScalarM2[2] = {0.0,0.0};
	double SumNormM = 0; //No sum of squares term, as it is the same as SumScalarM2
	Vec3 SumM[2]; SumM[0].V3Zeros(); SumM[1].V3Zeros();
	Mat33 SumM2(0.0); //No sum of squares term
	Vec3 SumCosn[2]; SumCosn[0].V3Zeros(); SumCosn[1].V3Zeros();
	
	//Accumulate sums of block averages (element 0) and their squares (element 1)
	for(unsigned int i = 0; i < Statistics.num_blocks; i++){
		
		//Energy and its square
		for (unsigned int j = 0; j < NUM_V_STORE; j++){
			SumU[0][j] += BlockU[i][j];
			SumU[1][j] += BlockU[i][j]*BlockU[i][j];
		}
		SumTotalU[0] += BlockTotalU[i];
		SumTotalU[1] += BlockTotalU[i]*BlockTotalU[i];
		SumTotalU2[0] += BlockTotalU2[i];
		SumTotalU2[1] += BlockTotalU4[i];
		
		//Volume
		SumV[0] += BlockV[i];
		SumV[1] += BlockV[i]*BlockV[i];
		
		//total dipole moment in E-direction
		SumME[0] += BlockME[i];
		SumME[1] += BlockME[i]*BlockME[i];
		
		//Dipole moment, its norm, scalar square, and matrix square
		SumM[0] += BlockM[i];
		Vec3Amult(BlockM[i], BlockM[i], MMTemp); //As a sidenote, the 'scalar math' header really needs cleanup - LEJ
		SumM[1] += MMTemp;
		ScalarM2Temp = BlockM2[i].M3Trace(); //Simply <M'*M>, used for dielectric
		SumScalarM2[0] += ScalarM2Temp;
		SumScalarM2[1] += ScalarM2Temp*ScalarM2Temp;
		SumNormM += sqrt(ScalarM2Temp);
		SumM2 += BlockM2[i]; //<M*M'> for covariance matrix
		
		//Order parameters
		SumCosn[0] += BlockCosn[i];
		Vec3Amult(BlockCosn[i], BlockCosn[i], CosnTemp);
		SumCosn[1] += CosnTemp;
	}
	
	//Denominators for overall averages (1/n_blocks, 1/NkT, etc.)
	double inv_n = 1.0/(double)(Statistics.num_blocks);
	double inv_n2 = inv_n*inv_n;
	double inv_sqrt_n = sqrt(inv_n);
	
	//Calculate averages and standard errors (standard error of the mean), convert units if needed, and store in Statistics
	//Energy components
	cout << "\tCalculating average thermodynamic properties...\n";
	for (unsigned int j = 0; j < NUM_V_STORE; j++){
		Statistics.MeanU_NkT[j] = SumU[0][j]*inv_n*inv_NkT;
		Statistics.MeanU_SI[j] = SumU[0][j]*inv_n*kJ_mol;
		StdUTemp = sqrt(SumU[1][j]*inv_n - (SumU[0][j]*SumU[0][j]*inv_n2))*inv_sqrt_n;
		Statistics.StdU_NkT[j] = StdUTemp*inv_NkT;
		Statistics.StdU_SI[j] = StdUTemp*kJ_mol;
	}
	
	//Total energy and its variance
	Statistics.TotalU_NkT[0] = SumTotalU[0]*inv_n*inv_NkT;
	Statistics.TotalU_SI[0] = SumTotalU[0]*inv_n*kJ_mol;
	VarTotalUTemp = SumTotalU[1]*inv_n - (SumTotalU[0]*SumTotalU[0])*inv_n2;
	StdUTemp = sqrt(VarTotalUTemp)*inv_sqrt_n; // stddev_mean = 1/sqrt(n)*stddev
	Statistics.TotalU_NkT[1] = StdUTemp*inv_NkT;
	Statistics.TotalU_SI[1] = StdUTemp*kJ_mol;
	VarTotalUTemp = SumTotalU2[0]*inv_n - (SumTotalU[0]*SumTotalU[0])*inv_n2;
	Statistics.VarTotalU_NkT[0] = VarTotalUTemp*inv_NkT*inv_NkT;
	Statistics.VarTotalU_SI[0] = VarTotalUTemp*kJ_mol*kJ_mol;
/* Gaussian error propagation of
 *
 * var(U) = <U^2> - <U>^2
 * => var(var(U)) = var(<U^2>) + (2*<U>)^2*var(<U>) = 1/N*var(U^2) + 1/N*(2*<U>)^2*var(U)
 *
 * --AT
 */
	double varU2 = (SumTotalU2[1]*inv_n - (SumTotalU2[0]*SumTotalU2[0])*inv_n2);
	StdUTemp = sqrt(varU2*inv_n2 + qqpwr(2.0*SumTotalU[0]*inv_n,2)*VarTotalUTemp*inv_n2); // use stddevs of mean with 1/n for variance
	Statistics.VarTotalU_NkT[1] = StdUTemp*inv_NkT*inv_NkT;
	Statistics.VarTotalU_SI[1] = StdUTemp*kJ_mol*kJ_mol;
	
/* Heat capacity from variance of energy. Note that for rigid-body calculations, vibrations are neglected entirely.
 * Start with component from intermolecular potentials
 *
 * Units matching:
 * U is in unit of picoergs, varU in pErg^2, kB in pErg/K, and T in K
 * => C is in units of pErg^2/(pErg/K*K*K) = pErg/K
 */
	double C = VarTotalUTemp/(kB*configuration->T*configuration->T); // contribution from potential energy variations
	// Add rotational and translational DOFs (could be precalculated in configuration but since it is only needed once here for the moment ...)
	// -- we are adding units of kB (pErg/K) to C, so everything's alright
	for(unsigned int i=0; i<configuration->num_element_types; i++) C += configuration->element_types[i]->number*configuration->element_types[i]->dof/2.0*kB; //elements
	for(unsigned int i=0; i<configuration->num_groups; i++) C += configuration->groups[i]->number*configuration->groups[i]->Type->dof/2.0*kB; //groups
	// normalize to one entity
	// -- C is in pErg/K, kJ_mol converts pErg->kJ/mol
	Statistics.Heat_Capacity[0] = C*kJ_mol*1000.0; // units then are J/(mol*K), which is what we want indeed
	// StdUTemp (see above) is in pErg^2
	// => pErg^2/(pErg/K*K*K) = pErg/K, using kJ_mol this becomes kJ/(mol*K)
	Statistics.Heat_Capacity[1] = StdUTemp/(kB*configuration->T*configuration->T)*kJ_mol*1000.0;
	
	//Volume and concentration
	Statistics.Volume[0] = SumV[0]*inv_n;
	vartemp=SumV[1]*inv_n -SumV[0]*SumV[0]*inv_n2;
	if(vartemp<EPS) vartemp=0.0;
	Statistics.Volume[1] = sqrt(vartemp);
	Statistics.Number_Density[0] = (double)(configuration->N)/Statistics.Volume[0]*1E24; //In particles/cc
	Statistics.Number_Density[1] = Statistics.Volume[1]*(double)(configuration->N)/(Statistics.Volume[0]*Statistics.Volume[0])*1E24*inv_sqrt_n;
	Statistics.Molarity[0] = Statistics.Number_Density[0]*1E3/NA; //In mol/L
	Statistics.Molarity[1] = Statistics.Number_Density[1]*1E3/NA;
	
#if CALC_ENTROPY_dG>0
	Statistics.S[0]+=kB*perg_to_kJ_per_mol*1000*((double)(!configuration->NpT)+46.6826192+log(Statistics.Mass*configuration->T/Statistics.Number_Density[0]))+Statistics.Heat_Capacity[0];
#endif
	
	//Dipole moment
	cout << "\tCalculating average electrostatic properties...\n";
	Statistics.ME[0] = SumME[0]*inv_n;
	Statistics.ME[1] = SumME[1]*inv_n-Statistics.ME[0]*Statistics.ME[0];
	Statistics.MeanM = SumM[0]*inv_n;
	Vec3Amult(SumM[0],SumM[0],MMTemp);
	StdMTemp = (SumM[1]*inv_n - MMTemp*inv_n2)*inv_n; //multiply by 1/n here since sqrt done elementwise on next line
	Statistics.StdM.vec[0] = sqrt(StdMTemp.vec[0]); Statistics.StdM.vec[1] = sqrt(StdMTemp.vec[1]); Statistics.StdM.vec[2] = sqrt(StdMTemp.vec[2]);
	Statistics.M2[0] = SumScalarM2[0]*inv_n;
	Statistics.M2[1] = sqrt(SumScalarM2[1]*inv_n - SumScalarM2[0]*SumScalarM2[0]*inv_n2)*inv_sqrt_n;
	Statistics.VarM[0] = Statistics.M2[0] - (Statistics.MeanM*Statistics.MeanM);
/* 
 * I would think of it more in terms of Gaussian error propagation (GEP):
 *
 * var(M) = <M^2> - <M>^2
 * => stddev(var(M))=sqrt((d(var(M))/d(<M^2>)*sigma(M^2))^2 + (d(var(M))/d(<M>)*sigma(M))^2)
 * => stddev(var(M))=sqrt(sigma(M^2)^2 + (2*<M>*sigma(<M>))^2)
 *
 * --AT
 */
	vartemp=qqpwr(Statistics.M2[1],2)+qqpwr(2.0*(Statistics.MeanM*Statistics.StdM),2);
	if(vartemp<EPS) vartemp=0.0;
	Statistics.VarM[1] = sqrt(vartemp);
	Statistics.NormM[0] = SumNormM*inv_n;
	Statistics.NormM[1] = sqrt(SumScalarM2[0]*inv_n - Statistics.NormM[0]*Statistics.NormM[0]);
	Statistics.CovarM = SumM2*inv_n - (Statistics.MeanM.V3TensProd(Statistics.MeanM));
	
	//Dielectric constant and fields
	average_dielectric(Statistics,configuration,edirection);
	
	//Order parameters
	cout << "\tCalculating average order parameters...\n";
	Statistics.MeanCos = SumCosn[0]*inv_n;
	Vec3Amult(SumCosn[0],SumCosn[0],CosnTemp);
	CosnTemp = (SumCosn[1]*inv_n - CosnTemp*inv_n2)*inv_n; //Use temporary vector since no sqrt function for vec3, mult by inv_n since inside sqrt
	Statistics.StdCos.vec[0] = sqrt(CosnTemp.vec[0]); Statistics.StdCos.vec[1] = sqrt(CosnTemp.vec[1]); Statistics.StdCos.vec[2] = sqrt(CosnTemp.vec[2]);
	Statistics.Var_Skew[0] = Statistics.MeanCos.vec[1] - Statistics.MeanCos.vec[0]*Statistics.MeanCos.vec[0]; //Var(cos) -- overall, not block average
	Statistics.Var_Skew[1] = Statistics.MeanCos.vec[2]-3.0*Statistics.MeanCos.vec[0]*Statistics.MeanCos.vec[1]+2.0*Statistics.MeanCos.vec[0]*Statistics.MeanCos.vec[0]*Statistics.MeanCos.vec[0]; //Skew(cos) -- overall, not block average
	Statistics.P_n_theta.vec[0] = Statistics.MeanCos.vec[0]; //<P_1(theta)> = <cos>
	Statistics.P_n_theta.vec[1] = (3.0*Statistics.MeanCos.vec[1] - 1.0)/2.0; //<P_2(theta)>, not LC <P_2>
	Statistics.P_n_theta.vec[2] = (5.0*Statistics.MeanCos.vec[2]-3.0*Statistics.MeanCos.vec[0])/2.0;
	Statistics.Std_Pn.vec[0] = Statistics.StdCos.vec[0]; //error same as <cos>
	Statistics.Std_Pn.vec[1] = 1.5*Statistics.StdCos.vec[1];
	Statistics.Std_Pn.vec[2] = sqrt(2.5*Statistics.StdCos.vec[2]*Statistics.StdCos.vec[2] + 1.5*Statistics.StdCos.vec[0]*Statistics.StdCos.vec[0]);
	
	//Acceptance statistics
	cout << "\tCalculating acceptance statistics...\n";
	Statistics.Acceptance=Vec3(0.0);
	if(tries[0]>0) Statistics.Acceptance.vec[0] = 100.0*(double)accepted[0]/(double)tries[0];
	if(tries[1]>0) Statistics.Acceptance.vec[1] = 100.0*(double)accepted[1]/(double)tries[1];
	if(tries[2]>0) Statistics.Acceptance.vec[2] = 100.0*(double)accepted[2]/(double)tries[2];
	Statistics.RMS_Move[0] = sqrt(msmoved[0]/(double)configuration->steps);
	Statistics.RMS_Move[1] = sqrt(msmoved[1]/(double)configuration->steps);
	
	// Output running averages
	if(Statistics.start_step<=configuration->randsteps){
		cout << "\tCalculating running averages for each block of +/- " << Statistics.running_average << " blocks...\n";
		ofstream statfn;
		cout << "\t\t-> Outputting to " << Statistics.stats_fn << ".dat\n";
		statfn.open((Statistics.stats_fn+".dat").c_str());
		if(statfn.fail()) {
			cout << "ERROR: Unable to open output file.\n";
			exit(1);
		}
		statfn.setf(ios::fixed,ios::floatfield);
		statfn.precision(6);
		statfn << "# Running average of +/- " << Statistics.running_average*Statistics.block_size << " steps\n";
		statfn << "#Step\tVolume\tStdErr\tMeanM_x\tMeanM_y\tMeanM_z\tStdErr_x\tStdErr_y\tStdErr_z\tM^2\tStdErr\tVarM\tStdErr\tEpsilon\tStdErr\tEpsilon_RF\tStdErr\tHeat capacity\tStdErr\tU_total\tStdErr";
		for(unsigned int i=0; i<NUM_V_STORE; i++) statfn << "\t" << Statistics.potentials[i] << "\tStdErr";
		statfn << "\tcos\tcos^2\tcos^3\tStdErr_cos\tStdErr_cos^2\tStdErr_cos^3\n";
#if DEBUG_LEVEL>1
		unsigned int avgnr=0;
		double eps_avg=0.0; double eps_s_avg=0.0;
#endif
		for(unsigned int i=Statistics.running_average; i<Statistics.num_blocks-Statistics.running_average; i++){
			stats_t BlockStatistics;
			BlockStatistics.num_blocks=Statistics.running_average*2+1;
			inv_n = 1.0/(double)(BlockStatistics.num_blocks);
			inv_n2 = inv_n*inv_n;
			inv_sqrt_n = sqrt(inv_n);
			for(unsigned int k=0; k<NUM_V_STORE; k++){
				SumU[0][k]=0.0;
				SumU[1][k]=0.0;
			}
			SumTotalU[0] = 0.0; SumTotalU[1] = 0.0;
			SumTotalU2[0] = 0.0; SumTotalU2[1] = 0.0;
			SumV[0] = 0.0; SumV[1] = 0.0;
			SumM[0].V3Zeros(); SumM[1].V3Zeros();
			SumScalarM2[0] = 0.0; SumScalarM2[1] = 0.0;
			SumCosn[0].V3Zeros(); SumCosn[1].V3Zeros();
			for(int j=-Statistics.running_average; j<=(int)Statistics.running_average; j++){
				for (unsigned int k = 0; k < NUM_V_STORE; k++){
					SumU[0][k] += BlockU[i+j][k];
					SumU[1][k] += BlockU[i+j][k]*BlockU[i+j][k];
				}
				SumTotalU[0] += BlockTotalU[i+j];
				SumTotalU[1] += BlockTotalU[i+j]*BlockTotalU[i+j];
				SumTotalU2[0] += BlockTotalU2[i+j];
				SumTotalU2[1] += BlockTotalU4[i+j];
				// epsilon depends on Volume (+stderr), StdM, M2 (+stderr), VarM (+stderr), and MeanM, so let's go ahead and calculate these per block
				// volume
				SumV[0] += BlockV[i+j];
				SumV[1] += BlockV[i+j]*BlockV[i+j];
				SumM[0] += BlockM[i+j];
				Vec3Amult(BlockM[i+j], BlockM[i+j], MMTemp); //As a sidenote, the 'scalar math' header really needs cleanup - LEJ
				SumM[1] += MMTemp;
				ScalarM2Temp = BlockM2[i+j].M3Trace(); //Simply <M'*M>, used for dielectric
				SumScalarM2[0] += ScalarM2Temp;
				SumScalarM2[1] += ScalarM2Temp*ScalarM2Temp;
				SumCosn[0] += BlockCosn[i];
				SumCosn[1] += BlockCosn[i]*BlockCosn[i];
			}
			// heat capacity
			VarTotalUTemp = SumTotalU2[0]*inv_n - (SumTotalU[0]*SumTotalU[0])*inv_n2;
			varU2 = SumTotalU2[1]*inv_n - (SumTotalU2[0]*SumTotalU2[0])*inv_n2;
			StdUTemp = sqrt(varU2 + 4.0*qqpwr(SumTotalU[0]*inv_n,2)*VarTotalUTemp)*inv_sqrt_n;
			C = VarTotalUTemp/(kB*configuration->T*configuration->T); // contribution from potential energy variations
			//Add rotational and translational DOFs (could be precalculated in configuration but since it is only needed once here for the moment ...)
			for(unsigned int k=0; k<configuration->num_element_types; k++) C += configuration->element_types[k]->number*configuration->element_types[k]->dof/2.0*kB; //elements
			for(unsigned int k=0; k<configuration->num_groups; k++) C += configuration->groups[k]->number*configuration->groups[k]->Type->dof/2.0*kB; //groups
			// normalize to one entity
			BlockStatistics.Heat_Capacity[0] = C*kJ_mol*1000.0;
			BlockStatistics.Heat_Capacity[1] = StdUTemp/(kB*configuration->T*configuration->T)*kJ_mol*1000.0;
			// volume
			BlockStatistics.Volume[0] = SumV[0]*inv_n;
			vartemp=SumV[1]*inv_n - SumV[0]*SumV[0]*inv_n2;
			if(vartemp<EPS) vartemp=0.0;
			BlockStatistics.Volume[1] = sqrt(vartemp);
			// average M
			BlockStatistics.MeanM = SumM[0]*inv_n;
			Vec3Amult(SumM[0],SumM[0],MMTemp);
			// stdM
			StdMTemp = (SumM[1]*inv_n - MMTemp*inv_n2)*inv_n; //multiply by 1/n here since sqrt done elementwise on next line
			BlockStatistics.StdM.vec[0] = sqrt(StdMTemp.vec[0]); BlockStatistics.StdM.vec[1] = sqrt(StdMTemp.vec[1]); BlockStatistics.StdM.vec[2] = sqrt(StdMTemp.vec[2]);
			// M^2
			BlockStatistics.M2[0] = SumScalarM2[0]*inv_n;
			BlockStatistics.M2[1] = sqrt(SumScalarM2[1]*inv_n - SumScalarM2[0]*SumScalarM2[0]*inv_n2)*inv_sqrt_n;
			// VarM
			BlockStatistics.VarM[0] = (BlockStatistics.M2[0] - (BlockStatistics.MeanM*BlockStatistics.MeanM));
			vartemp=BlockStatistics.M2[1]+qqpwr(2.0*(BlockStatistics.MeanM*BlockStatistics.StdM),2);
			if(vartemp<EPS) vartemp=0.0;
			BlockStatistics.VarM[1] = sqrt(vartemp);
			// calculate dielectric constant
			average_dielectric(BlockStatistics,configuration,edirection);
#if DEBUG_LEVEL>1
			eps_avg+=BlockStatistics.Epsilon[0];
			eps_s_avg+=BlockStatistics.Epsilon_s[0];
			avgnr++;
#endif
			BlockStatistics.MeanCos = SumCosn[0]*inv_n;
			Vec3Amult(SumCosn[0],SumCosn[0],MMTemp);
			// stdM
			StdMTemp = (SumCosn[1]*inv_n - MMTemp*inv_n2)*inv_n; //multiply by 1/n here since sqrt done elementwise on next line
			BlockStatistics.StdCos.vec[0] = sqrt(StdMTemp.vec[0]); BlockStatistics.StdCos.vec[1] = sqrt(StdMTemp.vec[1]); BlockStatistics.StdCos.vec[2] = sqrt(StdMTemp.vec[2]);
			statfn << Statistics.start_step+(i+1)*Statistics.block_size << "\t" << BlockStatistics.Volume[0] << "\t" << BlockStatistics.Volume[1] << "\t" << BlockStatistics.MeanM.V3Str('\t') << "\t" << BlockStatistics.StdM.V3Str('\t') << "\t" << BlockStatistics.M2[0] << "\t" << BlockStatistics.M2[1] << "\t" << BlockStatistics.VarM[0] << "\t" << BlockStatistics.VarM[1] << "\t" << BlockStatistics.Epsilon[0] << "\t" << BlockStatistics.Epsilon[1] << "\t" << BlockStatistics.Epsilon_s[0] << "\t" << BlockStatistics.Epsilon_s[1] << "\t" << BlockStatistics.Heat_Capacity[0] << "\t" << BlockStatistics.Heat_Capacity[1] << "\t" << SumTotalU[0]*inv_n << "\t" << sqrt(SumTotalU[1]*inv_n-SumTotalU[0]*SumTotalU[0]*inv_n2)*inv_sqrt_n;
			for(unsigned int k=0; k<NUM_V_STORE; k++) statfn << "\t" << SumU[0][k]*inv_n << "\t" << sqrt(SumU[1][k]*inv_n-SumU[0][k]*SumU[0][k]*inv_n2)*inv_sqrt_n;
			statfn << "\t" << BlockStatistics.MeanCos.V3Str('\t') << "\t" << BlockStatistics.StdCos.V3Str('\t') << "\n";
		}
		statfn.close();
		// Calculate total dipole moment autocorrelation function,
		// uses R(k) = 1/((n-k)*var(x))*sum_i=1^(n-k)(x_i-x_avg)(x_(i+k)-x_avg)
		if(calc_dipole_autocorr){
			cout << "\tCalculating total dipole moment autocorrelation function (may take a while)\n";
			cout << "\t\t-> Progress: 0%";
			cout.flush();
			unsigned int percentage=0;
			double* R=new double[Statistics.nr_steps];
			double count=0.0;
//			double invmaxcount=8.0/((double)Statistics.nr_steps*(3.0*Statistics.nr_steps-2.0)); // sum_{i=1}^{N/2} N-i = N*(3*N-2)/8
			double invmaxcount=2.0/((double)Statistics.nr_steps*((double)Statistics.nr_steps-1.0))/AUTOCORR_CUTOFF; // sum_{i=1}^N N-i = N*(N-1)/2
			for(unsigned int k=0; k<Statistics.nr_steps*AUTOCORR_CUTOFF; k++){
				Vec3 xbar_i(0.0);
				Vec3 xbar_ik(0.0);
				double xbar2_i=0.0;
				double xbar2_ik=0.0;
				unsigned int Nk=Statistics.nr_steps-k;
				for(unsigned int i=0; i<Nk; i++){
					Vec3 x_i=total_dipole[i];
					Vec3 x_ik=total_dipole[i+k];
					xbar_i+=x_i;
					xbar2_i+=x_i*x_i;
					xbar_ik+=x_ik;
					xbar2_ik+=x_ik*x_ik;
				}
				xbar_i/=Nk;
				xbar_ik/=Nk;
				double xvar_i=xbar2_i/Nk-xbar_i*xbar_i;
				double xvar_ik=xbar2_ik/Nk-xbar_ik*xbar_ik;
				R[k]=0.0;
				for(unsigned int i=0; i<Nk; i++){
					Vec3 x_i=total_dipole[i];
					Vec3 x_ik=total_dipole[i+k];
					R[k]+=(x_i-xbar_i)*(x_ik-xbar_ik);
					count+=1.0;
				}
//				R[k]=(x_iik-(xbar*x_i_ik)+(double)(Statistics.nr_steps-k)*(xbar*xbar))/((double)Statistics.nr_steps*xvar);
//				R[k]=(x_iik-(xbar*x_i_ik))/(Nk*xvar)+(xbar*xbar)/xvar;
				R[k]/=Nk*sqrt(xvar_i*xvar_ik);
				
				double percent=100.0*count*invmaxcount;
				if(percent>=percentage+5){
					percentage=(unsigned int)floor(percent/5.0)*5;
					if(percentage%20==0) cout << percentage << "%"; else cout << ".";
					cout.flush();
				}
			}
			cout << "\n\t\t-> Outputting to " << Statistics.stats_fn << ".autocorr.dat\n";
			ofstream autocorr;
			autocorr.open((Statistics.stats_fn+".autocorr.dat").c_str());
			if(autocorr.fail()) {
				cout << "ERROR: Unable to open output file.\n";
				exit(1);
			}
			autocorr.setf(ios::fixed,ios::floatfield);
			autocorr.precision(6);
			for(unsigned int k=0; k<Statistics.nr_steps*AUTOCORR_CUTOFF; k++) autocorr << k << "\t" << R[k] << "\n";
			autocorr.close();
			delete[] R;
		}
#if DEBUG_LEVEL>2
		eps_avg/=avgnr;
		eps_s_avg/=avgnr;
		cout << "<eps> = " << eps_avg << "\n";
		cout << "<eps_s> = " << eps_s_avg << "\n";
#endif
	} else print_statistics(Statistics,configuration);
	
	//Clean up
	delete[] BlockU;
	delete[] BlockTotalU;
	delete[] BlockTotalU2;
	delete[] BlockTotalU4;
	delete[] BlockV;
	delete[] BlockM;
	delete[] BlockM2;
	delete[] BlockCosn;
	
	cout << "<- Done.\n";
}


///Entry point for stats code
int main (int argc, char * const argv[]) {
	
	cout << "\nRobinson Group statistics analysis tool for MCfig\n";
	cout << "Compiled " << __DATE__ << " (Build " << BUILD << ", " << VERSION << ")\n";
	
	MC_Config config;
	Config_Data* configuration = &config.parameters;
	
	string conffile="";
	
	// Check to make sure there are enough command line parameters
	if(argc<2){
		cout << "Syntax: " << argv[0] << " <trajectory file> <optional: initial step> <optional w/ initial step: ending step>\n\n";
		cout << "Note: If initial step is within the randomization steps raw statistics is outputted to a .stat.dat file.\n\n";
		exit(1);
	} else{
		//Trajectory file name
		conffile=argv[1];
	}
	
	//Load configuration from file
	cout << "-> Reading configuration file...\n";
#ifndef USE_CMWC4096
	configuration->idum = new __int32_t; // so we don't get segfaults ...
#endif
	config.GetFromFile(conffile.c_str());
	phys_configuration(configuration);
	
	//Check starting step to make sure it is valid
	unsigned int start_step;
	bool calc_dipole_autocorr=true;
	int id;
	if(argc>2){ // second parameter: requested starting step
		if(compare_strings(argv[2],"raw") || compare_strings(argv[2],"all")){
			start_step=configuration->randsteps;
			if(start_step>=configuration->last_step){
				cout << "ERROR: Not enough data after initial randomization.\n";
				exit(2);
			}
		} else{
			bool fail=from_string(id,argv[2]);
			if(!fail) if(id<0) fail=true;
			if(!fail){
				start_step=(unsigned int)id;
			} else{
				cout << "ERROR: Second command line parameter (starting step) needs to be an integer greater or equal to zero.\n";
				exit(1);
			}
		}
		if(argc>3){
			if(compare_strings(argv[3],"no_dipole_autocorr")) calc_dipole_autocorr=false;
			if(argc>4){
				if(compare_strings(argv[4],"no_dipole_autocorr")) calc_dipole_autocorr=false;
			}
		}
	} else{
		if(configuration->last_step>configuration->laststep){
			start_step = configuration->last_step-configuration->laststep;
		} else{
			start_step=configuration->randsteps;
			if(start_step>=configuration->last_step){
				cout << "ERROR: Not enough data after initial randomization.\n";
				exit(2);
			}
		}
	}
	
	//Check ending step to make sure it is valid
	unsigned int end_step = configuration->last_step;
	if(argc>3){ // third parameter: requested ending step
		bool fail=from_string(id,argv[3]);
		if(!fail) if(id<0) fail=true;
		if(!fail){
			end_step=(unsigned int)id;
		} else{
			cout << "ERROR: Third command line parameter (ending step) needs to be an integer greater or equal to zero.\n";
			exit(1);
		}
		if(end_step < start_step){
			cout << "ERROR: Ending step cannot be before starting step.\n";
			exit(1);
		} else{
			if(end_step > configuration->last_step) {
				cout << "ERROR: Requested ending step is beyond last step in trajectory file.\n";
				exit(1);
			}
		}
	}
	
	//Allocate memory for reading statistics
	Vstore* energy_components = new Vstore[end_step-start_step];
	Vec3* total_dipole = new Vec3[end_step-start_step];
	Vec3* cos_n = new Vec3[end_step-start_step];
	double* V = new double[end_step-start_step];
	double msmoved[2];
	unsigned int accepted[3], tries[3];
	
	//Read in and process statistics data
	if(config.GetStatistics(start_step,end_step, energy_components, total_dipole, cos_n, V, msmoved, accepted, tries)){
		cout << "-> Calculating statistics from step " << start_step << " to " << end_step << "...\n";
		stats_t Statistics;
		Statistics.NpT=configuration->NpT;
		Statistics.start_step = start_step;
		Statistics.end_step = end_step;
		Statistics.nr_steps = end_step-start_step;
		Statistics.potentials[0] = "UPol"; Statistics.potentials[1] = "UES"; Statistics.potentials[2] = "ULJ"; Statistics.potentials[3] = "UGrp"; //THIS SHOULD GO IN CONFIG.H!!!
		Statistics.Mass=configuration->Mass;
		
		cout << "\t-> Initial step is ";
		if(Statistics.start_step>configuration->randsteps) cout << "not ";
		cout << "within the first " << configuration->randsteps << " randomization steps.\n";
		
		//Set file name and block size
		Statistics.stats_fn = name_statfile(configuration->trajectoryfile);
		if(configuration->correlation_length!=0) Statistics.block_size = configuration->correlation_length; else Statistics.block_size = configuration->grfreq;
		Statistics.running_average = configuration->running_average;
		if((Statistics.end_step-Statistics.start_step)%Statistics.block_size!=0){
			cout << "WARNING: Unable to evenly divide trajectory into blocks of " << Statistics.block_size << ". Truncating data (not using last " << (Statistics.end_step-Statistics.start_step)%Statistics.block_size << " steps).\n";
			Statistics.nr_steps=(Statistics.nr_steps/Statistics.block_size)*Statistics.block_size;
		}
		Statistics.num_blocks = Statistics.nr_steps/Statistics.block_size;
		
		if(Statistics.num_blocks < 1){
			cout << "Insufficient data to calculate average.\n";
			exit(1);
		}
		
		Vec3 edirection=Vec3(0.0,0.0,1.0);
		if(configuration->Efield.V3Norm()>EPS) edirection=configuration->Efield;
		if((configuration->rotate_Efield_steps>0) && configuration->use_trajectory){
			if(configuration->last_step>=configuration->rotate_Efield_steps){
				Vstore potentials;
				Vec3 cosmeans;
				double Vs, msmoved[2];
				unsigned int accepted[3], tries[3];
				cout << "Electric field was rotated into total dipole moment direction at " << configuration->rotate_Efield_steps << " steps, adjusting...\n";
				config.GetStatistics(configuration->rotate_Efield_steps,configuration->rotate_Efield_steps+1,&potentials,&edirection,&cosmeans,&Vs,msmoved,accepted,tries);
				cout << "<- Done.\n";
			}
		}
		edirection/=edirection.V3Norm();
		calc_means(configuration, Statistics, energy_components, total_dipole, cos_n, V, msmoved, accepted, tries, calc_dipole_autocorr, edirection);
	} else cout << "ERROR: Cannot read statistics.\n";
	
	//Clean up
	delete[] energy_components;
	delete[] total_dipole;
	delete[] cos_n;
	delete[] V;
	
	return 0;
}
