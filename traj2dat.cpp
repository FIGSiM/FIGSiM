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

#include "traj2dat.h"
#include "x3d_functions.h"

inline void correlation(cfparams_t &cf_params, unsigned int &analysis_nr)
{
	bool has_gammaalpha=false;
	bool has_sdf=false;
	bool has_Q=false;
	bool is_single=false;
	unsigned int sdf_nr=cf_params.nr_correlations;
	unsigned int any_nr_output=22; // number of output fields for Q-tensor
	cout << "\t\t\t-> Calculating ";
	for(unsigned int i=0; i<cf_params.nr_correlations; i++){
		cout << output_corr_types[cf_params.types[i]];
		if(cf_params.typenr[i]>1) cout << "^" << cf_params.typenr[i];
		if((cf_params.types[i]==gr) || (cf_params.types[i]==gmu)) cout << "(r)";
		if(cf_params.types[i]==gammaalpha) has_gammaalpha=true;
		if(cf_params.types[i]==Q) has_Q=true;
		if(cf_params.types[i]==sdf){
			if(sdf_nr<cf_params.nr_correlations){
				cout << "\nERROR: Only one spatial distribution function (SDF) allowed per analysis section\n";
				exit(1);
			}
			has_sdf=true;
			sdf_nr=i;
		}
		if(cf_params.types[i]==alphabetagamma){
			if(!is_single && (i>0)){
				cout << "\nERROR: Single element analyses cannot be mixed with multiple element analyses.\n";
				exit(2);
			}
			is_single=true;
		} else{
			if(is_single){
				cout << "\nERROR: Single element analyses cannot be mixed with multiple element analyses.\n";
				exit(2);
			}
		}
		if(cf_params.nr_correlations>1){
			if(i+1<cf_params.nr_correlations) cout << ", ";
			if(i+2==cf_params.nr_correlations) cout << "and ";
		}
	}
	cout << " ...\n";
	unsigned int nr_results_total=0;
	double rmax=0.0;
	double averageV=0.0;
	for(unsigned int frame=0; frame<cf_params.n_frames; frame++) averageV+=cf_params.V[frame+cf_params.startframe-cf_params.start_frame_in_data];
	averageV/=cf_params.n_frames;
	if(!has_Q && !is_single){ // for the moment, limit ourselves to use rmax only for (pair-)correlations (no property-type calculations, i.e. Q-tensor)
		// find rmax first
		cout << "\t\t\t\t-> Determining rmax: ";
		if(cf_params.rmax<0.0){
			double sbl=cf_params.configuration->boxlength[0];
			if(sbl>cf_params.configuration->boxlength[1]) sbl=cf_params.configuration->boxlength[1];
			if(sbl>cf_params.configuration->boxlength[2]) sbl=cf_params.configuration->boxlength[2];
			rmax=sbl;
			for(unsigned int frame=0; frame<cf_params.n_frames; frame++){
				double newV=cf_params.V[frame+cf_params.startframe-cf_params.start_frame_in_data];
				// adjust volume for boundary conditions
				if(cf_params.configuration->LJwall_calc && cf_params.configuration->LJwall_fixed){
					double newXm=cf_params.Xm[frame+cf_params.startframe-cf_params.start_frame_in_data];
					double scale=sqrt((newV*cf_params.configuration->boxlength[0])/(cf_params.configuration->V*2.0*newXm));
					cf_params.configuration->LJwall_xm=newXm;
					cf_params.configuration->boxlength[0]=2.0*newXm;
					cf_params.configuration->boxlength[1]*=scale; cf_params.configuration->boxlength[2]*=scale;
					cf_params.configuration->nndist*=scale;
					cf_params.configuration->V=cf_params.configuration->boxlength[0]*cf_params.configuration->boxlength[1]*cf_params.configuration->boxlength[2];
				} else update_volume(cf_params.configuration,newV);
				sbl=cf_params.configuration->boxlength[0];
				if(sbl>cf_params.configuration->boxlength[1]) sbl=cf_params.configuration->boxlength[1];
				if(sbl>cf_params.configuration->boxlength[2]) sbl=cf_params.configuration->boxlength[2];
				if(rmax>sbl) rmax=sbl;
			}
			
			rmax/=2.0; // half the boxsize
		} else rmax=cf_params.rmax;
		cout << rmax << " Angström\n";
		if(rmax<EPS){
			cout << "ERROR: rmax should not be zero.\n";
			exit(3);
		}
	}
	for(unsigned int i=0; i<cf_params.items_involved; i++){
		if(!GetPosDipole(cf_params,i,0,cf_params.positions[i],cf_params.vectors[i],cf_params.max_vectors[i],cf_params.rots[i],cf_params.total_nr[i])){
			cout << "\nERROR: Could not obtain position and dipole from trajectory file.\n";
			exit(5);
		}
	}
	if(!has_Q && !is_single) cout << "\t\t\t\t-> Calculating histograms (" << cf_params.nr_bins << " bins): 0%"; else cout << "\t\t\t\t-> Calculating values: 0%";
	cout.flush();
	double dr=rmax/(double)cf_params.nr_bins;
	if((dr<EPS) && !has_Q && !is_single){
		cout << "\nERROR: Bin width is zero. Impossible to calculate histograms this way.\n";
		exit(6);
	}
	unsigned int ga_bins=(unsigned int)round(sqrt(4.0*cf_params.nr_bins));
	double ga_dr=rmax/(double)ga_bins;
	unsigned int a=0;
	if((cf_params.items_involved==2) && !cf_params.any_number_of_items){
		if(cf_params.total_nr[0]<cf_params.total_nr[1]) a=1;
	}
	unsigned int percentage=0;
	unsigned int results_per_frame=cf_params.total_nr[a];
	if(cf_params.any_number_of_items && !has_Q) results_per_frame=cf_params.total_nr[1];
	bool exclude_idx=false;
	unsigned int multiplier=1;
	if(cf_params.items_involved==2){
		exclude_idx=((cf_params.items[0].group_idx==cf_params.items[1].group_idx) || ((cf_params.items[0].is_element && cf_params.items[1].is_element) && (cf_params.items[0].element_idx==cf_params.items[1].element_idx)));
		multiplier=cf_params.total_nr[(a+1)%2]-(unsigned int)exclude_idx;
		if(multiplier==0) multiplier=1; // safety net
	}
	if(!cf_params.group_internal) results_per_frame*=multiplier;
	cf_params.results=new double**[cf_params.n_frames/cf_params.average_frames+1];
	unsigned int** gammaalpha_nr=new unsigned int*[cf_params.n_frames/cf_params.average_frames+1];
	for(unsigned int i=0; i<cf_params.n_frames/cf_params.average_frames+1; i++) gammaalpha_nr[i]=new unsigned int[ga_bins*ga_bins];
	CVec3 cev;
	unsigned int mu_nr, ii;
	Vec3 distance, u, musum, usum, eigenvalues, P1s, cos3s, avec;
	Vec3* max_dipole = new Vec3[cf_params.items_involved];
	Mat33 Qtensor, half_I, eigenvectors;
	half_I.M3Eye();
	half_I/=2.0;
	double r,x,cosg,cosa,theta,phi,mu_max, cost, p;
	
	unsigned int x2, xn;
	double* average=new double[cf_params.nr_correlations];
	double* average2=new double[cf_params.nr_correlations];
	Vec3* Vec3average=new Vec3[cf_params.nr_correlations];
	Vec3* Vec3average2=new Vec3[cf_params.nr_correlations];
	unsigned int* Qorder=new unsigned int[3*cf_params.nr_correlations];
	for(unsigned int i=0; i<cf_params.nr_correlations; i++){
		average[i]=0.0;
		average2[i]=0.0;
		Vec3average[i]=Vec3(0.0);
		Vec3average2[i]=Vec3(0.0);
		Qorder[3*i]=0; Qorder[3*i+1]=0; Qorder[3*i+2]=0;
	}
	for(unsigned int frame=0; frame<=cf_params.n_frames; frame++){
		for(unsigned int i=0; i<cf_params.items_involved; i++) max_dipole[i]=Vec3(0.0);
		if(frame%cf_params.average_frames==0){
			cf_params.results[frame/cf_params.average_frames]=new double*[cf_params.nr_correlations];
			for(unsigned int i=0; i<cf_params.nr_correlations; i++){
				switch(cf_params.types[i]){
					case sdf:
						cf_params.results[frame/cf_params.average_frames][i]=new double[ga_bins*ga_bins*cf_params.nr_bins+2*ga_bins*ga_bins]; // holds voxel sdf+additional space for most likely positions
						break;
					case gammaalpha:
						cf_params.results[frame/cf_params.average_frames][i]=new double[2*ga_bins*ga_bins];
						break;
					case Q:
						cf_params.results[frame/cf_params.average_frames][i]=new double[2*any_nr_output];
						break;
					case alphabetagamma:
						cf_params.results[frame/cf_params.average_frames][i]=new double[3*cf_params.nr_bins];
						break;
					default:
						cf_params.results[frame/cf_params.average_frames][i]=new double[cf_params.nr_bins];
						break;
				}
				for(unsigned int j=0; j<ga_bins*ga_bins; j++){
					switch(cf_params.types[i]){
						case sdf:
							for(unsigned int k=0; k<cf_params.nr_bins; k++) cf_params.results[frame/cf_params.average_frames][i][j*cf_params.nr_bins+k]=0.0;
							cf_params.results[frame/cf_params.average_frames][i][ga_bins*ga_bins*cf_params.nr_bins+2*j]=0.0;
							cf_params.results[frame/cf_params.average_frames][i][ga_bins*ga_bins*cf_params.nr_bins+2*j+1]=0.0;
							break;
						case gammaalpha:
							gammaalpha_nr[frame/cf_params.average_frames][j]=0;
							cf_params.results[frame/cf_params.average_frames][i][j]=0.0;
							cf_params.results[frame/cf_params.average_frames][i][j+ga_bins*ga_bins]=0.0; // used for holding standard deviation
							break;
						case Q:
							if(j<2*any_nr_output){
								cf_params.results[frame/cf_params.average_frames][i][j]=0.0;
							}
							break;
						case alphabetagamma:
							if(j<3*cf_params.nr_bins){
								cf_params.results[frame/cf_params.average_frames][i][j]=0.0;
							}
							break;
						default:
							if(j<cf_params.nr_bins){
								cf_params.results[frame/cf_params.average_frames][i][j]=0.0;
							}
							break;
					}
				}
			}
		}
		unsigned int result_nr=0;
		Qtensor.M3Zeros(); // Q-tensor is only valid per frame (one system state)
		musum=Vec3(0.0);
		usum=musum; // = 0
		// adjust volume for boundary conditions
		double newV=cf_params.V[frame+cf_params.startframe-cf_params.start_frame_in_data];
		if(cf_params.configuration->LJwall_calc && cf_params.configuration->LJwall_fixed){
			double newXm=cf_params.Xm[frame+cf_params.startframe-cf_params.start_frame_in_data];
			double scale=sqrt((newV*cf_params.configuration->boxlength[0])/(cf_params.configuration->V*2.0*newXm));
			cf_params.configuration->LJwall_xm=newXm;
			cf_params.configuration->boxlength[0]=2.0*newXm;
			cf_params.configuration->boxlength[1]*=scale; cf_params.configuration->boxlength[2]*=scale;
			cf_params.configuration->nndist*=scale;
			cf_params.configuration->V=cf_params.configuration->boxlength[0]*cf_params.configuration->boxlength[1]*cf_params.configuration->boxlength[2];
		} else update_volume(cf_params.configuration,newV);
		
		for(unsigned int i=0; i<cf_params.items_involved; i++){
			if(!GetPosDipole(cf_params,i,frame+cf_params.startframe-cf_params.start_frame_in_data,cf_params.positions[i],cf_params.vectors[i],cf_params.max_vectors[i],cf_params.rots[i],cf_params.total_nr[i])){
				cout << "\nERROR: Could not obtain position and dipole from trajectory file.\n";
				exit(5);
			}
		}
		for(unsigned int i=0; i<cf_params.total_nr[a]; i++){ // go over each individual entity (group or element) in system
			if(cf_params.group_internal || has_Q || is_single){ // or single-element case (need only go once over each element)
				// flag down analysis with any number of elements (i.e. Q-tensor) here as well
				if((cf_params.total_nr[0]==cf_params.total_nr[1%cf_params.items_involved]) || has_Q || is_single){ // element-in-group case
					distance=cf_params.positions[1%cf_params.items_involved][i]-cf_params.positions[0][i];
					apply_PBCs(distance,cf_params.configuration);
					r=distance.V3Norm();
					x=floor(r/dr);
					unsigned int k;
					if((unsigned int)x>=cf_params.nr_bins) x=cf_params.nr_bins-1;
					for(unsigned int j=0; j<cf_params.nr_correlations; j++){
						switch(cf_params.types[j]){
							default:
							case gr:
							case distances:
								xn=floor(qqpwr(r,cf_params.typenr[j])/(qqpwr(rmax,cf_params.typenr[j])/(double)cf_params.nr_bins));
								if((unsigned int)xn>=cf_params.nr_bins) xn=cf_params.nr_bins-1;
								if(r<rmax) cf_params.results[frame/cf_params.average_frames][j][(unsigned int)xn]+=1.0;
								break;
							case alphabetagamma:
								u=Rot2AlphaBetaGamma(cf_params.rots[0][i]);
								Vec3average[j]+=u;
								Vec3average2[j].vec[0]+=u.vec[0]*u.vec[0]; Vec3average2[j].vec[1]+=u.vec[1]*u.vec[1]; Vec3average2[j].vec[2]+=u.vec[2]*u.vec[2];
								xn=floor((pi+u.vec[0])/(2.0*pi)*(double)cf_params.nr_bins);
								if(xn<0) xn=0;
								if((unsigned int)xn>=cf_params.nr_bins) xn=cf_params.nr_bins-1;
								cf_params.results[frame/cf_params.average_frames][j][(unsigned int)3*xn]+=1.0;
								xn=floor(u.vec[1]/pi*(double)cf_params.nr_bins);
								if(xn<0) xn=0;
								if((unsigned int)xn>=cf_params.nr_bins) xn=cf_params.nr_bins-1;
								cf_params.results[frame/cf_params.average_frames][j][(unsigned int)3*xn+1]+=1.0;
								xn=floor((pi+u.vec[2])/(2.0*pi)*(double)cf_params.nr_bins);
								if(xn<0) xn=0;
								if((unsigned int)xn>=cf_params.nr_bins) xn=cf_params.nr_bins-1;
								cf_params.results[frame/cf_params.average_frames][j][(unsigned int)3*xn+2]+=1.0;
								break;
							case sdf:
								if(r<rmax){
									avec=cf_params.rots[0][i].ColumnVec3(2); // z-direction of source element
									cosa=VectorCos(avec,distance); // [-1;+1]
									x2=(unsigned int)fastfloor((cosa+1.0)/2.0*ga_bins);
									if(x2>=ga_bins) x2=ga_bins-1;
									phi=atan2(distance*cf_params.rots[0][i].ColumnVec3(1),distance*cf_params.rots[0][i].ColumnVec3(0)); // [-pi;+pi]
									xn=(unsigned int)fastfloor((phi/pi+1.0)/2.0*ga_bins);
									if(xn>=ga_bins) xn=ga_bins-1;
									cf_params.results[frame/cf_params.average_frames][j][(x2*ga_bins+xn)*cf_params.nr_bins+(unsigned int)x]+=1.0;
									// determine (and keep up-to-date) maximum likely position data
									if(cf_params.results[frame/cf_params.average_frames][j][(x2*ga_bins+xn)*cf_params.nr_bins+(unsigned int)x]>cf_params.results[frame/cf_params.average_frames][j][ga_bins*ga_bins*cf_params.nr_bins+2*(x2*ga_bins+xn)]){
										cf_params.results[frame/cf_params.average_frames][j][ga_bins*ga_bins*cf_params.nr_bins+2*(x2*ga_bins+xn)]=cf_params.results[frame/cf_params.average_frames][j][(x2*ga_bins+xn)*cf_params.nr_bins+(unsigned int)x];
										cf_params.results[frame/cf_params.average_frames][j][ga_bins*ga_bins*cf_params.nr_bins+2*(x2*ga_bins+xn)+1]=(unsigned int)x;
									}
								}
								break;
							case angle:
								theta=VectorAngle(cf_params.vectors[0][i],cf_params.vectors[1][i]);
								x2=(unsigned int)fastfloor(theta/pi*cf_params.nr_bins);
								if(x2>=cf_params.nr_bins) x2=cf_params.nr_bins-1;
								cf_params.results[frame/cf_params.average_frames][j][x2]+=1.0;
								break;
							case gammaalpha:
								cosa=VectorCos(cf_params.vectors[0][i],distance);
								x2=(unsigned int)fastfloor(0.5*(cosa+1.0)*ga_bins);
								if(x2>=ga_bins) x2=ga_bins-1;
								cosg=VectorCos(cf_params.vectors[0][i],cf_params.vectors[1][i]);
								if(r<rmax){
									gammaalpha_nr[frame/cf_params.average_frames][(unsigned int)(floor(r/ga_dr)*ga_bins+x2)]+=1;
									cf_params.results[frame/cf_params.average_frames][j][(unsigned int)(floor(r/ga_dr)*ga_bins+x2)]+=cosg;
									cf_params.results[frame/cf_params.average_frames][j][(unsigned int)(floor(r/ga_dr)*ga_bins+x2)+ga_bins*ga_bins]+=cosg*cosg;
								}
								break;
							case Q:
								// Q-tensor definitions according to:
								//   - A. Majumdar, A. Zarnescu: "Landau-de Gennes theory of nematic liquid crystals: the Oseen-Frank limit and beyond", OCCAM report 09/10 (2008)
								//   - D. H. Laidlaw, J. Weickert: "Visualization and Processing of Tensor Fields: Advances and Perspectives", Springer (2009), p 216ff
								//   - J. M. Ball, A. Majumdar: "Nematic liquid crystals: from Maier-Saupe to a continuum theory", OCCAM report 09/33 (2009)
								//   - R. J. Low: "Measuring order and biaxiality", Eur. J. Phys. 23, 111 (2002)
								k=0;
								mu_nr=0;
								mu_max=0.0;
								while(k<cf_params.items_involved){ // allow Q-tensor to be comprised of multiple items' dipole moments
									u=cf_params.vectors[k][i];
									max_dipole[k]+=cf_params.max_vectors[k][i];
									musum+=u;
									u/=u.V3Norm();
									usum+=u*cf_params.weights[k];
									Qtensor+=(u.V3TensProd(u)*1.5-half_I)*cf_params.weights[k];
									mu_nr+=cf_params.weights[k]*cf_params.total_nr[k];
									if(i+1==cf_params.total_nr[0]) mu_max+=max_dipole[k].V3Norm();
									k++;
								}
								if(i+1==cf_params.total_nr[0]){
									Qtensor/=mu_nr;
									cev=Qtensor.Eigenvalues();
									eigenvalues=cev.Im(); // use eigenvalues vector to temporarily store imaginary parts
									if(eigenvalues*eigenvalues>EPS*EPS){
										cout << "ERROR: Q-tensor eigenvalues have imaginary parts.\n";
										exit(7);
									}
									eigenvalues=cev.Re();
									eigenvectors=Qtensor.Eigenvectors(eigenvalues,true);
									// sort by magnitude and keep initial order
									if((Qorder[3*j]==0) && (Qorder[3*j+1]==0) && (Qorder[3*j+2]==0)){
#if DEBUG_LEVEL>3
										cout << "\n" << eigenvalues.V3Str(',') << "\n";
#endif
										Qorder[3*j]=0; Qorder[3*j+1]=1; Qorder[3*j+2]=2;
										unsigned int l;
										if(fabs(eigenvalues.vec[Qorder[3*j+1]]*eigenvectors.ColumnVec3(Qorder[3*j+1]).vec[2])<fabs(eigenvalues.vec[Qorder[3*j+2]]*eigenvectors.ColumnVec3(Qorder[3*j+2]).vec[2])){ // if 1<2, switch
//										if(eigenvalues.vec[Qorder[3*j+1]]<eigenvalues.vec[Qorder[3*j+2]]){ // if 1<2, switch
											l=Qorder[3*j+1];
											Qorder[3*j+1]=Qorder[3*j+2];
											Qorder[3*j+2]=l;
#if DEBUG_LEVEL>3
											cout << "1->2\n";
#endif
										}
										if(fabs(eigenvalues.vec[Qorder[3*j]]*eigenvectors.ColumnVec3(Qorder[3*j]).vec[2])<fabs(eigenvalues.vec[Qorder[3*j+1]]*eigenvectors.ColumnVec3(Qorder[3*j+1]).vec[2])){ // if 0<1, switch (now 0 is definitely the largest)
//										if(eigenvalues.vec[Qorder[3*j]]<eigenvalues.vec[Qorder[3*j+1]]){ // if 0<1, switch (now 0 is definitely the largest)
											l=Qorder[3*j];
											Qorder[3*j]=Qorder[3*j+1];
											Qorder[3*j+1]=l;
#if DEBUG_LEVEL>3
											cout << "0->1\n";
#endif
										}
										if(fabs(eigenvalues.vec[Qorder[3*j+1]]*eigenvectors.ColumnVec3(Qorder[3*j+1]).vec[2])<fabs(eigenvalues.vec[Qorder[3*j+2]]*eigenvectors.ColumnVec3(Qorder[3*j+2]).vec[2])){ // if 1<2, switch (now 1 is second largest and 2 smallest)
//										if(eigenvalues.vec[Qorder[3*j+1]]<eigenvalues.vec[Qorder[3*j+2]]){ // if 1<2, switch (now 1 is second largest and 2 smallest)
											l=Qorder[3*j+1];
											Qorder[3*j+1]=Qorder[3*j+2];
											Qorder[3*j+2]=l;
#if DEBUG_LEVEL>3
											cout << "1->2\n";
#endif
										}
#if DEBUG_LEVEL>3
										cout << eigenvalues.vec[Qorder[3*j]] << ", " << eigenvalues.vec[Qorder[3*j+1]] << ", " << eigenvalues.vec[Qorder[3*j+2]] << "\n";
#endif
									}
									Vec3 switchorder=eigenvalues;
									eigenvalues.vec[0]=switchorder.vec[Qorder[3*j]];
									eigenvalues.vec[1]=switchorder.vec[Qorder[3*j+1]];
									eigenvalues.vec[2]=switchorder.vec[Qorder[3*j+2]];
									cf_params.results[frame/cf_params.average_frames][j][0]+=eigenvalues.vec[0];
									cf_params.results[frame/cf_params.average_frames][j][1]+=eigenvalues.vec[1];
									cf_params.results[frame/cf_params.average_frames][j][2]+=eigenvalues.vec[2];
									cf_params.results[frame/cf_params.average_frames][j][0+any_nr_output]+=eigenvalues.vec[0]*eigenvalues.vec[0];
									cf_params.results[frame/cf_params.average_frames][j][1+any_nr_output]+=eigenvalues.vec[1]*eigenvalues.vec[1];
									cf_params.results[frame/cf_params.average_frames][j][2+any_nr_output]+=eigenvalues.vec[2]*eigenvalues.vec[2];
#if DEBUG_LEVEL>3
									cout << "\nQ-tensor:\n" << Qtensor.M3Str() << "\n";
									cout << "Tr(Q) = " << Qtensor.M3Trace() << "\n";
#endif
									eigenvectors=Qtensor.Eigenvectors(eigenvalues,true);
									cf_params.results[frame/cf_params.average_frames][j][3]+=eigenvectors.ColumnVec3(0).vec[0];
									cf_params.results[frame/cf_params.average_frames][j][4]+=eigenvectors.ColumnVec3(0).vec[1];
									cf_params.results[frame/cf_params.average_frames][j][5]+=eigenvectors.ColumnVec3(0).vec[2];
									cf_params.results[frame/cf_params.average_frames][j][3+any_nr_output]+=eigenvectors.ColumnVec3(0).vec[0]*eigenvectors.ColumnVec3(0).vec[0];
									cf_params.results[frame/cf_params.average_frames][j][4+any_nr_output]+=eigenvectors.ColumnVec3(0).vec[1]*eigenvectors.ColumnVec3(0).vec[1];
									cf_params.results[frame/cf_params.average_frames][j][5+any_nr_output]+=eigenvectors.ColumnVec3(0).vec[2]*eigenvectors.ColumnVec3(0).vec[2];
									cf_params.results[frame/cf_params.average_frames][j][6]+=eigenvectors.ColumnVec3(1).vec[0];
									cf_params.results[frame/cf_params.average_frames][j][7]+=eigenvectors.ColumnVec3(1).vec[1];
									cf_params.results[frame/cf_params.average_frames][j][8]+=eigenvectors.ColumnVec3(1).vec[2];
									cf_params.results[frame/cf_params.average_frames][j][6+any_nr_output]+=eigenvectors.ColumnVec3(1).vec[0]*eigenvectors.ColumnVec3(1).vec[0];
									cf_params.results[frame/cf_params.average_frames][j][7+any_nr_output]+=eigenvectors.ColumnVec3(1).vec[1]*eigenvectors.ColumnVec3(1).vec[1];
									cf_params.results[frame/cf_params.average_frames][j][8+any_nr_output]+=eigenvectors.ColumnVec3(1).vec[2]*eigenvectors.ColumnVec3(1).vec[2];
									cf_params.results[frame/cf_params.average_frames][j][9]+=eigenvectors.ColumnVec3(2).vec[0];
									cf_params.results[frame/cf_params.average_frames][j][10]+=eigenvectors.ColumnVec3(2).vec[1];
									cf_params.results[frame/cf_params.average_frames][j][11]+=eigenvectors.ColumnVec3(2).vec[2];
									cf_params.results[frame/cf_params.average_frames][j][9+any_nr_output]+=eigenvectors.ColumnVec3(2).vec[0]*eigenvectors.ColumnVec3(2).vec[0];
									cf_params.results[frame/cf_params.average_frames][j][10+any_nr_output]+=eigenvectors.ColumnVec3(2).vec[1]*eigenvectors.ColumnVec3(2).vec[1];
									cf_params.results[frame/cf_params.average_frames][j][11+any_nr_output]+=eigenvectors.ColumnVec3(2).vec[2]*eigenvectors.ColumnVec3(2).vec[2];
									P1s.vec[0]=fabs(usum*eigenvectors.ColumnVec3(0));
									P1s.vec[1]=fabs(usum*eigenvectors.ColumnVec3(1));
									P1s.vec[2]=fabs(usum*eigenvectors.ColumnVec3(2));
									P1s/=mu_nr;
									// Go over data again for cos^3 since we just now got the directors
									ii=0;
									cos3s=Vec3(0.0);
									while(ii<cf_params.total_nr[0]){
										k=0;
										mu_nr=0;
										while(k<cf_params.items_involved){
											u=cf_params.vectors[k][ii];
											u/=u.V3Norm();
											cos3s.vec[0]+=qqpwr((u*eigenvectors.ColumnVec3(0)),3)*cf_params.weights[k];
											cos3s.vec[1]+=qqpwr((u*eigenvectors.ColumnVec3(1)),3)*cf_params.weights[k];
											cos3s.vec[2]+=qqpwr((u*eigenvectors.ColumnVec3(2)),3)*cf_params.weights[k];
											mu_nr+=cf_params.weights[k]*cf_params.total_nr[k];
											k++;
										}
										ii++;
									}
									cos3s/=mu_nr;
									cf_params.results[frame/cf_params.average_frames][j][12]+=P1s.vec[0];
									cf_params.results[frame/cf_params.average_frames][j][13]+=P1s.vec[1];
									cf_params.results[frame/cf_params.average_frames][j][14]+=P1s.vec[2];
									cf_params.results[frame/cf_params.average_frames][j][12+any_nr_output]+=P1s.vec[0]*P1s.vec[0];
									cf_params.results[frame/cf_params.average_frames][j][13+any_nr_output]+=P1s.vec[1]*P1s.vec[1];
									cf_params.results[frame/cf_params.average_frames][j][14+any_nr_output]+=P1s.vec[2]*P1s.vec[2];
									cf_params.results[frame/cf_params.average_frames][j][15]+=1.5*eigenvalues.vec[0];
									cf_params.results[frame/cf_params.average_frames][j][16]+=sgn(eigenvalues.vec[0])*(-0.5*(2.0*eigenvalues.vec[1]+eigenvalues.vec[0]));
									cf_params.results[frame/cf_params.average_frames][j][15+any_nr_output]+=(1.5*eigenvalues.vec[0])*(1.5*eigenvalues.vec[0]);
									cf_params.results[frame/cf_params.average_frames][j][16+any_nr_output]+=(sgn(eigenvalues.vec[0])*(-0.5*(2.0*eigenvalues.vec[1]+eigenvalues.vec[0])))*(sgn(eigenvalues.vec[0])*(-0.5*(2.0*eigenvalues.vec[1]+eigenvalues.vec[0])));
									cf_params.results[frame/cf_params.average_frames][j][17]+=musum.V3Norm()/mu_max;
									cf_params.results[frame/cf_params.average_frames][j][18]+=fabs(usum*(musum/musum.V3Norm()))/mu_nr;
									cf_params.results[frame/cf_params.average_frames][j][17+any_nr_output]+=(musum.V3Norm()/mu_max)*(musum.V3Norm()/mu_max);
									cf_params.results[frame/cf_params.average_frames][j][18+any_nr_output]+=(fabs(usum*(musum/musum.V3Norm()))/mu_nr)*(fabs(usum*(musum/musum.V3Norm()))/mu_nr);
									cf_params.results[frame/cf_params.average_frames][j][19]+=cos3s.vec[0];
									cf_params.results[frame/cf_params.average_frames][j][20]+=cos3s.vec[1];
									cf_params.results[frame/cf_params.average_frames][j][21]+=cos3s.vec[2];
									cf_params.results[frame/cf_params.average_frames][j][19+any_nr_output]+=cos3s.vec[0]*cos3s.vec[0];
									cf_params.results[frame/cf_params.average_frames][j][20+any_nr_output]+=cos3s.vec[1]*cos3s.vec[1];
									cf_params.results[frame/cf_params.average_frames][j][21+any_nr_output]+=cos3s.vec[2]*cos3s.vec[2];
#if DEBUG_LEVEL>3
									cout << "sorted eigenvalues = " << eigenvalues.V3Str(',') << "\n";
									cout << "S = " << 1.5*eigenvalues.vec[0] << "    ;    R = " << sgn(eigenvalues.vec[0])*(-0.5*(2.0*eigenvalues.vec[1]+eigenvalues.vec[0])) << "\n";
									cout << "normalized eigenvectors (director axis column vectors):\n" << eigenvectors.M3Str() << "\n";
									cout << "Polarization in director axis directions: " <<  P1s.V3Str(',') << "\n";
									cout << "P_max = |M|/|M_max| = " << musum.V3Norm()/mu_max << ", P_dipole = |<e_mu*e_M>| = " << fabs(usum*(musum/musum.V3Norm()))/mu_nr << ", |P_1| = " << P1s.V3Norm() << "\n";
#endif
								}
								break;
							case gmu:
								cosg=VectorCos(cf_params.vectors[0][i],cf_params.vectors[1][i]);
								if(r<rmax) cf_params.results[frame/cf_params.average_frames][j][(unsigned int)x]+=qqpwr(cosg,cf_params.typenr[j]);
								break;
							case galpha:
								cosa=VectorCos(cf_params.vectors[0][i],distance);
								if(r<rmax) cf_params.results[frame/cf_params.average_frames][j][(unsigned int)x]+=qqpwr(cosa,cf_params.typenr[j]);
								break;
							case cosine:
								k=1;
								while(k<cf_params.items_involved){
									cosg=qqpwr(VectorCos(cf_params.vectors[0][i],cf_params.vectors[k][i]),cf_params.typenr[j]);
									average[j]+=cosg/(cf_params.items_involved-1);
									average2[j]+=cosg*cosg/(cf_params.items_involved-1);
									if(cf_params.typenr[j]%2==0){
										x2=(unsigned int)fastfloor(cosg*cf_params.nr_bins);
									} else x2=(unsigned int)fastfloor(0.5*(cosg+1.0)*cf_params.nr_bins);
									if(x2>=cf_params.nr_bins) x2=cf_params.nr_bins-1;
									cf_params.results[frame/cf_params.average_frames][j][x2]+=1.0;
									k++;
								}
								break;
							case dihedral:
								// distances between center and other two elements
								Vec3 a=(cf_params.positions[0][i]-cf_params.positions[1][i]); // vector A_1->A_2
								Vec3 b=(cf_params.positions[2][i]-cf_params.positions[3][i]); // vector B_1->B_2
								Vec3 delta=(cf_params.positions[1][i]-cf_params.positions[2][i]); // vector A_2->B_1 (rotation axis)
								b.V3Cross(delta); // b is normal on plane spanned by B_1->B_2 and A_2->B_1 vectors
								double sine=(a*b)*delta.V3Norm(); // sin(theta) = [(B_1->B_2) x (A_2->B_1)]*(A_1->A_2) * |A_2->B_1| / (|(A_1->A_2) x (A_2->B_1)|*|(B_1->B_2) x (A_2->B_1)|)
								a.V3Cross(delta); // a is normal on plane spanned by A_1->A_2 and A_2->B_1 vectors
								double cosine=a*b; // cos(theta) = [(A_1->A_2) x (A_2->B_1)]*[(B_1->B_2) x (A_2->B_1)] / (|(A_1->A_2) x (A_2->B_1)|*|(B_1->B_2) x (A_2->B_1)|)
								// denominator is omitted for sine and cosine:
								// theta=atan(sine/cosine) <- it cancels, right?
								double theta=0.0;
								if(abs(cosine)>1E-12) theta=atan(sine/cosine);
								average[j]+=theta;
								average2[j]+=theta*theta;
								x2=(unsigned int)fastfloor((theta+pi/2.0)/pi*cf_params.nr_bins);
								if(x2>=cf_params.nr_bins) x2=cf_params.nr_bins-1;
								cf_params.results[frame/cf_params.average_frames][j][x2]+=1.0;
								break;
						}
					}
				} else{ // element-in-group vs. additional type (e.g. z-Axis)
					if(cf_params.total_nr[(a+1)%2]>1){
						cout << "Something went wrong here, there can only be one additional type.\n";
						exit(42);
					}
					distance=cf_params.positions[(a+1)%2][0]-cf_params.positions[a][i];
					apply_PBCs(distance,cf_params.configuration);
					r=distance.V3Norm();
					x=floor(r/dr);
					if((unsigned int)x>=cf_params.nr_bins) x=cf_params.nr_bins-1;
					for(unsigned int j=0; j<cf_params.nr_correlations; j++){
						switch(cf_params.types[j]){
							default:
							case gr:
							case distances:
								xn=floor(qqpwr(r,cf_params.typenr[j])/(qqpwr(rmax,cf_params.typenr[j])/(double)cf_params.nr_bins));
								if((unsigned int)xn>=cf_params.nr_bins) xn=cf_params.nr_bins-1;
								if(r<rmax) cf_params.results[frame/cf_params.average_frames][j][(unsigned int)xn]+=1.0;
								break;
							case sdf:
								if(r<rmax){
									avec=cf_params.rots[a][i].ColumnVec3(2); // z-direction of source element
									cosa=VectorCos(avec,distance); // [-1;+1]
									x2=(unsigned int)fastfloor((cosa+1.0)/2.0*ga_bins);
									if(x2>=ga_bins) x2=ga_bins-1;
									phi=atan2(distance*cf_params.rots[a][i].ColumnVec3(1),distance*cf_params.rots[a][i].ColumnVec3(0)); // [-pi;+pi]
									xn=(unsigned int)fastfloor((phi/pi+1.0)/2.0*ga_bins);
									if(xn>=ga_bins) xn=ga_bins-1;
									cf_params.results[frame/cf_params.average_frames][j][(x2*ga_bins+xn)*cf_params.nr_bins+(unsigned int)x]+=1.0;
									// determine (and keep up-to-date) maximum likely position data
									if(cf_params.results[frame/cf_params.average_frames][j][(x2*ga_bins+xn)*cf_params.nr_bins+(unsigned int)x]>cf_params.results[frame/cf_params.average_frames][j][ga_bins*ga_bins*cf_params.nr_bins+2*(x2*ga_bins+xn)]){
										cf_params.results[frame/cf_params.average_frames][j][ga_bins*ga_bins*cf_params.nr_bins+2*(x2*ga_bins+xn)]=cf_params.results[frame/cf_params.average_frames][j][(x2*ga_bins+xn)*cf_params.nr_bins+(unsigned int)x];
										cf_params.results[frame/cf_params.average_frames][j][ga_bins*ga_bins*cf_params.nr_bins+2*(x2*ga_bins+xn)+1]=(unsigned int)x;
									}
								}
								break;
							case angle:
								theta=VectorAngle(cf_params.vectors[a][i],cf_params.vectors[(a+1)%2][0]);
								x2=(unsigned int)fastfloor(theta/pi*cf_params.nr_bins);
								if(x2>=cf_params.nr_bins) x2=cf_params.nr_bins-1;
								cf_params.results[frame/cf_params.average_frames][j][x2]+=1.0;
								break;
							case gammaalpha:
								cosa=VectorCos(cf_params.vectors[a][i],distance);
								x2=(unsigned int)fastfloor(0.5*(cosa+1.0)*ga_bins);
								if(x2>=ga_bins) x2=ga_bins-1;
								cosg=VectorCos(cf_params.vectors[a][i],cf_params.vectors[(a+1)%2][0]);
								if(r<rmax){
									gammaalpha_nr[frame/cf_params.average_frames][(unsigned int)floor(r/ga_dr)*ga_bins+x2]++;
									cf_params.results[frame/cf_params.average_frames][j][(unsigned int)floor(r/ga_dr)*ga_bins+x2]+=cosg;
									cf_params.results[frame/cf_params.average_frames][j][(unsigned int)floor(r/ga_dr)*ga_bins+x2+ga_bins*ga_bins]+=cosg*cosg;
								}
								break;
							case gmu:
								cosg=VectorCos(cf_params.vectors[a][i],cf_params.vectors[(a+1)%2][0]);
								if(r<rmax) cf_params.results[frame/cf_params.average_frames][j][(unsigned int)x]+=qqpwr(cosg,cf_params.typenr[j]);
								break;
							case galpha:
								cosa=VectorCos(cf_params.vectors[a][i],distance);
								if(r<rmax) cf_params.results[frame/cf_params.average_frames][j][(unsigned int)x]+=qqpwr(cosa,cf_params.typenr[j]);
								break;
							case cosine:
								cosg=qqpwr(VectorCos(cf_params.vectors[a][i],cf_params.vectors[(a+1)%2][0]),cf_params.typenr[j]);
								average[j]+=cosg;
								average2[j]+=cosg*cosg;
								if(cf_params.typenr[j]%2==0){
									x2=(unsigned int)fastfloor(cosg*cf_params.nr_bins);
								} else x2=(unsigned int)fastfloor(0.5*(cosg+1.0)*cf_params.nr_bins);
								if(x2>=cf_params.nr_bins) x2=cf_params.nr_bins-1;
								cf_params.results[frame/cf_params.average_frames][j][x2]+=1.0;
								break;
						}
					}
				}
				result_nr++;
				nr_results_total++;
			} else{
				unsigned int l;
				for(unsigned int j=0; j<cf_params.total_nr[(a+1)%2]; j++){
					if(!(exclude_idx && (i==j))){
						distance=cf_params.positions[(a+1)%2][j]-cf_params.positions[a][i];
						apply_PBCs(distance,cf_params.configuration);
						r=distance.V3Norm();
						x=floor(r/dr);
						if((unsigned int)x>=cf_params.nr_bins) x=cf_params.nr_bins-1;
						for(unsigned int k=0; k<cf_params.nr_correlations; k++){
							switch(cf_params.types[k]){
								default:
								case gr:
									if(r<rmax) cf_params.results[frame/cf_params.average_frames][k][(unsigned int)x]+=cf_params.V[frame+cf_params.startframe-cf_params.start_frame_in_data]/Vshell(x*dr,dr);
									break;
								case distances:
									if(r<rmax){
										xn=floor(qqpwr(r,cf_params.typenr[k])/(qqpwr(rmax,cf_params.typenr[k])/(double)cf_params.nr_bins));
										if((unsigned int)xn>=cf_params.nr_bins) xn=cf_params.nr_bins-1;
										cf_params.results[frame/cf_params.average_frames][k][(unsigned int)xn]+=1.0;
									}
									break;
								case sdf:
									if(r<rmax){
										avec=cf_params.rots[a][i].ColumnVec3(2); // z-direction of source element
										cosa=VectorCos(avec,distance); // [-1;+1]
										x2=(unsigned int)fastfloor((cosa+1.0)/2.0*ga_bins);
										if(x2>=ga_bins) x2=ga_bins-1;
										phi=atan2(distance*cf_params.rots[a][i].ColumnVec3(1),distance*cf_params.rots[a][i].ColumnVec3(0)); // [-pi;+pi]
										xn=(unsigned int)fastfloor((phi/pi+1.0)/2.0*ga_bins);
										if(xn>=ga_bins) xn=ga_bins-1;
										// rho(r,theta,phi) = <rho(theta,phi)>*g(r,theta,phi) => g(r) = rho(r,theta,phi)/<rho(theta,phi)>
										// local density:
										// 		rho(r,theta,phi) = N(r,theta,phi)/Vvoxel(r,theta,phi)
										// average number density per slab (there are ga_bins^2 slabs of identical size):
										// 		<rho(theta,phi)> = (N/ga_bins^2)/(Vtotal/ga_bins^2) = N/Vtotal
										// Note: since we divide by cos(theta) and phi N(theta,phi) = N/ga_bins^2 and Vslab(theta,phi) = Vtotal/ga_bins^2
										// => g(r,theta,phi) = N(r,theta,phi)*Vtotal/Vvoxel(r,theta,phi)/N
										cf_params.results[frame/cf_params.average_frames][k][(x2*ga_bins+xn)*cf_params.nr_bins+(unsigned int)x]+=cf_params.V[frame+cf_params.startframe-cf_params.start_frame_in_data]/Vvoxel(x*dr,dr,2.0/ga_bins,2.0*pi/ga_bins);
										// determine (and keep up-to-date) maximum likely position data
										if(cf_params.results[frame/cf_params.average_frames][k][(x2*ga_bins+xn)*cf_params.nr_bins+(unsigned int)x]>cf_params.results[frame/cf_params.average_frames][k][ga_bins*ga_bins*cf_params.nr_bins+2*(x2*ga_bins+xn)]){
											cf_params.results[frame/cf_params.average_frames][k][ga_bins*ga_bins*cf_params.nr_bins+2*(x2*ga_bins+xn)]=cf_params.results[frame/cf_params.average_frames][k][(x2*ga_bins+xn)*cf_params.nr_bins+(unsigned int)x];
											cf_params.results[frame/cf_params.average_frames][k][ga_bins*ga_bins*cf_params.nr_bins+2*(x2*ga_bins+xn)+1]=(unsigned int)x;
										}
									}
									break;
								case angle:
									theta=VectorAngle(cf_params.vectors[a][i],cf_params.vectors[(a+1)%2][j]);
									x2=(unsigned int)fastfloor(theta/pi*cf_params.nr_bins);
									if(x2>=cf_params.nr_bins) x2=cf_params.nr_bins-1;
									cf_params.results[frame/cf_params.average_frames][k][x2]+=1.0;
									break;
								case gammaalpha:
									cosa=VectorCos(cf_params.vectors[a][i],distance);
									x2=(unsigned int)fastfloor(0.5*(cosa+1.0)*ga_bins);
									if(x2>=ga_bins) x2=ga_bins-1;
									cosg=VectorCos(cf_params.vectors[a][i],cf_params.vectors[(a+1)%2][j]);
									if(r<rmax){
										gammaalpha_nr[frame/cf_params.average_frames][(unsigned int)floor(r/ga_dr)*ga_bins+x2]++;
										cf_params.results[frame/cf_params.average_frames][k][(unsigned int)floor(r/ga_dr)*ga_bins+x2]+=cosg;
										cf_params.results[frame/cf_params.average_frames][k][(unsigned int)floor(r/ga_dr)*ga_bins+x2+ga_bins*ga_bins]+=cosg*cosg;
									}
									break;
								case gmu:
									cosg=VectorCos(cf_params.vectors[a][i],cf_params.vectors[(a+1)%2][j]);
									if(r<rmax) cf_params.results[frame/cf_params.average_frames][k][(unsigned int)x]+=qqpwr(cosg,cf_params.typenr[k])*cf_params.V[frame+cf_params.startframe-cf_params.start_frame_in_data]/Vshell(x*dr,dr);
									break;
								case galpha:
									cosa=VectorCos(cf_params.vectors[a][i],distance);
									if(r<rmax) cf_params.results[frame/cf_params.average_frames][k][(unsigned int)x]+=qqpwr(cosa,cf_params.typenr[k])*cf_params.V[frame+cf_params.startframe-cf_params.start_frame_in_data]/Vshell(x*dr,dr);
									break;
								case cosine:
									l=1;
									cosg=0.0;
									while(l<cf_params.items_involved){
										cosa=qqpwr(VectorCos(cf_params.vectors[0][i],cf_params.vectors[l][j]),cf_params.typenr[k]);
										if(!cf_params.find_max_item){
											cosg+=cosa/(cf_params.items_involved-1);
										} else{
											if(cf_params.find_max_abs_item){
												if((fabs(cosa)>fabs(cosg)) || (l==1)) cosg=cosa;
											} else{
												if((cosa>cosg) || (l==1)) cosg=cosa;
											}
										}
										l++;
									}
									average[k]+=cosg;
									average2[k]+=cosg*cosg;
#if DEBUG_LEVEL>3
									cout << (cf_params.startframe+frame)*cf_params.configuration->grfreq << ": " << cosg << " (" << cf_params.vectors[a][i].V3Str(',') << "), (" << cf_params.vectors[(a+1)%2][j].V3Str(',') << ")\n";
#endif
									if(cf_params.typenr[k]%2==0){
										x2=(unsigned int)fastfloor(cosg*cf_params.nr_bins);
									} else x2=(unsigned int)fastfloor(0.5*(cosg+1.0)*cf_params.nr_bins);
									if(x2>=cf_params.nr_bins) x2=cf_params.nr_bins-1;
									cf_params.results[frame/cf_params.average_frames][k][x2]+=1.0;
									break;
							}
						}
						result_nr++;
						nr_results_total++;
					}
				}
			}
		}
		double percent=100.0*double(frame+1)/double(cf_params.n_frames);
		if(percent>=percentage+5){
			percentage=(unsigned int)floor(percent/5.0)*5;
			if(percentage%20==0) cout << percentage << "%"; else cout << ".";
			cout.flush();
		}
	}
	cout << "\n";
	
	if(cf_params.individual_output){
		fstream output;
		string filename = cf_params.configuration->fileout+"_analysis_"+int2str(analysis_nr+1)+".dat";
		if(cf_params.report) filename = cf_params.configuration->fileout+"_Analysis_"+int2str(analysis_nr+1)+".dat";
		output.open(filename.c_str(), fstream::out | fstream::trunc);
		output << "#Title:\t";
		if(cf_params.group_internal && !has_Q){
			output << "Group-internal ";
		} else output << "System-wide ";
		if(!cf_params.any_number_of_items){
			if(cf_params.items_involved>1) output << "correlation between: " << cf_params.item_names[0] << " and " << cf_params.item_names[1] << "\n"; else output << "autocorrelation for: " << cf_params.item_names[0] << "\n";
		} else{
			output << "analysis for: ";
			for(unsigned int i=0; i<cf_params.items_involved; i++){
				if(cf_params.weights[i]!=1.0) output << cf_params.weights[i] << "*";
				output << cf_params.item_names[i];
				if(i+1<cf_params.items_involved) output << ", ";
			}
			output << "\n";
		}
		string* x_label=new string[cf_params.nr_correlations];
		string* y_label=new string[cf_params.nr_correlations];
		double* x_min=new double[cf_params.nr_correlations];
		double* x_max=new double[cf_params.nr_correlations];
		double* dx=new double[cf_params.nr_correlations];
		string to_the_n;
		unsigned int additional_corrs=0;
		for(unsigned int i=0; i<cf_params.nr_correlations; i++){
			x_label[i]="";
			y_label[i]="";
			x_min[i]=0.0;
			x_max[i]=0.0;
			dx[i]=0.0;
			switch(cf_params.types[i]){
				default:
				case gr:
				case gmu:
				case galpha:
					x_label[i]="Distance [Angström]";
					x_min[i]=0.0;
					x_max[i]=rmax;
					y_label[i]=output_corr_types[cf_params.types[i]];
					if((cf_params.types[i]==gmu) && (cf_params.typenr[i]>1)) y_label[i]+="^"+int2str(cf_params.typenr[i]);
					y_label[i]+="(r)";
					dx[i]=dr;
					break;
				case sdf:
					x_label[i]="cos(theta)\tphi\tr_max [Angström]\tsdf(r)\tDistance [Angström]\tcos(theta)\tphi";
					x_min[i]=0.0;
					x_max[i]=rmax;
					y_label[i]=output_corr_types[cf_params.types[i]];
					y_label[i]+="(r)";
					dx[i]=dr;
					break;
				case distances:
					x_min[i]=0.0;
					x_max[i]=qqpwr(rmax,cf_params.typenr[i]);
					to_the_n="";
					if(cf_params.typenr[i]>1) to_the_n="^"+int2str(cf_params.typenr[i]);
					average[i]/=nr_results_total;
					average2[i]/=nr_results_total;
					y_label[i]="Normalized counts";
					dx[i]=(qqpwr(rmax,cf_params.typenr[i])/(double)cf_params.nr_bins);
					x_label[i]="Distance"+to_the_n+" [Angström"+to_the_n+"]";
					break;
				case angle:
					x_label[i]="Angle [degree]";
					average[i]/=nr_results_total;
					average2[i]/=nr_results_total;
					x_min[i]=0.0;
					x_max[i]=180.0;
					dx[i]=180.0/cf_params.nr_bins;
					y_label[i]="Normalized counts";
					break;
				case alphabetagamma:
					x_label[i]="Alpha\tBeta\tGamma";
					average[i]/=nr_results_total;
					average2[i]/=nr_results_total;
					x_min[i]=-pi;
					x_max[i]=pi;
					dx[i]=2.0*pi/cf_params.nr_bins;
					y_label[i]="Normalized counts\tNormalized counts\tNormalized counts";
					Vec3average[i]/=nr_results_total;
					Vec3average2[i]/=nr_results_total;
					additional_corrs+=2;
					break;
				case gammaalpha:
					x_label[i]="cos(alpha)\tDistance [Angström]";
					x_min[i]=0.0;
					x_max[i]=rmax;
					dx[i]=ga_dr;
					y_label[i]="<cos(gamma)>\t+/-";
					break;
				case Q:
					x_label[i]="lambda_1\t+/-\tlambda_2\t+/-\tlambda_3\t+/-\tv_1_x\t+/-\tv_1_y\t+/-\tv_1_z\t+/-\tv_2_x\t+/-\tv_2_y\t+/-\tv_2_z\t+/-\tv_3_x\t+/-\tv_3_y\t+/-\tv_3_z\t+/-\tP_v_1\t+/-\tP_v_2\t+/-\tP_v_3\t+/-\tS\t+/-\tR\t+/-\tP_max\t+/-\tP_dipole\t+/-\tcos^3(v_1)\t+/-\tcos^3(v_2)\t+/-\tcos^3(v_3)\t+/-";
					y_label[i]="";
					x_min[i]=-1.0;
					x_max[i]=1.0;
					dx[i]=0.0;
					break;
				case cosine:
					x_label[i]="cos";
					if(cf_params.typenr[i]>1) x_label[i]+="^"+int2str(cf_params.typenr[i]);
					x_label[i]+="(gamma)";
					average[i]/=nr_results_total;
					average2[i]/=nr_results_total;
					dx[i]=1.0/cf_params.nr_bins;
					if(cf_params.typenr[i]%2==0){
						dx[i]=1.0/cf_params.nr_bins;
						x_min[i]=0.0;
					} else{
						dx[i]=2.0/cf_params.nr_bins;
						x_min[i]=-1.0;
					}
					x_max[i]=1.0;
					y_label[i]="Normalized counts";
					break;
				case dihedral:
					x_label[i]="dihedral";
					average[i]/=nr_results_total;
					average2[i]/=nr_results_total;
					dx[i]=pi/cf_params.nr_bins;
					x_min[i]=-pi/2.0;
					x_max[i]=pi/2.0;
					y_label[i]="Normalized counts";
					break;
			}
		}
		output << "#N correlations:\t" << cf_params.nr_correlations+additional_corrs << "\n";
		output << "#N per frame:\t" << results_per_frame << "\n";
		output << "#N total:\t" << nr_results_total << "\n";
		output << "#N frames averaged:\t" << cf_params.average_frames << "\n";
		output << "#N histograms:\t" << cf_params.n_frames/cf_params.average_frames << "\n";
		output.precision(6);
		output.setf(ios::fixed, ios::floatfield);
		output << "#averages:";
		for(unsigned int i=0; i<cf_params.nr_correlations; i++){
			if(cf_params.types[i]==alphabetagamma) output << "\t" << Vec3average[i].V3Str('\t'); else{
				if(compare_strings(y_label[i].c_str(),"Normalized counts")) output << "\t" << average[i]; else output << "\t-";
			}
		}
		output << "\n";
		output << "#stddevs:";
		for(unsigned int i=0; i<cf_params.nr_correlations; i++){
			if(cf_params.types[i]==alphabetagamma) output << "\t" << sqrt(Vec3average2[i].vec[0]-Vec3average[i].vec[0]*Vec3average[i].vec[0]) << "\t" << sqrt(Vec3average2[i].vec[1]-Vec3average[i].vec[1]*Vec3average[i].vec[1]) << "\t" << sqrt(Vec3average2[i].vec[2]-Vec3average[i].vec[2]*Vec3average[i].vec[2]); else{
				if(compare_strings(y_label[i].c_str(),"Normalized counts")) output << "\t" << sqrt(average2[i]-average[i]*average[i]); else output << "\t-";
			}
		}
		output << "\n";
		output << "#x_label:";
		for(unsigned int i=0; i<cf_params.nr_correlations; i++) output << "\t" << x_label[i];
		output << "\n";
		output << "#y_label:";
		for(unsigned int i=0; i<cf_params.nr_correlations; i++) if(y_label[i]!="") output << "\t" << y_label[i];
		output << "\n";
		output << "#x_min:";
		for(unsigned int i=0; i<cf_params.nr_correlations; i++){
			output << "\t" << x_min[i];
			if(cf_params.types[i]==alphabetagamma) output << "\t" << 0.0 << "\t" << x_min[i];
		}
		output << "\n";
		output << "#x_max:";
		for(unsigned int i=0; i<cf_params.nr_correlations; i++){
			output << "\t" << x_max[i];
			if(cf_params.types[i]==alphabetagamma) output << "\t" << x_max[i] << "\t" << x_max[i];
		}
		output << "\n";
		output << "#binwidth:";
		for(unsigned int i=0; i<cf_params.nr_correlations; i++){
			output << "\t" << dx[i];
			if(cf_params.types[i]==alphabetagamma) output << "\t" << dx[i]/2.0 << "\t" << dx[i];
		}
		output << "\n";
		output << "#Step";
		for(unsigned int i=0; i<cf_params.nr_correlations; i++){
			output << "\t" << x_label[i];
			if(y_label[i]!="") output << "\t" << y_label[i];
		}
		output << "\n";
		cout.flush();
		percentage=0;
		double dcos=1.0/cf_params.nr_bins;
		unsigned int output_rows=cf_params.nr_bins;
		if(has_gammaalpha) output_rows=ga_bins*ga_bins;
		ofstream vizout;
		double sdf_min=1.0;
		double sdf_max=0.0;
		double sdf_avg=0.0;
		string sdf_coords="";
		string sdf_colors="";
		double** sdf_plot=NULL;
		if(has_sdf){
			sdf_plot=new double*[cf_params.n_frames/cf_params.average_frames+2];
			for(unsigned int frame=0; frame<cf_params.n_frames+(cf_params.average_frames==1); frame++){
				if(frame%cf_params.average_frames==0){
					sdf_plot[frame/cf_params.average_frames]=new double[2*ga_bins*ga_bins];
					for(unsigned int i=0; i<ga_bins*ga_bins; i++){
						sdf_plot[frame/cf_params.average_frames][2*i]=0.0;
						sdf_plot[frame/cf_params.average_frames][2*i+1]=0.0;
					}
				}
			}
			output_rows=ga_bins*ga_bins*cf_params.nr_bins;
			cout << "\t\t\t\t-> Generating SDF X3D file\n";
			// Create file
			string vizname = cf_params.configuration->fileout+"_analysis_"+int2str(analysis_nr+1)+".x3d";;
			vizout.open(vizname.c_str());
			if (vizout.fail()){
				cout << "ERROR: Unable to open output file.\n";
				exit(1);
			}
			// Configure viewpoint -- want to have left-handed coordinate system (x to the right, y to the back, z up)
			double position_scale=cbrt(averageV/cf_params.configuration->V);
			viewpoint_t viewpoint;
			viewpoint.Center = Vec3(0,0,0);
			double lookdownangle=-pi/7.854;
			viewpoint.Orientation = Vec4(1,0,0,pi*0.5+lookdownangle);
			viewpoint.Position = Vec3(0.0,position_scale*0.75*cf_params.configuration->boxlength[1]/lookdownangle,position_scale*0.75*cf_params.configuration->boxlength[2]);
			viewpoint.Description = "Default viewpoint";
			// Output file header info
			output_header(vizout,vizname,viewpoint,false,0);
			// Define elements/groups used in SDF analysis
			vizout << "\t\t<!-- Define element types -->\n";
			if(cf_params.items[0].group_idx<0){ // independent element
				Element_Type* type=cf_params.configuration->element_types[cf_params.items[0].element_idx];
				define_element(vizout,type,NULL,false,0,cf_params.configuration);
				vizout << "\t\t<!-- done -->\n";
				// Output the element centered at the SDF source position
				vizout << "\t\t<Transform>\n";
				vizout <<  "\t\t\t<ProtoInstance name='" << type->name << "' DEF='" << type->name << "'>\n";
				vizout <<  "\t\t\t\t<fieldValue name='position' value='0 0 0'/>\n";
				vizout <<  "\t\t\t\t<fieldValue name='rotation' value='1 0 0 0'/>\n";
				vizout <<  "\t\t\t\t<fieldValue name='transparency' value='" << type->transparency << "'/>\n";
				vizout <<  "\t\t\t</ProtoInstance>\n";
				vizout <<  "\t\t</Transform>\n";
			} else{
				if(cf_params.items[0].element_idx>=0){ // element in group
					Element_Group* group=cf_params.configuration->groups[cf_params.items[0].traj_group];
					Element_Group* sdf_group=cf_params.configuration->groups[cf_params.items[0].group_idx];
					if((sdf_group!=group) && (group->levelofdetail>0)) group->Type->LOD->visual_original[group->levelofdetail-1]=true; // want to see underlying fully-atomistic structure if requested sdf is of lower LOD
					// Output the group centered at the SDF source position
					Vec3 sdf_center=cf_params.configuration->group_elements[sdf_group->elements[cf_params.items[0].element_idx]].center;
					Mat33 sdf_rot=cf_params.configuration->group_elements[sdf_group->elements[cf_params.items[0].element_idx]].rot.M3Transpose(); // rotate everthing into frame of sdf center element
					for(unsigned int i=0; i<group->nr_elements; i++){
						Element* element=&cf_params.configuration->group_elements[group->elements[i]];
						Vec3 position=sdf_rot*(element->center-sdf_center);
						double transparency=element->MyType->transparency;
						if(group->Type->transparency>=0.0) transparency=group->Type->transparency;
						if(group->Type->LOD){
							if((group->Type->LOD->nr_components>0) && (group->Type->LOD->element_in_component[i]>=0)){
								if(group->Type->LOD->component_transparency[group->Type->LOD->element_in_component[i]]>=0.0){
									transparency=group->Type->LOD->component_transparency[group->Type->LOD->element_in_component[i]];
								}
							}
						}
						define_element(vizout,cf_params.configuration->element_types[cf_params.configuration->group_elements[group->elements[i]].mytype],group,false,0,cf_params.configuration,true,false,position,Rot2AxisAngle(sdf_rot*element->rot),transparency);
					}
				} else{ // whole group (and since there is no frame of reference, error)
					cout << "ERROR: Whole group has no static frame of reference for SDF visualization.\n";
					exit(2);
				}
			}
		}
		unsigned int x_nr, y_nr, k;
		if(has_Q) output_rows=1;
		cout << "\t\t\t\t-> Creating individual data output: 0%";
		for(unsigned int frame=0; frame<cf_params.n_frames+(cf_params.average_frames==1); frame++){
			if(frame%cf_params.average_frames==0){
				for(unsigned int i=0; i<output_rows; i++){
					output << (cf_params.startframe+frame)*cf_params.configuration->grfreq;
					for(unsigned int j=0; j<cf_params.nr_correlations; j++){
						output << "\t";
						switch(cf_params.types[j]){
							default:
							case gr:
							case distances:
							case gmu:
							case galpha:
								if(i<cf_params.nr_bins){
									output << (i+0.5)*dx[j] << "\t" << cf_params.results[frame/cf_params.average_frames][j][i]/(cf_params.average_frames*results_per_frame);
								} else output << "-\t-";
								break;
							case sdf:
								x_nr=(i/cf_params.nr_bins)/ga_bins;
								y_nr=(i/cf_params.nr_bins)%ga_bins;
								if(i<ga_bins*ga_bins){
									cost=(double)(i/ga_bins+0.5)/ga_bins*2.0-1.0;
									phi=((double)(i%ga_bins+0.5)/ga_bins*2.0-1.0)*pi;
									r=(cf_params.results[frame/cf_params.average_frames][j][ga_bins*ga_bins*cf_params.nr_bins+2*i+1]+0.5)*dx[j];
									p=cf_params.results[frame/cf_params.average_frames][j][ga_bins*ga_bins*cf_params.nr_bins+2*i]/(cf_params.average_frames*results_per_frame);
									// filter out values which are above or below
									output << cost << "\t" << phi << "\t" << r << "\t" << p << "\t";
								} else output << "-\t-\t-\t-\t";
								r=(i%cf_params.nr_bins+0.5)*dx[j];
								output << r << "\t" << (double)(x_nr+0.5)/ga_bins*2.0-1.0 << "\t" << ((double)(y_nr+0.5)/ga_bins*2.0-1.0)*pi << "\t";
								// only output values which have counts (otherwise output "-" instead)
								p=cf_params.results[frame/cf_params.average_frames][j][i]/(cf_params.average_frames*results_per_frame);
								if(p>EPS) output << p; else output << "-";
								// determine (and keep up-to-date) maximum likely position data
								if(!((p<cf_params.additional_parameters[j][2]) || ((cf_params.additional_parameters[j][3]>0.0) && (p>cf_params.additional_parameters[j][3])))){
									if(p>sdf_plot[frame/cf_params.average_frames][2*(x_nr*ga_bins+y_nr)]){
										sdf_plot[frame/cf_params.average_frames][2*(x_nr*ga_bins+y_nr)]=p;
										sdf_plot[frame/cf_params.average_frames][2*(x_nr*ga_bins+y_nr)+1]=r;
									}
								}
								break;
							case angle:
								if(i<cf_params.nr_bins){
									output << (i+0.5)*180.0/cf_params.nr_bins << "\t" << cf_params.results[frame/cf_params.average_frames][j][i]/(cf_params.average_frames*results_per_frame);
								} else output << "-\t-";
								break;
							case alphabetagamma:
								if(i<cf_params.nr_bins){
									output << (i+0.5)*2.0*pi/cf_params.nr_bins-pi << "\t" << cf_params.results[frame/cf_params.average_frames][j][3*i]/(cf_params.average_frames*results_per_frame) << "\t" << (i+0.5)*pi/cf_params.nr_bins << "\t" << cf_params.results[frame/cf_params.average_frames][j][3*i+1]/(cf_params.average_frames*results_per_frame) << "\t" << (i+0.5)*2.0*pi/cf_params.nr_bins-pi << "\t" << cf_params.results[frame/cf_params.average_frames][j][3*i+2]/(cf_params.average_frames*results_per_frame);
								} else output << "-\t-";
								break;
							case gammaalpha:
								if(i<ga_bins*ga_bins){
									x_nr=i%ga_bins;
									y_nr=i/ga_bins;
									output << 2.0*(double)x_nr/ga_bins-1.0+1.0/ga_bins << "\t" << (y_nr+0.5)*ga_dr << "\t";
									if(gammaalpha_nr[frame/cf_params.average_frames][i]>0){
										cf_params.results[frame/cf_params.average_frames][j][i]/=gammaalpha_nr[frame/cf_params.average_frames][i];
										cf_params.results[frame/cf_params.average_frames][j][i+ga_bins*ga_bins]/=gammaalpha_nr[frame/cf_params.average_frames][i];
										cf_params.results[frame/cf_params.average_frames][j][i+ga_bins*ga_bins]-=cf_params.results[frame/cf_params.average_frames][j][i]*cf_params.results[frame/cf_params.average_frames][j][i];
										cf_params.results[frame/cf_params.average_frames][j][i+ga_bins*ga_bins]=sqrt(cf_params.results[frame/cf_params.average_frames][j][i+ga_bins*ga_bins]);
										output << cf_params.results[frame/cf_params.average_frames][j][i] << "\t" << cf_params.results[frame/cf_params.average_frames][j][i+ga_bins*ga_bins];
									} output << "-\t-";
								} else output << "-\t-";
								break;
							case Q:
								k=0;
								while((k<any_nr_output) && (i==0)){
									double avg=cf_params.results[frame/cf_params.average_frames][j][k]/(cf_params.average_frames);
									double avg2=cf_params.results[frame/cf_params.average_frames][j][k+any_nr_output]/(cf_params.average_frames);
									double var=avg2-avg*avg;
									if(var<EPS) var=0.0;
									output << avg << "\t" << sqrt(var);
									if(k+1<any_nr_output) output << "\t";
									k++;
								}
								break;
							case cosine:
								if(i<cf_params.nr_bins){
									if(cf_params.typenr[j]%2==0) output << i*dcos; else output << 2.0*i*dcos-1.0+dcos;
									output << "\t" << cf_params.results[frame/cf_params.average_frames][j][i]/(cf_params.average_frames*results_per_frame);
								} else output << "-\t-";
								break;
							case dihedral:
								if(i<cf_params.nr_bins){
									output << pi*i/cf_params.nr_bins-pi/2.0+pi/(2.0*cf_params.nr_bins) << "\t" << cf_params.results[frame/cf_params.average_frames][j][i]/(cf_params.average_frames*results_per_frame);
								} else output << "-\t-";
								break;
						}
					}
					output << "\n";
				}
			}
			double percent=100.0*double(frame+1)/double(cf_params.n_frames);
			if(percent>=percentage+5){
				percentage=(unsigned int)floor(percent/5.0)*5;
				if(percentage%20==0) cout << percentage << "%"; else cout << ".";
				cout.flush();
			}
		}
		output.close();
		cout << "\n";
		if(has_sdf){
			unsigned int sdf_plot_nr=cf_params.n_frames/cf_params.average_frames+1;
			sdf_plot[sdf_plot_nr]=new double[2*ga_bins*ga_bins];
			memset(sdf_plot[sdf_plot_nr],0,sizeof(double)*2*ga_bins*ga_bins);
			// Draw sdf around molecule
			for(unsigned int frame=0; frame<cf_params.n_frames+(cf_params.average_frames==1); frame++){
				if(frame%cf_params.average_frames==0){
					for(unsigned int i=0; i<ga_bins*ga_bins; i++){
						sdf_plot[sdf_plot_nr][2*i]+=sdf_plot[frame/cf_params.average_frames][2*i];
						sdf_plot[sdf_plot_nr][2*i+1]+=sdf_plot[frame/cf_params.average_frames][2*i+1];
					}
				}
			}
			for(unsigned int i=0; i<ga_bins*ga_bins; i++){
				sdf_plot[sdf_plot_nr][2*i]/=(cf_params.n_frames/cf_params.average_frames);
				sdf_plot[sdf_plot_nr][2*i+1]/=(cf_params.n_frames/cf_params.average_frames);
				p=sdf_plot[sdf_plot_nr][2*i];
				if(p<sdf_min) sdf_min=p;
				if(p>sdf_max) sdf_max=p;
				sdf_avg+=p;
			}
			sdf_avg/=ga_bins*ga_bins;
			double r_top=0.0;
			double p_top=0.0;
			double r_bottom=0.0;
			double p_bottom=0.0;
//			double r_s, ct_s, st_s, p_s;
			Vec3* points = new Vec3[ga_bins*ga_bins];
			Vec4* rgba = new Vec4[ga_bins*ga_bins];
			unsigned int top_count=0;
			unsigned int bottom_count=0;
			for(unsigned int i=0; i<ga_bins*ga_bins; i++){
				cost=(double)(i/ga_bins+0.5)/ga_bins*2.0-1.0;
				double sint=sqrt(1.0-cost*cost);
				phi=((double)(i%ga_bins+0.5)/ga_bins*2.0-1.0)*pi;
				r=sdf_plot[sdf_plot_nr][2*i+1];
/*				if(i%ga_bins==0){
					r_s=r;
					ct_s=cost;
					st_s=sint;
					p_s=phi;
					if(i>0) sdf_coords+=", "+double2str(r_s*st_s*cos(p_s))+" "+double2str(r_s*st_s*sin(p_s))+" "+double2str(ct_s*r_s);
				}*/
				p=sdf_plot[sdf_plot_nr][2*i];
				double alpha=(p-sdf_min)/(sdf_max-sdf_min);
				if(alpha<cf_params.additional_parameters[sdf_nr][0]) alpha=0.0;
				if(alpha>cf_params.additional_parameters[sdf_nr][1]) alpha=1.0;
				rgba[i]=Vec4(ColorMappingVec(p,sdf_min,sdf_max,sdf_avg),alpha); // alpha-channel (opacity) is proportional to probability of finding element at vertex
				if(i>0){
					sdf_coords+=", ";
					sdf_colors+=", ";
				}
				points[i].vec[0]=r*sint*cos(phi);
				points[i].vec[1]=r*sint*sin(phi);
				points[i].vec[2]=cost*r;
				sdf_coords+=points[i].V3Str(' ');
//				if(i==ga_bins*ga_bins-1) sdf_coords+=", "+double2str(r_s*st_s*cos(p_s))+" "+double2str(r_s*st_s*sin(p_s))+" "+double2str(ct_s*r_s);
				sdf_colors+=rgba[i].V4Str(' ');
				// calculate bottom point on sphere
				if(i<ga_bins){
					if(p>EPS){
						r_bottom+=r;
						p_bottom+=p;
						bottom_count++;
					}
				}
				if(i==ga_bins){
					if(bottom_count>0){
						r_bottom/=bottom_count;
						p_bottom/=bottom_count;
					}
					alpha=(p_bottom-sdf_min)/(sdf_max-sdf_min);
					if(alpha<cf_params.additional_parameters[sdf_nr][0]) alpha=0.0;
					if(alpha>cf_params.additional_parameters[sdf_nr][1]) alpha=1.0;
					// use point as first point in list
					sdf_coords="0 0 "+double2str(r_bottom*(-1.0))+", "+sdf_coords;
//					for(unsigned int j=0; j<ga_bins+1; j++) sdf_coords="0 0 -"+double2str(r_bottom)+", "+sdf_coords;
					sdf_colors=(Vec4(ColorMappingVec(p_bottom,sdf_min,sdf_max,sdf_avg),alpha)).V4Str(' ')+", "+sdf_colors;
				}
				// calculate top point on sphere
				if(i>=(ga_bins-1)*ga_bins){
					if(p>EPS){
						r_top+=r;
						p_top+=p;
						top_count++;
					}
				}
				if(i+1==ga_bins*ga_bins){
					if(top_count>0){
						r_top/=top_count;
						p_top/=top_count;
					}
					alpha=(p_top-sdf_min)/(sdf_max-sdf_min);
					if(alpha<cf_params.additional_parameters[sdf_nr][0]) alpha=0.0;
					if(alpha>cf_params.additional_parameters[sdf_nr][1]) alpha=1.0;
					// use point as last point in list
					sdf_coords+=", 0 0 "+double2str(r_top);
//					for(unsigned int j=0; j<ga_bins+1; j++) sdf_coords+=", 0 0 "+double2str(r_top);
					sdf_colors+=", "+(Vec4(ColorMappingVec(p_top,sdf_min,sdf_max,sdf_avg),alpha)).V4Str(' ');
				}
			}
			// calculate geometric centers of quadrilaterals
			for(unsigned int i=0; i<ga_bins*ga_bins-ga_bins; i++){
				// points of the quadrilateral:
				// - (B) i		(C) (i+1)%ga_bins+(i/ga_bins)*ga_bins (<- which is a fancy way of say to loop back to beginning of phi space)
				//             (E)
				// - (A) i+ga_bins	(D) (i+1)%ga_bins+(i/ga_bins+1)*ga_bins (like above, but moved down one slice)
				double norm=((double)(((double)sdf_plot[sdf_plot_nr][2*((i+1)%ga_bins+(i/ga_bins+1)*ga_bins)]>EPS)+((double)sdf_plot[sdf_plot_nr][2*((i+1)%ga_bins+(i/ga_bins)*ga_bins)]>EPS)+((double)sdf_plot[sdf_plot_nr][2*(i+ga_bins)]>EPS)+((double)sdf_plot[sdf_plot_nr][2*i]>EPS)));
				if(norm<EPS) norm=1.0; // to avoid 0/0 aka nan
				sdf_coords+=", "+((points[i+ga_bins]*((double)sdf_plot[sdf_plot_nr][2*(i+ga_bins)]>EPS)+points[i]*((double)sdf_plot[sdf_plot_nr][2*(i)]>EPS)+points[(i+1)%ga_bins+(i/ga_bins)*ga_bins]*((double)sdf_plot[sdf_plot_nr][2*((i+1)%ga_bins+(i/ga_bins)*ga_bins)]>EPS)+points[(i+1)%ga_bins+(i/ga_bins+1)*ga_bins]*((double)sdf_plot[sdf_plot_nr][2*((i+1)%ga_bins+(i/ga_bins+1)*ga_bins)]>EPS))/norm).V3Str(' ');
				sdf_colors+=", "+((rgba[i+ga_bins]*((double)sdf_plot[sdf_plot_nr][2*(i+ga_bins)]>EPS)+rgba[i]*((double)sdf_plot[sdf_plot_nr][2*i]>EPS)+rgba[(i+1)%ga_bins+(i/ga_bins)*ga_bins]*((double)sdf_plot[sdf_plot_nr][2*((i+1)%ga_bins+(i/ga_bins)*ga_bins)]>EPS)+rgba[(i+1)%ga_bins+(i/ga_bins+1)*ga_bins]*((double)sdf_plot[sdf_plot_nr][2*((i+1)%ga_bins+(i/ga_bins+1)*ga_bins)]>EPS))/norm).V4Str(' ');
			}
			delete[] points;
			delete[] rgba;
			// we got data coordinates+top and bottom average points and colors+alpha
			// -> now onto creating indices of vertices and also midpoints of each four adjacent points (forming a quadrilateral) which we'll break up into four triangles
			// -> first thing to do draw bottom piece (triangles between bottom point and each neighboring two points)
			string sdf_indices="";
			for(unsigned int i=0; i<ga_bins; i++){
				if((sdf_plot[sdf_plot_nr][2*i]>EPS) && (sdf_plot[sdf_plot_nr][2*((i+1)%ga_bins)]>EPS)){
					if(i>0) sdf_indices+=" ";
					sdf_indices+="0 "+int2str(i+1)+" "+int2str((i+1)%ga_bins+1); // data is arranged as slices (varying phi) with varying cos(theta)
				}
			}
			// -> draw everything in between (quadrilateral composed of four triangles)
			for(unsigned int i=0; i<ga_bins*ga_bins-ga_bins; i++){
				// points of the quadrilateral:
				// - (B) i		(C) (i+1)%ga_bins+(i/ga_bins)*ga_bins
				//             (E) ga_bins*ga_bins+i+1
				// - (A) i+ga_bins	(D) (i+1)%ga_bins+(i/ga_bins+1)*ga_bins
				// triangle ABE
				if((sdf_plot[sdf_plot_nr][2*(i+ga_bins)]>EPS) && (sdf_plot[sdf_plot_nr][2*(i)]>EPS))
					sdf_indices+=" "+int2str(i+ga_bins+1)+" "+int2str(i+1)+" "+int2str(ga_bins*ga_bins+i+2);
				// BCE
				if((sdf_plot[sdf_plot_nr][2*(i)]>EPS) && (sdf_plot[sdf_plot_nr][2*((i+1)%ga_bins+(i/ga_bins)*ga_bins)]>EPS))
					sdf_indices+=" "+int2str(i+1)+" "+int2str((i+1)%ga_bins+(i/ga_bins)*ga_bins+1)+" "+int2str(ga_bins*ga_bins+i+2);
				// CDE
				if((sdf_plot[sdf_plot_nr][2*((i+1)%ga_bins+(i/ga_bins)*ga_bins)]>EPS) && (sdf_plot[sdf_plot_nr][2*((i+1)%ga_bins+(i/ga_bins+1)*ga_bins)]>EPS))
					sdf_indices+=" "+int2str((i+1)%ga_bins+(i/ga_bins)*ga_bins+1)+" "+int2str((i+1)%ga_bins+(i/ga_bins+1)*ga_bins+1)+" "+int2str(ga_bins*ga_bins+i+2);
				// DAE
				if((sdf_plot[sdf_plot_nr][2*((i+1)%ga_bins+(i/ga_bins+1)*ga_bins)]>EPS) && (sdf_plot[sdf_plot_nr][2*(i+ga_bins)]>EPS))
					sdf_indices+=" "+int2str((i+1)%ga_bins+(i/ga_bins+1)*ga_bins+1)+" "+int2str(i+ga_bins+1)+" "+int2str(ga_bins*ga_bins+i+2);
			}
			// -> draw top piece
			for(unsigned int i=ga_bins*ga_bins-ga_bins; i<ga_bins*ga_bins; i++){
				if((sdf_plot[sdf_plot_nr][2*i]>EPS) && (sdf_plot[sdf_plot_nr][2*((i+1)%ga_bins+ga_bins*ga_bins-ga_bins)]>EPS))
					sdf_indices+=" "+int2str(i+1)+" "+int2str((i+1)%ga_bins+ga_bins*ga_bins-ga_bins+1)+" "+int2str(ga_bins*ga_bins+1); // data is arranged as slices (varying phi) with varying cos(theta)
			}
			vizout << "\t\t<Transform>\n";
			vizout << "\t\t\t<Shape>\n";
			vizout << "\t\t\t\t<IndexedTriangleSet index='" << sdf_indices << "' solid='false'>\n";
//			vizout << "\t\t\t\t<NurbsPatchSurface uDimension='" << ga_bins+1 << "' vDimension='" << ga_bins+2 << "' uClosed='true' vClosed='false' solid='false'>\n";
			vizout << "\t\t\t\t\t<Coordinate point='" << sdf_coords << "'/>\n";
//			vizout << "\t\t\t\t\t<Coordinate containerField='controlPoint' point='" << sdf_coords << "'/>\n";
			vizout << "\t\t\t\t\t<ColorRGBA color='" << sdf_colors << "'/>\n";
			vizout << "\t\t\t\t</IndexedTriangleSet>\n";
//			vizout << "\t\t\t\t</NurbsPatchSurface>\n";
//			vizout << "\t\t\t\t<Appearance>\n";
//			vizout << "\t\t\t\t\t<Material diffuseColor='0.5 0.5 0.5' shininess='0.8' transparency='0'/>\n";
//			vizout << "\t\t\t\t</Appearance>\n";
			vizout << "\t\t\t</Shape>\n";
			vizout << "\t\t</Transform>\n";
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
			cout << "\t\t\t\t<- SDF X3D file succesfully created.\n";
			for(unsigned int frame=0; frame<cf_params.n_frames+(cf_params.average_frames==1); frame++){
				if(frame%cf_params.average_frames==0){
					delete[] sdf_plot[frame/cf_params.average_frames];
				}
			}
			delete[] sdf_plot;
		}
	}
	analysis_nr++;
	delete[] gammaalpha_nr;
	delete[] average;
	delete[] average2;
	delete[] Vec3average;
	delete[] Vec3average2;
	delete[] max_dipole;
	delete[] Qorder;
}

// Entry point for correlation function code
int main(int argc, char * const argv[])
{
	cout << "\nRobinson Group trajectory file analysis tool\n";
	cout << "Compiled " << __DATE__ << " (Build " << BUILD << ", " << VERSION << ")\n";
	
	MC_Config config;
	Config_Data* configuration = &config.parameters;
	
	string conffilename="";
	
	// Check to make sure there are enough command line parameters
	if(argc<2){
		cout << "Syntax: " << argv[0] << " <trajectory file> <optional: initial step>\n\n";
		exit(1);
	} else{
		//Trajectory file name
		conffilename=argv[1];
	}
	
	// load configuration from file
#ifndef USE_CMWC4096
	configuration->idum = new __int32_t; // so we don't get segfaults ...
#endif
	config.GetFromFile(conffilename.c_str());
	phys_configuration(configuration);
	configuration->fileout=stringNcopy(configuration->trajectoryfile.c_str(),configuration->trajectoryfile.length()-5);
	
	//Load starting step and make sure it is valid
	readfile conffile;
	conffile.filename=configuration->configfile;
	conffile.directory=configuration->configdir;
	GetDirectory(conffile);
	conffile.file.open((conffile.directory+conffile.filename).c_str(),ios::binary);
	
	unsigned int start_step=configuration->steps-configuration->laststep;
	int step=start_step;
	unsigned int position=0;
	unsigned int nr_analyses=config.GetNrSections(&conffile,"Analysis",position,false);
	position=0;
	bool report=false;
	int report_nr=-1;
	for(unsigned int i=0; i<nr_analyses; i++){
		string type="";
		char* analysis=config.GetExclusiveSection(&conffile,"Analysis",position,type); // subname is <group name, analysis type> *xor* <analysis type>
#if DEBUG_LEVEL>3
		cout << analysis << "\n--- (" << position << ")\n";
#endif
		config.SetParam(start_step,"startstep",analysis,configuration->steps);
		if((int)start_step<step) step=start_step;
		if(!report){
			config.SetParam(report,"report",analysis,false);
			if(report) report_nr=i;
		}
		delete[] analysis;
	}
	start_step=(unsigned int)step;
	unsigned int initialstep=configuration->steps-configuration->laststep;
	if(configuration->use_trajectory){
		if(argc>2){ // optional second parameter: starting step (overrides configuration file starting steps)
			int id;
			bool fail=from_string(id,argv[2]);
			if(!fail) if(id<0) fail=true;
			if(!fail){
				initialstep=(unsigned int)id;
			} else{
				cout << "Starting step needs to be an integer greater or equal to zero.\n";
				exit(1);
			}
		}
		if(initialstep<start_step) start_step=initialstep;
		if(start_step%configuration->grfreq!=0){
			cout << "WARNING: Trajectory file does not contain data for requested starting step. Rounding up to next available step...\n";
			step = start_step+(start_step%configuration->grfreq>1)*configuration->grfreq-(start_step%configuration->grfreq);
		} else step = start_step;
		if(step>(int)configuration->last_step){
			cout << "ERROR: Requested starting step is beyond last step in trajectory file.\n";
			exit(5);
		}
	} else cout << "\n**********************************\n* Program is in calculation mode *\n**********************************\n";
	
	Vec3 rotatedE_direction=Vec3(0.0,0.0,1.0);
	if(configuration->Efield.V3Norm()>EPS) rotatedE_direction=configuration->Efield;
	if((configuration->rotate_Efield_steps>0) && configuration->use_trajectory){
		if(configuration->last_step>=configuration->rotate_Efield_steps){
			Vstore potentials;
			Vec3 cosmeans;
			double Vs, msmoved[2];
			unsigned int accepted[3], tries[3];
			config.GetStatistics(configuration->rotate_Efield_steps,configuration->rotate_Efield_steps+1,&potentials,&rotatedE_direction,&cosmeans,&Vs,msmoved,accepted,tries);
		}
	}
	rotatedE_direction/=rotatedE_direction.V3Norm();
	
	// get step data
	Traj_EP* elements = NULL;
	unsigned int* individual_steps=NULL;
	double* V=NULL;
	double* Xm=NULL;
	double* time=NULL;
	Traj_EP** movements=NULL;
	if((configuration->n_oids>0) && configuration->use_trajectory) elements = new Traj_EP[configuration->n_oids];
	if(elements || !configuration->use_trajectory){
		cout << "\n";
		bool goahead=true;
		if(configuration->use_trajectory) goahead=config.GetElementProperties(step,elements);
		if(goahead){
			unsigned int steps=configuration->last_step/configuration->grfreq+(configuration->last_step%configuration->grfreq>0)+1-step/configuration->grfreq;
			
			char** type_names=new char*[nr_corr_types];
			unsigned int* max_type_nr= new unsigned int[nr_corr_types];
			for(unsigned int i=0; i<nr_corr_types; i++){
				string name=corr_types[i];
				type_names[i]=stringNcopy(name.c_str(),name.length());
				if(i<nr_corr_base_types) max_type_nr[i]=0; else max_type_nr[i]=8;
			}
			// Check for inconsistencies and set item names (including group names, LODs, and Components)
			unsigned int n_element_oids=(configuration->n_oids-configuration->n_group_oids);
			if(!configuration->use_trajectory) n_element_oids=configuration->num_element_types;
			unsigned int item_names_indices=n_element_oids*(nr_additional_subtypes+1);
			if(config.max_nr_components>0) item_names_indices+=nr_additional_subtypes+1;
			if(config.max_levelofdetail>0) item_names_indices+=config.max_levelofdetail*(nr_additional_subtypes+1);
			unsigned int group_indices_start=item_names_indices;
			char** item_names=(char**)malloc(sizeof(char*)*item_names_indices);
			unsigned int* item_nr_elements=(unsigned int*)malloc(sizeof(unsigned int)*item_names_indices);
			unsigned int* group_mapping=NULL;
			string name;
			unsigned int max_group_nr_elements=0;
			for(unsigned int i=0; i<configuration->num_groups; i++){
				if(configuration->groups[i]->number>0 || !configuration->use_trajectory){ // use all groups in calculation mode *or* groups which do exist in trajectory file
					if(configuration->groups[i]->nr_elements>max_group_nr_elements) max_group_nr_elements=configuration->groups[i]->nr_elements;
					item_names_indices+=nr_additional_subtypes+1;
					item_names=(char**)realloc(item_names,sizeof(char*)*item_names_indices);
					if(!item_names){
						cout << "Not enough memory.\n";
						exit(2);
					}
					item_nr_elements=(unsigned int*)realloc(item_nr_elements,sizeof(unsigned int)*item_names_indices);
					if(!item_nr_elements){
						cout << "Not enough memory.\n";
						exit(2);
					}
					group_mapping=(unsigned int*)realloc(group_mapping,sizeof(unsigned int)*(item_names_indices-group_indices_start));
					if(!group_mapping){
						cout << "Not enough memory.\n";
						exit(2);
					}
					group_mapping[item_names_indices-nr_additional_subtypes-1-group_indices_start]=i;
					name=configuration->groups[i]->Type->name;
					item_names[item_names_indices-nr_additional_subtypes-1]=stringNcopy(name.c_str(),name.length());
					item_nr_elements[item_names_indices-nr_additional_subtypes-1]=configuration->groups[i]->nr_elements;
					for(unsigned int j=0; j<nr_additional_subtypes; j++){
						group_mapping[item_names_indices-nr_additional_subtypes-group_indices_start+j]=i;
						name=configuration->groups[i]->Type->name+"."+additional_subtypes[j];
						item_names[item_names_indices-nr_additional_subtypes+j]=stringNcopy(name.c_str(),name.length());
						item_nr_elements[item_names_indices-nr_additional_subtypes+j]=configuration->groups[i]->nr_elements;
					}
//					if((configuration->groups[i]->levelofdetail>0) && (configuration->groups[i]->Type->LOD->groups[0]->number==0)){
					if(configuration->groups[i]->Type->LOD){
						item_names_indices+=nr_additional_subtypes+1;
						if(configuration->groups[i]->Type->LOD->nr_components>0) item_names_indices+=nr_additional_subtypes+1;
						item_names=(char**)realloc(item_names,sizeof(char*)*item_names_indices);
						if(!item_names){
							cout << "Not enough memory.\n";
							exit(2);
						}
						item_nr_elements=(unsigned int*)realloc(item_nr_elements,sizeof(unsigned int)*item_names_indices);
						if(!item_nr_elements){
							cout << "Not enough memory.\n";
							exit(2);
						}
						group_mapping=(unsigned int*)realloc(group_mapping,sizeof(unsigned int)*(item_names_indices-group_indices_start));
						if(!group_mapping){
							cout << "Not enough memory.\n";
							exit(2);
						}
						unsigned int base_nr=0;
						name=configuration->groups[i]->Type->LOD->groups[0]->Type->name;
						while(base_nr<configuration->num_groups){
							if(configuration->groups[base_nr]==configuration->groups[i]->Type->LOD->groups[0]) break;
							base_nr++;
						}
						group_mapping[item_names_indices-(configuration->groups[i]->Type->LOD->nr_components>0)*(nr_additional_subtypes+1)-nr_additional_subtypes-1-group_indices_start]=base_nr;
						item_names[item_names_indices-(configuration->groups[i]->Type->LOD->nr_components>0)*(nr_additional_subtypes+1)-nr_additional_subtypes-1]=stringNcopy(name.c_str(),name.length());
						item_nr_elements[item_names_indices-(configuration->groups[i]->Type->LOD->nr_components>0)*(nr_additional_subtypes+1)-nr_additional_subtypes-1]=configuration->groups[i]->Type->LOD->groups[0]->nr_elements;
						if(configuration->groups[i]->Type->LOD->groups[0]->nr_elements>max_group_nr_elements) max_group_nr_elements=configuration->groups[i]->Type->LOD->groups[0]->nr_elements;
						for(unsigned int j=0; j<nr_additional_subtypes; j++){
							group_mapping[item_names_indices-(configuration->groups[i]->Type->LOD->nr_components>0)*(nr_additional_subtypes+1)-nr_additional_subtypes-group_indices_start+j]=i;
							name=configuration->groups[i]->Type->LOD->groups[0]->Type->name+"."+additional_subtypes[j];
							item_names[item_names_indices-(configuration->groups[i]->Type->LOD->nr_components>0)*(nr_additional_subtypes+1)-nr_additional_subtypes+j]=stringNcopy(name.c_str(),name.length());
							item_nr_elements[item_names_indices-(configuration->groups[i]->Type->LOD->nr_components>0)*(nr_additional_subtypes+1)-nr_additional_subtypes+j]=configuration->groups[i]->Type->LOD->groups[0]->nr_elements;
						}
						if(configuration->groups[i]->Type->LOD->nr_components>0){
							group_mapping[item_names_indices-nr_additional_subtypes-1-group_indices_start]=base_nr;
							name=configuration->groups[i]->Type->LOD->groups[0]->Type->name+".Component";
							item_names[item_names_indices-nr_additional_subtypes-1]=stringNcopy(name.c_str(),name.length());
							item_nr_elements[item_names_indices-nr_additional_subtypes-1]=configuration->groups[i]->Type->LOD->nr_components;
							for(unsigned int j=0; j<nr_additional_subtypes; j++){
								group_mapping[item_names_indices-nr_additional_subtypes-group_indices_start+j]=base_nr;
								name=configuration->groups[i]->Type->LOD->groups[0]->Type->name+".Component"+"."+additional_subtypes[j];
								item_names[item_names_indices-nr_additional_subtypes+j]=stringNcopy(name.c_str(),name.length());
								item_nr_elements[item_names_indices-nr_additional_subtypes+j]=configuration->groups[i]->Type->LOD->nr_components;
							}
						}
					}
				}
			}
			unsigned int additional_types_indices_start=item_names_indices;
			if(nr_additional_types+nr_additional_subtypes>0){
				item_names_indices+=nr_additional_types+nr_additional_subtypes;
				item_names=(char**)realloc(item_names,sizeof(char*)*item_names_indices);
				if(!item_names){
					cout << "Not enough memory.\n";
					exit(2);
				}
				item_nr_elements=(unsigned int*)realloc(item_nr_elements,sizeof(unsigned int)*item_names_indices);
				if(!item_nr_elements){
					cout << "Not enough memory.\n";
					exit(2);
				}
				for(unsigned int j=0; j<nr_additional_subtypes; j++){
					item_names[additional_types_indices_start+j]=stringNcopy(additional_subtypes[j].c_str(),additional_subtypes[j].length());
					item_nr_elements[additional_types_indices_start+j]=max_group_nr_elements;
				}
				for(unsigned int j=0; j<nr_additional_types; j++){
					item_names[additional_types_indices_start+nr_additional_subtypes+j]=stringNcopy(additional_types[j].c_str(),additional_types[j].length());
					item_nr_elements[additional_types_indices_start+nr_additional_subtypes+j]=0;
				}
			}
			if(config.max_nr_components>0){
				name="Component";
				item_names[n_element_oids*(nr_additional_subtypes+1)]=stringNcopy(name.c_str(),name.length());
				item_nr_elements[n_element_oids*(nr_additional_subtypes+1)]=config.max_nr_components;
				for(unsigned int j=0; j<nr_additional_subtypes; j++){
					name="Component."+additional_subtypes[j];
					item_names[n_element_oids*(nr_additional_subtypes+1)+j+1]=stringNcopy(name.c_str(),name.length());
					item_nr_elements[n_element_oids*(nr_additional_subtypes+1)+j+1]=config.max_nr_components;
				}
			}
			if(config.max_levelofdetail>0){
				for(unsigned int i=0; i<config.max_levelofdetail; i++){
					name="LOD"+int2str(i+1);
					item_names[n_element_oids*(nr_additional_subtypes+1)+(config.max_nr_components>0)*(nr_additional_subtypes+1)+i*(nr_additional_subtypes+1)]=stringNcopy(name.c_str(),name.length());
					item_nr_elements[n_element_oids*(nr_additional_subtypes+1)+(config.max_nr_components>0)*(nr_additional_subtypes+1)+i*(nr_additional_subtypes+1)]=config.max_LOD_elements;
					for(unsigned int j=0; j<nr_additional_subtypes; j++){
						name="LOD"+int2str(i+1)+"."+additional_subtypes[j];
						item_names[n_element_oids*(nr_additional_subtypes+1)+(config.max_nr_components>0)*(nr_additional_subtypes+1)+i*(nr_additional_subtypes+1)+j+1]=stringNcopy(name.c_str(),name.length());
						item_nr_elements[n_element_oids*(nr_additional_subtypes+1)+(config.max_nr_components>0)*(nr_additional_subtypes+1)+i*(nr_additional_subtypes+1)+j+1]=config.max_nr_components;
					}
				}
			}
			if(configuration->use_trajectory){
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
						name=configuration->element_types[elements[i].element_type]->name;
						item_names[(i-configuration->n_group_oids)*(nr_additional_subtypes+1)]=stringNcopy(name.c_str(),name.length());
						item_nr_elements[(i-configuration->n_group_oids)*(nr_additional_subtypes+1)]=1;
						for(unsigned int j=0; j<nr_additional_subtypes; j++){
							name=configuration->element_types[elements[i].element_type]->name+"."+additional_subtypes[j];
							item_names[(i-configuration->n_group_oids)*(nr_additional_subtypes+1)+j+1]=stringNcopy(name.c_str(),name.length());
							item_nr_elements[(i-configuration->n_group_oids)*(nr_additional_subtypes+1)+j+1]=1;
						}
					}
				}
			} else{
				for(unsigned int i=0; i<configuration->num_element_types; i++){
					name=configuration->element_types[i]->name;
					item_names[i*(nr_additional_subtypes+1)]=stringNcopy(name.c_str(),name.length());
					item_nr_elements[i*(nr_additional_subtypes+1)]=1;
					for(unsigned int j=0; j<nr_additional_subtypes; j++){
						name=configuration->element_types[i]->name+"."+additional_subtypes[j];
						item_names[i*(nr_additional_subtypes+1)+j+1]=stringNcopy(name.c_str(),name.length());
						item_nr_elements[i*(nr_additional_subtypes+1)+j+1]=1;
					}
				}
			}
#if DEBUG_LEVEL>3
			cout << "-> Recognized item names (" << item_names_indices << ":" << group_indices_start << "):\n";
			for(unsigned int i=0; i<item_names_indices; i++){
				cout << i+1 << "\t" << item_names[i] << "\n";
				if(i+1==group_indices_start) cout << "---\n";
			}
#endif
			if(configuration->use_trajectory){
				// Load entire requested trajectory in one pass
				individual_steps=new unsigned int[steps];
				V=new double[steps];
				if(configuration->LJwall_calc) Xm=new double[steps];
				time=new double[steps];
				movements=new Traj_EP*[configuration->n_oids];
				if(!config.GetAllElementProperties(movements,individual_steps,V,Xm,time,steps,step)) exit(1);
			}
			position=0;
			cout << "-> " << nr_analyses << " data ";
			if(nr_analyses==1) cout << "analysis"; else cout << "analyses";
			cout << " found in configuration file:\n";
			cfparams_t* cf_params=new cfparams_t[nr_analyses];
			string type, groupname;
			unsigned int analysis_nr=0;
			int groupnr;
			for(unsigned int i=0; i<nr_analyses; i++){
				cf_params[i].rotatedE_direction=rotatedE_direction;
				cf_params[i].startstep = step;
				cf_params[i].endstep = configuration->last_step;
				cf_params[i].start_frame_in_data = step/configuration->grfreq;
				cf_params[i].configuration=configuration;
				cf_params[i].individual_output = true; 
				cf_params[i].trajectory=movements;
				cf_params[i].V=V;
				cf_params[i].Xm=Xm;
				cf_params[i].time=time;
				cf_params[i].startframe=cf_params[i].startstep/configuration->grfreq;
				cf_params[i].endframe = configuration->last_step/configuration->grfreq+(configuration->last_step%configuration->grfreq>0)+1;
				cf_params[i].n_frames = cf_params[i].endframe-cf_params[i].startframe;
				cf_params[i].nr_bins = configuration->grinc;
				cf_params[i].report = ((int)i==report_nr);
				type="";
				groupname="";
				groupnr=-1;
				char* analysis=config.GetExclusiveSection(&conffile,"Analysis",position,type); // subname is <group name, analysis type> *xor* <analysis type>
				const char* comma=strrchr(type.c_str(),','); // find last occurence of ","
				if(comma){ // since no commas are allowed in analysis types we found <group name, analysis type>
					groupname=stringNcopy(type.c_str(),comma-type.c_str());
					type=stringNcopy(comma+1,strlen(type.c_str())-(comma-type.c_str())-1);
					// Take care of whitespace (can only be at end of groupname and/or beginning of type)
					while((groupname.length()>0) && ((groupname[groupname.length()-1]==' ') || (groupname[groupname.length()-1]=='\t'))) groupname.erase(groupname.length()-1,1);
					while((type.length()>0) && ((type[0]==' ') || (type[0]=='\t'))) type.erase(0,1);
				}
				// Check if type string is a type or group
				bool found=false;
				cftypes_t basetype;
				basetype=gr;
				if(type!=""){
					for(unsigned int j=0; j<nr_corr_types; j++){
						if(compare_strings(type.c_str(),corr_types[j].c_str())){
							found=true;
							basetype=(cftypes_t)j;
							break;
						}
					}
					if(!found) groupname=type;
				}
				if(groupname!=""){
					found=false;
					for(unsigned int j=0; j<configuration->num_groups; j++){
						if(compare_strings(groupname.c_str(),configuration->groups[j]->Type->name.c_str())){
							groupnr=j;
							found=true;
							break;
						}
					}
					if(!found){
						cout << "ERROR: Specified group <" << groupname << "> could not be found.\n";
						exit(3);
					}
				}
				cout << "\t-> Calculating ";
				if(groupnr>=0) cout << "group-internal"; else cout << "system-wide";
				cout << " correlation function(s)";
				if(groupnr>=0) cout << " for group <" << configuration->groups[groupnr]->Type->name << ">";
				cout << "\n";
				config.SetParam(cf_params[i].nr_bins,"nr_bins",analysis,100);
				config.SetParam(cf_params[i].rmax,"rmax",analysis,-1.0);
				if(configuration->use_trajectory){
					config.SetParam(cf_params[i].startstep,"startstep",analysis,initialstep);
					if(cf_params[i].startstep%configuration->grfreq!=0){
						cout << "WARNING: Trajectory file does not contain data for requested starting step. Rounding up to next available step...\n";
						cf_params[i].startstep = cf_params[i].startstep+(cf_params[i].startstep%configuration->grfreq>1)*configuration->grfreq-(cf_params[i].startstep%configuration->grfreq);
					}
					if((int)cf_params[i].startstep<step){
						cout << "WARNING: Specified starting step " << cf_params[i].startstep << " is too low. Adjusting to " << step << ".\n";
						cf_params[i].startstep=(unsigned int)step;
					}
					config.SetParam(cf_params[i].endstep,"endstep",analysis,configuration->last_step);
					if(cf_params[i].endstep%configuration->grfreq!=0){
						cout << "WARNING: Trajectory file does not contain data for requested end step. Rounding up to next available step...\n";
						cf_params[i].endstep = cf_params[i].endstep+(cf_params[i].endstep%configuration->grfreq>1)*configuration->grfreq-(cf_params[i].endstep%configuration->grfreq);
					}
					if(cf_params[i].endstep>configuration->last_step){
						cout << "WARNING: Specified endstep in configuration is beyond trajectory file stored steps. Adjusting to " << configuration->laststep << ".\n";
						cf_params[i].endstep=configuration->last_step;
					}
					if(cf_params[i].startstep>=cf_params[i].endstep){
						cout << "ERROR: Starting step (" << cf_params[i].startstep << ") is larger or equal to endstep (" << cf_params[i].endstep << ").\n";
						exit(5);
					}
					cf_params[i].startframe=cf_params[i].startstep/configuration->grfreq;
					cf_params[i].endframe=cf_params[i].endstep/configuration->grfreq;
					cf_params[i].n_frames=cf_params[i].endframe-cf_params[i].startframe;
					config.SetParam(cf_params[i].average_frames,"average_frames",analysis,cf_params[i].n_frames);
					config.SetParam(cf_params[i].find_max_item,"find_max_item",analysis,false);
					config.SetParam(cf_params[i].find_max_abs_item,"find_max_abs_item",analysis,false);
					if(configuration->noEfield) cf_params[i].find_max_item=false;
				}
				double* user_types=config.get_flex_tupel("type",analysis,type_names,max_type_nr,nr_corr_types);
				unsigned int count=0;
				unsigned int idx, subnumber;
				cf_params[i].nr_correlations=0;
				cf_params[i].items_involved=2;
				cf_params[i].types=NULL;
				cf_params[i].typenr=NULL;
				while(user_types && user_types[count]>0){
					if(user_types[count]>1){
						cout << "ERROR: Only one correlation function type allowed per analysis, e.g. \"type = {gr|gmu}\".\n";
						exit(6);
					}
					if(user_types[count+1]>0){
						cout << "ERROR: Only allowed correlation function types are: ";
						for(unsigned int j=0; j<nr_corr_types; j++){
							cout << corr_types[j];
							if(j+1<nr_corr_types) cout << ", ";
						}
						cout << ".\n";
						exit(6);
					}
					idx=(unsigned int)floor(-user_types[count+1]);
					if((int)qround(user_types[count+1])!=user_types[count+1]) subnumber=(unsigned int)qround((max_type_nr[idx]+1)*(-user_types[count+1]-(double)idx)); else subnumber=1;
					cf_params[i].nr_correlations++;
					cf_params[i].types=(cftypes_t*)realloc(cf_params[i].types,cf_params[i].nr_correlations*sizeof(cftypes_t));
					cf_params[i].typenr=(unsigned int*)realloc(cf_params[i].typenr,cf_params[i].nr_correlations*sizeof(unsigned int));
					if(!cf_params[i].types || !cf_params[i].typenr){
						cout << "Not enough memory to store correlation function types.\n";
						exit(3);
					}
					cf_params[i].types[cf_params[i].nr_correlations-1]=(cftypes_t)idx;
					if(((cftypes_t)idx==dihedral) && (groupnr<0)){
						cout << "ERROR: Dihedral only allowed for group-internal analyses.\n";
						exit(4);
					}
					cf_params[i].typenr[cf_params[i].nr_correlations-1]=subnumber;
					count+=2;
				}
				if(cf_params[i].nr_correlations==0){
					cf_params[i].nr_correlations=1;
					cf_params[i].types=(cftypes_t*)malloc(sizeof(cftypes_t));
					cf_params[i].types[0]=basetype;
					if(((cftypes_t)basetype==dihedral) && (groupnr<0)){
						cout << "ERROR: Dihedral only allowed for group-internal analyses.\n";
						exit(4);
					}
					cf_params[i].typenr=(unsigned int*)malloc(sizeof(unsigned int));
					cf_params[i].typenr[0]=1;
				}
				// Now that the types are set let's make sure the number of items involved is consistent for all of them
				for(unsigned int j=1; j<cf_params[i].nr_correlations; j++){
					if(corr_types_items_involved[cf_params[i].types[j]]!=corr_types_items_involved[cf_params[i].types[0]]){
						cout << "ERROR: Number of items involved per analysis type needs to be constant:\n";
						for(unsigned int k=0; k<cf_params[i].nr_correlations; k++){
							cout << "\t" << output_corr_types[cf_params[i].types[k]] << ": ";
							if(corr_types_items_involved[cf_params[i].types[k]]>0) cout << corr_types_items_involved[cf_params[i].types[k]]; else cout << "any number of";
							cout << " items";
							if(corr_types_items_involved[cf_params[i].types[k]]<0) cout << " (cannot be mixed with fixed number analyses)";
							cout << "\n";
						}
						exit(3);
					}
				}
				// Load analysis type specific parameters
				cf_params[i].additional_parameters=new double*[cf_params[i].nr_correlations];
				for(unsigned int j=0; j<cf_params[i].nr_correlations; j++){
					switch(cf_params[i].types[j]){
						case sdf:
							cf_params[i].additional_parameters[j]=new double[4];
							config.SetParam(cf_params[i].additional_parameters[j][0],"sdf_rel_min",analysis,0.0);
							config.SetParam(cf_params[i].additional_parameters[j][1],"sdf_rel_max",analysis,1.0);
							config.SetParam(cf_params[i].additional_parameters[j][2],"sdf_abs_min",analysis,0.0);
							config.SetParam(cf_params[i].additional_parameters[j][3],"sdf_abs_max",analysis,-1.0);
							break;
						default:
							cf_params[i].additional_parameters[j]=NULL;
							break;
					}
				}
				cf_params[i].any_number_of_items=(corr_types_items_involved[cf_params[i].types[0]]<0);
				if(cf_params[i].any_number_of_items) cf_params[i].items_involved=0; else cf_params[i].items_involved=corr_types_items_involved[cf_params[i].types[0]];
				double* corr_elements=config.get_flex_tupel("elements",analysis,item_names,item_nr_elements,item_names_indices);
				if(!corr_elements) cout << "WARNING: No elements specified -> Nothing to do here.\n";
				double* weights=config.get_flex_tupel("weights",analysis);
				unsigned int w_count=0;
				count=0;
				while(corr_elements && corr_elements[count]>0){
					if((corr_elements[count]==cf_params[i].items_involved) || (cf_params[i].any_number_of_items)){
						if(!cf_params[i].any_number_of_items){
							if(cf_params[i].items_involved>1) cout << "\t\t-> Correlation between items: "; else cout << "\t\t-> Correlation for item: ";
						} else{
							cout << "\t\t-> Analysis for items: ";
							cf_params[i].items_involved=(unsigned int)corr_elements[count];
						}
						cf_params[i].items=new Corr_Item[cf_params[i].items_involved];
						cf_params[i].weights=new double[cf_params[i].items_involved];
						cf_params[i].item_names=new string[cf_params[i].items_involved];
						cf_params[i].positions=new Vec3*[cf_params[i].items_involved];
						cf_params[i].vectors=new Vec3*[cf_params[i].items_involved];
						cf_params[i].rots=new Mat33*[cf_params[i].items_involved];
						cf_params[i].max_vectors=new Vec3*[cf_params[i].items_involved];
						cf_params[i].total_nr=new unsigned int[cf_params[i].items_involved];
						cf_params[i].group_internal=(groupnr>=0);
						cf_params[i].groupnr=groupnr;
						for(unsigned int j=0; j<(unsigned int)corr_elements[count]; j++){
							cf_params[i].weights[j]=1.0; // default is 1.0
							if(weights && (weights[w_count]>0) && ((int)j<weights[w_count])) cf_params[i].weights[j]=weights[w_count+j+1];
							cf_params[i].items[j].element_idx=-1; // stores subnumber but is -1 if none specified
							cf_params[i].items[j].group_idx=groupnr;
							cf_params[i].items[j].type_idx=-1;
							cf_params[i].items[j].LOD_nr=0;
							cf_params[i].items[j].is_element=false;
							cf_params[i].items[j].is_group=(groupnr>=0);
							cf_params[i].items[j].is_component=false;
							cf_params[i].items[j].is_LOD=false;
							if(corr_elements[count+j+1]>0){
								if(cf_params[i].group_internal){
									if((unsigned int)corr_elements[count+j+1]>configuration->groups[groupnr]->nr_elements){
										cout << "ERROR: Requested element is not in group, group has " << configuration->groups[groupnr]->nr_elements << " elements.\n";
										exit(4);
									}
								} else{
									if((unsigned int)corr_elements[count+j+1]>configuration->n_oids){
										cout << "ERROR: Requested element does not exist, only " << configuration->n_oids << " elements are created.\n";
										exit(4);
									}
								}
								cf_params[i].items[j].element_idx=(unsigned int)corr_elements[count+j+1]-1;
							} else{ // found special item (element, group names, component, LOD, etc. ...)
								bool sub=((int)qround(corr_elements[count+j+1])!=corr_elements[count+j+1]);
								idx=(unsigned int)floor(-corr_elements[count+j+1]);
								subnumber=0;
								if(idx>=item_names_indices){ // this should not happen
									cout << "A weird thing just happened: index limits have been exceeded. Contact your optometrist.\n";
									exit(42);
								}
								cf_params[i].items[j].type_idx=idx%(nr_additional_subtypes+1)-1;
								if(idx<n_element_oids){ // found item is an element type
									cf_params[i].items[j].is_element=true;
									cf_params[i].items[j].element_idx=idx/(nr_additional_subtypes+1);
									if(configuration->use_trajectory) cf_params[i].items[j].element_idx=elements[idx/(nr_additional_subtypes+1)+configuration->n_group_oids].element_type; // element type
								} else{
									if(idx<group_indices_start){
										if((config.max_nr_components>0) && ((idx>=n_element_oids*(nr_additional_subtypes+1)) && (idx<n_element_oids*(nr_additional_subtypes+1)+nr_additional_subtypes+1))){
											cf_params[i].items[j].is_component=true;
										} else{
											cf_params[i].items[j].is_LOD=true;
											cf_params[i].items[j].LOD_nr=config.max_levelofdetail-(group_indices_start-1-idx)/(nr_additional_subtypes+1);
										}
									} else{ // found item belongs to group in trajectory file
										if(idx<additional_types_indices_start){
											cf_params[i].items[j].is_group=true;
											cf_params[i].items[j].group_idx=group_mapping[idx-group_indices_start];
										} else{
											cf_params[i].items[j].type_idx=idx-additional_types_indices_start;
											if(!cf_params[i].group_internal && (idx<additional_types_indices_start+nr_additional_subtypes)){
												cout << "\nERROR: There is no global " << additional_subtypes[idx-additional_types_indices_start] << ".\n";
												exit(4);
											}
											if(cf_params[i].group_internal && (idx>=additional_types_indices_start+nr_additional_subtypes)) cf_params[i].items[j].is_group=false; // additional types are not in a group even for group-internal analysis
										}
									}
								}
								if(sub){
									if(cf_params[i].items[j].type_idx>=(int)nr_additional_subtypes){
										cout << additional_types[cf_params[i].items[j].type_idx-nr_additional_subtypes] << " does not have elements.\n";
										exit(4);
									}
									subnumber=(unsigned int)qround((item_nr_elements[idx]+1)*(-corr_elements[count+j+1]-(double)idx));
								}
								if(cf_params[i].group_internal){ // sanity check (Is component and LOD number consistent with current group (if it is one)?)
									if(((cf_params[i].items[j].type_idx>=0) && ((unsigned int)cf_params[i].items[j].type_idx<nr_additional_subtypes)) && (subnumber>configuration->groups[(unsigned int)groupnr]->nr_elements)){
										cout << "\nERROR: Requested element number out of range. Group <" << configuration->groups[(unsigned int)groupnr]->Type->name << "> only has " << configuration->groups[(unsigned int)groupnr]->nr_elements << " elements.\n";
										exit(4);
									}
									if(cf_params[i].items[j].is_element || (cf_params[i].items[j].is_group && (cf_params[i].items[j].group_idx!=groupnr))){
										cout << "\nERROR: Group-internal correlation function needs to stay within same group.\n";
										exit(4);
									}
									if(cf_params[i].items[j].is_component){
										if(subnumber>configuration->groups[(unsigned int)groupnr]->Type->LOD->nr_components){
											cout << "\nERROR: Requested component number out of range. Group <" << configuration->groups[(unsigned int)groupnr]->Type->name << "> only has " << configuration->groups[(unsigned int)groupnr]->Type->LOD->nr_components << " components.\n";
											exit(4);
										}
									} else{
										if(cf_params[i].items[j].is_LOD){
											if(!configuration->groups[(unsigned int)groupnr]->Type->LOD){
												cout << "\nERROR: Group <" << configuration->groups[(unsigned int)groupnr]->Type->name << "> does not have any LOD levels.";
											}
											if(cf_params[i].items[j].LOD_nr>configuration->groups[(unsigned int)groupnr]->Type->LOD->levels){
												cout << "\nERROR: Requested LOD level " << cf_params[i].items[j].LOD_nr << " does not exist for group <" << configuration->groups[(unsigned int)groupnr]->Type->name << ">.";
												cout << " Only " << configuration->groups[(unsigned int)groupnr]->Type->LOD->levels << " LOD level is specified.\n";
												exit(4);
											}
											if(subnumber>configuration->groups[(unsigned int)groupnr]->Type->LOD->groups[cf_params[i].items[j].LOD_nr]->nr_elements){
												cout << "\nERROR: Requested LOD" << cf_params[i].items[j].LOD_nr << " elements for group <" << configuration->groups[(unsigned int)groupnr]->Type->name << "> does not exist. There are only " << configuration->groups[(unsigned int)groupnr]->Type->LOD->groups[cf_params[i].items[j].LOD_nr]->nr_elements << " ellispoids in this LOD level.\n";
												exit(4);
											}
											// adjust group nr to reflect LOD level
											cf_params[i].items[j].group_idx=configuration->groups[(unsigned int)groupnr]->Type->LOD->groups[cf_params[i].items[j].LOD_nr]->type;
										}
									}
								} else{
									if(cf_params[i].items[j].is_component || cf_params[i].items[j].is_LOD){
										cout << "\nERROR: There is no global ";
										if(cf_params[i].items[j].is_component) cout << "component";
										if(cf_params[i].items[j].is_LOD) cout << "level of detail";
										cout << ".\n";
										exit(4);
									}
									if(cf_params[i].items[j].is_group){ // If global entry is a group, check if it's a component or LOD
										char* pos=find_string(item_names[idx],".Component");
										if(pos) cf_params[i].items[j].is_component=true;
										pos=find_string(item_names[idx],".LOD");
										if(pos){
											cf_params[i].items[j].is_LOD=true;
											char* nr=stringNcopy(pos+4,strlen(item_names[idx])-(unsigned int)(pos-item_names[idx])-4);
											if(from_string(cf_params[i].items[j].LOD_nr,nr)){ // This should not happen (since code generate item_names)
												cout << "\nERROR: cannot determine LOD nr for <" << item_names[idx] << ">.\n";
												exit(4);
											}
											delete[] nr;
										}
									}
								}
								if(sub) if(!cf_params[i].items[j].is_element) cf_params[i].items[j].element_idx=subnumber-1;
								if(!sub && cf_params[i].items[j].is_component){
									cout << "\nERROR: No component nr specified.\n";
									exit(4);
								}
							}
							cf_params[i].item_names[j]="";
							if(((cf_params[i].items[j].type_idx>=0) && ((unsigned int)cf_params[i].items[j].type_idx<nr_additional_subtypes)) || (cf_params[i].items[j].type_idx<0)){
								if(cf_params[i].group_internal){
									if(cf_params[i].items[j].is_component){
										cf_params[i].item_names[j]+="Component:"+int2str(cf_params[i].items[j].element_idx+1);
									} else{
										if(cf_params[i].items[j].is_LOD){
											cf_params[i].item_names[j]+="LOD"+int2str(cf_params[i].items[j].LOD_nr)+":"+int2str(cf_params[i].items[j].element_idx+1);
										} else{
											cf_params[i].item_names[j]+=configuration->group_elements[configuration->groups[groupnr]->elements[cf_params[i].items[j].element_idx]].MyType->name;
										}
									}
								} else{
									if(cf_params[i].items[j].is_element){
										cf_params[i].item_names[j]+=configuration->element_types[cf_params[i].items[j].element_idx]->name;
									} else{
										if(cf_params[i].items[j].is_group){
											if(cf_params[i].items[j].is_LOD){
												cf_params[i].item_names[j]+=configuration->groups[cf_params[i].items[j].group_idx]->Type->LOD->groups[0]->Type->name;
												cf_params[i].item_names[j]+=".LOD"+int2str(cf_params[i].items[j].LOD_nr);
											} else cf_params[i].item_names[j]+=configuration->groups[cf_params[i].items[j].group_idx]->Type->name;
											if(cf_params[i].items[j].is_component) cf_params[i].item_names[j]+=".Component";
											if(cf_params[i].items[j].element_idx>=0) cf_params[i].item_names[j]+=":"+int2str(cf_params[i].items[j].element_idx+1);
										} else{
											cf_params[i].item_names[j]+=configuration->element_types[cf_params[i].items[j].element_idx]->name;
											if(configuration->use_trajectory) cf_params[i].item_names[j]+=configuration->element_types[elements[cf_params[i].items[j].element_idx].element_type]->name;
										}
									}
								}
							}
							if(cf_params[i].items[j].type_idx>=0){
								if((unsigned int)cf_params[i].items[j].type_idx<nr_additional_subtypes){
									cf_params[i].item_names[j]+=" "+additional_subtypes[cf_params[i].items[j].type_idx];
								} else{
									cf_params[i].item_names[j]=additional_types[cf_params[i].items[j].type_idx-nr_additional_subtypes];
								}
							}
							cf_params[i].positions[j]=NULL;
							cf_params[i].vectors[j]=NULL;
							cf_params[i].rots[j]=NULL;
							cf_params[i].max_vectors[j]=NULL;
							cf_params[i].total_nr[j]=0;
							
							cout << cf_params[i].item_names[j];
							if(j+1<(unsigned int)corr_elements[count]) cout << ", ";
						}
						cout << "\n";
						
						// Calculate correlation function
						if(configuration->use_trajectory){
							correlation(cf_params[i],analysis_nr);
						} else{
							bool data_analysis=true;
							for(unsigned int j=0; j<cf_params[i].nr_correlations; j++){
								switch(cf_params[i].types[j]){
									case VLJ:
										data_analysis=false;
										break;
									default: break;
								}
							}
							if(data_analysis) cout << "\t\t\t-> Ignoring trajectory file data analysis.\n";
						}
						
						delete[] cf_params[i].items;
						delete[] cf_params[i].item_names;
						for(unsigned int j=0; j<(unsigned int)corr_elements[count]; j++){
							if(cf_params[i].positions[j]) free(cf_params[i].positions[j]);
							if(cf_params[i].vectors[j]) free(cf_params[i].vectors[j]);
							if(cf_params[i].rots[j]) free(cf_params[i].rots[j]);
							if(cf_params[i].max_vectors[j]) free(cf_params[i].max_vectors[j]);
						}
						delete[] cf_params[i].positions;
						delete[] cf_params[i].vectors;
						delete[] cf_params[i].max_vectors;
						delete[] cf_params[i].total_nr;
						count+=(unsigned int)corr_elements[count]+1;
						if(weights && weights[w_count]>0) w_count+=(unsigned int)weights[w_count]+1;
					} else{
						cout << "Requested correlation function is between " << cf_params[i].items_involved << " items but " << (unsigned int)corr_elements[count] << " are provided.\n";
						exit(4);
					}
				}
				delete[] analysis;
				delete[] corr_elements;
				delete[] weights;
			}
			conffile.file.close();
			free(item_names);
			free(item_nr_elements);
			free(group_mapping);
			
			// Clean up
			if(movements){
				for(unsigned int j=0; j<configuration->n_oids; j++) delete[] movements[j];
				delete[] movements;
			}
			if(individual_steps) delete[] individual_steps;
			if(V) delete[] V;
			if(Xm) delete[] Xm;
			if(time) delete[] time;
			for(unsigned int i=0; i<nr_analyses; i++){
				delete[] cf_params[i].types;
				delete[] cf_params[i].typenr;
				for(unsigned int j=0; j<cf_params[i].nr_correlations; j++) if(cf_params[i].additional_parameters[j]) delete[] cf_params[i].additional_parameters[j];
				delete[] cf_params[i].additional_parameters;
				if(configuration->use_trajectory) delete[] cf_params[i].results;
			}
			delete[] cf_params;
		} else{
			cout << "Error getting element properties.\n";
			exit(1);
		}
		if(elements) delete[] elements;
	} else{
		cout << "No elements specified in configuration file => No output.\n";
		exit(42);
	}
	
	cout << "<- Done.\n";
	return 0;
}
