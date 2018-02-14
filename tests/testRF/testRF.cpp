/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*
testfmm.cpp
	Test fast multipole method implementation
	written by Andreas Tillack Jan 22, 2011
*/

#define n2 1.7*1.7
#define kT 300*1.38e-4

#include <iostream>
#include <string.h>
#include <sstream>
#include <fstream>
#include <cstdlib>

#include "../../ScalarMat.h"
#include "../../VecMat.h"

template <class T> bool from_string(T &t, const char* s);
template <class T> bool from_string(T &t, const char* s){
	std::stringstream stream(s);
	return (stream >> t).fail();
}

using namespace std;

inline void apply_PBCs(double bl, double ibl, Vec3 &rvec)
{
	// Rectangular box PBCs
	rvec.vec[0] -= bl*specialround(rvec.vec[0]*ibl);
	rvec.vec[1] -= bl*specialround(rvec.vec[1]*ibl);
	rvec.vec[2] -= bl*specialround(rvec.vec[2]*ibl);
}

inline double binomial(unsigned int n, unsigned int k)
{
	if(k>n) return 0; // sanity check
	unsigned int l=n-k;
	unsigned int m=k;
	if(k>n-k){
		l=k;
		m=n-k;
	}
	double num=1.0;
	double den=1.0;
	//                             { n-k>=k: (n-k+1)*..*n/(1*..*k) -> n-(n-k+1) = k-1 terms
	// binomial = n!/(k!*(n-k)!) = {
	//                             { k>n-k: (k+1)*..*n/(1*..*(n-k)) -> n-(k+1) = n-k-1 terms
	for(unsigned int i=l+1; i<=n; i++){
		num*=i; // l+1*..*n
		        //                    { n-k>=k: k*..*(k-(n-(n-k)-1)) = k*..*1
		den*=m; // m*..*(m-(n-l-1)) = {
		m--;    //                    { k>n-k: (n-k)*..*((n-k)-(n-k-1)) = (n-k)*..*1
	}
	return num/den;
}

double sphericalK(unsigned int l, int m)
{
	if(abs(m)<=l){
		double result=(2*l+1)/(4.0*pi);
		unsigned int a=l-m;
		unsigned int b=l+m;
		// calculate (l-m)!/(l+m)! = a!/b!
		if(a>=b){
			for(unsigned int i=b+1; i<=a; i++) result*=(double)i;
		} else for(unsigned int i=a+1; i<=b; i++) result/=(double)i;
		return sqrt(result);
	} else return 0.0;
}

int main(int argc, char* argv[])
{
	cout << "*** testRF.cpp ***\n";
#ifdef USE_CMWC4096
	init_CMWC4096(8);
	zigset();
#else
	__int32_t idum=-8;
#endif
	unsigned int N=10000;
	// Check if there is a command line parameter
	if(argc>1){
		int id;
		bool fail=from_string(id,argv[1]);
		if(!fail) if(id<1) fail=true;
		if(!fail){
			N=(unsigned int)id;
		} else{
			cout << "Command line parameter (N) needs to be an integer greater than zero.\n";
			exit(1);
		}
	}
	if(N>1E6){ // don't allow more than 1 million charges
		cout << "WARNING: Only upto one million charges are allowed, changing N accordingly.\n";
		N=1E6;
	}
	// boxlength and inverse
	double bl=(double)N;
	double ibl=1.0/bl;
	// limit number of loops (to increase timing accuracy) to no more than one billion evaluations
	unsigned int loops=1E9/(N*N);
	if(loops==0) loops=1;
	cout << "-> Creating " << N << " charges between -0.5e and 0.5e randomly placed\n";
	Vec4* charges=new Vec4[N];
	double qsum=0.0;
	// Create N charges in a box of N^3 Angstrom size
	for(unsigned int i=0; i<N; i++){
#ifdef USE_CMWC4096
		charges[i].vec[0]=(ranQ()-0.5)*bl;
		charges[i].vec[1]=(ranQ()-0.5)*bl;
		charges[i].vec[2]=(ranQ()-0.5)*bl;
		charges[i].vec[3]=(ranQ()-0.5)*e_in_esu; // charges is between -0.5 and 0.5 e
#else
		charges[i].vec[0]=(ran2(idum)-0.5)*bl;
		charges[i].vec[1]=(ran2(idum)-0.5)*bl;
		charges[i].vec[2]=(ran2(idum)-0.5)*bl;
		charges[i].vec[3]=(ran2(idum)-0.5)*e_in_esu; // charges is between -0.5 and 0.5 e
#endif
		qsum+=charges[i].vec[3];
	}
	qsum*=-1.0;
	qsum/=(double)N;
	cout << "\t-> Enforcing charge neutrality in box by adding " << qsum/e_in_esu << "e to each charge\n";
	Vec3 M(0.0); // calculate M while there
	for(unsigned int i=0; i<N; i++){
		charges[i].vec[3]+=qsum;
		M.vec[0]+=charges[i].vec[3]*charges[i].vec[0];
		M.vec[1]+=charges[i].vec[3]*charges[i].vec[1];
		M.vec[2]+=charges[i].vec[3]*charges[i].vec[2];
	}
	cout << "<- Done.\n";
	cout << "-> Calculating electrostatic potential over whole box using direct summation.\n";
	double tstart = clock();
	double direct_sum;
	double q;
	Vec3 rmu;
	double phiq;
	for(unsigned int i=0; i<loops; i++){
		direct_sum=0.0;
		for(unsigned int j=0; j<N-1; j++){
			phiq=0.0;
			for(unsigned int k=j+1; k<N; k++){
				//Calculate distances
				rmu=Vec3(charges[k].vec[0]-charges[j].vec[0],charges[k].vec[1]-charges[j].vec[1],charges[k].vec[2]-charges[j].vec[2]); // rmu is vector from element j to element k
				apply_PBCs(bl,ibl,rmu);
				q=charges[k].vec[3];
				phiq+=q/rmu.V3Norm();
			}
			// direct potential sum
			direct_sum += charges[j].vec[3]*phiq;
		}
		break;
	}
	direct_sum/=n2;
	double tend = clock();
	// 1/4pie0 e^2/10^-10 m = 1/4pie0 2.5664 10^-28 J = 2.5664 10^-28/1.1126 10^-10 J = 2.3066 10^-18 J = 23.066 perg
	// (4.803 10^-10 esu/e*e)^2/10^-10 m = 1/4pie0 (4.803 10^-10/2997924580)^2/10^-10 m = 1/4pie0 e^2/10^-10 m = 23.066 perg
	// (4.803)^2/1 = 23.066 perg
	cout << "\t-> Electrostatic potential over whole box: " << direct_sum << " perg\n";
	cout << "<- Done, took " << (tend-tstart)/CLOCKS_PER_SEC*1000/loops << " ms on average to calculate.\n\n";
	double epsRF=n2*1.1;
/*	double M2 = M*M;
	cout << "-> Calculating epsRF from total dipole moment (" << M.V3Str(',') << ") D: ";
	// Calculate q = yg
	q = 4*pi*M2/(9*N*N*N*kT);
	// eps_RF = 1/4*(n^2+9*q +/- 3*sqrt(n^4+2*n^2*q+9*q^2))
	// - out of +/- sqrt solution only + is solution b/c for q->0 : eps_RF=1/4*(n^2+3*sqrt(n^4)) = 1/4*(4*n^2)
	epsRF = 0.25*n2+2.25*q+0.75*sqrt(n2*n2+2.0*n2*q+9.0*q*q);
	cout << epsRF << "\n";*/
	cout << "-> Calculating electrostatic potential and reaction field with half boxlength cutoff using direct summation.\n";
	// evaluate potential energy directly
	double rcut=bl/2.0;
	double rfcut3=rcut*rcut*rcut;
	double rf_correction = 2.0*(epsRF-n2)/(rfcut3*(2.0*epsRF+n2));
	double cutoff2=rcut*rcut; // half the boxsize is cutoff
	double direct_potential;
	double reaction_field;
	// variable definitions
	double rmu2;
	double RFp;
	
	tstart = clock();
	for(unsigned int i=0; i<loops; i++){
		direct_potential=0.0;
		reaction_field=0.0;
		for(unsigned int j=0; j<N; j++){
			phiq=0.0;
			RFp=0.0;
			double sigma=charges[j].vec[3];
			for(unsigned int k=0; k<N; k++){
				if(j!=k){
					//Calculate distances
					rmu=Vec3(charges[k].vec[0]-charges[j].vec[0],charges[k].vec[1]-charges[j].vec[1],charges[k].vec[2]-charges[j].vec[2]); // rmu is vector from element j to element k
					apply_PBCs(bl,ibl,rmu);
					rmu2=rmu*rmu;
					if(rmu2<cutoff2){
						q=charges[k].vec[3];
						sigma+=q;
						// Add interacting charge to reaction field
						RFp += rmu2*q*0.5;
						phiq+=q/sqrt(rmu2);
					}
				}
			}
			// correct for residual charge
			phiq-=sigma/rcut;
			RFp-=0.5*sigma*rcut*rcut;
			// direct potential sum
			direct_potential += charges[j].vec[3]*phiq;
			// Calculate reaction field contribution
			reaction_field += RFp*charges[j].vec[3];
		}
		break;
	}
	direct_potential/=n2*2.0;
	reaction_field*=rf_correction*0.5;
	tend = clock();
	cout << "\t-> Nearfield potential sum: " << direct_potential << " perg\n";
	cout << "\t-> Reaction field direct sum: " << reaction_field << " perg\n";
	cout << "\t-> Electrostatics potential overall: " << direct_potential+reaction_field << " perg\n";
	cout << "<- Done, took " << (tend-tstart)/CLOCKS_PER_SEC*1000/loops << " ms on average to calculate.\n\n";
	delete[] charges;
	cout << "*** Finished. ***\n\n";
	return 0;
}
