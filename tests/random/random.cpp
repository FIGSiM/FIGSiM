/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*
testmatrix.cpp
	Test matrix implementation
	written by Andreas Tillack Jan 22, 2011
*/

#define N 1E9
#define bins 1000

#include <iostream>
#include <string.h>
#include "../../ScalarMat.h"

inline string Byte2Hex(unsigned short b)
{
	string result="";
	char hex[16]={'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};
	
	result+=hex[(b>>4)&0xF];
	result+=hex[b&0xF];
	
	return result;
}

template <class T> bool from_string(T &t, const char* s){
	std::stringstream stream(s);
	return (stream >> t).fail();
}

int main(int argc, char * const argv[]) {
	cout << "*** random.cpp ***\n";
	int millions=-1;
	if(argc<2){
		cout << "-> No binary random number output.\n\n";
	} else{
		bool fail=from_string(millions,argv[1]);
		if(!fail) if(millions<0) fail=true;
		if(fail){
			cout << "ERROR: Second command line parameter (millions of random numbers to be outputted) needs to be an integer greater or equal to zero.\n";
			exit(1);
		}
	}
	cout << "-> Creating " << N << " uniformly distributed random numbers and calculating distribution.\n";
	__int32_t idum=-8;
	double tstart = clock();
	unsigned int* distribution=new unsigned int[bins+1];
	double tend = clock();
	double tmem=tend-tstart;
	tstart=clock();
	for(unsigned i=0; i<N; i++){
		double nr=ran2(idum);
		distribution[(unsigned int)floor(nr*bins)]++;
	}
	tend=clock();
	double timeran=tend-tstart-tmem;
	double sum=0.0;
	double sum2=0.0;
	cout << "-> Outputting distribution.\n";
	fstream* outfile;
	string filename="ran2.dat";
	outfile=new fstream(filename.c_str(),ios::out|ios::trunc);
	string outstring="";
	if(outfile->fail()==true){
		cout << "Could not open output file.\n";
		exit(1);
	}
	for(unsigned int i=0; i<bins; i++){
		sum+=(double)distribution[i];
		sum2+=(double)distribution[i]*distribution[i];
		outstring=double2str((i+0.5)/(double)bins)+"\t"+double2str(distribution[i])+"\n";
		outfile->write(outstring.c_str(),outstring.length());
	}
	outfile->close();
	double avg=sum/bins;
	double stddev=sqrt(1.0/bins*(sum2-(sum*sum/bins)));
	cout << "-> Finished, took " << timeran/CLOCKS_PER_SEC << " ns per random value.\n";
	cout << "\t-> Average count per bin: " << avg << " +/- " << stddev << "\n\n";
	
	cout << "-> Creating " << N << " uniformly distributed random numbers (using CMWC4096) and calculating distribution.\n";
	init_CMWC4096(8);
	for(unsigned int i=0; i<bins; i++) distribution[i]=0;
	tstart=clock();
	for(unsigned i=0; i<N; i++){
		double nr=ranQ();
		distribution[(unsigned int)floor(nr*bins)]++;
	}
	tend=clock();
	timeran=tend-tstart-tmem;
	sum=0.0;
	sum2=0.0;
	cout << "-> Outputting distribution.\n";
	filename="CMWC4096.dat";
	outfile=new fstream(filename.c_str(),ios::out|ios::trunc);
	if(outfile->fail()==true){
		cout << "Could not open output file.\n";
		exit(1);
	}
	for(unsigned int i=0; i<bins; i++){
		sum+=(double)distribution[i];
		sum2+=(double)distribution[i]*distribution[i];
		outstring=double2str((i+0.5)/(double)bins)+"\t"+double2str(distribution[i])+"\n";
		outfile->write(outstring.c_str(),outstring.length());
	}
	outfile->close();
	avg=sum/bins;
	stddev=sqrt(1.0/bins*(sum2-(sum*sum/bins)));
	cout << "-> Finished, took " << timeran/CLOCKS_PER_SEC << " ns per random value.\n";
	cout << "\t-> Average count per bin: " << avg << " +/- " << stddev << "\n\n";
	
	cout << "-> Creating " << N << " normal distributed random numbers (Ziggurat using CMWC4096, with tail) and calculating distribution.\n";
	zigset();
	for(unsigned int i=0; i<bins; i++) distribution[i]=0;
	double inv_r=1.0/(2.0*3.442619855899);
	tstart=clock();
	for(unsigned i=0; i<N; i++){
		double nr=0.5+ran_n()*inv_r;
		if(nr<0.0) nr=0.0;
		if(nr>1.0-EPS) nr=1.0-EPS;
		distribution[(unsigned int)floor(nr*bins)]++;
	}
	tend=clock();
	timeran=tend-tstart-tmem;
	cout << "-> Outputting distribution.\n";
	filename="gauss.dat";
	outfile=new fstream(filename.c_str(),ios::out|ios::trunc);
	if(outfile->fail()==true){
		cout << "Could not open output file.\n";
		exit(1);
	}
	for(unsigned int i=0; i<bins; i++){
		outstring=double2str((i+0.5)/(double)bins)+"\t"+double2str(distribution[i])+"\n";
		outfile->write(outstring.c_str(),outstring.length());
	}
	outfile->close();
	cout << "-> Finished, took " << timeran/CLOCKS_PER_SEC << " ns per random value.\n\n";
	
	cout << "-> Creating " << N << " normal distributed random numbers (Ziggurat using CMWC4096, without tail) and calculating distribution.\n";
	for(unsigned int i=0; i<bins; i++) distribution[i]=0;
	tstart=clock();
	for(unsigned i=0; i<N; i++){
		double nr=0.5+ran_n2()*inv_r;
		if(nr<0.0) nr=0.0;
		if(nr>1.0-EPS) nr=1.0-EPS;
		distribution[(unsigned int)floor(nr*bins)]++;
	}
	tend=clock();
	timeran=tend-tstart-tmem;
	cout << "-> Outputting distribution.\n";
	filename="gauss_notail.dat";
	outfile=new fstream(filename.c_str(),ios::out|ios::trunc);
	if(outfile->fail()==true){
		cout << "Could not open output file.\n";
		exit(1);
	}
	for(unsigned int i=0; i<bins; i++){
		outstring=double2str((i+0.5)/(double)bins)+"\t"+double2str(distribution[i])+"\n";
		outfile->write(outstring.c_str(),outstring.length());
	}
	outfile->close();
	cout << "-> Finished, took " << timeran/CLOCKS_PER_SEC << " ns per random value.\n\n";
	
	delete[] distribution;
	
	if(millions>0){
		cout << "Outputting " << millions << " million ran2 random numbers in binary format for statistical analysis.\n";
		filename="ran2.bin";
		outfile=new fstream(filename.c_str(),ios::out|ios::trunc);
		if (outfile->fail()==true){
			cout << "Could not open output file.\n";
			exit(1);
		}
		char* dout=new char[4];
		__uint32_t u;
		for(unsigned int i=0; i<(unsigned int)millions*1E6; i++){
			u=(__uint32_t)(ran2(idum)*4294967296);
			memcpy(dout,&u,4);
			outfile->put(dout[3]);
			outfile->put(dout[2]);
			outfile->put(dout[1]);
			outfile->put(dout[0]);
		}
		outfile->close();
		
		cout << "\nOutputting " << millions << " million CMWC4096 random numbers in binary format for statistical analysis.\n";
		filename="CMWC4096.bin";
		outfile=new fstream(filename.c_str(),ios::out|ios::trunc);
		if(outfile->fail()==true){
			cout << "Could not open output file.\n";
			exit(1);
		}
		for(unsigned int i=0; i<(unsigned int)millions*1E6; i++){
			u=CMWC4096();
			memcpy(dout,&u,4);
			outfile->put(dout[3]);
			outfile->put(dout[2]);
			outfile->put(dout[1]);
			outfile->put(dout[0]);
		}
		outfile->close();
		
		delete[] dout;
	}
	cout << "*** Finished. ***\n\n";
	return 0;
}
