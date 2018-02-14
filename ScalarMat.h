/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * ScalarMat.h : include file for standard system include files,
 * or project specific include files that are used frequently, but
 * are changed infrequently
 *
 * updated by Andreas Tillack on Jan 28, 2011
 * - moved GetCnsts() and global_Cnsts to Cnsts.h (more logical place for them to be)
 * - added "specialround(A)" function optimized for |A|<1.5 (falls back to qround otherwise)
 * - added "qqpwr" which is slightly faster with certain exponentials (and not slower at others ...)
 */

#ifndef INCLUDED_SCALARMAT
#define INCLUDED_SCALARMAT

//Standard library
#include <cstdio>
#include <ctime>
#include <string>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <immintrin.h>

//Custom Classes
#ifndef INCLUDED_CONFIG
#include "Config.h"
#endif
#ifndef INCLUDED_VECTORMAT
#include "VecMat.h"
#endif

union float_union{
	__uint32_t i32;
	float f32;
};

union double_union{
	__uint64_t i64;
	double d64;
};

// Global functions
// Placeholder linear algebra functions have been removed; their replacements are in VecMat.h LEJ 01/08/10

// The Robin Block
double ran2(__int32_t &cnst_idum); ///< Random number generator

// The Lewis block
double average(double *A, int numel); ///< Calculate the mean of an array of length numel
double average_range(double *A, int first, int last); ///< Calculate the mean of part of an array
void average2x(double **A, double *B, int numelx, int numely); ///< Same as above, but 2D
void average_range2x(double **A, double *B, int firstx, int lastx, int firsty, int lasty); ///< Same as above, but 2D
int totalcounti(int *A, int numel); ///< Sum of an integer array
double totalcountd(double *A, int numel); ///< Sum of a double array
void TensorProduct3(Vec3 &A, Vec3 &B, Mat33 &C); ///< REPLACE - see VecMat.h for details
void Vec3Amult(Vec3 &A, Vec3 &B, Vec3 &C); ///< Multiply two Vec3s as if they are normal arrays, instead of doing a dot product
void Vec3Adiv(Vec3 &A, Vec3 &B, Vec3 &C); ///< Divide two Vec3s as if they are normal arrays
void advanceline(fstream &file, int nline); ///< File read utility for advancing a line
int gettdint(fstream &file); ///< Extracts tab-delimited integers from a line of a file
double gettddoub(fstream &file); ///< Extracts tab-delimited doubles from a line of a file

// The Andreas block
template <typename T> int sgn(T val);
void init_CMWC4096(__uint32_t x);
__uint32_t CMWC4096();
double ranQ();
double ran_n();
double ran_n2();
void zigset();
__uint64_t ran_count();
__int32_t ran2_int(__int32_t &idum);

// The inlined Lewis+Andreas block

template <typename T> int sgn(T val)
{
	return (T(0)<val) - (val<T(0));
}

inline __uint32_t OneBitSet(__uint32_t* bitfield, unsigned int nr)
{
	__uint32_t result=0;
	
	for(unsigned int iii=0; iii<nr; nr++) result|=bitfield[iii];
	
	return result;
}

inline unsigned int BitCount(__uint32_t i)
{
	__uint32_t temp=i-((i>>1) & 0x55555555);
	temp=(temp & 0x33333333) + ((temp >> 2) & 0x33333333);
	return (((temp + (temp >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24; // 0xF0F0F0F == 0x0F0F0F0F ...
}

inline string int2str(int i)
{
	stringstream converter;
	converter << i;
	return converter.str();
}

inline string int2str_fixed(int i, int length)
{
	string s1=int2str(i);
	string s2="";
	for(unsigned int j=0; j<length-s1.length(); j++) s2+=" ";
	return s2+s1;
}

inline string double2str(double d)
{
	stringstream converter;
	converter << d;
	return converter.str();
}

inline string double2str(double d, int precision)
{
	stringstream converter;
	converter.precision(precision);
	converter << d;
	return converter.str();
}

#ifdef __SSE2__
inline double sse_sqrt(double d)
{
	double result;
	__m128d value = _mm_load_sd(&d);
	_mm_store_sd(&result,_mm_sqrt_sd(value,value));
	return result;
}

inline double sse_invsqrt(double d)
{
	double result;
	__m128d value = _mm_load_sd(&d);
	const __m128d one = {1.0,1.0};
	_mm_store_sd(&result,_mm_div_sd(one,_mm_sqrt_sd(value,value)));
	return result;
}
#endif

inline double taylor_sqrt(double d)
{
	double a=d-0.5;
//     1/sqrt(2) * [ 1       +       a - 1/2*a^2 + 1/2*a^3 - 5/8*a^4 + 7/8*a^5
	return 1.0/sqrt(2.0)*(1.0+a*(1.0-0.5*a*(1.0-a*(1.0-0.25*a*(5.0-7.0*a)))));
}

inline double custom_sqrt(double d)
{
#ifdef __SSE2__
	return sse_sqrt(d);
#else
	return sqrt(d);
#endif
}

inline double custom_invsqrt(double d)
{
#ifdef __SSE2__
	return sse_invsqrt(d);
#else
	return 1.0/sqrt(d);
#endif
}

/// the Quake III "magic" invsqr code with Lormont's "magic" numbers
#ifdef __x86_64__
inline double magic_invsqr(double d) // on 64bit machines it is faster to use 64bit doubles
{
	double_union x;
	x.d64 = d;
	double half=0.5*x.d64;
	x.i64 = 0x5fe6ec85e7de30da - (x.i64 >> 1);
	x.d64 = x.d64*(1.5-half*x.d64*x.d64); // Newton-Raphson for 1/x^2=a
	x.d64 = x.d64*(1.5-half*x.d64*x.d64);
	x.d64 = x.d64*(1.5-half*x.d64*x.d64);
	x.d64 = x.d64*(1.5-half*x.d64*x.d64);
	return x.d64;
}
#else
inline float magic_invsqr(double d)
{
	float_union x;
	x.f32 = (float)d; // on 32bit machines it is faster to use 32bit floats
	float half=0.5*x.f32;
	x.i32 = 0x5f375a86 - (x.i32 >> 1);
	x.f32 = x.f32*(1.5-half*x.f32*x.f32); // Newton-Raphson for 1/x^2=a
	x.f32 = x.f32*(1.5-half*x.f32*x.f32);
	x.f32 = x.f32*(1.5-half*x.f32*x.f32);
	x.f32 = x.f32*(1.5-half*x.f32*x.f32);
	return x.f32;
}
#endif

/// Optimized floor algorithm suggested by DTrace
inline double fastfloor(double f)
{
	register double twoTo52 = 4503599627370496.0;
	double c = ( f >= 0.0 ? -twoTo52 : twoTo52 );
	double result = (f - c) + c;
	if( f < result ) result -= 1.0;
	return result;
}

/// Fast rounding algorithm using fastfloor
inline double qround(double A)
{
	double rA = fastfloor(A+0.5);
	return rA;
}

/// Special round for |A|<1.5, fallback to fast rounding otherwise -- AT
inline double specialround(double A){
	if(fabs(A)<1.5){
		if(A>=0.5)  return 1.0;
			else if(A<=-0.5) return -1.0;
	} else return qround(A);
	return 0.0;
}

/// Seems very obvious, right?
inline double pwr3(double a)
{
	return a*a*a;
}

/// Fast power algorithm that only works for positive integers
inline double qpwr(double a, int b)
{
	double c = 1;
	for (int i = 0; i < b; i++) c *= a;
	return c;
}

/// Slighlty faster recursive a^b (b being positive integer) -- AT
/* inline double qqpwr(const double a, const int b)
{
	if(b==0) return 1.0; // anything to the zeroth power (including zero) is one (end of story)
	if (fabs(a)<EPS) return 0.0;
	switch(b){
		case 1: return a;
		case 2: return a*a;
		case 3: return a*a*a;
	}
	double result=qqpwr(a*a,b>>1); // bitshift right by one is integer/2
	if(b & 1) result *= a;
	return result;
}*/

/// Slighlty faster non-recursive a^b (b being positive integer) -- AT
inline double qqpwr(double a, int b)
{
	switch(b){
		case 0: return 1.0;
		case 1: return a;
		case 2: return a*a;
		case 3: return a*a*a;
	}
	double result=1.0;
	while(b){
		if(b & 1) result *= a;
		b >>= 1; // bitshift right by one is integer/2
		a *= a;
	}
	return result;
}

inline double SphereOverlapVolume(double r1, double r2, double d)
{
	if(d<r1+r2){
		// see Weisstein, Eric W. "Sphere-Sphere Intersection.", From Mathworld-A Wolfram Web Resource. http://mathworld.wolfram.com/Sphere-Sphere-Intersection.html
		return pi*qqpwr(r2+r1-d,2)*(d*d+2.0*d*r1-3.0*r1*r1+2.0*d*r2+6*r1*r2-3*r2*r2)/(12.0*d);
	} else{
		if(fabs(d)<=EPS){
			if(r1<r2) return 4.0*pi/3.0*qqpwr(r1,3); else return 4.0*pi/3.0*qqpwr(r2,3);
		} else return 0.0;
	}
}

inline double EllipsoidSurface(double a, double b, double c)
{
// Knud Thomsen's formula
	double p=1.6075; // ~1% error
	double k=0.0942;
	return 4*pi*pow((pow(a*b,p)+pow(a*c,p)+pow(b*c,p))/(3.0-k*(1.0-27.0*a*b*c/qqpwr(a+b+c,3))),1.0/p);
}

/// Maximum value for an array of doubles
inline double maxd(double *A, int l)
{
	double Max = A[0];
	for (int i=1; i<l; i++) if(A[i]>Max) Max = A[i];
	return Max;
}

/// Minimum value for an array of doubles
inline double mind(double *A, int l)
{
	double Min = A[0];
	for (int i=1; i<l; i++) if(A[i]<Min) Min = A[i];
	return Min;
}

#endif

