/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 *  VecMat.h
 *  Created by ljohnson on 7/24/09.
 *
 *  Lightweight linear algebra library intended for rapid three-dimensional computation. This library cannot handle
 *	vectors and matrices of any other dimensionality; a more general library should be used if variable or higher
 *	dimensions are required. The library consists of two classes: Vec3 (vectors) and Mat33 (matrices). Vectors are not
 *	tracked as row or column vectors, with behavior instead determined by the function acting on the vector.
 *	Overloaded operators will always assume that vectors are oriented such that their inner dimensions agree.
 *	All double-precision arithmatic operations are overloaded for vector-scalar and matrix-scalar operations, but are not currently
 *	commutative, meaning that the scalar must be placed *AFTER* the vector or matrix in order to avoid an error, as the
 *	opposite arrangement is not yet defined. All legal vector-vector, vector-matrix, and matrix-matrix BLAS operations are overloaded,
 *	and should not have a problem with commutivity. Only the == and != comparison operators are defined for vectors and matrices. All
 *	elementary math functions are inlined and vector-scalar and vector-vector operations were written unrolled. Other functions may be
 *	partially unrolled. The inverse, linear solve, and string conversion functions are not inlined.
 *
 *	IMPORTANT NOTE: Arguments for operators are defined as constant. Expressions of the form B = A*B will result in a segfault.
 *	This may be avoided by using an accumulating operator; e.g. B *= A instead.
 *
 *  Updated on Jan 22, 2011 by Andreas Tillack:
 *  - added complex vector class
 *  - added solver for third order reduced polynomial (using method of Cardano)
 *  - use third order polynomial solver to get eigenvalues of 3x3 matrix
 *  - defined 1/2*cubicroot(3) as constant for third order polynomial solving
 *  - added Vector subtraction from 3x3 Matrix diagonal (to help solve for eigenvektors)
 *  another update on Feb 5, 2011 (AT)
 *  - forward declaration of Mat33 fixes formerly "uncompilable" code, which can now be used
 *
 * Updated on May 17, 2011 by LEJ
 * - added output functions with user-specified delimiters
 * - added Vec4 class; no BLAS functions but has minimal I/O needed for X3D output
 * - added converter between rotation matrices and axis-angle rotations (Vec4), in addition to AT's vec2rot converter
 *
 * Updated on June 10, 2011 by AT
 * - fixed bad logical fallacy for "!=" operator (all classes in here):
 *   - code was 1:1 inversion of "==" operator's code
 *   - returned false when *any* one of the components were equal (which is wrong in most cases ...)
 *
 * Updated on Oct 24, 2011 by AT
 * - added double4 struct for OpenCL
 * - use double4 for quaternions
 *
 * Updated on Feb 28, 2012 by AT
 * - fixed bad bug in V3TensProd: return matrix A started out as unit matrix (needs to be zero)
 *
 */

#ifndef INCLUDED_VECTORMAT
#define INCLUDED_VECTORMAT

#include <complex>
#include <cstdio>
#include <ctime>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>

#ifndef INCLUDED_CONFIG
#include "Config.h"
#endif
#ifndef INCLUDED_MC_OPENCL
#include "MC_OpenCL.h"
#endif

using namespace std;

const double sqrt3 = sqrt(3.0);

class Mat33; // Forward declaration of class Mat33 - AT Feb 5, 2011

/// Vec3 class - dedicated class for 3-dimensional vectors
class Vec3
{
	public:
		double vec[3];
		friend class Mat33;
		
	// Constructors and destructors
		Vec3();
		Vec3(const double, const double, const double);
		Vec3(const double);
		Vec3(const double*); // construct from array
		Vec3(const Vec3 &);
		~Vec3() {};
	
	// Utility functions
		Vec3& operator = (const Vec3 &);
		void V3Zeros();
		void V3Swap(Vec3 &);
		
	//Comparison operators
		bool operator== (const Vec3 &);
		bool operator!= (const Vec3 &);

	//Level 1 BLAS
		Vec3 operator+ (const double);
		Vec3& operator+= (const double);
		Vec3 operator- (const double);
		Vec3& operator-= (const double);
		Vec3 operator* (const double);
		Vec3& operator*= (const double);
		Vec3 operator/ (const double);
		Vec3& operator/= (const double);
		
		Vec3 operator+ (const Vec3 &);
		Vec3& operator+= (const Vec3 &);
		Vec3 operator- (const Vec3 &);
		Vec3& operator-= (const Vec3 &);
		double operator* (const Vec3 &) const;
	
	//Level 2 BLAS
		//These do not compile properly. Ask someone who knows C++ better than I do. - LEJ 01/07/10
		// -- solved by AT Feb 5, 2011 (needed forward declaration of Mat33 class)
		Vec3 operator* (const Mat33 &);
		Vec3& operator*= (const Mat33 &);
		
	//Special functions
		Vec3 V3Cross(const Vec3 &);
		Mat33 V3TensProd(const Vec3 &); //Does not compile properly. Ask someone who knows C++ better than I do. - LEJ 01/07/10 - does now (AT, Feb 5, 2011)
		double V3Norm() const;
		double V3Sum() const;
		double V3Prod() const;
		std::string V3Str() const;
		std::string V3Str(const char &) const; 
};

/// complex class - 3d complex vector class, only minimal function set defined needed for eigenvalue problem
class CVec3{
	public:
		complex<double> cvec[3];
	// Constructors and destructor
		CVec3();
		CVec3(const complex<double>, const complex<double>, const complex<double>);
		CVec3(const Vec3 &);
		CVec3(const CVec3 &);
		~CVec3(){};
	// Operators
		bool operator== (const CVec3 &);
		bool operator!= (const CVec3 &);
		bool operator== (const Vec3 &);
		bool operator!= (const Vec3 &);
		CVec3& operator = (const CVec3 &);
		CVec3 operator- (const double); 
		CVec3& operator-= (const double);
	// Functions
		Vec3 Re();
		Vec3 Im();
		Vec3 Abs();
		std::string CV3Str() const;
		std::string CV3Str(const char &) const;
};

///Four-element vector (needed for X3D rotations)
class Vec4 {
	public:
		double vec[4];
		friend class Vec3;
		friend class Mat33;
	//Constructors and destructors
		Vec4();
		Vec4(const double, const double, const double, const double);
		Vec4(const Vec3 &, const double);
		Vec4(const Vec4 &);
		~Vec4() {};
	//Operators
		bool operator== (const Vec4 &);
		bool operator!= (const Vec4 &);
		Vec4 operator+ (const Vec4 &);
		Vec4 operator* (const double);
		Vec4 operator/ (const double);
	//Functions
		std::string V4Str() const;
		std::string V4Str(const char &) const; 
};

/// Mat33 class - for when you need to do transformations on 3-vectors…
class Mat33
{
	public:
		double mat[3][3];
		friend class Vec3;
		
		//Constructors and destructors
		inline Mat33();
		Mat33(const Mat33 &);
		Mat33(const Vec3 &, const Vec3 &);
		Mat33(const double);
		Mat33(const double a, const double b, const double c); // diagonal constructor
		Mat33(const double*);
		~Mat33(){};
		
		//Utility functions
		Mat33& operator= (const Mat33 &);
		void M3Eye();
		inline void M3Zeros();
		
		//Comparison operators
		bool operator== (const Mat33 &);
		bool operator!= (const Mat33 &);
		
		//Matrix-scalar operations
		Mat33 operator+ (const double);
		Mat33& operator+= (const double);
		Mat33 operator- (const double);
		Mat33& operator-= (const double);
		Mat33 operator* (const double);
		Mat33& operator*= (const double);
		Mat33 operator/ (const double);
		Mat33& operator/= (const double);
		
		//Level 2 BLAS
		Vec3 operator* (const Vec3 &);
		//No *= operator due to datatype mismatch.
		Mat33 operator- (const Vec3 &); // subtract vector from diagonal
		Mat33& operator-= (const Vec3 &); // subtract vector from diagonal
		
		//Level 3 BLAS
		Mat33 operator+ (const Mat33 &);
		Mat33& operator+= (const Mat33 &);
		Mat33 operator- (const Mat33 &);
		Mat33& operator -= (const Mat33 &);
		Mat33 operator* (const Mat33 &);
		Mat33& operator*= (const Mat33 &);
		//These do not compile properly. Ask someone who knows C++ better than I do. - LEJ 01/07/10
		Mat33 operator/ (const Mat33 &);
		Mat33& operator/= (const Mat33 &);
		
		//Special functions
		CVec3 Eigenvalues();
		Mat33 Eigenvectors(Vec3 &ew){ bool* multiples=new bool[3]; Mat33 A=Eigenvectors(ew,multiples,false); delete[] multiples; return A; };
		Mat33 Eigenvectors(Vec3 &ew, bool normalize){ bool* multiples=new bool[3]; Mat33 A=Eigenvectors(ew,multiples,normalize); delete[] multiples; return A; };
		Mat33 Eigenvectors(Vec3 &ew, bool* multiples, bool normalize);
		inline Mat33 M3MulDiag(const double a, const double b, const double c);
		Mat33 M3Transpose();
		Mat33 M3RowSwap(const int, const int);
		inline Vec3 ColumnVec3(const int);
		inline Mat33 M3Inv(bool &sing);
		inline Mat33 MulSymM3(const Mat33 &C);
		inline Mat33 TransMulM3(const Mat33 &C);
		inline Vec3 TransMulVec(const Vec3 &u);
		inline Mat33 SymMulM3(const Mat33 &C);
		inline Mat33 SymMulSymM3(const Mat33 &C);
		inline Mat33 SymM3Inv(bool &sing);
		inline Vec3 SymM3InvMult(bool &sing, Vec3 &z);
		double M3Trace();
		inline double M3Det();
		Vec3 M3Diag();
		void M3GEPP(Vec3 &x);
		Vec3 M3LinSolve(Vec3 &u,bool &multiples);
		void LUDecomposition();
		bool M3BackSub(Vec3 &x); // returns true if multiples of vector are also solutions
		std::string M3Str();
		std::string M3RowStr(const int);
		std::string M3RowStr(const int, const char &);
};

/*	Function implementations are below. Math functions and constructors are inlined for speed.
 *	Non-math functions are kept in .cpp file. The vast majority of the class is contained in
 *	This file, organized in the same order that the prototypes are listed
 */

/// get solutions to third order polynomial in reduced form: z^3 + p*z + q = 0
CVec3 SolvePolynomial3(const double p, const double q);


/*	***Vec3***	*/

//Constructors and Destructors

/// Default constructor (zeros)
inline Vec3::Vec3()
{
	vec[0] = 0.0; vec[1] = 0.0; vec[2] = 0.0;
}

/// Construct from doubles
inline Vec3::Vec3(const double a, const double b, const double c)
{
	vec[0] = a;
	vec[1] = b;
	vec[2] = c;
}

inline Vec3::Vec3(const double d)
{
	vec[0] = d; vec[1] = d; vec[2] = d;
}
/// Construct from an array of 3 doubles
inline Vec3::Vec3(const double* A)
{
	vec[0] = A[0];
	vec[1] = A[1];
	vec[2] = A[2];
}

/// Copy constructor
inline Vec3::Vec3(const Vec3 &u)
{
	vec[0] = u.vec[0];
	vec[1] = u.vec[1];
	vec[2] = u.vec[2];
}


//Non-BLAS utility functions

/// Zero vector v = [0 0 0]
inline void Vec3::V3Zeros()
{
	vec[0] = 0.0;
	vec[1] = 0.0;
	vec[2] = 0.0;
}

/// Assignment operator
inline Vec3& Vec3::operator=(const Vec3 &u)
{
	vec[0] = u.vec[0];
	vec[1] = u.vec[1];
	vec[2] = u.vec[2];
	return *this;
}

/// Swap values in two vectors u <-> v
inline void Vec3::V3Swap(Vec3 &u)
{
	double temp;
	for (int i = 0; i < 3; i ++)
	{
		temp = vec[i];
		vec[i] = u.vec[i];
		u.vec[i] = temp;
	}
}


// Comparison operators
inline bool Vec3::operator==(const Vec3 &u)
{
	if (fabs(u.vec[0]-vec[0])>EPS) return false;
	if (fabs(u.vec[1]-vec[1])>EPS) return false;
	if (fabs(u.vec[2]-vec[2])>EPS) return false;
	return true;
}

inline bool Vec3::operator!=(const Vec3 &u)
{
	if (fabs(u.vec[0]-vec[0])>EPS) return true;
	if (fabs(u.vec[1]-vec[1])>EPS) return true;
	if (fabs(u.vec[2]-vec[2])>EPS) return true;
	return false;
}


//Level 1 BLAS

/// Vector scalar addition v = u + a
inline Vec3 Vec3::operator+ (const double a)
{
	Vec3 v(a+vec[0], a+vec[1], a+vec[2]);
	return v;
}

/// Vector scalar addition u = u + a
inline Vec3& Vec3::operator += (const double a)
{
	vec[0] +=a;
	vec[1] +=a;
	vec[2] +=a;
	return *this;
}

/// Vector scalar subtraction v = u - a
inline Vec3 Vec3::operator- (const double a)
{
	Vec3 v(vec[0]-a, vec[1]-a, vec[2]-a);
	return v;
}

/// Vector scalar subtraction u = u - a
inline Vec3& Vec3::operator -= (const double a)
{
	vec[0] -=a;
	vec[1] -=a;
	vec[2] -=a;
	return *this;
}

/// Vector scalar multiplication v = u*a
inline Vec3 Vec3::operator* (const double a)
{
	Vec3 v(a*vec[0], a*vec[1], a*vec[2]);
	return v;
}


/// Vector scalar multiplication u = u*a
inline Vec3& Vec3::operator*= (const double a)
{
	vec[0] *=a;
	vec[1] *=a;
	vec[2] *=a;
	return *this;
}

/// Vector scalar division v = u/a
inline Vec3 Vec3::operator/ (const double a)
{
	double b = 1.0/a;
	Vec3 v(vec[0]*b, vec[1]*b, vec[2]*b);
	return v;
}

/// Vector scalar division u = u/a
inline Vec3& Vec3::operator/= (const double a)
{
	double b = 1.0/a;
	vec[0] *=b;
	vec[1] *=b;
	vec[2] *=b;
	return *this;
}

/// Vector addition v = u + w
inline Vec3 Vec3::operator+ (const Vec3 &w)
{
	Vec3 v(vec[0]+w.vec[0], vec[1]+w.vec[1], vec[2] +w.vec[2]);
	return v;
}

/// Vector addition u = u + w
inline Vec3& Vec3::operator+= (const Vec3 &w)
{
	vec[0] += w.vec[0];
	vec[1] += w.vec[1];
	vec[2] += w.vec[2];
	return *this;
}

/// Vector subtraction v = u - w
inline Vec3 Vec3::operator- (const Vec3 &w)
{
	Vec3 v(vec[0]-w.vec[0], vec[1]-w.vec[1], vec[2] -w.vec[2]);
	return v;
}

/// Vector subtraction u = u - w
inline Vec3& Vec3::operator-= (const Vec3 &w)
{
	vec[0] -= w.vec[0];
	vec[1] -= w.vec[1];
	vec[2] -= w.vec[2];
	return *this;
}

/// Dot product a = v'*u
inline double Vec3::operator* (const Vec3 &B) const
{
	double a = vec[0]*B.vec[0] + vec[1]*B.vec[1] + vec[2]*B.vec[2];
	return a;
}

//Level 2 BLAS

/// Vector-Matrix Multiply v = u'*A
inline Vec3 Vec3::operator* (const Mat33 &A)
{
	cout << "here\n";
	Vec3 v(	vec[0]*A.mat[0][0]+vec[1]*A.mat[0][1]+vec[2]*A.mat[0][2],
		vec[0]*A.mat[1][0]+vec[1]*A.mat[1][1]+vec[2]*A.mat[1][2],
		vec[0]*A.mat[2][0]+vec[1]*A.mat[2][1]+vec[2]*A.mat[2][2]);
	return v;
}

/// Vector-Matrix multiply u' = u'*A
inline Vec3& Vec3::operator*= (const Mat33 &A)
{
	Vec3 u(*this);
	vec[0] = u.vec[0]*A.mat[0][0]+u.vec[1]*A.mat[0][1]+u.vec[2]*A.mat[0][2];
	vec[1] = u.vec[0]*A.mat[1][0]+u.vec[1]*A.mat[1][1]+u.vec[2]*A.mat[1][2];
	vec[2] = u.vec[0]*A.mat[2][0]+u.vec[1]*A.mat[2][1]+u.vec[2]*A.mat[2][2];
	return *this;
}


//Special functions

/// Vector norm |v|
inline double Vec3::V3Norm() const
{
	return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}

/// Cross product v = u x w
inline Vec3 Vec3::V3Cross(const Vec3 &w)
{
	Vec3 v(	vec[1]*w.vec[2] - vec[2]*w.vec[1],
		vec[2]*w.vec[0] - vec[0]*w.vec[2],
		vec[0]*w.vec[1] - vec[1]*w.vec[0]);
	vec[0] = v.vec[0];
	vec[1] = v.vec[1];
	vec[2] = v.vec[2];
	return v;
}

/// Tensor product A = u*w'
inline Mat33 Vec3::V3TensProd(const Vec3 &w)
{
	Mat33 A;
	A.M3Zeros();
	
	for(int i = 0; i < 3; i++){
		A.mat[i][0] += vec[i]*w.vec[0];
		A.mat[i][1] += vec[i]*w.vec[1];
		A.mat[i][2] += vec[i]*w.vec[2];
	}
	
	return A;
}

/// Sum of vector elements
inline double Vec3::V3Sum() const
{
	double a = vec[0]+vec[1]+vec[2];
	return a;
}

/// Product of vector elements
inline double Vec3::V3Prod() const
{
	double a = vec[0]*vec[1]*vec[2];
	return a;
}

/*	***CVec3***	*/

// Constructors

/// Default constructor (everything is zero)
inline CVec3::CVec3(){
	cvec[0] = complex<double>(0,0);
	cvec[1] = complex<double>(0,0);
	cvec[2] = complex<double>(0,0);
}

/// construct from complex<double>
inline CVec3::CVec3(const complex<double> a, const complex<double> b, const complex<double> c){
	cvec[0] = a;
	cvec[1] = b;
	cvec[2] = c;
}

/// construct from Vec3
inline CVec3::CVec3(const Vec3 &u){
	cvec[0] = u.vec[0];
	cvec[1] = u.vec[1];
	cvec[2] = u.vec[2];
}

/// Copy constructor
inline CVec3::CVec3(const CVec3 &u){
	cvec[0] = u.cvec[0];
	cvec[1] = u.cvec[1];
	cvec[2] = u.cvec[2];
}

// Comparison operators
inline bool CVec3::operator==(const CVec3 &u)
{
	if (u.cvec[0] != cvec[0]) return false;
	if (u.cvec[1] != cvec[1]) return false;
	if (u.cvec[2] != cvec[2]) return false;
	return true;
}

inline bool CVec3::operator!=(const CVec3 &u)
{
	if (u.cvec[0] != cvec[0]) return true;
	if (u.cvec[1] != cvec[1]) return true;
	if (u.cvec[2] != cvec[2]) return true;
	return false;
}
inline bool CVec3::operator==(const Vec3 &u)
{
	if (cvec[0] != complex<double>(u.vec[0])) return false;
	if (cvec[1] != complex<double>(u.vec[1])) return false;
	if (cvec[2] != complex<double>(u.vec[2])) return false;
	return true;
}

inline bool CVec3::operator!=(const Vec3 &u)
{
	if (cvec[0] != complex<double>(u.vec[0])) return true;
	if (cvec[1] != complex<double>(u.vec[1])) return true;
	if (cvec[2] != complex<double>(u.vec[2])) return true;
	return true;
}

/// Assignment operator
inline CVec3& CVec3::operator=(const CVec3 &u){
	cvec[0] = u.cvec[0];
	cvec[1] = u.cvec[1];
	cvec[2] = u.cvec[2];
	return *this;
}

/// Scalar subtraction z = u - d
inline CVec3 CVec3::operator- (const double d)
{
	CVec3 cv(cvec[0]-d, cvec[1]-d, cvec[2]-d);
	return cv;
}

/// Scalar subtraction z = u - d
inline CVec3& CVec3::operator -= (const double d)
{
	cvec[0] -=d;
	cvec[1] -=d;
	cvec[2] -=d;
	return *this;
}

// Functions

inline Vec3 CVec3::Re(){
	Vec3 v(real(cvec[0]),real(cvec[1]),real(cvec[2]));
	return v;
}

inline Vec3 CVec3::Im(){
	Vec3 v(imag(cvec[0]),imag(cvec[1]),imag(cvec[2]));
	return v;
}

inline Vec3 CVec3::Abs(){
	Vec3 v(abs(cvec[0]),abs(cvec[1]),abs(cvec[2]));
	return v;
}

/*	***Mat33***	*/

// Constructors and destructors

/// Default constructor (identity matrix)
inline Mat33::Mat33()
{
	mat[0][0] = 1.0; mat[0][1] = 0.0; mat[0][2] = 0.0;
	mat[1][0] = 0.0; mat[1][1] = 1.0; mat[1][2] = 0.0;
	mat[2][0] = 0.0; mat[2][1] = 0.0; mat[2][2] = 1.0;
}

/// Copy constructor
inline Mat33::Mat33(const Mat33 &A)
{
	for (int i = 0; i < 3; i++)
	{
		mat[i][0] = A.mat[i][0];
		mat[i][1] = A.mat[i][1];
		mat[i][2] = A.mat[i][2];
	}
}

/// Tensor constructor
inline Mat33::Mat33(const Vec3 &u, const Vec3 &w)
{
	for(int i = 0; i < 3; i++){
		mat[i][0] = u.vec[i]*w.vec[0];
		mat[i][1] = u.vec[i]*w.vec[1];
		mat[i][2] = u.vec[i]*w.vec[2];
	}
}

/// Constant constructor
inline Mat33::Mat33(const double a)
{
	mat[0][0] = a;
	mat[0][1] = a;
	mat[0][2] = a;
	
	mat[1][0] = a;
	mat[1][1] = a;
	mat[1][2] = a;
	
	mat[2][0] = a;
	mat[2][1] = a;
	mat[2][2] = a;
}

/// Diagonal constructor
inline Mat33::Mat33(const double a, const double b, const double c)
{
	mat[0][0] = a;
	mat[0][1] = 0.0;
	mat[0][2] = 0.0;
	
	mat[1][0] = 0.0;
	mat[1][1] = b;
	mat[1][2] = 0.0;
	
	mat[2][0] = 0.0;
	mat[2][1] = 0.0;
	mat[2][2] = c;
}

/// Constant constructor from array
inline Mat33::Mat33(const double* a)
{
	mat[0][0] = a[0];
	mat[0][1] = a[1];
	mat[0][2] = a[2];
	
	mat[1][0] = a[3];
	mat[1][1] = a[4];
	mat[1][2] = a[5];
	
	mat[2][0] = a[6];
	mat[2][1] = a[7];
	mat[2][2] = a[8];
}

//Utility functions

/// Assignment operator operator B = A
inline Mat33& Mat33::operator= (const Mat33 &A)
{
	for (int i = 0; i < 3; i++){
		mat[i][0] = A.mat[i][0];
		mat[i][1] = A.mat[i][1];
		mat[i][2] = A.mat[i][2];
	}
	return *this;
}

/// Identity matrix
inline void Mat33::M3Eye()
{
	mat[0][0] = 1.0; mat[0][1] = 0.0; mat[0][2] = 0.0;
	mat[1][0] = 0.0; mat[1][1] = 1.0; mat[1][2] = 0.0;
	mat[2][0] = 0.0; mat[2][1] = 0.0; mat[2][2] = 1.0;
}

/// Zero matrix
inline void Mat33::M3Zeros()
{
	mat[0][0] = 0.0; mat[0][1] = 0.0; mat[0][2] = 0.0;
	mat[1][0] = 0.0; mat[1][1] = 0.0; mat[1][2] = 0.0;
	mat[2][0] = 0.0; mat[2][1] = 0.0; mat[2][2] = 0.0;
}


/// Comparison operators
inline bool Mat33::operator== (const Mat33 &A)
{
	for(int i = 0; i < 3; i++){
		if(fabs(A.mat[i][0]-mat[i][0])>EPS) return false;
		if(fabs(A.mat[i][1]-mat[i][1])>EPS) return false;
		if(fabs(A.mat[i][2]-mat[i][2])>EPS) return false;
	}
	return true;
}

inline bool Mat33::operator!= (const Mat33 &A)
{
	for(int i = 0; i < 3; i++){
		if(fabs(A.mat[i][0]-mat[i][0])>EPS) return true;
		if(fabs(A.mat[i][1]-mat[i][1])>EPS) return true;
		if(fabs(A.mat[i][2]-mat[i][2])>EPS) return true;
	}
	return false;
}


//Matrix-scalar operations

/// Matrix-scalar addition B = A + a
inline Mat33 Mat33::operator+ (const double a)
{
	Mat33 B;
	B.mat[0][0] = mat[0][0]+a;
	B.mat[0][1] = mat[0][1]+a;
	B.mat[0][2] = mat[0][2]+a;
	B.mat[1][0] = mat[1][0]+a;
	B.mat[1][1] = mat[1][1]+a;
	B.mat[1][2] = mat[1][2]+a;
	B.mat[2][0] = mat[2][0]+a;
	B.mat[2][1] = mat[2][1]+a;
	B.mat[2][2] = mat[2][2]+a;
	return B;
}

/// Matrix-scalar addition A = A + a
inline Mat33& Mat33::operator+= (const double a)
{
	mat[0][0] += a;
	mat[0][1] += a;
	mat[0][2] += a;
	mat[1][0] += a;
	mat[1][1] += a;
	mat[1][2] += a;
	mat[2][0] += a;
	mat[2][1] += a;
	mat[2][2] += a;
	return *this;
}

/// Matrix-scalar subtraction B = A - a
inline Mat33 Mat33::operator- (const double a)
{
	Mat33 B;
	
	B.mat[0][0] = mat[0][0]-a;
	B.mat[0][1] = mat[0][1]-a;
	B.mat[0][2] = mat[0][2]-a;
	B.mat[1][0] = mat[1][0]-a;
	B.mat[1][1] = mat[1][1]-a;
	B.mat[1][2] = mat[1][2]-a;
	B.mat[2][0] = mat[2][0]-a;
	B.mat[2][1] = mat[2][1]-a;
	B.mat[2][2] = mat[2][2]-a;
	
	return B;
}

/// Matrix-scalar subtraction A = A - a
inline Mat33& Mat33::operator-= (const double a)
{
	mat[0][0] -= a;
	mat[0][1] -= a;
	mat[0][2] -= a;
	mat[1][0] -= a;
	mat[1][1] -= a;
	mat[1][2] -= a;
	mat[2][0] -= a;
	mat[2][1] -= a;
	mat[2][2] -= a;
	return *this;
}

/// Matrix-vector subtraction from diagonal
inline Mat33 Mat33::operator- (const Vec3 &u)
{
	Mat33 B;
	B.mat[0][0] = mat[0][0]-u.vec[0];
	B.mat[1][1] = mat[1][1]-u.vec[1];
	B.mat[2][2] = mat[2][2]-u.vec[2];
	return B;
}

/// Matrix-vector subtraction from diagonal
inline Mat33& Mat33::operator-= (const Vec3 &u)
{
	mat[0][0] -= u.vec[0];
	mat[1][1] -= u.vec[1];
	mat[2][2] -= u.vec[2];
	return *this;
}

/// Matrix-scalar multiplication B = A*a
inline Mat33 Mat33::operator* (const double a)
{
	Mat33 B;
	
	B.mat[0][0] = mat[0][0]*a;
	B.mat[0][1] = mat[0][1]*a;
	B.mat[0][2] = mat[0][2]*a;
	
	B.mat[1][0] = mat[1][0]*a;
	B.mat[1][1] = mat[1][1]*a;
	B.mat[1][2] = mat[1][2]*a;
	
	B.mat[2][0] = mat[2][0]*a;
	B.mat[2][1] = mat[2][1]*a;
	B.mat[2][2] = mat[2][2]*a;
	
	return B;
}

/// Matrix-scalar multiplication A = A*a
inline Mat33& Mat33::operator*= (const double a)
{
	mat[0][0] *= a;
	mat[0][1] *= a;
	mat[0][2] *= a;
	mat[1][0] *= a;
	mat[1][1] *= a;
	mat[1][2] *= a;
	mat[2][0] *= a;
	mat[2][1] *= a;
	mat[2][2] *= a;
	return *this;
}

/// Matrix-scalar division B = A/a
inline Mat33 Mat33::operator/ (const double a)
{
	Mat33 B(*this);
	double inva = 1.0/a;
	B.mat[0][0] = mat[0][0]*inva;
	B.mat[0][1] = mat[0][1]*inva;
	B.mat[0][2] = mat[0][2]*inva;
	
	B.mat[1][0] = mat[1][0]*inva;
	B.mat[1][1] = mat[1][1]*inva;
	B.mat[1][2] = mat[1][2]*inva;
	
	B.mat[2][0] = mat[2][0]*inva;
	B.mat[2][1] = mat[2][1]*inva;
	B.mat[2][2] = mat[2][2]*inva;
	return B;
}

/// Matrix-scalar division A = A/a
inline Mat33& Mat33::operator/= (const double a)
{
	double inva = 1.0/a;
	mat[0][0] *= inva;
	mat[0][1] *= inva;
	mat[0][2] *= inva;
	mat[1][0] *= inva;
	mat[1][1] *= inva;
	mat[1][2] *= inva;
	mat[2][0] *= inva;
	mat[2][1] *= inva;
	mat[2][2] *= inva;
	return *this;
}


//Level 2 BLAS

/// Matrix-vector multiplication v = A*u
inline Vec3 Mat33::operator* (const Vec3 &u)
{
	return Vec3(mat[0][0]*u.vec[0]+mat[0][1]*u.vec[1]+mat[0][2]*u.vec[2],
		mat[1][0]*u.vec[0]+mat[1][1]*u.vec[1]+mat[1][2]*u.vec[2],
		mat[2][0]*u.vec[0]+mat[2][1]*u.vec[1]+mat[2][2]*u.vec[2]);
}


//Level 3 BLAS

/// Matrix addition B = A + C
inline Mat33 Mat33::operator+ (const Mat33 &C)
{
	Mat33 B(*this);
	for(int i = 0; i < 3; i++){
		B.mat[i][0] += C.mat[i][0];
		B.mat[i][1] += C.mat[i][1];
		B.mat[i][2] += C.mat[i][2];
	}
	return B;
}

/// Matrix addition C = A + C
inline Mat33& Mat33::operator+= (const Mat33 &C)
{
	for(int i = 0; i < 3; i++){
		mat[i][0] += C.mat[i][0];
		mat[i][1] += C.mat[i][1];
		mat[i][2] += C.mat[i][2];
	}
	return *this;
}

/// Matrix subtraction B = A - C
inline Mat33 Mat33::operator- (const Mat33 &C)
{
	Mat33 B(*this);
	for(int i = 0; i < 3; i++){
		B.mat[i][0] -= C.mat[i][0];
		B.mat[i][1] -= C.mat[i][1];
		B.mat[i][2] -= C.mat[i][2];
	}
	return B;
}

/// Matrix subtraction A = A - C
inline Mat33& Mat33::operator -= (const Mat33 &C)
{
	for(int i = 0; i < 3; i++){
		mat[i][0] -= C.mat[i][0];
		mat[i][1] -= C.mat[i][1];
		mat[i][2] -= C.mat[i][2];
	}
	return *this;
}

/// Matrix multiplication B = A*C
inline Mat33 Mat33::operator* (const Mat33 &C)
{
	Mat33 B;
	
	B.mat[0][0] = mat[0][0]*C.mat[0][0] + mat[0][1]*C.mat[1][0] + mat[0][2]*C.mat[2][0];
	B.mat[0][1] = mat[0][0]*C.mat[0][1] + mat[0][1]*C.mat[1][1] + mat[0][2]*C.mat[2][1];
	B.mat[0][2] = mat[0][0]*C.mat[0][2] + mat[0][1]*C.mat[1][2] + mat[0][2]*C.mat[2][2];
	
	B.mat[1][0] = mat[1][0]*C.mat[0][0] + mat[1][1]*C.mat[1][0] + mat[1][2]*C.mat[2][0];
	B.mat[1][1] = mat[1][0]*C.mat[0][1] + mat[1][1]*C.mat[1][1] + mat[1][2]*C.mat[2][1];
	B.mat[1][2] = mat[1][0]*C.mat[0][2] + mat[1][1]*C.mat[1][2] + mat[1][2]*C.mat[2][2];
	
	B.mat[2][0] = mat[2][0]*C.mat[0][0] + mat[2][1]*C.mat[1][0] + mat[2][2]*C.mat[2][0];
	B.mat[2][1] = mat[2][0]*C.mat[0][1] + mat[2][1]*C.mat[1][1] + mat[2][2]*C.mat[2][1];
	B.mat[2][2] = mat[2][0]*C.mat[0][2] + mat[2][1]*C.mat[1][2] + mat[2][2]*C.mat[2][2];
	
	return B;
}

/// Matrix multiplication A = A*C
inline Mat33& Mat33::operator*= (const Mat33 &C)
{
	Mat33 A(*this);
	
	mat[0][0] = A.mat[0][0]*C.mat[0][0] + A.mat[0][1]*C.mat[1][0] + A.mat[0][2]*C.mat[2][0];
	mat[0][1] = A.mat[0][0]*C.mat[0][1] + A.mat[0][1]*C.mat[1][1] + A.mat[0][2]*C.mat[2][1];
	mat[0][2] = A.mat[0][0]*C.mat[0][2] + A.mat[0][1]*C.mat[1][2] + A.mat[0][2]*C.mat[2][2];
	mat[1][0] = A.mat[1][0]*C.mat[0][0] + A.mat[1][1]*C.mat[1][0] + A.mat[1][2]*C.mat[2][0];
	mat[1][1] = A.mat[1][0]*C.mat[0][1] + A.mat[1][1]*C.mat[1][1] + A.mat[1][2]*C.mat[2][1];
	mat[1][2] = A.mat[1][0]*C.mat[0][2] + A.mat[1][1]*C.mat[1][2] + A.mat[1][2]*C.mat[2][2];
	mat[2][0] = A.mat[2][0]*C.mat[0][0] + A.mat[2][1]*C.mat[1][0] + A.mat[2][2]*C.mat[2][0];
	mat[2][1] = A.mat[2][0]*C.mat[0][1] + A.mat[2][1]*C.mat[1][1] + A.mat[2][2]*C.mat[2][1];
	mat[2][2] = A.mat[2][0]*C.mat[0][2] + A.mat[2][1]*C.mat[1][2] + A.mat[2][2]*C.mat[2][2];
	
	return *this;
}

//Matrix right division B = A/C
//In need of attention by a C++ expert LEJ 01/07/2010
/*inline Mat33 Mat33::operator/ (const Mat33 &C)
{
	Mat33 B(0.0);
	Mat33 D = C.M3Inv(NULL);
	for(int i = 0; i < 3; i++)
	{
		B.mat[i][0] += mat[i][0]*D.mat[0][0] + mat[i][1]*D.mat[1][0] + mat[i][2]*D.mat[2][0];
		B.mat[i][1] += mat[i][0]*D.mat[0][1] + mat[i][1]*D.mat[1][1] + mat[i][2]*D.mat[2][1];
		B.mat[i][2] += mat[i][0]*D.mat[0][2] + mat[i][1]*D.mat[1][2] + mat[i][2]*D.mat[2][2];
	}
	return B;
}

//Matrix right division A = A/C
inline Mat33& Mat33::operator/= (const Mat33 &C)
{
	Mat33 A(*this);
	Mat33 D = C.M3Inv(NULL);
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			mat[i][j] = A.mat[i][0]*C.mat[0][j] + A.mat[i][1]*C.mat[1][j] + A.mat[i][2]*C.mat[2][j];
		}
	}
	return *this;
}*/
								

//Special functions

inline Mat33 Mat33::TransMulM3(const Mat33 &C)
{
	Mat33 B;
	
	B.mat[0][0] = mat[0][0]*C.mat[0][0] + mat[1][0]*C.mat[1][0] + mat[2][0]*C.mat[2][0];
	B.mat[0][1] = mat[0][0]*C.mat[0][1] + mat[1][0]*C.mat[1][1] + mat[2][0]*C.mat[2][1];
	B.mat[0][2] = mat[0][0]*C.mat[0][2] + mat[1][0]*C.mat[1][2] + mat[2][0]*C.mat[2][2];
	
	B.mat[1][0] = mat[0][1]*C.mat[0][0] + mat[1][1]*C.mat[1][0] + mat[2][1]*C.mat[2][0];
	B.mat[1][1] = mat[0][1]*C.mat[0][1] + mat[1][1]*C.mat[1][1] + mat[2][1]*C.mat[2][1];
	B.mat[1][2] = mat[0][1]*C.mat[0][2] + mat[1][1]*C.mat[1][2] + mat[2][1]*C.mat[2][2];
	
	B.mat[2][0] = mat[0][2]*C.mat[0][0] + mat[1][2]*C.mat[1][0] + mat[2][2]*C.mat[2][0];
	B.mat[2][1] = mat[0][2]*C.mat[0][1] + mat[1][2]*C.mat[1][1] + mat[2][2]*C.mat[2][1];
	B.mat[2][2] = mat[0][2]*C.mat[0][2] + mat[1][2]*C.mat[1][2] + mat[2][2]*C.mat[2][2];
	
	return B;
}

/// Matrix-vector multiplication v = A*u
inline Vec3 Mat33::TransMulVec(const Vec3 &u)
{
	return Vec3(mat[0][0]*u.vec[0]+mat[1][0]*u.vec[1]+mat[2][0]*u.vec[2],
		mat[0][1]*u.vec[0]+mat[1][1]*u.vec[1]+mat[2][1]*u.vec[2],
		mat[0][2]*u.vec[0]+mat[1][2]*u.vec[1]+mat[2][2]*u.vec[2]);
}

inline Vec3 Mat33::ColumnVec3(const int col)
{
	Vec3 v(mat[0][col],mat[1][col],mat[2][col]);
	return v;
}

/// Matrix determinant d = Det(A)
inline double Mat33::M3Det()
{
	double det = mat[0][0]*(mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])
			+ mat[1][0]*(mat[2][1]*mat[0][2]-mat[0][1]*mat[2][2])
			+ mat[2][0]*(mat[0][1]*mat[1][2]-mat[0][2]*mat[1][1]);
	return det;
}

// Multiplication of non-symmetric matrix with symmetric matrix C
inline Mat33 Mat33::MulSymM3(const Mat33 &C)
{
	Mat33 B;
	
	B.mat[0][0] = mat[0][0]*C.mat[0][0] + mat[0][1]*C.mat[0][1] + mat[0][2]*C.mat[0][2];
	B.mat[0][1] = mat[0][0]*C.mat[0][1] + mat[0][1]*C.mat[1][1] + mat[0][2]*C.mat[1][2];
	B.mat[0][2] = mat[0][0]*C.mat[0][2] + mat[0][1]*C.mat[1][2] + mat[0][2]*C.mat[2][2];
	
	B.mat[1][0] = mat[1][0]*C.mat[0][0] + mat[1][1]*C.mat[0][1] + mat[1][2]*C.mat[0][2];
	B.mat[1][1] = mat[1][0]*C.mat[0][1] + mat[1][1]*C.mat[1][1] + mat[1][2]*C.mat[1][2];
	B.mat[1][2] = mat[1][0]*C.mat[0][2] + mat[1][1]*C.mat[1][2] + mat[1][2]*C.mat[2][2];
	
	B.mat[2][0] = mat[2][0]*C.mat[0][0] + mat[2][1]*C.mat[0][1] + mat[2][2]*C.mat[0][2];
	B.mat[2][1] = mat[2][0]*C.mat[0][1] + mat[2][1]*C.mat[1][1] + mat[2][2]*C.mat[1][2];
	B.mat[2][2] = mat[2][0]*C.mat[0][2] + mat[2][1]*C.mat[1][2] + mat[2][2]*C.mat[2][2];
	
	return B;
}

// Multiplication of symmetric matrix with non-symmetric matrix C
inline Mat33 Mat33::SymMulM3(const Mat33 &C)
{
	Mat33 B;
	
	B.mat[0][0] = mat[0][0]*C.mat[0][0] + mat[0][1]*C.mat[1][0] + mat[0][2]*C.mat[2][0];
	B.mat[0][1] = mat[0][0]*C.mat[0][1] + mat[0][1]*C.mat[1][1] + mat[0][2]*C.mat[2][1];
	B.mat[0][2] = mat[0][0]*C.mat[0][2] + mat[0][1]*C.mat[1][2] + mat[0][2]*C.mat[2][2];
	
	B.mat[1][0] = mat[0][1]*C.mat[0][0] + mat[1][1]*C.mat[1][0] + mat[1][2]*C.mat[2][0];
	B.mat[1][1] = mat[0][1]*C.mat[0][1] + mat[1][1]*C.mat[1][1] + mat[1][2]*C.mat[2][1];
	B.mat[1][2] = mat[0][1]*C.mat[0][2] + mat[1][1]*C.mat[1][2] + mat[1][2]*C.mat[2][2];
	
	B.mat[2][0] = mat[0][2]*C.mat[0][0] + mat[1][2]*C.mat[1][0] + mat[2][2]*C.mat[2][0];
	B.mat[2][1] = mat[0][2]*C.mat[0][1] + mat[1][2]*C.mat[1][1] + mat[2][2]*C.mat[2][1];
	B.mat[2][2] = mat[0][2]*C.mat[0][2] + mat[1][2]*C.mat[1][2] + mat[2][2]*C.mat[2][2];
	
	return B;
}

// Multiplication of two symmetric matrices
inline Mat33 Mat33::SymMulSymM3(const Mat33 &C)
{
	Mat33 B;
	
	B.mat[0][0] = mat[0][0]*C.mat[0][0] + mat[0][1]*C.mat[0][1] + mat[0][2]*C.mat[0][2];
	B.mat[0][1] = mat[0][0]*C.mat[0][1] + mat[0][1]*C.mat[1][1] + mat[0][2]*C.mat[1][2];
	B.mat[0][2] = mat[0][0]*C.mat[0][2] + mat[0][1]*C.mat[1][2] + mat[0][2]*C.mat[2][2];
	
	B.mat[1][0] = mat[0][1]*C.mat[0][0] + mat[1][1]*C.mat[0][1] + mat[1][2]*C.mat[0][2];
	B.mat[1][1] = mat[0][1]*C.mat[0][1] + mat[1][1]*C.mat[1][1] + mat[1][2]*C.mat[1][2];
	B.mat[1][2] = mat[0][1]*C.mat[0][2] + mat[1][1]*C.mat[1][2] + mat[1][2]*C.mat[2][2];
	
	B.mat[2][0] = mat[0][2]*C.mat[0][0] + mat[1][2]*C.mat[0][1] + mat[2][2]*C.mat[0][2];
	B.mat[2][1] = mat[0][2]*C.mat[0][1] + mat[1][2]*C.mat[1][1] + mat[2][2]*C.mat[1][2];
	B.mat[2][2] = mat[0][2]*C.mat[0][2] + mat[1][2]*C.mat[1][2] + mat[2][2]*C.mat[2][2];
	
	return B;
}

/*!
 * Cramer's rule inversion for 3x3 matrix
 */
inline Mat33 Mat33::M3Inv(bool &sing)
{
	Mat33 A(*this);
	
	A.mat[0][0] = mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1];
	A.mat[0][1] = mat[2][1]*mat[0][2]-mat[0][1]*mat[2][2];
	A.mat[0][2] = mat[0][1]*mat[1][2]-mat[0][2]*mat[1][1];
	
	double det = mat[0][0]*A.mat[0][0] + mat[1][0]*A.mat[0][1] + mat[2][0]*A.mat[0][2];
	
	if(fabs(det) < EPS){
		cout << "WARNING: Matrix is close to singular. det = " << det << "\n";
		sing=true;
	}
	
	double invdet = 1.0/det;
	
	A.mat[0][0] *= invdet;
	A.mat[0][1] *= invdet;
	A.mat[0][2] *= invdet;
	
	A.mat[1][0] = (mat[2][0]*mat[1][2]-mat[1][0]*mat[2][2])*invdet;
	A.mat[1][1] = (mat[0][0]*mat[2][2]-mat[2][0]*mat[0][2])*invdet;
	A.mat[1][2] = (mat[1][0]*mat[0][2]-mat[0][0]*mat[1][2])*invdet;
	
	A.mat[2][0] = (mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0])*invdet;
	A.mat[2][1] = (mat[2][0]*mat[0][1]-mat[0][0]*mat[2][1])*invdet;
	A.mat[2][2] = (mat[0][0]*mat[1][1]-mat[1][0]*mat[0][1])*invdet;
	
	return A;
}

/*!
 * Cramer's rule inversion for *symmetric* 3x3 matrix (expects upper half to be filled, WARNING: no checking if matrix is actually symmetric)
 */
inline Mat33 Mat33::SymM3Inv(bool &sing)
{
	Mat33 A;
	
	A.mat[0][0] = mat[1][1]*mat[2][2]-mat[1][2]*mat[1][2];
	A.mat[0][1] = mat[1][2]*mat[0][2]-mat[0][1]*mat[2][2];
	A.mat[0][2] = mat[0][1]*mat[1][2]-mat[0][2]*mat[1][1];
	double det = mat[0][0]*A.mat[0][0] + mat[0][1]*A.mat[0][1] + mat[0][2]*A.mat[0][2];
	
	if(fabs(det) < EPS){
		cout << "WARNING: Matrix is close to singular. det = " << det << "\n";
		sing=true;
	}
	
	double invdet = 1.0/det;
	
	A.mat[0][0] *= invdet;
	A.mat[0][1] *= invdet;
	A.mat[0][2] *= invdet;
	
	A.mat[1][0] = A.mat[0][1];
	A.mat[1][1] = (mat[0][0]*mat[2][2]-mat[0][2]*mat[0][2])*invdet;
	A.mat[1][2] = (mat[0][1]*mat[0][2]-mat[0][0]*mat[1][2])*invdet;
	
	A.mat[2][0] = A.mat[0][2];
	A.mat[2][1] = A.mat[1][2];
	A.mat[2][2] = (mat[0][0]*mat[1][1]-mat[0][1]*mat[0][1])*invdet;
	
	return A;
}

/*!
 * Cramer's rule inversion for *symmetric* 3x3 matrix (expects upper half to be filled, WARNING: no checking if matrix is actually symmetric)
 * followed by multiplication with vector
 */
inline Vec3 Mat33::SymM3InvMult(bool &sing, Vec3 &z)
{
	Vec3 v, r;
	
	v.vec[0]=mat[1][1]*mat[2][2]-mat[1][2]*mat[1][2];
	v.vec[1]=mat[1][2]*mat[0][2]-mat[0][1]*mat[2][2];
	v.vec[2]=mat[0][1]*mat[1][2]-mat[0][2]*mat[1][1];
	
	double det = mat[0][0]*v.vec[0] + mat[0][1]*v.vec[1] + mat[0][2]*v.vec[2];
	if(fabs(det) < EPS){
		cout << "WARNING: Matrix is close to singular. det = " << det << "\n";
		sing=true;
	}
	
	r.vec[0] = v.vec[0]*z.vec[0]+v.vec[1]*z.vec[1]+v.vec[2]*z.vec[2];
	v.vec[0] = (mat[0][1]*mat[0][2]-mat[0][0]*mat[1][2]); // a_12
	r.vec[1] = v.vec[1]*z.vec[0]+(mat[0][0]*mat[2][2]-mat[0][2]*mat[0][2])*z.vec[1]+v.vec[0]*z.vec[2];
	r.vec[2] = v.vec[2]*z.vec[0]+v.vec[0]*z.vec[1]+(mat[0][0]*mat[1][1]-mat[0][1]*mat[0][1])*z.vec[2];
	r/=det;
	
	return r;
}

/// Multiply matrix with diagonal matrix specified through a,b,c C = A*[a,0,0 ; 0,b,0 ; 0,0,c ]
inline Mat33 Mat33::M3MulDiag(const double a, const double b, const double c)
{
	Mat33 B;
	
	B.mat[0][0] = mat[0][0]*a;
	B.mat[0][1] = mat[0][1]*b;
	B.mat[0][2] = mat[0][2]*c;
	
	B.mat[1][0] = mat[1][0]*a;
	B.mat[1][1] = mat[1][1]*b;
	B.mat[1][2] = mat[1][2]*c;
	
	B.mat[2][0] = mat[2][0]*a;
	B.mat[2][1] = mat[2][1]*b;
	B.mat[2][2] = mat[2][2]*c;
	
	return B;
}

/// Matrix transpose C = A'
inline Mat33 Mat33::M3Transpose()
{
	Mat33 C;
	C.mat[0][0] = mat[0][0];
	C.mat[0][1] = mat[1][0];
	C.mat[0][2] = mat[2][0];
	C.mat[1][0] = mat[0][1];
	C.mat[1][1] = mat[1][1];
	C.mat[1][2] = mat[2][1];
	C.mat[2][0] = mat[0][2];
	C.mat[2][1] = mat[1][2];
	C.mat[2][2] = mat[2][2];
	return C;
}

/// Matrix row swap
inline Mat33 Mat33::M3RowSwap(const int a, const int b)
{
	Mat33 A;
	double temp;
	temp = mat[b][0];
	A.mat[b][0] = mat[a][0];
	A.mat[a][0] = temp;
	temp = mat[b][1];
	A.mat[b][1] = mat[a][1];
	A.mat[a][1] = temp;
	temp = mat[b][2];
	A.mat[b][2] = mat[a][2];
	A.mat[a][2] = temp;
	return A;
}

/// Extract diagonal elements of matrix v = Aii -- goes in Mat33 class
inline Vec3 Mat33::M3Diag()
{
	Vec3 v(mat[0][0],mat[1][1],mat[2][2]);
	return v;
}

/// Trace of a matrix a = Tr(A)
inline double Mat33::M3Trace()
{
	double a = mat[0][0]+mat[1][1]+mat[2][2];
	return a;
}

/// Back-substitution on 3x3 triangular system
inline bool Mat33::M3BackSub(Vec3 &x)
{
	bool multiples=false;
	if (fabs(mat[2][2])>EPS){
		x.vec[2] = x.vec[2]/mat[2][2];
	} else{
		x.vec[2]=1.0;
		multiples=true;
	}
	if (fabs(mat[1][1])>EPS){
		x.vec[1] = (x.vec[1] - mat[1][2]*x.vec[2])/mat[1][1];
	} else{
		x.vec[1]=1.0;
		multiples=true;
	}
	if (fabs(mat[0][0])>EPS){
		x.vec[0] = (x.vec[0] - mat[0][1]*x.vec[1] - mat[0][2]*x.vec[2])/mat[0][0];
	} else{
		x.vec[0]=1.0;
		multiples=true;
	}
#if DEBUG_LEVEL>2
	cout << M3Str() << "\n---\n";
#endif
	return multiples;
}

/// Linear Solver for 3x3 systems using GEPP + back-substitution
inline Vec3 Mat33::M3LinSolve(Vec3 &u, bool &multiples)
{
	Vec3 v(u);
	Mat33 A(*this);
	A.M3GEPP(v);
	multiples=A.M3BackSub(v);
	return v;
}

inline Mat33 Vec2Rot(Vec3 &rot_vec, const unsigned int a0) // a0 is reference axis {0,1,2} = {x,y,z}
{
	Mat33 a,b;
	rot_vec/=rot_vec.V3Norm(); // normalize (just in case)
	unsigned int a1=(a0+1)%3;
	unsigned int a2=(a0+2)%3;
	double sqr_argument=rot_vec.vec[a1]*rot_vec.vec[a1]+rot_vec.vec[a2]*rot_vec.vec[a2];
	if(sqr_argument>EPS){
		double invhyp=1.0/sqrt(sqr_argument);
		double cp=rot_vec.vec[a2]*invhyp;
		double sp=rot_vec.vec[a1]*invhyp;
		double st=sqrt(1.0-rot_vec.vec[a0]*rot_vec.vec[a0]);
		
		// we're doing axis rotation here, that's why the sign of the sin term sp is opposite
		a.mat[a0][a0]=1.0;	a.mat[a0][a1]=0.0;	a.mat[a0][a2]=0.0;
		a.mat[a1][a0]=0.0;	a.mat[a1][a1]=cp;	a.mat[a1][a2]=sp;
		a.mat[a2][a0]=0.0;	a.mat[a2][a1]=-1.0*sp;	a.mat[a2][a2]=cp;
		
		b.mat[a0][a0]=rot_vec.vec[a0];	b.mat[a0][a1]=0.0;	b.mat[a0][a2]=-1.0*st;
		b.mat[a1][a0]=0.0;		b.mat[a1][a1]=1.0;	b.mat[a1][a2]=0.0;
		b.mat[a2][a0]=st;		b.mat[a2][a1]=0.0;	b.mat[a2][a2]=rot_vec.vec[a0];
		
		return a*b; // return rotation matrix if vector not pointing in z-direction
	}
	// if we get here then the vector is in axis-direction
	if(rot_vec.vec[a0]<0.0){ // rotate 180 degrees around next axis
		a.mat[a0][a0]=-1.0;
		a.mat[a2][a2]=-1.0;
	}
	return a; // return unit matrix if no rotation necessary (also applies to zero-vector)
}

/*	***Vec4***	*/

/// Default constructor (zeros)
inline Vec4::Vec4()
{
	vec[0] = 0.0; vec[1] = 0.0; vec[2] = 0.0; vec[3] = 0.0;
}

/// Construct from doubles
inline Vec4::Vec4(const double a, const double b, const double c, const double d)
{
	vec[0] = a;
	vec[1] = b;
	vec[2] = c;
	vec[3] = d;
}

/// Construct from Vec3 and double
inline Vec4::Vec4(const Vec3 &u, const double d)
{
	vec[0] = u.vec[0];
	vec[1] = u.vec[1];
	vec[2] = u.vec[2];
	vec[3] = d;
}

/// Copy constructor
inline Vec4::Vec4(const Vec4 &u)
{
	vec[0] = u.vec[0];
	vec[1] = u.vec[1];
	vec[2] = u.vec[2];
	vec[3] = u.vec[3];
}

/// Vector addition v = u + w
inline Vec4 Vec4::operator+ (const Vec4 &w)
{
	Vec4 v(vec[0]+w.vec[0], vec[1]+w.vec[1], vec[2] +w.vec[2], vec[3] +w.vec[3]);
	return v;
}

/// Vector scalar multiplication v = u*a
inline Vec4 Vec4::operator* (const double a)
{
	Vec4 v(a*vec[0], a*vec[1], a*vec[2], a*vec[3]);
	return v;
}

/// Vector scalar division v = u/a
inline Vec4 Vec4::operator/ (const double a)
{
	double b = 1.0/a;
	Vec4 v(b*vec[0], b*vec[1], b*vec[2], b*vec[3]);
	return v;
}

/// Comparison operators
inline bool Vec4::operator==(const Vec4 &u)
{
	if(fabs(u.vec[0]-vec[0])>EPS) return false;
	if(fabs(u.vec[1]-vec[1])>EPS) return false;
	if(fabs(u.vec[2]-vec[2])>EPS) return false;
	if(fabs(u.vec[3]-vec[3])>EPS) return false;
	return true;
}

inline bool Vec4::operator!=(const Vec4 &u)
{
	if(fabs(u.vec[0]-vec[0])>EPS) return true;
	if(fabs(u.vec[1]-vec[1])>EPS) return true;
	if(fabs(u.vec[2]-vec[2])>EPS) return true;
	if(fabs(u.vec[3]-vec[3])>EPS) return true;
	return false;
}

/// Rotation axis and angle from rotation matrix
inline Vec4 Rot2AxisAngle(Mat33 R)
{
	Vec4 AaA(1.0,0.0,0.0,0.0); // using axis=(1,0,0), theta=0 as failsafe b/c X3D animation would fail with (0,0,0,0) ...
	double cost=(R.M3Trace()-1.0)/2.0;
	if(cost>1.0) cost=1.0; // safety first
	if(cost<-1.0) cost=-1.0;
	double theta = acos(cost);
	double sintheta = sin(theta);
	double n_norm=1.0;
	if(fabs(sintheta)>EPS*EPS){
		n_norm = 1.0/(2.0*sintheta);
	}
	AaA.vec[0] = n_norm*(R.mat[2][1]-R.mat[1][2]);
	AaA.vec[1] = n_norm*(R.mat[0][2]-R.mat[2][0]);
	AaA.vec[2] = n_norm*(R.mat[1][0]-R.mat[0][1]);
	if(AaA.vec[0]*AaA.vec[0]+AaA.vec[1]*AaA.vec[1]+AaA.vec[2]*AaA.vec[2]<EPS) return Vec4(1.0,0.0,0.0,0.0);
	AaA.vec[3] = theta;
	return AaA;
}


/// Rotation matrix from axis and angle
inline Mat33 AxisAngle2Rot(Vec4 AxisAngle)
{
	// normalize axis part, just in case
	double norm=AxisAngle.vec[0]*AxisAngle.vec[0]+AxisAngle.vec[1]*AxisAngle.vec[1]+AxisAngle.vec[2]*AxisAngle.vec[2];
	Mat33 rot;
	norm=sqrt(norm);
	if(norm>=EPS){ // safety first
		AxisAngle.vec[0]/=norm; AxisAngle.vec[1]/=norm; AxisAngle.vec[2]/=norm;
		
		double oneminuscos=1.0-cos(AxisAngle.vec[3]);
		double sine=sin(AxisAngle.vec[3]);
		double x=AxisAngle.vec[0]; double y=AxisAngle.vec[1]; double z=AxisAngle.vec[2];
		
		rot.mat[0][0]+=oneminuscos*(x*x-1.0);	rot.mat[0][1]=-z*sine+oneminuscos*x*y;	rot.mat[0][2]=y*sine+oneminuscos*x*z;
		rot.mat[1][0]=z*sine+oneminuscos*x*y;	rot.mat[1][1]+=oneminuscos*(y*y-1.0);	rot.mat[1][2]=-x*sine+oneminuscos*y*z;
		rot.mat[2][0]=-y*sine+oneminuscos*x*z;	rot.mat[2][1]=x*sine+oneminuscos*y*z;	rot.mat[2][2]+=oneminuscos*(z*z-1.0);
	}
	return rot;
}

/// Rotation matrix from axis and angle
inline Mat33 AxisAngle2Rot(Vec3 &axis, double &angle)
{
	Vec4 aangle(axis.vec[0],axis.vec[1],axis.vec[2],angle);
	return AxisAngle2Rot(aangle);
}

inline bool Point_in_Ellipsoid(Vec3 &point, Vec3 &saxes, Vec3 &center, Mat33 &rot)
{
	// first rotate point-to-center distance vector so that ellipsoid is aligned with coordinate system
	Vec3 dist=rot.M3Transpose()*(point-center);
	dist.vec[0]/=saxes.vec[0]; // x/a
	dist.vec[1]/=saxes.vec[1]; // y/b
	dist.vec[2]/=saxes.vec[2]; // z/c
	
	if(dist*dist<=1.0) return true; else return false; // point is inside if (x^2/a^2 + y^2/b^2 + z^2/c^2 <= 1)
}

#define atan_a -0.012299380859105
#define atan_b 0.054082655552459
#define atan_c -0.11769677376706
#define atan_d 0.19402227554937
#define atan_e -0.33269718723178
#define atan_f 0.99998657415361

inline double fastatan(double x)
{
	double arg=fabs(x);
	double arg2=x*x;
	if(arg<=1.0){
		return copysign((((((atan_a*arg2+atan_b)*arg2+atan_c)*arg2+atan_d)*arg2+atan_e)*arg2+atan_f)*arg,x);
	} else{
		arg=1.0/arg;
		arg2=arg*arg;
		return copysign(pi/2-(((((atan_a*arg2+atan_b)*arg2+atan_c)*arg2+atan_d)*arg2+atan_e)*arg2+atan_f)*arg,x);
	}
}

inline double fastatan2(double y, double x)
{
	if(x>0.0) return fastatan(y/x);
	if(y>=0.0) return fastatan(y/x)+pi;
	return fastatan(y/x)-pi;
}

inline void UnitVec2ThetaPhi(Vec3 &v, double &theta, double &phi)
{
	if(v.vec[2]>1.0) v.vec[2]=1.0; // safety first
	if(v.vec[2]<-1.0) v.vec[2]=-1.0;
	theta=acos(v.vec[2]);
	phi=atan(v.vec[1]/v.vec[0]);
}

inline void Vec2ThetaPhi(Vec3 &v, double &theta, double &phi)
{
	Vec3 r=v/v.V3Norm();
	if(r.vec[2]>1.0) r.vec[2]=1.0; // safety first
	if(r.vec[2]<-1.0) r.vec[2]=-1.0;
	theta=acos(r.vec[2]);
	phi=atan(r.vec[1]/r.vec[0]); // atan is 1/0 safe ;-)
}

inline void Vec2ThetaPhi(Vec3 &v, double &theta, double &phi, double &r)
{
	r=v.V3Norm();
	Vec3 rvec=v/r;
	if(rvec.vec[2]>1.0) rvec.vec[2]=1.0; // safety first
	if(rvec.vec[2]<-1.0) rvec.vec[2]=-1.0;
	theta=acos(rvec.vec[2]);
	phi=atan(rvec.vec[1]/rvec.vec[0]);
}

inline double VecDist2ThetaPhi(Vec3 &v, double r, double &theta, double &phi)
{
	Vec3 rvec=v;
	rvec.vec[2]/=r;
	if(rvec.vec[2]>1.0) rvec.vec[2]=1.0; // safety first
	if(rvec.vec[2]<-1.0) rvec.vec[2]=-1.0;
	theta=acos(rvec.vec[2]);
	phi=atan(rvec.vec[1]/rvec.vec[0]);
	return rvec.vec[2];
}

inline double VecDist2Phi(Vec3 &v, double r, double &phi)
{
	double cost=v.vec[2]/r;
	if(cost>1.0) cost=1.0; // safety first
	if(cost<-1.0) cost=-1.0;
	phi=fastatan2(v.vec[1],v.vec[0]);
	return cost;
}

inline double EllipsoidRmin(double &theta, double &phi, Vec3 &saxes)
{
	// first rotate point-to-center distance vector so that ellipsoid is aligned with coordinate system
	double sint=sin(theta);
	Vec3 dist(sint*cos(phi),sint*sin(phi),cos(theta));
	dist.vec[0]/=saxes.vec[0]; // x/a
	dist.vec[1]/=saxes.vec[1]; // y/b
	dist.vec[2]/=saxes.vec[2]; // z/c
	// Now solve x^2/a^2+y^2/b^2+z^2/c^2=1/r^2 => r = sqrt(1/(x^2/a^2+y^2/b^2+z^2/c^2))
	return 1.0/sqrt(dist*dist);
}

inline double EllipsoidRmin(double &theta, double &phi, Vec3 &saxes, Mat33 &rot)
{
	// first rotate point-to-center distance vector so that ellipsoid is aligned with coordinate system
	double sint=sin(theta);
	Vec3 direction(sint*cos(phi),sint*sin(phi),cos(theta));
	Vec3 dist=rot.M3Transpose()*direction;
	dist.vec[0]/=saxes.vec[0]; // x/a
	dist.vec[1]/=saxes.vec[1]; // y/b
	dist.vec[2]/=saxes.vec[2]; // z/c
	// Now solve x^2/a^2+y^2/b^2+z^2/c^2=1/r^2 => r = sqrt(1/(x^2/a^2+y^2/b^2+z^2/c^2))
	return 1.0/sqrt(dist*dist);
}

inline double EllipsoidRmin(Vec3 &direction, Vec3 &saxes, Mat33 rot)
{
	// first rotate point-to-center distance vector so that ellipsoid is aligned with coordinate system
	Vec3 dist=rot.M3Transpose()*direction/direction.V3Norm();
	dist.vec[0]/=saxes.vec[0]; // x/a
	dist.vec[1]/=saxes.vec[1]; // y/b
	dist.vec[2]/=saxes.vec[2]; // z/c
	// Now solve x^2/a^2+y^2/b^2+z^2/c^2=1/r^2 => r = sqrt(1/(x^2/a^2+y^2/b^2+z^2/c^2))
	return 1.0/sqrt(dist*dist);
}

// Based on P.P. Klein, "On the Ellipsoid and Plane Intersection Equation", Applied Mathematics 3, 1634 (2012)
// see page 1639
inline double Ellipsoid_Cross_Section(Vec3 &invsaxes2, Vec3 direction)
{
	double a=direction*direction;
	if(a>EPS*EPS){
		a=1.0/a;
		double c=direction.vec[0]*direction.vec[0]*(invsaxes2.vec[1]*invsaxes2.vec[2])+direction.vec[1]*direction.vec[1]*(invsaxes2.vec[0]*invsaxes2.vec[2])+direction.vec[2]*direction.vec[2]*(invsaxes2.vec[0]*invsaxes2.vec[1]);
		// a*beta^2 + b*beta + c = 0
		// beta_+/- = -b/2a +/- sqrt(b^2 - 4*a*c)/2a
		// beta+ * beta_ = (-b + sqrt(b^2 - 4*a*c)) * (-b - sqrt(b^2 - 4*a*c))/(2a)^2 = (b^2 - (b^2 - 4ac))/4a^2 = 4ac/4a^2 = c/a
		return pi/sqrt(c*a);
	} else return 0.0;
}

inline double VectorAngle(Vec3 &a, Vec3 &b)
{
	double theta;
	double cosine=(a*b)/sqrt((a*a)*(b*b));
	if(cosine>1.0){ // the Cauchy-Schwartz inequality can be broken by computers ...
		theta=0.0;
	} else{
		if(cosine<-1.0){
			theta=pi;
		} else theta=acos(cosine);
	}
	return theta;
}

/// Obtain rotation matrix mapping vector a to b direction
inline Mat33 RotAtoB(Vec3 &a, Vec3 &b)
{
	Vec3 axis=a;
	axis.V3Cross(b);
	double angle=VectorAngle(a,b);
	if(axis*axis<EPS){ // vectors are colinear, two solutions here: everything stays as is (angle=0) or rotate pi around vector perpendicular to a (and b)
		if(fabs(fabs(angle)-pi)<EPS){
			// find non zero component of a
			if(fabs(a.vec[0])>EPS){
				axis.vec[1]=1.0;
				axis.vec[2]=1.0;
				// solve for axis*a = 0
				// => a_x*axis_x+a_y*axis_y+a_z*axis_z
				// => axis_x = -(a_y+a_z)/a_x
				axis.vec[0]=-(a.vec[1]+a.vec[2])/a.vec[0];
			} else{
				if(fabs(a.vec[1])>EPS){
					axis.vec[0]=1.0;
					axis.vec[2]=1.0;
					axis.vec[1]=-(a.vec[0]+a.vec[2])/a.vec[1];
				} else{
					if(fabs(a.vec[2])>EPS){
						axis.vec[0]=1.0;
						axis.vec[1]=1.0;
						axis.vec[2]=-(a.vec[0]+a.vec[1])/a.vec[2];
					}
				}
			}
		}
	}
	return AxisAngle2Rot(axis,angle);
}

inline double VectorCos(Vec3 &a, Vec3 &b)
{
	double cosine=(a*b)/sqrt((a*a)*(b*b));
	if(cosine>1.0){ // the Cauchy-Schwartz inequality can be broken by computers ...
		cosine=1.0;
	} else{
		if(cosine<-1.0) cosine=-1.0;
	}
	return cosine;
}

inline Vec3 Rot2AlphaBetaGamma(Mat33 &rot)
{
	Vec3 result;
	/*    Z(a)            Y(b)           Z(c)
	 * ( ca  sa  0) * ( cb  0  sb) * ( cc  sc  0)   ( ca  sa  0) * ( cb*cc  cb*sc sb)   ( ca*cb*cc+sa*sc  ca*cb*sc+sa*cc  ca*sb)
	 * (-sa  ca  0) * ( 0   1   0) * (-sc  cc  0) = (-sa  ca  0) * (  -sc     cc  0 ) = (-sa*cb*cc-ca*sc -sa*cb*sc+ca*cc -sa*sb)
	 * ( 0   0   1) * (-sb  0  cb) * ( 0   0   1)   ( 0   0   1) * (-sb*cc -sb*sc cb)   (    -sb*cc          -sb*sc         cb )
	 */
	
	double cb=rot.mat[2][2];
	if(cb>1.0) cb=1.0;
	if(cb<-1.0) cb=-1.0;
	double y=-rot.mat[1][2];
	double x=rot.mat[0][2];
	
	result.vec[0] = atan2(y,x);
	result.vec[1] = acos(cb);
	
	y=-rot.mat[2][1];
	x=-rot.mat[2][0];
	
	result.vec[2] = atan2(y,x);
	return result;
}

#ifndef OpenCL

struct double4
{
	double x,y,z,w;
	double4(){}
	double4(double v){ x = y = z = w = v; }
	
	double4 operator*(const double4& other)
	{
		double4 tmp;
		tmp.x = x*other.x;
		tmp.y = y*other.y;
		tmp.z = z*other.z;
		tmp.w = w*other.w;
		return tmp;
	}
	
	double4 operator*(const double& other)
	{
		double4 tmp;
		tmp.x = x*other;
		tmp.y = y*other;
		tmp.z = z*other;
		tmp.w = w*other;
		return tmp;
	}
	
	double4& operator+=(const double4& other)
	{
		x += other.x;
		y += other.y;
		z += other.z;
		w += other.w;
		return *this;
	}
	
	double4& operator-=(const double4& other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;
		w -= other.w;
		return *this;
	}
	
	double4& operator*=(double scalar)
	{
		x *= scalar;
		y *= scalar;
		z *= scalar;
		w *= scalar;
		return (*this);
	}
	
	double4& operator/=(double scalar)
	{
		x /= scalar;
		y /= scalar;
		z /= scalar;
		w /= scalar;
		return (*this);
	}
};

inline double4 fabs(const double4& a)
{
	double4 tmp;
	tmp.x = a.x < 0.f ? 0.f  : a.x;
	tmp.y = a.y < 0.f ? 0.f  : a.y;
	tmp.z = a.z < 0.f ? 0.f  : a.z;
	tmp.w = a.w < 0.f ? 0.f  : a.w;
	return tmp;
}

inline double4 operator+(const double4& a,const double4& b)
{
	double4 tmp;
	tmp.x = a.x + b.x;
	tmp.y = a.y + b.y;
	tmp.z = a.z + b.z;
	tmp.w = a.w + b.w;
	return tmp;
}

inline double4 operator-(const double4& a,const double4& b)
{
	double4 tmp;
	tmp.x = a.x - b.x;
	tmp.y = a.y - b.y;
	tmp.z = a.z - b.z;
	tmp.w = a.w - b.w;
	return tmp;
}

inline double4 operator*(const double4& a,const double& s)
{
	double4 tmp;
	tmp.x = a.x*s;
	tmp.y = a.y*s;
	tmp.z = a.z*s;
	tmp.w = a.w*s;
	return tmp;
}

inline double4 operator/(const double4& a,const double& s)
{
	double4 tmp;
	tmp.x = a.x/s;
	tmp.y = a.y/s;
	tmp.z = a.z/s;
	tmp.w = a.w/s;
	return tmp;
}

inline double4 cross(const double4& p0, const double4& p1)
{
	double4 result;
	result.w=0.0;
	
	result.x=p0.y*p1.z-p0.z*p1.y;
	result.y=p0.z*p1.x-p0.x*p1.z;
	result.z=p0.x*p1.y-p0.y*p1.x;
	
	return result;
}

inline double dot(const double4& p0, const double4& p1)
{
	return p0.x*p1.x+p0.y*p1.y+p0.z*p1.z+p0.w*p1.w;
}

inline double4 normalize(const double4& a)
{
	double norm=dot(a,a);
	norm=1.0/sqrt(norm);
	return a*norm;
}

typedef double double16[16];

#else // ifndef OpenCL
typedef cl_double4 double4;
typedef cl_double16 double16;
#endif

inline double4 create_double4(double a, double b, double c, double d)
{
	double4 result;
	result.x = a;
	result.y = b;
	result.z = c;
	result.w = d;
	return result;
}

inline double4 create_double4(Vec3 v)
{
	double4 result;
	result.x = v.vec[0];
	result.y = v.vec[1];
	result.z = v.vec[2];
	result.w = 0.0;
	return result;
}

inline double4 create_double4(double d)
{
	double4 result;
	result.x = d;
	result.y = d;
	result.z = d;
	result.w = d;
	return result;
}

inline double4 AxisAngle2Quaternion(Vec4& aa)
{
	double4 result;
	result.w=cos(aa.vec[3]/2.0);
	double norm=sin(aa.vec[3]/2.0)/sqrt(aa.vec[0]*aa.vec[0]+aa.vec[1]*aa.vec[1]+aa.vec[2]*aa.vec[2]);
	result.x=aa.vec[0]*norm;
	result.y=aa.vec[1]*norm;
	result.z=aa.vec[2]*norm;
	
	return result;
}

inline Mat33 RotFromEigenvectors(Mat33 ev)
{
	Mat33 rot;
	// check if eigenvectors are orthogonal
#if DEBUG_LEVEL>2
	cout << "Eigenvectors:\n#1 = (" << ev.ColumnVec3(0).V3Str(',') << ")\n#2 = (" << ev.ColumnVec3(1).V3Str(',') << ")\n#3 = (" << ev.ColumnVec3(2).V3Str(',') << ")\n";
#endif
	bool no12=(fabs(ev.ColumnVec3(0)*ev.ColumnVec3(1))>EPS);
	bool no13=(fabs(ev.ColumnVec3(0)*ev.ColumnVec3(2))>EPS);
	bool no23=(fabs(ev.ColumnVec3(1)*ev.ColumnVec3(2))>EPS);
	bool all_equal=false;
	if(no12 || no13 || no23){
#if DEBUG_LEVEL>2
		cout << "first eigenvector times second: " << ev.ColumnVec3(0)*ev.ColumnVec3(1) << " (" << (ev.ColumnVec3(0)-ev.ColumnVec3(1)).V3Norm() << ")\n";
		cout << "first eigenvector times third: " << ev.ColumnVec3(0)*ev.ColumnVec3(2) << " (" << (ev.ColumnVec3(0)-ev.ColumnVec3(2)).V3Norm() << ")\n";
		cout << "second eigenvector times third: " << ev.ColumnVec3(1)*ev.ColumnVec3(2) << " (" << (ev.ColumnVec3(1)-ev.ColumnVec3(2)).V3Norm() << ")\n";
#endif
		bool diff12=((ev.ColumnVec3(0)-ev.ColumnVec3(1)).V3Norm()<EPS);
		bool diff13=((ev.ColumnVec3(0)-ev.ColumnVec3(2)).V3Norm()<EPS);
		bool diff23=((ev.ColumnVec3(1)-ev.ColumnVec3(2)).V3Norm()<EPS);
		if(!(diff12 || diff13 || diff23)){ // numerical issue (calculating the respective crossproduct again should get solved)
			if(no12) diff12=true;
			if(no13) diff13=true;
			if(no23) diff23=true;
		}
		if((unsigned int)(diff12+diff13+diff23)<=1){
			Vec3 a,b;
			unsigned int col=0;
#if DEBUG_LEVEL>1
			cout << "WARNING: Eigenvectors ";
#endif
			if(diff12){ // #1=#2
				col=1;
				a=ev.ColumnVec3(0);
				b=ev.ColumnVec3(2);
#if DEBUG_LEVEL>1
				cout << "#1 and #2";
#endif
			} else{
				if(diff13){ // #1=#3
					col=2;
					a=ev.ColumnVec3(0);
					b=ev.ColumnVec3(1);
#if DEBUG_LEVEL>1
					cout << "#1 and #3";
#endif
				} else{ // diff23 (#2=#3)
					col=2;
					a=ev.ColumnVec3(0);
					b=ev.ColumnVec3(1);
#if DEBUG_LEVEL>1
					cout << "#2 and #3";
#endif
				}
			}
			Vec3 new_ev=a.V3Cross(b);
			ev.mat[0][col]=new_ev.vec[0];
			ev.mat[1][col]=new_ev.vec[1];
			ev.mat[2][col]=new_ev.vec[2];
#if DEBUG_LEVEL>1
			cout << " were not orthogonal";
#endif
#if DEBUG_LEVEL>2
			cout << ", new eigenvector #" << col+1 << " is: " << new_ev.V3Str(',') << "\n";
#else
			cout << ".\n";
#endif
		} else all_equal=true;
	}
	// now determine rotation matrix
	rot.M3Eye(); // start with identity matrix
	if(!all_equal){
		bool xz=false; bool yz=false; bool zz=false;
		unsigned int zerocount=0;
		Vec3 x=ev.ColumnVec3(0);
		double r=x.V3Norm();
		if(r>EPS) x/=r; else{ xz=true; zerocount++; }
		Vec3 y=ev.ColumnVec3(1);
		r=y.V3Norm();
		if(r>EPS) y/=r; else{ yz=true; zerocount++; }
		Vec3 z=ev.ColumnVec3(2);
		r=z.V3Norm();
		if(r>EPS) z/=r; else{ zz=true; zerocount++; }
		if(zerocount<=1){
			if(xz){
				x=y;
				x.V3Cross(z);
			} else{
				if(yz){
					y=z;
					y.V3Cross(x);
				} else{
					if(zz){
						z=x;
						x.V3Cross(y);
					}
				}
			}
			rot.mat[0][0]=x.vec[0]; rot.mat[0][1]=y.vec[0]; rot.mat[0][2]=z.vec[0];
			rot.mat[1][0]=x.vec[1]; rot.mat[1][1]=y.vec[1]; rot.mat[1][2]=z.vec[1];
			rot.mat[2][0]=x.vec[2]; rot.mat[2][1]=y.vec[2]; rot.mat[2][2]=z.vec[2];
		} else{
			if(zerocount==2){
				if(!xz){
					rot=Vec2Rot(x,0);
				} else{
					if(!yz){
						rot=Vec2Rot(y,1);
					} else{
						if(!zz){
							rot=Vec2Rot(z,2);
						}
					}
				}
			}
		}
	}
	// Need to check determinat being +1, if -1 coordinate system ended up being right-handed ...
	double determinant=rot.M3Det();
	if(fabs(fabs(determinant)-1.0)>EPS){ // looks weird, but is correct b/c det(rot)=-1 is also acceptable here (recoverable)
		cout << "Could not determine proper rotation matrix, det(rot)=" << determinant << ", which is not +1\n";
		exit(2);
	}
	if(determinant<0.0){ // right-handed coordinate system needs to be changed to left-handed
		rot.mat[0][0]*=-1.0; // do this by flipping x around
		rot.mat[1][0]*=-1.0;
		rot.mat[2][0]*=-1.0;
#if DEBUG_LEVEL>2
		cout << "Fixed improper rotation matrix\n";
#endif
	}
#if DEBUG_LEVEL>2
	cout << "rotation matrix:\n" << rot.M3Str() << "\n";
	cout << "det(rot) = " << rot.M3Det() << "\n";
	bool sing=false;
	cout << "inverse rotation matrix:\n" << (rot.M3Inv(sing)).M3Str() << "\n";
	cout << "rotated x: " << (rot*Vec3(1.0,0.0,0.0)).V3Str(',') << "\n";
	cout << "rotated y: " << (rot*Vec3(0.0,1.0,0.0)).V3Str(',') << "\n";
	cout << "rotated z: " << (rot*Vec3(0.0,0.0,1.0)).V3Str(',') << "\n";
#endif
	return rot;
}

inline void BackSubstitute(Mat33 &A, Vec3& solution, unsigned l)
{
	bool diagzero[3], columnzero[3];
	for(unsigned int j=0; j<3; j++){
		diagzero[j]=(fabs(A.mat[j][j])<=EPS);
		columnzero[j]=false;
		for(unsigned int k=0; k<3; k++) columnzero[j]&=(fabs(A.mat[k][j])<=EPS);
	}
	bool goty=false;
	if(diagzero[2]){
		if(!columnzero[2] && (l==2)){
			if(diagzero[1]){
				if(fabs(A.mat[1][2])>EPS){
					solution.vec[2] = solution.vec[2]/A.mat[1][2]; // z-solution is uniquely defined
				} else{ // element 02 is non-zero
					if(fabs(A.mat[0][1])>EPS){ // y-solution is defined here and can compensate any z-solution (may as well be 1 then)
						solution.vec[2]=1.0;
					} else{ // y-solution is arbitrary
						if(diagzero[0]){ // x-solution is arbitary, but z-solution is defined
							solution.vec[2] = solution.vec[2]/A.mat[0][2];
						} else solution.vec[2] = 1.0; // since x-solution exists it can compensate any z-solution (may as well be 1)
					}
				}
			} else{ // unique y-solution exists
				if(fabs(A.mat[1][2])>EPS){ // we care for z-solution existing but not really for y-solution, so let's try to eliminate the y-solution
					// element 01 *must* exist (otherwise Gaussian elimination is not sorted)
					// no matter if x-solution exists or not, y and z-solutions need to be the
					// same for the first and the second row (and we can set y-solution to zero)
					solution.vec[2] = solution.vec[1]/A.mat[1][2];
					solution.vec[1] = 0.0;
					goty=true;
				} else{ // second row has nothing to do with z-solution which defines y-solution (and means z-solution is in first row)
					solution.vec[1] = solution.vec[1]/A.mat[1][1];
					goty=true;
					if(fabs(A.mat[0][1])>EPS){
						if(diagzero[0]){ // x-solution is arbitrary (has nothing to do with first row), hence y-solution solves unique z-solution
							solution.vec[2] = (solution.vec[0] - solution.vec[1]*A.mat[0][1])/A.mat[0][2];
						} else solution.vec[2] = 1.0; // x-solution can compensate any y and z-solution (which may as well be 1 then)
					} else solution.vec[2] = 1.0; // z-solution is independent of y-solution, and no matter what the x-solution is the z-solution is either arbitrary or compensate (aka may as well be 1)
				}
			}
		} else solution.vec[2]=0.0; // it does not matter what number we choose for z-solution, may as well get rid of it ...
	} else solution.vec[2] = solution.vec[2]/A.mat[2][2]; // z-solution is uniquely defined
	if(!goty){
		if(diagzero[1]){
			if(!columnzero[1] && (l==1)){ // y-solution is defined in first row and we care (b/c it's corresponding to the largest eigenvalue)
				if(diagzero[0]){ // x-solution is arbitrary and has nothing to do with anything, hence y-solution is defined
					solution.vec[1] = (solution.vec[0] - solution.vec[2]*A.mat[0][2])/A.mat[0][1];
				} else{ // x-solution exists and can compensate for y-solution (which may as well be 1 then)
					solution.vec[1] = 1.0;
				}
			} else solution.vec[1] = 0.0; // we can't be bothered
		} else solution.vec[1] = (solution.vec[1] - A.mat[1][2]*solution.vec[2])/A.mat[1][1]; // y-solution is uniquely defined
	}
	if(diagzero[0]){ // x-solution can be anything
		if(l==0) solution.vec[0]=1.0; else solution.vec[0]=0.0;
	} else solution.vec[0] = (solution.vec[0] - A.mat[0][1]*solution.vec[1] - A.mat[0][2]*solution.vec[2])/A.mat[0][0];
}

inline bool M3equal(Mat33 A, Mat33 B, double error)
{
	for(int i = 0; i < 3; i++){
		if(fabs(A.mat[i][0]-B.mat[i][0])>error) return false;
		if(fabs(A.mat[i][1]-B.mat[i][1])>error) return false;
		if(fabs(A.mat[i][2]-B.mat[i][2])>error) return false;
	}
	return true;
}

inline double touch_sphere_sigma(Vec3 &r, Vec3 &saxes, Mat33 &rot, double rT)
{
	if(rT<EPS) return EllipsoidRmin(r,saxes,rot);
	// Minimization vectors
	Vec3 V, CI, saxes2;
	double t;
	
	// Create lab frame version of A and B matrices (both are symmetric matrices R * L_A/B * R^T)
	// A/B_ij=sum_k L_A/B_k*R_ik*R_jk
	// move into A ellipsoid frame of reference
	// doing so:
	// - replaces 36 multiplication and 12 additions with 36 multiplications and 24 additions once
	// - saves 4 multiplications and 6 additions in loop
	Vec3 z; // multply with kk's transposed (inverse) rotation matrix
	z.vec[0] = r.vec[0]*rot.mat[0][0]+r.vec[1]*rot.mat[1][0]+r.vec[2]*rot.mat[2][0];
	z.vec[1] = r.vec[0]*rot.mat[0][1]+r.vec[1]*rot.mat[1][1]+r.vec[2]*rot.mat[2][1];
	z.vec[2] = r.vec[0]*rot.mat[0][2]+r.vec[1]*rot.mat[1][2]+r.vec[2]*rot.mat[2][2];
	double r2=r*r;
	double d=r2/(saxes.vec[0]*saxes.vec[1]*saxes.vec[2]);
	double e=d*d;
	double rT2 = rT*rT*e;
	saxes2.vec[0] = saxes.vec[0]*saxes.vec[0]*e;
	saxes2.vec[1] = saxes.vec[1]*saxes.vec[1]*e;
	saxes2.vec[2] = saxes.vec[2]*saxes.vec[2]*e;
	z*=d;
	
	double lambda=0.5;
	double xlx=1.0; // x=lambda/(1-lambda) -> do manual calculation with value above
	double Var = 1.0; // trying different forms found 6 loops gives 7 figs for Fx
	double VAV, VBV, det, Vz, V22;
	while(Var > 1E-6){ // loop until variance in distance is sufficiently small
		Var = lambda; // keep lambda around but do CI matrix in terms of x (scale CI by 1/(1-lambda))
		//Populate CI matrix
		CI.vec[0] = xlx*rT2+saxes2.vec[0];
		CI.vec[1] = xlx*rT2+saxes2.vec[1];
		CI.vec[2] = xlx*rT2+saxes2.vec[2];
		
		// Solve z = CI*V for V using inverse => V=CI^-1*z
		V.vec[0] = CI.vec[1]*CI.vec[2]*z.vec[0];
		V.vec[1] = CI.vec[0]*CI.vec[2]*z.vec[1];
		V.vec[2] = CI.vec[0]*CI.vec[1]*z.vec[2];
		
		// VAV=V*A_LF*V (uses fact that A_LF is symmetric)
		V22=V.vec[2]*V.vec[2];
		VAV = V.vec[0]*V.vec[0]*saxes2.vec[0]+V.vec[1]*V.vec[1]*saxes2.vec[1]+V22*saxes2.vec[2];
		// denominator=V*B_LF*V (uses fact that B_LF is symmetric)
		VBV = (V.vec[0]*V.vec[0]+V.vec[1]*V.vec[1]+V22)*rT2;
		//Calculate minimization parameter lambda
/*		if(VBV < EPS*EPS){
			cout << "ERROR: Denominator between oids in touch is too close to zero (" << VBV << ").\n";
			exit(3);
		}*/
		xlx = sqrt(VAV/VBV); // independent of z-scaling (and determinant) -> also, interesting note: the sqrt is better than anything else in terms of speed and convergence
		lambda = xlx/(1.0+xlx);
		Var -= lambda;
		Var *= Var;
	}
	
	//Reconstruct CI and run a final iteration once converged
	CI.vec[0] = xlx*rT2+saxes2.vec[0];
	CI.vec[1] = xlx*rT2+saxes2.vec[1];
	CI.vec[2] = xlx*rT2+saxes2.vec[2];
	
	t=CI.vec[1]*CI.vec[2];
	det = CI.vec[0]*t;
	V22=2.0*z.vec[2];
	Vz=z.vec[0]*z.vec[0]*t+CI.vec[0]*CI.vec[2]*z.vec[1]*z.vec[1]+CI.vec[0]*CI.vec[1]*z.vec[2]*z.vec[2];
	// return sigma
	return sqrt(r2*det/(lambda*Vz));
}

#endif

