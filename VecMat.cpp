/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 *  VecMat.cpp
 *  Created by ljohnson on 7/24/09. Last modified on 01/06/10
 *  Three-dimensional vector and matrix classes. Most functions are contained in the header; only non-inlined functions
 *	(mostly helper functions such as string conversion operators) are located here.
 */

#include "VecMat.h"

/*!
 * String conversion functions
 * Stringify a vector, with tabs between elements
 */
std::string Vec3::V3Str() const
{
	std::stringstream converter;
	std::cout.precision(6);
	
	converter << std::fixed <<vec[0]<<"\t"<<vec[1]<<"\t"<<vec[2];
	
	std::string vecstring = converter.str();
	return vecstring;
}

/*!
 * String conversion functions
 * Stringify a vector, with user-selected delimiters between elements
 */
std::string Vec3::V3Str(const char &dlm) const
{
	//Check to make sure delimiter is valid (tab, space, comma, or line break). Default to tab
	char delim;
	if (!((dlm == '\t') || (dlm == ' ') || (dlm == ',') || (dlm == '\n'))) {
		delim = '\t';
	}
	else delim = dlm;
	
	std::stringstream converter;
	std::cout.precision(6);
	
	converter << std::fixed <<vec[0]<< delim <<vec[1]<< delim <<vec[2];
	
	std::string vecstring = converter.str();
	return vecstring;
}

/// Stringify a complex vector, with tabs between elements
std::string CVec3::CV3Str() const
{
	std::stringstream converter;
	std::cout.precision(6);
	
	converter << std::fixed <<cvec[0]<<"\t"<<cvec[1]<<"\t"<<cvec[2];
	
	return converter.str();
}

/// Stringify a complex vector, with user-specified delimiters between elements
std::string CVec3::CV3Str(const char &dlm) const
{
	//Check to make sure delimiter is valid (tab, space, comma, or line break). Default to tab
	char delim;
	
	if (!((dlm == '\t') || (dlm == ' ') || (dlm == ',') || (dlm == '\n'))){
		delim = '\t';
	} else delim = dlm;
	
	std::stringstream converter;
	std::cout.precision(6);
	
	converter << std::fixed << cvec[0] << delim << cvec[1] << delim <<cvec[2];
	
	return converter.str();
}

/// Stringify a matrix, with tabs between elements
std::string Mat33::M3Str()
{
	std::stringstream converter;
	std::cout.precision(6);
	
	converter << std::fixed <<mat[0][0]<<"\t"<<mat[0][1]<<"\t"<<mat[0][2]<<"\n";
	converter << std::fixed <<mat[1][0]<<"\t"<<mat[1][1]<<"\t"<<mat[1][2]<<"\n";
	converter << std::fixed <<mat[2][0]<<"\t"<<mat[2][1]<<"\t"<<mat[2][2];
	
	std::string matstring = converter.str();
	return matstring;
}

/// Stringify a row of a matrix, with tabs between elements
std::string Mat33::M3RowStr(int i)
{
	std::stringstream converter;
	std::cout.precision(6);
	
	converter << std::fixed <<mat[i][0]<<"\t"<<mat[i][1]<<"\t"<<mat[i][2];
	
	std::string matstring = converter.str();
	return matstring;
}

/// Stringify a row of a matrix, with tabs between elements
std::string Mat33::M3RowStr(int i, const char &dlm)
{
	//Check to make sure delimiter is valid (tab, space, comma, or line break). Default to tab
	char delim;
	
	if (!((dlm == '\t') || (dlm == ' ') || (dlm == ',') || (dlm == '\n'))) {
		delim = '\t';
	} else delim = dlm;
	
	std::stringstream converter;
	std::cout.precision(6);
	
	converter << std::fixed << mat[i][0] << delim << mat[i][1] << delim <<mat[i][2];
	
	std::string matstring = converter.str();
	return matstring;
}

///Stringify a Vec4, with tabs between elements
std::string Vec4::V4Str() const
{
	std::stringstream converter;
	std::cout.precision(6);
	
	converter << std::fixed <<vec[0]<<"\t"<<vec[1]<<"\t"<<vec[2]<<"\t"<<vec[3];
	
	std::string vecstring = converter.str();
	return vecstring;
}

///Stringify a Vec4, with user-selected delimiters between elements
std::string Vec4::V4Str(const char &dlm) const
{
	//Check to make sure delimiter is valid (tab, space, comma, or line break). Default to tab
	char delim;
	if (!((dlm == '\t') || (dlm == ' ') || (dlm == ',') || (dlm == '\n'))) {
		delim = '\t';
	}
	else delim = dlm;
	
	std::stringstream converter;
	std::cout.precision(6);
	
	converter << std::fixed <<vec[0]<< delim <<vec[1]<< delim <<vec[2] << delim << vec[3];
	
	std::string vecstring = converter.str();
	return vecstring;
}

/*!
 * Gaussian elimination with partial pivoting for solving Ax = b
 * Used within M3LinSolve, but not inlined due to size
 */
void Mat33::M3GEPP(Vec3 &x)
{
	// Perform elimination
	for(unsigned int i=0; i<3; i++){
		// Find largest pivot (to be on diagonal - AT), with index pidx
		unsigned int pidx = i;
		for(unsigned int j=i+1; j<3; j++){
			if(fabs(mat[j][i]) > fabs(mat[pidx][i])) pidx = j; // was compared to mat[pidx][j], but no need to care for the diagonal down the road ...
		}
		
		// Swap rows -- replace when figure out how to use class function within each other LEJ 01/07/10
		if(pidx != i){
			for(unsigned int j=0; j<3; j++){
				double tempij = mat[pidx][j];
				mat[pidx][j] = mat[i][j];
				mat[i][j] = tempij;
			}
			double tempij =  x.vec[pidx];
			x.vec[pidx] = x.vec[i];
			x.vec[i] = tempij;
		}
		// Eliminate
		for(unsigned k=i+1; k<3; k++){
			x.vec[k] -= mat[k][i]/mat[i][i]*x.vec[i];
			for(int j=2; j>(int)i; j--) mat[k][(unsigned int)j] -= (mat[k][i]/mat[i][i])*mat[i][(unsigned int)j];
			mat[k][i]=0.0;
		}
	}
}

/// LU decomposition routine. WARNING: No pivoting at the moment -- AT
void Mat33::LUDecomposition()
{
	int i,j,k;
	for(i=0; i<3; i++){
		for(j=i; j<3; j++){
			for(k=0; k<i-1; k++){
				mat[i][j] -= mat[i][k]*mat[k][j];
			}
		}
		for(j=i+1; j<3; j++){
			for(k=0; k<i-1; k++){
				mat[j][i] -= mat[j][k]*mat[k][i];
			}
#if DEBUG_LEVEL>2
			cout << "a[" << j << "][" << i << "] = " << mat[j][i] << ", a[" << i << "][" << i << "] = " << mat[i][i]<< "\n";
#endif
			mat[j][i] /= mat[i][i];
		}
	}
}

/// get solutions to third order polynomial in reduced form: z^3 + p*z + q = 0 -- AT
CVec3 SolvePolynomial3(const double p, const double q)
{
	double D=q*q/4.0+p*p*p/27.0;
#if DEBUG_LEVEL>2
	cout << "Start SolvePolynomial3\n";
	cout << "p = " << p << ", q = " << q << " => D = " << D << "\n";
#endif
	double u, v;
	double minusqhalf=-0.5*q;
	CVec3 cv;
	
	if(D>EPS*EPS){ // one real, two complex solutions
		double sqrtD=sqrt(D);
		u=cbrt(sqrtD+minusqhalf);
		v=cbrt(-sqrtD+minusqhalf);
#if DEBUG_LEVEL>2
	cout << "u = " << u << ", v = " << v << "\n";
#endif
		cv.cvec[0]=u+v;
		cv.cvec[1]=complex<double>(-0.5*(u+v),0.5*sqrt3*(u-v));
		cv.cvec[2]=conj(cv.cvec[1]);
	} else{
		if(D<-EPS*EPS){ // three different real solutions, now things get complex to calculate ;-)
			// for any complex number z = Re+Im*i = A*e^i*theta ; A = sqrt(Re^2+Im^2) ; theta = atan(Im/Re)
			// sqrt(D) = i sqrt(-D)
			double sqrtD=sqrt(-D); // is imaginary now ...
			// u = cbrt(minusqhalf+i*sqrt(-D)) = cbrt(sqrt(minusqhalf^2-D)*e^(i*atan(sqrt(-D)/minusqhalf)))
			// => u = cbrt(sqrt(minusqhalf^2-D)*e^(i*1/3*atan(sqrt(-D)/minusqhalf))
			// similar for v: v = cbrt(sqrt(minusqhalf^2-D))*e^(-i*1/3*atan(sqrt(-D)/minusqhalf));
			double uv_length=cbrt(sqrt(minusqhalf*minusqhalf-D)); // cubicroot(sqrt(Re^2 + Im^2) (Im^2 = sqrt(-D)^2)
			double uv_phase=atan2(sqrtD,minusqhalf)/3.0;
#if DEBUG_LEVEL>2
			cout << "magnitude = " << uv_length << ", phase = " << uv_phase << "\n";
#endif
			u=uv_length*cos(uv_phase); // here, u is real part
			v=uv_length*sin(uv_phase); // here, v is imaginary part
#if DEBUG_LEVEL>2
			cout << "u = v* = (" << u << ", " << v << ")\n";
#endif
			cv.cvec[0]=u+u;
			cv.cvec[1]=sqrt3*v-u;
			cv.cvec[2]=-sqrt3*v-u;
		} else{ // D=0: either threefold real solution, or real and two-fold real solution (two distinct solutions)
#if DEBUG_LEVEL>2
			cout << "u = v = " << cbrt(minusqhalf) << "\n";
#endif
			cv.cvec[0]=0.0;
			if(fabs(p)>EPS) cv.cvec[0]=3*q/p;
			cv.cvec[1]=cv.cvec[0];
			cv.cvec[1]*=-0.5;
			cv.cvec[2]=cv.cvec[1];
		}
	}
#if DEBUG_LEVEL>2
	cout << "Finished SolvePolynomial3.\n";
#endif
	return cv;
}

/*!
 * Determine eigenvalues of 3x3 matrix -- AT
 * characteristic equation: det(A-lambda*E)=0 (lambda ... eigenvalues; E=3x3 unit diagonal matrix)
 *
 * |	a11-lambda	a12	a13	|
 * |	a21	a22-lambda	a23	| = 0 = lambda^3 + alpha*lambda^2 + beta*lambda + gamma
 * |	a31	a32	a33-lambda	|
 *
 * alpha = -a11 - a22 - a33 ;
 * beta = a11*(a22+a33) - alphax - a13*a31 - a12*a21
 * gamma = a11*alphax - a12*(a23*a31-a21*a33) - a13*(a21*a32-a31*a22)
 * alphax = a23*a32 - a22*a33
 *
 * reduced form after substition with lambda = z - alpha/3
 * z^3 + p*z + q = 0
 *
 * p = beta - alpha*alpha/3
 * q = (2*alpha^3 - 9*alpha*beta + 27*gamma)/27 = gamma + alpha/3*(2*(alpha/3)^2 - beta)
 */
CVec3 Mat33::Eigenvalues()
{
#if DEBUG_LEVEL>2
	cout << "Start Mat33::Eigenvalue\n";
#endif
	double alpha = -mat[0][0]-mat[1][1]-mat[2][2];
	double alphax = mat[1][2]*mat[2][1]-mat[1][1]*mat[2][2];
	double beta = mat[0][0]*(mat[1][1]+mat[2][2])-alphax-mat[0][2]*mat[2][0]-mat[0][1]*mat[1][0];
	double gamma = mat[0][0]*alphax-mat[0][1]*(mat[1][2]*mat[2][0]-mat[1][0]*mat[2][2])-mat[0][2]*(mat[1][0]*mat[2][1]-mat[2][0]*mat[1][1]);
#if DEBUG_LEVEL>2
	cout << "alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << ", alphax = " << alphax << "\n";
#endif
	CVec3 cv;
	alphax=alpha/3.0; // alphax redefined
	cv=SolvePolynomial3(beta-alpha*alphax,gamma+alphax*(2.0*alphax*alphax-beta));
	cv-=alphax; // lambda = z - alpha/3
#if DEBUG_LEVEL>2
	cout << "Finished Mat33::Eigenvalues\n";
#endif
	return cv;
}

Mat33 Mat33::Eigenvectors(Vec3 &ew, bool* multiples, bool normalize)
{
	Mat33 result;
	for(unsigned int i=0; i<3; i++){
		Mat33 A(*this);
		A.mat[0][0]-=ew.vec[i]; A.mat[1][1]-=ew.vec[i]; A.mat[2][2]-=ew.vec[i];
		Vec3 zero(0.0);
		Vec3 v=A.M3LinSolve(zero,multiples[i]);
		if(normalize) v/=v.V3Norm();
		result.mat[0][i]=v.vec[0]; result.mat[1][i]=v.vec[1]; result.mat[2][i]=v.vec[2];
	}
	return result;
}

