/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

CL_TO_STRING(

inline void M3mult(double* A, double* B, double* C);
inline void M3multT(double* A, double* B, double* C);
inline double4 M3Vec_mult(double* A, double4 u);
inline double4 M3TVec_mult(double* A, double4 u);
inline double M3det(double* A);
inline void M3Inv(double* A, unsigned int* sing);
inline double qvec_norm(double4 q);
inline double4 vector_angles(double4 v); // stores theta, phi, cos(theta)
inline double4 quaternion_rotate(double4 v, double4 rot);
inline double4 quaternion_transpose(double4 rot);
short Point_in_Ellipsoid(double4 point, double4 saxes, double4 center, double* rot);
double ran_cl(int* idum, int* iv, int* iy, int* idum2);
double qpwr_cl(double a, unsigned int b);
double EllipsoidRmin(double theta, double phi, double4 saxes, double* rot);
double EllipsoidRminVec(double4 direction, double4 saxes, double* rot);
inline double specialround_CL(double A);
inline void apply_PBCs(double4* rvec, __constant Simulation_Attribs* configuration);

/* random function
 * needs to be called with negative idum value to initialize iv[32], iy, and idum2 
 */
double ran_cl(int* idum, int* iv, int* iy, int* idum2)
{
	int j;
	int k;
	
	if(*idum<=0){
		if(-(*idum)<1) *idum=1; else *idum=-(*idum);
		(*idum2)=*idum;
		for(j=39;j>=0;j--){
			k=(*idum)/53668;
			*idum=40014*(*idum)-k*(53668*40014+12211);
			if(*idum < 0) *idum += 2147483563;
			if(j < 32) iv[j] = *idum;
		};
		(*iy)=iv[0];
	}
	k=(*idum)/53668;
	*idum=40014*(*idum)-k*(53668*40014+12211);
	if(*idum < 0) *idum += 2147483563;
	k=(*idum2)/52774;
	(*idum2)=40692*(*idum2)-k*(52774*40692+3791);
	if((*idum2) < 0) (*idum2) += 2147483399;
	j=(*iy)/67108862;
	(*iy)=iv[j]-(*idum2);
	iv[j] = (*idum);
	if ((*iy) < 1) (*iy) += 2147483562;
	double temp=(double)(*iy)/(double)2147483563;
	if(temp > 1.0-EPS_CL) return 1.0-EPS_CL; else return temp;
}

/* Fast power algorithm that only works for positive integers */
inline double qpwr_cl(double a, unsigned int b)
{
	double c = 1;
	for(unsigned int i = 0; i < b; i++) c *= a;
	return c;
}

/// Matrix multiplication A*B=C
inline void M3mult(double* A, double* B, double* C)
{
	double result[9];
	result[0] = A[0]*B[0] + A[1]*B[3] + A[2]*B[6];
	result[1] = A[0]*B[1] + A[1]*B[4] + A[2]*B[7];
	result[2] = A[0]*B[2] + A[1]*B[5] + A[2]*B[8];
	result[3] = A[3]*B[0] + A[4]*B[3] + A[5]*B[6];
	result[4] = A[3]*B[1] + A[4]*B[4] + A[5]*B[7];
	result[5] = A[3]*B[2] + A[4]*B[5] + A[5]*B[8];
	result[6] = A[6]*B[0] + A[7]*B[3] + A[8]*B[6];
	result[7] = A[6]*B[1] + A[7]*B[4] + A[8]*B[7];
	result[8] = A[6]*B[2] + A[7]*B[5] + A[8]*B[8];
	C[0]=result[0];
	C[1]=result[1];
	C[2]=result[2];
	C[3]=result[3];
	C[4]=result[4];
	C[5]=result[5];
	C[6]=result[6];
	C[7]=result[7];
	C[8]=result[8];
}

/// Matrix multiplication A*B^T=C
inline void M3multT(double* A, double* B, double* C)
{
	double result[9];
	result[0] = A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
	result[1] = A[0]*B[3] + A[1]*B[4] + A[2]*B[5];
	result[2] = A[0]*B[6] + A[1]*B[7] + A[2]*B[8];
	result[3] = A[3]*B[0] + A[4]*B[1] + A[5]*B[2];
	result[4] = A[3]*B[3] + A[4]*B[4] + A[5]*B[5];
	result[5] = A[3]*B[6] + A[4]*B[7] + A[5]*B[8];
	result[6] = A[6]*B[0] + A[7]*B[1] + A[8]*B[2];
	result[7] = A[6]*B[3] + A[7]*B[4] + A[8]*B[5];
	result[8] = A[6]*B[6] + A[7]*B[7] + A[8]*B[8];
	C[0]=result[0];
	C[1]=result[1];
	C[2]=result[2];
	C[3]=result[3];
	C[4]=result[4];
	C[5]=result[5];
	C[6]=result[6];
	C[7]=result[7];
	C[8]=result[8];
}

/// Matrix-vector multiplication A*u=v
inline double4 M3Vec_mult(double* A, double4 u)
{
	double4 result;
	result.w=u.w;
	
	result.x=A[0]*u.x+A[1]*u.y+A[2]*u.z;
	result.y=A[3]*u.x+A[4]*u.y+A[5]*u.z;
	result.z=A[6]*u.x+A[7]*u.y+A[8]*u.z;
	
	return result;
}

/// transpose Matrix-vector multiplication A^T*u=v
inline double4 M3TVec_mult(double* A, double4 u)
{
	double4 result;
	result.w=u.w;
	
	result.x=A[0]*u.x+A[3]*u.y+A[6]*u.z;
	result.y=A[1]*u.x+A[4]*u.y+A[7]*u.z;
	result.z=A[2]*u.x+A[5]*u.y+A[8]*u.z;
	
	return result;
}

/// Matrix determinant d = Det(A)
inline double M3det(double* A)
{
	double det = A[0]*(A[4]*A[8]-A[5]*A[7])
		   + A[3]*(A[7]*A[2]-A[1]*A[8])
		   + A[6]*(A[1]*A[5]-A[2]*A[4]);
	return det;
}

/*!
 * Cramer's rule inversion for 3x3 matrix
 */
inline void M3Inv(double* A, unsigned int* sing)
{
	*sing=0;
	double det = M3det(A);
	
	if(fabs(det)<EPS_CL) *sing=1;
	
	double invdet = 1.0/det;
	
	double inv[9];
	inv[0] = (A[4]*A[8]-A[7]*A[5])*invdet;
	inv[1] = (A[7]*A[2]-A[1]*A[8])*invdet;
	inv[2] = (A[1]*A[5]-A[4]*A[2])*invdet;
	
	inv[3] = (A[6]*A[5]-A[3]*A[8])*invdet;
	inv[4] = (A[0]*A[8]-A[6]*A[2])*invdet;
	inv[5] = (A[3]*A[2]-A[0]*A[5])*invdet;
	
	inv[6] = (A[3]*A[7]-A[4]*A[6])*invdet;
	inv[7] = (A[6]*A[1]-A[0]*A[7])*invdet;
	inv[8] = (A[0]*A[4]-A[3]*A[1])*invdet;
	
	A[0]=inv[0];
	A[1]=inv[1];
	A[2]=inv[2];
	A[3]=inv[3];
	A[4]=inv[4];
	A[5]=inv[5];
	A[6]=inv[6];
	A[7]=inv[7];
	A[8]=inv[8];
}

inline double qvec_norm(double4 q)
{
	return sqrt(q.x*q.x+q.y*q.y+q.z*q.z);
}

inline double4 quaternion_rotate(double4 v, double4 rot)
{
	double4 result;
	
	result.w=0.0;
	double4 z=cross(rot,v)*2.0;
	result = v + z*rot.w + cross(rot,z);
	
	return result;
}

inline double4 quaternion_transpose(double4 rot)
{
	double4 result;
	result.w=rot.w;
	result.x=-rot.x;
	result.y=-rot.y;
	result.z=-rot.z;
	return result;
}

inline double4 vector_angles(double4 v)
{
	double4 result;
	result.w=0.0;
	
	// gives x=sin(theta)*cos(phi) ; y=sin(theta)*sin(phi) ; z=cos(theta)
	double4 r=normalize(v);
	if(r.z>1.0) r.z=1.0; // safety first
	if(r.z<-1.0) r.z=-1.0;
	result.z=r.z; // = cos(theta)
	result.x=acos(r.z); // = theta
	r.y/=r.x;
	result.y=atan(r.y);
	
	return result;
}

short Point_in_Ellipsoid(double4 point, double4 saxes, double4 center, double* rot)
{
//	/* first rotate point-to-center distance vector so that ellipsoid is aligned with coordinate system */
//	double4 dist=quaternion_rotate(point-center,quaternion_transpose(rot));
	double4 dist=M3TVec_mult(rot,point-center);
	dist.x=dist.x/saxes.x; /* x/a */
	dist.y=dist.y/saxes.y; /* y/b */
	dist.z=dist.z/saxes.z; /* z/c */
	
	// point is inside if (x^2/a^2 + y^2/b^2 + z^2/c^2 <= 1)
	if(dot(dist,dist)<=1.0) return 1;
	return 0;
}

double EllipsoidRmin(double theta, double phi, double4 saxes, double* rot)
{
	/* first rotate point-to-center distance vector so that ellipsoid is aligned with coordinate system */
	double sint=sin(theta);
	double4 direction;
	direction.w=0.0;
	direction.x=sint*cos(phi);
	direction.y=sint*sin(phi);
	direction.z=cos(theta);
//	double4 dist=quaternion_rotate(direction,quaternion_transpose(rot));
	double4 dist=direction;
	dist=M3TVec_mult(rot,direction);
	dist.x/=saxes.x; // x/a
	dist.y/=saxes.y; // y/b
	dist.z/=saxes.z; // z/c
	// Now solve x^2/a^2+y^2/b^2+z^2/c^2=1/r^2 => r = sqrt(1/(x^2/a^2+y^2/b^2+z^2/c^2))
	return 1.0/sqrt(dot(dist,dist));
}

double EllipsoidRminVec(double4 direction, double4 saxes, double* rot)
{
//	double4 dist=quaternion_rotate(normalize(direction),quaternion_transpose(rot));
	double4 dist=M3TVec_mult(rot,normalize(direction));
	dist.x/=saxes.x; // x/a
	dist.y/=saxes.y; // y/b
	dist.z/=saxes.z; // z/c
	// Now solve x^2/a^2+y^2/b^2+z^2/c^2=1/r^2 => r = sqrt(1/(x^2/a^2+y^2/b^2+z^2/c^2))
	return 1.0/sqrt(dot(dist,dist));
}

double IA(double4 dir_area, double4 invsaxes2)
{
	double a=1.0/(dir_area.x*dir_area.x+dir_area.y*dir_area.y+dir_area.z*dir_area.z);
	double c=dir_area.x*dir_area.x*(invsaxes2.y*invsaxes2.z)+dir_area.y*dir_area.y*(invsaxes2.x*invsaxes2.z)+dir_area.z*dir_area.z*(invsaxes2.x*invsaxes2.y);
	// a*beta^2 + b*beta + c = 0
	// beta_+/- = -b/2a +/- sqrt(b^2 - 4*a*c)/2a
	// beta+ * beta_ = (-b + sqrt(b^2 - 4*a*c)) * (-b - sqrt(b^2 - 4*a*c))/(2a)^2 = (b^2 - (b^2 - 4ac))/4a^2 = 4ac/4a^2 = c/a
	return dir_area.w/sqrt(c*a);
}

/// Special round for |A|<1.5, integer fallback otherwise -- AT
inline double specialround_CL(double A)
{
	double d=0.0;
	if ((A>=1.5) || (A<=-1.5)){
		double absd=fabs(d);
		double sgn=d/absd;
		d=sgn*((unsigned int)(absd+0.5));
	} else{
		if (A>=0.5) d=1.0;
			else if (A<=-0.5) d=-1.0;
	}
	return d;
}

inline void apply_PBCs(double4* rvec, __constant Simulation_Attribs* configuration)
{
	switch(configuration->latticetype){
		case 1:
		case 2: // Rectangular box (PBCs on arbitrary axes)
			// Apply PBC correction if PBCs are used along that axis
			if(configuration->PBCs[0]>0) (*rvec).x -= configuration->boxlength[0]*specialround_CL((*rvec).x/configuration->boxlength[0]);
			if(configuration->PBCs[1]>0) (*rvec).y -= configuration->boxlength[1]*specialround_CL((*rvec).y/configuration->boxlength[1]);
			if(configuration->PBCs[2]>0) (*rvec).z -= configuration->boxlength[2]*specialround_CL((*rvec).z/configuration->boxlength[2]);
			break;
		case 4: // Cylindrical box (optional PBCs on z-axis only)
			if(configuration->PBCs[2]>0) (*rvec).z -= configuration->boxlength[2]*specialround_CL((*rvec).z/configuration->boxlength[2]);
			break;
		case 3: // Spherical box (PBCs not yet supported)
		default: // Default condition: do nothing
			break;
	}
}

)

