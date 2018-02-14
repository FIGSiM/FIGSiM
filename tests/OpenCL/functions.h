/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

CL_TO_STRING(

inline double4 quaternion_rotate(double4 v, double4 rot);
inline double4 quaternion_transpose(double4 rot);
inline short Point_in_Ellipsoid(double4 point, double4 saxes, double4 center, double4 rot);
inline double ran_cl(int* idum, int* iv, int* iy, int* idum2);
inline double qpwr_cl(double a, unsigned int b);
inline double qqpwr_cl(double a, int b);
double EllipsoidRmin(double theta, double phi, double4 saxes, double4 rot);
double EllipsoidRminVec(double4 point, double4 saxes, double4 center, double4 rot);

// random function
// needs to be called with negative idum value to initialize iv[32], iy, and idum2 
inline double ran_cl(int* idum, int* iv, int* iy, int* idum2)
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
	if(temp > 1.0-1.2E-7) return 1.0-1.2E-7; else return temp;
}

/// Fast power algorithm that only works for positive integers
inline double qpwr_cl(double a, unsigned int b)
{
	double c = 1;
	for(unsigned int i = 0; i < b; i++) c *= a;
	return c;
}

/// Slighlty faster recursive a^b (b being positive integer) -- AT
inline double qqpwr_cl(double a, int b)
{
	if (a==0.0) return 0;
	switch(b){
		case 0: return 1.0;
		case 1: return a;
		case 2: return a*a;
		case 3: return a*a*a;
	}
	if (b%2==0) return qqpwr_cl(a*a,b>>1); else return a*qqpwr_cl(a*a,b>>1); // bitshift right by one is integer/2
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

inline short Point_in_Ellipsoid(double4 point, double4 saxes, double4 center, double4 rot)
{
	// first rotate point-to-center distance vector so that ellipsoid is aligned with coordinate system
	double4 dist=quaternion_rotate(point-center,quaternion_transpose(rot));
	dist.x=dist.x/saxes.x; // x/a
	dist.y=dist.y/saxes.y; // y/b
	dist.z=dist.z/saxes.z; // z/c
	
	// point is inside if (x^2/a^2 + y^2/b^2 + z^2/c^2 <= 1)
	if(dot(dist,dist)<=1.0) return 1;
	return 0;
}

double EllipsoidRmin(double theta, double phi, double4 saxes, double4 rot)
{
	// first rotate point-to-center distance vector so that ellipsoid is aligned with coordinate system
	double sint=sin(theta);
	double4 direction;
	direction.w=0.0;
	direction.x=sint*cos(phi);
	direction.y=sint*sin(phi);
	direction.z=cos(theta);
	double4 dist=quaternion_rotate(direction,quaternion_transpose(rot));
	dist.x/=saxes.x; // x/a
	dist.y/=saxes.y; // y/b
	dist.z/=saxes.z; // z/c
	// Now solve x^2/a^2+y^2/b^2+z^2/c^2=1/r^2 => r = sqrt(1/(x^2/a^2+y^2/b^2+z^2/c^2))
	return 1.0/sqrt(dot(dist,dist));
}

double EllipsoidRminVec(double4 point, double4 saxes, double4 center, double4 rot)
{
	double4 direction=point-center;
	direction=normalize(direction);
	double4 dist=quaternion_rotate(direction,quaternion_transpose(rot));
	dist.x/=saxes.x; // x/a
	dist.y/=saxes.y; // y/b
	dist.z/=saxes.z; // z/c
	// Now solve x^2/a^2+y^2/b^2+z^2/c^2=1/r^2 => r = sqrt(1/(x^2/a^2+y^2/b^2+z^2/c^2))
	return 1.0/sqrt(dot(dist,dist));
}


)

