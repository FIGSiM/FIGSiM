/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

CL_TO_STRING(

inline double touch(Element_Dynamic* element_kk, __global Element_Dynamic* element_i, __constant Element_Static* element_stat, double4* thermu, double distance);
inline double VSS_notouch(Element_Dynamic* element_kk, __global Element_Dynamic* element_i, __constant Simulation_Attribs* configuration, __global double* pre_params, double distance);
inline double simpleVLJ_notouch(Element_Dynamic* element_kk, __global Element_Dynamic* element_i, __constant Simulation_Attribs* configuration, __global double* pre_params, double distance);
inline double VLJ_notouch(Element_Dynamic* element_kk, __global Element_Dynamic* element_i, __constant Simulation_Attribs* configuration, __global double* pre_params, double distance);
inline double simpleVLJ_touch(Element_Dynamic* element_kk, __global Element_Dynamic* element_i, __constant Element_Static* element_stat, __constant Simulation_Attribs* configuration, __global double* pre_params, double4* thermu, double distance);
inline double VSE_touch(Element_Dynamic* element_kk, __global Element_Dynamic* element_i, __constant Element_Static* element_stat, __constant Simulation_Attribs* configuration, __global double* pre_params, double4* thermu, double distance);

/*!
 * Touch, an algorithm designed for determining the closest contact distance between two ellipsoids at arbitrary angles,
 * reducing the three-dimensional problem to a one-dimensional problem in terms of a scalar variable lambda, then iteratively
 * solving a system of linear equations to optimize lambda. Initially written in Matlab by BHR, ported to C++ by RSB in 02/09, then
 * harmonized with current codebase by LEJ in 04/09, and optimized by LEJ in 01/10, finally fixed by AT in 02/11, and ported to
 * OpenCL by AT in 07/12
 *
 * Return value in distance is the effective LJ sigma for the molecules at their current positions
 */
inline double touch(Element_Dynamic* element_kk, __global Element_Dynamic* element_i, __constant Element_Static* element_stat, double4* thermu, double distance)
{
	//Normalize distance vector
	//t is Rmu, nt is distance
	double invdist = 1.0/distance;
	double invdist2 = invdist*invdist;
	unsigned int sing=0;
	
	//Minimization vectors
	double4 z=*thermu*invdist;
	z.w=0.0;
	double4 V;
	
	//Create A and B matrices, 
	double A[9]; //3x3 with saxes for the kkth oid on the diagonal.
	double B[9]; //3x3 with saxes for the ith oid on the diagonal.
	
	A[0]=(element_stat[element_kk->type].saxes_vdw.x*element_stat[element_kk->type].saxes_vdw.x*invdist2);
	A[1]=0.0; A[2]=0.0; A[3]=0.0;
	A[4]=(element_stat[element_kk->type].saxes_vdw.y*element_stat[element_kk->type].saxes_vdw.y*invdist2);
	A[5]=0.0; A[6]=0.0; A[7]=0.0;
	A[8]=(element_stat[element_kk->type].saxes_vdw.z*element_stat[element_kk->type].saxes_vdw.z*invdist2);
	
	B[0]=(element_stat[element_i->type].saxes_vdw.x*element_stat[element_i->type].saxes_vdw.x*invdist2);
	B[1]=0.0; B[2]=0.0; B[3]=0.0;
	B[4]=(element_stat[element_i->type].saxes_vdw.y*element_stat[element_i->type].saxes_vdw.y*invdist2);
	B[5]=0.0; B[6]=0.0; B[7]=0.0;
	B[8]=(element_stat[element_i->type].saxes_vdw.z*element_stat[element_i->type].saxes_vdw.z*invdist2);
	
	
	//Create lab frame version of A and B matrices.
	double A_LF[9];
	double B_LF[9];
	double rot[9];
	for(unsigned m=0; m<9; m++) rot[m]=element_kk->rot[m];
	M3mult(rot,A,A_LF);
	M3multT(A_LF,rot,A_LF);
	for(unsigned m=0; m<9; m++) rot[m]=element_i->rot[m];
	M3mult(rot,B,B_LF);
	M3multT(B_LF,rot,B_LF);
	
	//Create CI matrix and its inverse.
	double CI[9];
	
	double d=1.0/M3det(B_LF);
	if(d>1E10) d=1E10;
	if(d>1.0){
		for(unsigned int i=0; i<9; i++){
			A_LF[i]*=d;
			B_LF[i]*=d;
		}
		z*=sqrt(d);
	}
	
	double lamx = 0.5;
	double lamoid = lamx;
	double Var = 1.0; // trying different forms found 6 loops gives 7 figs for Fx
	double numerator;
	double denominator;
	double xlx;
	double Var_temp;
	while((Var > 0.000001) && (sing==0)){ //loop until variance in distance is sufficiently small
		lamoid = lamx;
		//Populate CI matrix
		for(unsigned int j=0; j<9; j++) CI[j] = lamx*(B_LF[j]-A_LF[j])+A_LF[j];
		
		//Invert CI matrix
		M3Inv(CI,&sing);
		
		// Solve z = CI*V for V using inverse
		V = M3Vec_mult(CI,z);
		// numerator=V*A_LF*V
		numerator = V.x*(V.x*A_LF[0]+V.y*A_LF[3]+V.z*A_LF[6])+V.y*(V.x*A_LF[1]+V.y*A_LF[4]+V.z*A_LF[7])+V.z*(V.x*A_LF[2]+V.y*A_LF[5]+V.z*A_LF[8]);
		// denominator=V*B_LF*V
		denominator = V.x*(V.x*B_LF[0]+V.y*B_LF[3]+V.z*B_LF[6])+V.y*(V.x*B_LF[1]+V.y*B_LF[4]+V.z*B_LF[7])+V.z*(V.x*B_LF[2]+V.y*B_LF[5]+V.z*B_LF[8]);
		
		//Calculate minimization parameter lambda
		xlx = sqrt(numerator/denominator);
		lamx = xlx/(1.0+xlx);
		Var_temp = lamx-lamoid;
		Var = Var_temp*Var_temp;
	}
	
	//Reconstruct CI and run a final iteration once converged
	for(unsigned int j=0; j<9; j++) CI[j] = lamx*(B_LF[j]-A_LF[j])+A_LF[j];
	
	//Invert CI matrix
	M3Inv(CI,&sing);
	if(sing==0){
		// Solve z = CI*V for V using inverse
		V = M3Vec_mult(CI,z);
		//Calculate inverse of contact distance.
		double Fx = lamx*(1.0-lamx)*dot(z,V);
		if(Fx>0.0) return 1.0/Fx; else return -0.0;
	} else return -0.0;
}

inline double touch_sphere(__global Simulation_Attribs* attrib, double4 thermu, double distance)
{
	//Normalize distance vector
	//t is Rmu, nt is distance
	double invdist = 1.0/distance;
	double invdist2 = invdist*invdist;
	unsigned int sing=0;
	
	//Minimization vectors
	double4 z=thermu*invdist;
	z.w=0.0;
	double4 V;
	
	//Create A and B matrices, 
	double A[9]; //3x3 with saxes for the kkth oid on the diagonal.
	double B_LF[9]; // test sphere is always in lab frame (spherical symmetry)
	
	A[0]=(attrib->saxes.x*attrib->saxes.x*invdist2);
	A[1]=0.0; A[2]=0.0; A[3]=0.0;
	A[4]=(attrib->saxes.y*attrib->saxes.y*invdist2);
	A[5]=0.0; A[6]=0.0; A[7]=0.0;
	A[8]=(attrib->saxes.z*attrib->saxes.z*invdist2);
	
	B_LF[0]=(attrib->rp*attrib->rp*invdist2);
	B_LF[1]=0.0; B_LF[2]=0.0; B_LF[3]=0.0;
	B_LF[4]=(attrib->rp*attrib->rp*invdist2);
	B_LF[5]=0.0; B_LF[6]=0.0; B_LF[7]=0.0;
	B_LF[8]=(attrib->rp*attrib->rp*invdist2);
	
	
	//Create lab frame version of A and B matrices.
	double A_LF[9];
	double rot[9];
	for(unsigned m=0; m<9; m++) rot[m]=attrib->rot[m];
	M3mult(rot,A,A_LF);
	M3multT(A_LF,rot,A_LF);
	
	//Create CI matrix and its inverse.
	double CI[9];
	
	double d=1.0/M3det(B_LF);
	if(d>1E10) d=1E10;
	if(d>1.0){
		for(unsigned int i=0; i<9; i++){
			A_LF[i]*=d;
			B_LF[i]*=d;
		}
		z*=sqrt(d);
	}
	
	double lamx = 0.5;
	double lamoid = lamx;
	double Var = 1.0; // trying different forms found 6 loops gives 7 figs for Fx
	double numerator;
	double denominator;
	double xlx;
	double Var_temp;
	while((Var > 0.000001) && (sing==0)){ //loop until variance in distance is sufficiently small
		lamoid = lamx;
		//Populate CI matrix
		for(unsigned int j=0; j<9; j++) CI[j] = lamx*(B_LF[j]-A_LF[j])+A_LF[j];
		
		//Invert CI matrix
		M3Inv(CI,&sing);
		
		// Solve z = CI*V for V using inverse
		V = M3Vec_mult(CI,z);
		// numerator=V*A_LF*V
		numerator = V.x*(V.x*A_LF[0]+V.y*A_LF[3]+V.z*A_LF[6])+V.y*(V.x*A_LF[1]+V.y*A_LF[4]+V.z*A_LF[7])+V.z*(V.x*A_LF[2]+V.y*A_LF[5]+V.z*A_LF[8]);
		// denominator=V*B_LF*V
		denominator = V.x*(V.x*B_LF[0]+V.y*B_LF[3]+V.z*B_LF[6])+V.y*(V.x*B_LF[1]+V.y*B_LF[4]+V.z*B_LF[7])+V.z*(V.x*B_LF[2]+V.y*B_LF[5]+V.z*B_LF[8]);
		
		//Calculate minimization parameter lambda
		xlx = sqrt(numerator/denominator);
		lamx = xlx/(1.0+xlx);
		Var_temp = lamx-lamoid;
		Var = Var_temp*Var_temp;
	}
	
	//Reconstruct CI and run a final iteration once converged
	for(unsigned int j=0; j<9; j++) CI[j] = lamx*(B_LF[j]-A_LF[j])+A_LF[j];
	
	//Invert CI matrix
	M3Inv(CI,&sing);
	if(sing==0){
		// Solve z = CI*V for V using inverse
		V = M3Vec_mult(CI,z);
		//Calculate inverse of contact distance.
		double Fx = lamx*(1.0-lamx)*dot(z,V);
		if(Fx>0.0) return 1.0/Fx; else return -0.0;
	} else return -0.0;
}

/*!
 * Calculation of Lennard-Jones interactions. Notouch variants can only handle spheres, Touch variants handle anisotropy
 * using the method of Perram and Wertheim (J. Comp. Phys, 1985). Several functions are available, providing tradeoffs
 * between generality and speed
 *
 * Return value is the LJ energy for one pair of ellipsoids
 *
 * Soft sphere - only calculate nuclear repulsion term of LJ potential
 * pre_params[ikk] ... pre_eps
 * pre_params[ikk+1] ... pre_sigma
 * pre_params[ikk+2] ... pre_touch
 * pre_params[ikk+3] ... tcut
 */
inline double VSS_notouch(Element_Dynamic* element_kk, __global Element_Dynamic* element_i, __constant Simulation_Attribs* configuration, __global double* pre_params, double distance)
{
	unsigned int ikk = (element_i->type*configuration->num_element_types+element_kk->type)<<2;
	//Calculate energy
	return pre_params[ikk]*configuration->r*qpwr_cl(pre_params[ikk+1]/distance,configuration->LJexp[0]<<1); // bitshift left by one is multiplication by 2;
}

/// Simple Lennard-Jones potential (nuclear repulsion power is twice that of dispersion, e.g. LJ 12-6 potential)
inline double simpleVLJ_notouch(Element_Dynamic* element_kk, __global Element_Dynamic* element_i, __constant Simulation_Attribs* configuration, __global double* pre_params, double distance)
{
	unsigned int ikk = (element_i->type*configuration->num_element_types+element_kk->type)<<2;
	double disp = qpwr_cl(pre_params[ikk+1]/distance,configuration->LJexp[0]);
	
	//Calculate energy
	return pre_params[ikk]*(configuration->r*disp*(disp-1.0));
}

/// General LJ potential including optional Bruce correction (attenuating Gaussian function in bottom of well)
inline double VLJ_notouch(Element_Dynamic* element_kk, __global Element_Dynamic* element_i, __constant Simulation_Attribs* configuration, __global double* pre_params, double distance)
{
	unsigned int ikk = (element_i->type*configuration->num_element_types+element_kk->type)<<2;
	
	//Calculate energy
	double rmrmin = distance - configuration->Solvent[2];
	double sx = pre_params[ikk+1]/distance;
	
	//LJ + Bruce correction around r=rmin. Tanh functions have been replaced with a Gaussian for speed
	return pre_params[ikk]*configuration->r*(qpwr_cl(sx,configuration->LJexp[1])-qpwr_cl(sx,configuration->LJexp[0])) + configuration->Solvent[1]*exp(configuration->Solvent[0]*rmrmin*rmrmin);
}

/// Simplified VLJ, but using Touch algorithm to handle non-spherical ellipsoids. This ONLY works for even powers!
/// -- extended to be able to use square root LJ epsilon texture if present for either party
inline double simpleVLJ_touch(Element_Dynamic* element_kk, __global Element_Dynamic* element_i, __constant Element_Static* element_stat, __constant Simulation_Attribs* configuration, __global double* pre_params, double4* thermu, double distance)
{
	double Rx; //reduced distance (r0/r)^2
	unsigned int ikk = (element_i->type*configuration->num_element_types+element_kk->type)<<2;
	double epsilon = pre_params[ikk];
	
	//Determine whether touch should be run or the average radius should be used
	if((pre_params[ikk+2]) && ((!configuration->touchtrunc) || (distance<pre_params[ikk+3]))){
		Rx = touch(element_kk, element_i, element_stat, thermu, distance); // Square of reduced distance corrected based on closest contact
	} else{
		Rx = pre_params[ikk+1]/distance; // (r0/r)
		Rx *= Rx;
	}
	
	//Calculate parameters
	double disp = qpwr_cl(Rx,configuration->LJexp[0]>>1); // Since Rx is already squared, the small term must be halved (bitshift right by one is integer/2 - AT)
	
/*	//Calculate energy
	__global Element_Type* type_i=Elements[i].MyType;
	__global Element_Type* type_kk=Elements[kk].MyType;
	if(type_i->eps_texture || type_kk->eps_texture){
		Vec3 rmu_i_frame=Elements[i].rot.M3Transpose()*(Rmu*(-1.0)); // Rmu is vector from element kk to element i
		Vec3 rmu_kk_frame=Elements[kk].rot.M3Transpose()*Rmu;
		double theta, phi, cost;
		unsigned int t,p;
		double eps_i=determine_epsilon(type_i,sqrt(Ellipsoid_Cross_Section(type_i->saxes,rmu_i_frame))/pi);
		if(type_i->eps_texture){
			cost=VecDist2ThetaPhi(rmu_i_frame,distance,theta,phi);
			t=(unsigned int)((double)(1.0-cost)/2.0*theta_res-0.5);
			p=(phi_res+(int)((double)phi/(2.0*pi)*phi_res))%phi_res;
			eps_i*=type_i->eps_texture[t*phi_res+p];
		}
		double eps_kk = determine_epsilon(type_kk,sqrt(Ellipsoid_Cross_Section(type_kk->saxes,rmu_kk_frame))/pi);
		if(type_i->eps_texture){
			cost=VecDist2ThetaPhi(rmu_kk_frame,distance,theta,phi);
			t=(unsigned int)((double)(1.0-cost)/2.0*theta_res-0.5);
			p=(phi_res+(int)((double)phi/(2.0*pi)*phi_res))%phi_res;
			eps_kk*=type_kk->eps_texture[t*phi_res+p];
		}
		
		epsilon=eps_i*eps_kk;
	}*/
	return epsilon*(configuration->r*disp*(disp-1.0));
}

/// Soft-sphere potential, but using Touch algorithm to handle non-spherical ellipsoids
inline double VSE_touch(Element_Dynamic* element_kk, __global Element_Dynamic* element_i, __constant Element_Static* element_stat, __constant Simulation_Attribs* configuration, __global double* pre_params, double4* thermu, double distance)
{
	double Rx; //reduced distance (r0/r)^2
	unsigned int ikk = element_i->type*configuration->num_element_types+element_kk->type;
	
	//Determine whether touch should be run or the average radius should be used
	if((pre_params[ikk+2]) && ((configuration->touchtrunc>0) || (distance<pre_params[ikk+3]))){
		Rx = touch(element_kk, element_i, element_stat, thermu, distance); //Reduced distance (r0/r)^2, corrected based on closest contact
	} else{
		Rx = pre_params[ikk+1]/distance; // (r0/r)
		Rx *= Rx;
	}
	
	//Calculate energy
	return pre_params[ikk]*configuration->r*qpwr_cl(Rx,configuration->LJexp[0]); // since Rx is already squared, this gives the larger (repulsive term)
}

)

