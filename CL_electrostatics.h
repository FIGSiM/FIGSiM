/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

CL_TO_STRING(
inline double phi_q_mu(double4* thermu, __global double4* dipole, double4* RF);
inline double phi_q(__global Element_Dynamic* element_dyn, __constant Element_Static* element_stat, __global Group_Type* groups, __global double4* charges, double4* thermu, double4 *RF);
inline double4 E_mu_q(__global Element_Dynamic* element_dyn, __constant Element_Static* element_stat, __global Group_Type* groups, __global double4* charges, double4* thermu, double4 *RF);
inline double4 E_mu(__global Element_Dynamic* element_dyn, __constant Element_Static* element_stat, __global Group_Type* groups, __global double4* charges, double4* thermu, double distance, double4 *RF);

inline double phi_q_mu(double4* thermu, __global double4* dipole, double4* RF)
{
	// Add interacting charge to reaction field
	*RF += *dipole;
	// if dipole points in direction of charge to dipole vector then negative dipole charge is closer to charge q ...
	return -dot(*dipole,*thermu)/qpwr_cl(qvec_norm(*thermu),3); // attractive potentials are negative, which this one should be with a positive charge q
}

inline double phi_q(__global Element_Dynamic* element_dyn, __constant Element_Static* element_stat, __global Group_Type* groups, __global double4* charges, double4* thermu, double4 *RF)
{
	double result=0.0;
	(*RF).x=0.0; (*RF).y=0.0; (*RF).z=0.0; (*RF).w=0.0;
	result+=phi_q_mu(thermu,&(element_dyn->dipole),RF);
	for(unsigned int j=0; j<element_stat->nr_charges; j++){
		// thermu points from m-th charge on i to center of k
		double4 rq=*thermu+charges[j+element_dyn->charges_start];
		// Add interacting charge to reaction field
		double4 rd=rq;
		if(element_dyn->group_nr>=0){
			if(groups[element_dyn->group_nr].group_dipole>0) rd=element_dyn->center+charges[j+element_dyn->charges_start]-groups[element_dyn->group_nr].center;
		}
		*RF += (rd*rd.w);
		result+=rq.w/qvec_norm(rq); // negative (attractive) for opposite charges
	}
	return result;
}

inline double4 E_mu_q(__global Element_Dynamic* element_dyn, __constant Element_Static* element_stat, __global Group_Type* groups, __global double4* charges, double4* thermu, double4 *RF)
{
	double4 result;
	result.x=0.0; result.y=0.0; result.z=0.0; result.w=0.0;
	// element k is the one with charges
	// rmu points from dipole on i to center of k
	for(unsigned int j=0; j<element_stat->nr_charges; j++){
		// thermu points from dipole on i to center of k
		double4 rq=*thermu+charges[j+element_dyn->charges_start]; // vector from dipole on i to charge on k ...
		// Add interacting charge to reaction field
		double4 rd=rq;
		if(element_dyn->group_nr>=0){
			if(groups[element_dyn->group_nr].group_dipole>0) rd=element_dyn->center+charges[j+element_dyn->charges_start]-groups[element_dyn->group_nr].center;
		}
		*RF += (rd*rd.w);
		// if dipole points in direction of dipole to charge vector then positive dipole charge is closer to charge q ...
		result+=rq*rq.w/qpwr_cl(qvec_norm(rq),3); // potential should be negative (attractive) when charge q is negative
	}
	return result;
}

inline double4 E_mu(__global Element_Dynamic* element_dyn, __constant Element_Static* element_stat, __global Group_Type* groups, __global double4* charges, double4* thermu, double distance, double4 *RF)
{
	double4 result;
	result.x=0.0; result.y=0.0; result.z=0.0; result.w=0.0;
	(*RF).x=0.0; (*RF).y=0.0; (*RF).z=0.0; (*RF).w=0.0;
	
	double4 r_hat = *thermu;
	double r=distance;
	r_hat/=r; // r_hat now really is r_hat
	
	// Add interacting dipole to reaction field
	*RF += element_dyn->dipole;
	
	// Calculate the dipole-dipole energy (in picoerg)
	result += (element_dyn->dipole-r_hat*3.0*dot(element_dyn->dipole,r_hat))/(r*r*r); // will be negative (attractive) if both dipoles point in same direction congruent to distance vector
	
	if(element_stat->nr_charges>0) result+=E_mu_q(element_dyn,element_stat,groups,charges,thermu,RF);
	return result;
}

)

