/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

#include "CL_datatypes.h"
#include "CL_helperfunctions.h"
#include "CL_lennardjones.h"

CL_TO_STRING(

__kernel void match_epsilon(__global Simulation_Attribs* attrib,
                            __global Element_Dynamic* element_dyn,
                            __constant Element_Static* element_stat,
                            __global double* r_diff,
                            __global double* output,
                            const unsigned int k,
                            const unsigned int element_count,
                            const unsigned int count_t,
                            const unsigned int count_p,
                            const unsigned int count_r GUID_ARG)
{
	unsigned int i=get_global_id(0); // = t*count_p+p
	if(i<count_t*count_p+2){
		unsigned int t=i/count_p;
		unsigned int p=i%count_p;
		unsigned int idx=(i*count_r+(unsigned int)(k/attrib->r_overscan))<<2;
		
		double weight;
		double rmin;
		
		double neweps;
		
		double VLJ; double disp;
		
		double phi; double theta; double r;
		
		double VLJatomistic; double4 point;
		
		double4 pv; double rc; double rminc;
		
		phi=2.0*pi*(double)(p+0.5)/(double)count_p;
		double cost=1.0-2.0*(double)(t+0.5)/(double)count_t;
		double sint=sqrt(1.0-cost*cost);
		theta=acos(cost);
		if(i>=count_t*count_p){
			if(i==count_t*count_p){ // perfect theta = 0 case
				theta=0.0;
				cost=1.0;
				sint=0.0;
				phi=0.0;
			} else{ // perfect theta = pi case
				theta=pi;
				cost=-1.0;
				sint=0.0;
				phi=0.0;
			}
		}
		
		double rk=qpwr_cl((double)(k+attrib->k0*attrib->r_overscan)/((double)((count_r+attrib->k0)*attrib->r_overscan)),3);
		double rkp=qpwr_cl((double)(k+1+attrib->k0*attrib->r_overscan)/((double)((count_r+attrib->k0)*attrib->r_overscan)),3);
		r=(rk+rkp)*attrib->rmax/2.0;
		double voxelvolume=attrib->dphi_3*attrib->dcostheta*(qpwr_cl(rkp,3)-qpwr_cl(rk,3))*qpwr_cl(attrib->rmax,3);
		
		point.x=r*sint*cos(phi);
		point.y=r*sint*sin(phi);
		point.z=r*cost;
/*		double ina=1.0;
		if(attrib->inv_avg_area>0.0){
			point.w=attrib->inv_avg_area;
			ina*=IA(point,attrib->invsaxes2);
		}*/
		point.w=0.0;
		double rot[9];
		for(unsigned m=0; m<9; m++) rot[m]=attrib->rot[m];
		point=M3Vec_mult(rot,point);
		if(attrib->rp>EPS_CL){
			rmin = r*sqrt(touch_sphere(attrib,point,r)); // touch_sphere places the ellipsoid at the center
		} else rmin = EllipsoidRmin(theta,phi,attrib->saxes,rot)+attrib->rp;
		
		point=point+attrib->center;
		
		if(attrib->avg_width>EPS_CL){ // sigma-delta = avg_width => delta = sigma-avg_width
			r -= rmin - (attrib->avg_width+attrib->rp);
			rmin = attrib->avg_width+attrib->rp;
		}
		r-=r_diff[i]; // shift ellipsoid LJ so that zero-crossing of ellipsoid and fully-atomistic LJ are identical (without changing potential shape - otherwise would end up in rmin=not good)
		if(r<EPS_CL) r=EPS_CL;
		
		switch(attrib->vdwtype){
			case 1:
			case 4:
			default:
				disp = qpwr_cl(rmin/r,attrib->LJexp[0]);
				VLJ=attrib->r*disp*(disp-1.0);
				break;
			case 2:
				disp = rmin/r;
				// LJ + Bruce correction around r=rmin. Tanh functions have been replaced with a Gaussian for speed
				VLJ=attrib->r*(qpwr_cl(disp,attrib->LJexp[1])-qpwr_cl(disp,attrib->LJexp[0])) + attrib->Solvent[1]*exp(attrib->Solvent[0]*(r - attrib->Solvent[2])*(r - attrib->Solvent[2]));
				break;
			case 3:
			case 5:
				VLJ=attrib->r*qpwr_cl(rmin/r,attrib->LJexp[0]<<1); // bitshift left by one is multiplication by 2
				break;
		}
		
		// Calculate Lennard-Jones potential of fully atomistic model
		VLJatomistic=0.0;
		for(unsigned int l=0; l<element_count; l++){
			pv=point-element_dyn[l].center;
			rc=sqrt(pv.x*pv.x+pv.y*pv.y+pv.z*pv.z);
			for(unsigned m=0; m<9; m++) rot[m]=element_dyn[l].rot[m];
			rminc=EllipsoidRminVec(pv,element_stat[l].saxes_vdw,rot)+attrib->rp;
			switch(attrib->vdwtype){
				case 1:
				case 4:
				default:
					disp = qpwr_cl(rminc/rc,attrib->LJexp[0]);
					VLJatomistic+=element_stat[l].saxes_vdw.w*attrib->r*disp*(disp-1.0);
					break;
				case 2:
					disp = rminc/rc;
					// LJ + Bruce correction around r=rmin. Tanh functions have been replaced with a Gaussian for speed
					VLJatomistic+=element_stat[l].saxes_vdw.w*attrib->r*(qpwr_cl(disp,attrib->LJexp[1])-qpwr_cl(disp,attrib->LJexp[0])) + attrib->Solvent[1]*exp(attrib->Solvent[0]*(rc - attrib->Solvent[2])*(rc - attrib->Solvent[2]));
					break;
				case 3:
				case 5:
					VLJatomistic+=element_stat[l].saxes_vdw.w*attrib->r*qpwr_cl(rminc/rc,attrib->LJexp[0]<<1); // bitshift left by one is multiplication by 2
					break;
			}
		}
		weight=exp(-VLJ*attrib->weight_factor)*voxelvolume;
		neweps=VLJatomistic/VLJ;
		if((neweps>EPS_CL) && (VLJ<-0.01)){ // start from either LOD or fully atomistic surfaces
			output[idx]+=weight*neweps;
			output[idx+1]+=weight*neweps*neweps;
			output[idx+2]+=weight;
		}
	}
}

)

