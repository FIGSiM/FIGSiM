/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

#include "CL_datatypes.h"
#include "CL_helperfunctions.h"
#include "CL_electrostatics.h"
#include "CL_lennardjones.h"

CL_TO_STRING(

// uncomment for element rotations in the lab frame
// #define LAB_FRAME_ROTATION

__kernel void calculate_distances(__constant Simulation_Attribs* configuration,
                                  __global Element_Dynamic* element_dyn,
                                  __global One_Move* moves,
                                  __global Group_Type* groups,
                                  __global double4* output,
                                  const unsigned int element_count GUID_ARG)
{
	unsigned int gi=get_global_id(0);
	unsigned int i=gi;
	gi<<=1;
	if(i<2*element_count){
		One_Move trialmove=moves[0];
		Element_Dynamic element_kk=element_dyn[trialmove.kk];
		if(i>=element_count){ // apply trial move (rotation is not needed for distances ...)
			i-=element_count;
			element_kk.center+=trialmove.deltatrans;
		}
		//Calculate distances
		double4 dist2center;
		double4 temp;
		double4 rmu;
		double rdist=0.0;
		double rdist2;
		double rdist2lj;
		unsigned int calc=(i!=trialmove.kk);
		if(element_kk.group_nr>=0) dist2center=groups[element_kk.group_nr].center; else dist2center=element_kk.center;
		if(element_dyn[i].group_nr>=0){
			calc*=(element_dyn[i].group_nr!=element_kk.group_nr);
			temp=groups[element_dyn[i].group_nr].center-dist2center;
			apply_PBCs(&temp,configuration);
			rmu=element_dyn[i].center-groups[element_dyn[i].group_nr].center+temp;
			rmu-=(element_kk.center-groups[element_kk.group_nr].center)*(element_kk.group_nr>=0);
			rdist2=temp.x*temp.x+temp.y*temp.y+temp.z*temp.z;
			rdist2lj=rmu.x*rmu.x+rmu.y*rmu.y+rmu.z*rmu.z;
		} else{
			if(element_kk.group_nr>=0){
				calc*=(element_dyn[i].group_nr!=element_kk.group_nr);
				temp=element_dyn[i].center-dist2center;
				apply_PBCs(&temp,configuration);
				rmu=temp-element_kk.center-groups[element_kk.group_nr].center;
				rdist2=temp.x*temp.x+temp.y*temp.y+temp.z*temp.z;
				rdist2lj=rmu.x*rmu.x+rmu.y*rmu.y+rmu.z*rmu.z;
			} else{
				rmu=element_dyn[i].center-element_kk.center; // points from i to j
				apply_PBCs(&rmu,configuration);
				rdist2lj=rmu.x*rmu.x+rmu.y*rmu.y+rmu.z*rmu.z;
				rdist2=rdist2lj;
			}
		}
		calc*=(rdist2>=EPS_CL);
		rdist=sqrt(rdist2lj);
		
		temp.x=rdist;
		temp.y=rdist2;
		temp.z=rdist2lj;
		temp.w=calc;
		
		output[gi]=rmu;
		output[gi+1]=temp;
	}
}

__kernel void calculate_lennardjones(__constant Simulation_Attribs* configuration,
                                     __global Element_Dynamic* element_dyn,
                                     __constant Element_Static* element_stat,
                                     __global One_Move* moves,
                                     __global double4* distances,
                                     __global double* pre_params,
                                     __global double* output,
                                     const unsigned int element_count GUID_ARG)
{
	unsigned int gi=get_global_id(0);
	unsigned int i=gi;
	gi<<=1;
	if(i<2*element_count){
		One_Move trialmove=moves[0];
		Element_Dynamic element_kk=element_dyn[trialmove.kk];
		if(i>=element_count){ // apply trial move
			i-=element_count;
			element_kk.center+=trialmove.deltatrans;
#ifdef LAB_FRAME_ROTATION
			M3mult(trialmove.deltarot,element_kk.rot,element_kk.rot); // rotation in lab frame
#else
			M3mult(element_kk.rot,trialmove.deltarot,element_kk.rot); // rotation in element frame
#endif
		}
		double4 rmu=distances[gi];
		double4 dists=distances[gi+1];
		double VLJ=0.0;
		if(((dists.z<configuration->ljcut2) || (configuration->cut==0)) && (dists.w>EPS_CL)){
			switch(configuration->vdwtype){
				default:
				case 1:
					VLJ=simpleVLJ_notouch(&element_kk,&element_dyn[i],configuration,pre_params,dists.x);
					break;
				case 2:
					VLJ=VLJ_notouch(&element_kk,&element_dyn[i],configuration,pre_params,dists.x);
					break;
				case 3:
					VLJ=VSS_notouch(&element_kk,&element_dyn[i],configuration,pre_params,dists.x);
					break;
				case 4:
					VLJ=simpleVLJ_touch(&element_kk,&element_dyn[i],element_stat,configuration,pre_params,&rmu,dists.x);
					break;
				case 5:
					VLJ=VSE_touch(&element_kk,&element_dyn[i],element_stat,configuration,pre_params,&rmu,dists.x);
					break;
			}
		}
		output[gi]=VLJ;
	}
}

__kernel void calculate_electrostatics(__constant Simulation_Attribs* configuration,
                                       __global Element_Dynamic* element_dyn,
                                       __constant Element_Static* element_stat,
                                       __global One_Move* moves,
                                       __global double4* distances,
                                       __global Group_Type* groups,
                                       __global double4* initial_charges,
                                       __global double4* charges,
                                       __global double* output,
                                       const unsigned int element_count GUID_ARG)
{
	unsigned int gi=get_global_id(0);
	unsigned int i=gi;
	gi<<=1;
	if(i<2*element_count){
		One_Move trialmove=moves[0];
		Element_Dynamic element_kk=element_dyn[trialmove.kk];
		if(i>=element_count){
			i-=element_count;
			element_kk.center+=trialmove.deltatrans;
#ifdef LAB_FRAME_ROTATION
			M3mult(trialmove.deltarot,element_kk.rot,element_kk.rot); // rotation in lab frame
#else
			M3mult(element_kk.rot,trialmove.deltarot,element_kk.rot); // rotation in element frame
#endif
			element_kk.dipole=M3Vec_mult(element_kk.rot,element_stat[element_kk.type].initial_dipole);
		}
		double4 rmu=distances[gi];
		double4 dists=distances[gi+1];
		
		double4 rq;
		double4 RF;
		
		double Ves=0.0;
		double VRF=0.0;
		
		unsigned int updateRF=(fabs(configuration->rf_correction)>EPS_CL);
		
		if((dists.w>EPS_CL) && ((dists.y<configuration->escut2) || (configuration->cut==0))){
			double4 pos=rmu; // no offcenter dipole currently
			Ves += dot(element_kk.dipole,E_mu(&element_dyn[i],&element_stat[element_dyn[i].type],groups,charges,&pos,dists.x,&RF));
			if(updateRF>0) VRF -= dot(RF,element_kk.dipole);
			for(unsigned int j=0; j<element_stat[element_kk.type].nr_charges; j++){
				rq=M3Vec_mult(element_kk.rot,initial_charges[j+element_stat[element_kk.type].charges_start]);
				pos=rmu-rq;
				Ves += rq.w*phi_q(&element_dyn[i],&element_stat[element_dyn[i].type],groups,charges,&pos,&RF);
				double4 rd=(element_kk.center+rq);
				RF.w=0.0;
				VRF -= updateRF*dot(RF,rd)*rq.w;
			}
		}
		Ves += VRF*(configuration->rf_correction);
		Ves *= configuration->in2;
		output[gi+1]=Ves;
	}
}

__kernel void metropolis_condition(__constant Simulation_Attribs* configuration,
                                   __global Element_Dynamic* element_dyn,
                                   __constant Element_Static* element_stat,
                                   __global One_Move* moves,
                                   __global double4* distances,
                                   __global Group_Type* groups,
                                   __global double4* initial_charges,
                                   __global double4* charges,
                                   __global double* output,
                                   const unsigned int nr_moves,
                                   const unsigned int element_count GUID_ARG)
{
	if(get_global_id(0)==0){
		One_Move trialmove=moves[0];
		double deltaES=0.0;
		double deltaLJ=0.0;
		for(unsigned int i=0; i<element_count; i++){
			deltaES-=output[i<<1];
			deltaLJ-=output[(i<<1)+1];
		}
		for(unsigned int i=element_count; i<2*element_count; i++){
			deltaES+=output[(i<<1)];
			deltaLJ+=output[(i<<1)+1];
		}
		double deltaV=deltaES+deltaLJ;
		// Metropolis condition
		if(exp(-deltaV*configuration->beta) > trialmove.rand){
			Element_Dynamic element_kk=element_dyn[trialmove.kk];
			element_dyn[trialmove.kk].center+=trialmove.deltatrans;
#ifdef LAB_FRAME_ROTATION
			M3mult(trialmove.deltarot,element_kk.rot,element_kk.rot); // rotation in lab frame
#else
			M3mult(element_kk.rot,trialmove.deltarot,element_kk.rot); // rotation in element frame (local copy)
#endif
			element_dyn[trialmove.kk].rot[0]=element_kk.rot[0];
			element_dyn[trialmove.kk].rot[1]=element_kk.rot[1];
			element_dyn[trialmove.kk].rot[2]=element_kk.rot[2];
			element_dyn[trialmove.kk].rot[3]=element_kk.rot[3];
			element_dyn[trialmove.kk].rot[4]=element_kk.rot[4];
			element_dyn[trialmove.kk].rot[5]=element_kk.rot[5];
			element_dyn[trialmove.kk].rot[6]=element_kk.rot[6];
			element_dyn[trialmove.kk].rot[7]=element_kk.rot[7];
			element_dyn[trialmove.kk].rot[8]=element_kk.rot[8];
			element_dyn[trialmove.kk].dipole=M3Vec_mult(element_kk.rot,element_stat[element_kk.type].initial_dipole);
			for(unsigned int j=0; j<element_stat[element_kk.type].nr_charges; j++){
				charges[j+element_kk.charges_start]=M3Vec_mult(element_kk.rot,initial_charges[j+element_stat[element_kk.type].charges_start]);
			}
			output[(element_count<<2)]+=trialmove.deltatrans.w; // tmovtemp
			output[(element_count<<2)+1]+=trialmove.deltarot[9]; // rmovtemp
			output[(element_count<<2)+2]=0;
			output[(element_count<<2)+3]+=deltaES;
			output[(element_count<<2)+4]+=deltaLJ;
			output[(element_count<<2)+5]=0;
			output[(element_count<<2)+6]+=1.0;
		}
		unsigned new_count=trialmove.counter+1;
		if(new_count<nr_moves){ // otherwise we're done here ...
			moves[0]=moves[new_count];
			moves[0].counter=new_count;
		}
	}
}

)

