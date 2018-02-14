/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

#include "datatypes.h"
#include "functions.h"

CL_TO_STRING(

__kernel void invsqrt_cl(__global double* input, __global double* output, const unsigned int count GUID_ARG)
{
	unsigned int i = get_global_id(0);
	if(i<count) output[i] = 1.0/sqrt(input[i]);
}

__kernel void volume_inside_count(__global Frame_Attrib* frame, __global double4* random_nrs, __global Element_Dynamic* element_dyn, __global Element_Static* element_stat, __global short* output_count, const unsigned int element_count, const unsigned int count GUID_ARG)
{
	unsigned int i = get_global_id(0);
	if(i<count){
		double4 point=frame->box_corner;
		point.x+=random_nrs[i].x*frame->box.x;
		point.y+=random_nrs[i].y*frame->box.y;
		point.z+=random_nrs[i].z*frame->box.z;
		point=quaternion_rotate(point,frame->rot)+frame->center;
		
		short inside=0;
		for(unsigned int j=0; ((j<element_count) && (inside==0)); j++) inside=Point_in_Ellipsoid(point,element_stat[j].saxes,element_dyn[j].center,element_dyn[j].rot);
		output_count[i]=inside;
	}
}

);


