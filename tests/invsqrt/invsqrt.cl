/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

// We have to hide OpenCL preprocessing directives in this one line here ;-)
CL_TO_STRING( #ifdef KHR_DP_EXTENSION\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n#else\n#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n#endif\n

__kernel void invsqrt(__global double* input, __global double* output, const unsigned int count)
{
	int i = get_global_id(0);
	if(i<count) output[i] = rsqrt(input[i]);
}

);
