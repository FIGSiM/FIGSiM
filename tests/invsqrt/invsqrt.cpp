/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*
invsqrt.cpp
	Test inverse sqrt implementations
	written by Andreas Tillack April 12, 2011
*/

#define Nr 1024*512
#define low 0.2
#define high 13.29

#include <iostream>
#include "../../Config.h"
#ifdef OpenCL
#include "../../MC_OpenCL.h"
#endif
#include "../../ScalarMat.h"
#include "../../VecMat.h"

/* Magic to include OpenCL source code
 * from external file into executable
 */
#ifdef OpenCL
#define CL_TO_STRING(CL) #CL
const char* invsqrt_cl =
#include "invsqrt.cl"
#endif

int main(){
	cout << "*** invsqrt.cpp ***\n";
	cout << "-> Creating random array of " << Nr << " doubles within range [" << low << "," << high << "].\n";
	unsigned int i,k;
	double* values = new double[Nr];
	__int32_t idum=-8;
	for(i=0; i<Nr; i++) values[i]=low+ran2(idum)*(high-low);
	cout << "<- Done.\n\n";
	
	cout << "-> Calculating inverse square root using reference 1.0/sqrt().\n";
	double* reference = new double[Nr];
	double tstart = clock();
	for(k=0; k<1000; k++)
		for(i=0; i<Nr; i++) reference[i]=1.0/sqrt(values[i]);
	double tend = clock();
	cout << "<- Done, took " << (tend-tstart)/CLOCKS_PER_SEC << " ms per " << Nr << " values.\n\n";
	
	cout << "-> Calculating inverse square root using magic_invsqr() (uses 64bit doubles on x64, 32bit floats on x32).\n";
	double* compare = new double[Nr];
	double maxerror=0.0;
	unsigned int error_max=0;
	tstart = clock();
	for(k=0; k<1000; k++)
		for(i=0; i<Nr; i++) compare[i]=magic_invsqr(values[i]);
	tend = clock();
	for(i=0; i<Nr; i++){
		if(maxerror<fabs(compare[i]-reference[i])/reference[i]){
			maxerror=fabs(compare[i]-reference[i])/reference[i];
			error_max=i;
		}
	}
	cout << "<- Done, took " << (tend-tstart)/CLOCKS_PER_SEC << " ms per " << Nr << " values (maximum relative error: " << maxerror*100 << "% at value " << values[error_max] << ").\n\n";
	
	cout << "-> Calculating inverse square root using custom_invsqrt() (uses SSE2 if available -- use compiler switch -msse2).\n";
	tstart = clock();
	for(k=0; k<1000; k++)
		for(i=0; i<Nr; i++) compare[i]=custom_invsqrt(values[i]);
	tend = clock();
	maxerror=0.0;
	error_max=0;
	for(i=0; i<Nr; i++){
		if(maxerror<fabs(compare[i]-reference[i])/reference[i]){
			maxerror=fabs(compare[i]-reference[i])/reference[i];
			error_max=i;
		}
	}
	cout << "<- Done, took " << (tend-tstart)/CLOCKS_PER_SEC << " ms per " << Nr << " values (maximum relative error: " << maxerror*100 << "% at value " << values[error_max] << ").\n\n";
#ifdef OpenCL
	cout << "-> Calculating inverse square root using OpenCL.\n";
	cl_context Context=NULL;
	cl_command_queue Queue=NULL;
	cl_device_id deviceID=GetDevice(CL_DEVICE_TYPE_ALL,Context,Queue,true);
	cl_program program = CreateProgram(deviceID,Context,&invsqrt_cl);
	int Error_Code=0;
	cl_kernel kernel = clCreateKernel(program, "invsqrt", &Error_Code);
	if(!kernel || (Error_Code!=CL_SUCCESS)){
		cout << "ERROR: Could not create compute kernel.\n";
		exit(1);
	}
	
	cl_mem input = clCreateBuffer(Context,CL_MEM_READ_ONLY,sizeof(double)*Nr,NULL,NULL);
	cl_mem output = clCreateBuffer(Context,CL_MEM_WRITE_ONLY,sizeof(double)*Nr,NULL,NULL);
	if(!input || !output){
		cout << "ERROR: Could not allocate device memory.\n";
		exit(1);
	}
	Error_Code = clEnqueueWriteBuffer(Queue,input,CL_TRUE,0,sizeof(double)*Nr,values,0,NULL,NULL);
	if(Error_Code!=CL_SUCCESS){
		cout << "ERROR: Could not copy values to device.\n";
		exit(1);
	}
	Error_Code = clSetKernelArg(kernel, 0, sizeof(cl_mem), &input);
	Error_Code|= clSetKernelArg(kernel, 1, sizeof(cl_mem), &output);
	size_t global=Nr;
	unsigned int count=Nr;
	Error_Code|= clSetKernelArg(kernel, 2, sizeof(unsigned int), &count);
	if(Error_Code!=CL_SUCCESS){
		cout << "ERROR: Failed to set kernel arguments.\n";
		exit(1);
	}
	size_t local=0;
	Error_Code = clGetKernelWorkGroupInfo(kernel, deviceID, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local), &local, NULL);
	if(Error_Code!=CL_SUCCESS){
		cout << "ERROR: Could not determine workgroup information.\n";
		exit(1);
	}
	
	tstart = clock();
	for(k=0; k<1000; k++){
		Error_Code = clEnqueueNDRangeKernel(Queue,kernel,1,NULL,&global,&local,0,NULL,NULL);
		if(Error_Code!=CL_SUCCESS){
			cout << "ERROR: Failed to execute kernel: " << Error_Code << "\n";
			exit(1);
		}
		clFinish(Queue);
	}
	Error_Code = clEnqueueReadBuffer(Queue,output,CL_TRUE,0,sizeof(double)*Nr,compare,0,NULL,NULL);
	if(Error_Code!=CL_SUCCESS){
		cout << "ERROR: Could not read results.\n";
		exit(1);
	}
	tend = clock();
	
	clReleaseMemObject(input);
	clReleaseMemObject(output);
	clReleaseProgram(program);
	clReleaseKernel(kernel);
	clReleaseCommandQueue(Queue);
	clReleaseContext(Context);
	
	maxerror=0.0;
	error_max=0;
	for(i=0; i<Nr; i++){
		if(maxerror<fabs(compare[i]-reference[i])/reference[i]){
			maxerror=fabs(compare[i]-reference[i])/reference[i];
			error_max=i;
		}
	}
	
	cout << "<- Done, took " << (tend-tstart)/CLOCKS_PER_SEC << " ms (CPU time) per " << Nr << " values (maximum relative error: " << maxerror*100 << "% at value " << values[error_max] << ").\n\n";
#endif
	
	delete[] values;
	delete[] reference;
	delete[] compare;
	
	cout << "-> Creating random array of " << Nr << " Vec3 within range [" << low << "," << high << "] (for each dimension).\n";
	Vec3* vecs = new Vec3[Nr];
	for(i=0; i<Nr; i++){
		vecs[i].vec[0]=low+ran2(idum)*(high-low);
		vecs[i].vec[1]=low+ran2(idum)*(high-low);
		vecs[i].vec[2]=low+ran2(idum)*(high-low);
	}
	cout << "<- Done.\n\n";
	
	cout << "-> Calculating normalized vectors using reference Vec3/Vec3.V3Norm().\n";
	Vec3* calc_vecs = new Vec3[Nr];
	tstart = clock();
	for(k=0; k<1000; k++)
		for(i=0; i<Nr; i++) calc_vecs[i]=vecs[i]/vecs[i].V3Norm();
	tend = clock();
	maxerror=0.0;
	for(i=0; i<Nr; i++){
		if(maxerror<fabs(calc_vecs[i].V3Norm()-1.0)) maxerror=fabs(calc_vecs[i].V3Norm()-1.0);
	}
	cout << "<- Done, took " << (tend-tstart)/CLOCKS_PER_SEC << " ms per " << Nr << " vectors (maximum deviation from unit length: " << maxerror*100 << "%).\n\n";
	
	cout << "-> Calculating normalized vectors using Vec3*magic_invsqr(Vec3*Vec3).\n";
	tstart = clock();
	for(k=0; k<1000; k++)
		for(i=0; i<Nr; i++) calc_vecs[i]=vecs[i]*magic_invsqr(vecs[i]*vecs[i]);
	tend = clock();
	maxerror=0.0;
	for(i=0; i<Nr; i++){
		if(maxerror<fabs(calc_vecs[i].V3Norm()-1.0)) maxerror=fabs(calc_vecs[i].V3Norm()-1.0);
	}
	cout << "<- Done, took " << (tend-tstart)/CLOCKS_PER_SEC << " ms per " << Nr << " vectors (maximum deviation from unit length: " << maxerror*100 << "%).\n\n";
	
	cout << "-> Calculating normalized vectors using Vec3*custom_invsqrt(Vec3*Vec3).\n";
	tstart = clock();
	for(k=0; k<1000; k++)
		for(i=0; i<Nr; i++) calc_vecs[i]=vecs[i]*custom_invsqrt(vecs[i]*vecs[i]);
	tend = clock();
	maxerror=0.0;
	for(i=0; i<Nr; i++){
		if(maxerror<fabs(calc_vecs[i].V3Norm()-1.0)) maxerror=fabs(calc_vecs[i].V3Norm()-1.0);
	}
	cout << "<- Done, took " << (tend-tstart)/CLOCKS_PER_SEC << " ms per " << Nr << " vectors (maximum deviation from unit length: " << maxerror*100 << "%).\n\n";
	
	delete[] vecs;
	delete[] calc_vecs;
	cout << "*** Finished. ***\n\n";
	return 0;
}
