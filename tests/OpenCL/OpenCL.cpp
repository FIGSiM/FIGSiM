/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/* 
OpenCL.cpp
	Test OpenCL implementation
	written by Andreas Tillack September 2011
*/

#define Nr 1024*1024*10
#define low 0.2
#define high 13.29

#define EPS_CL EPS

#include <iostream>
#include "../../Config.h"
#include "../../MC_OpenCL.h"
#include "../../ScalarMat.h"
#include "../../VecMat.h"
#include <string.h>

/* Magic to include OpenCL source code
 * from external file into executable
 */
#ifdef OpenCL
#define CL_TO_STRING(CL) #CL
const char* OpenCL_cl = "#ifdef KHR_DP_EXTENSION\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n#else\n#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n#endif\n#define EPS_CL 1.2E-7\n#define pi M_PI\n#define GUID_ARG\n"
#include "OpenCL.cl"
#else
#define CL_TO_STRING(CL) CL
const char* OpenCL_cl = "";
#include "MiniCL/cl_MiniCL_Defs.h"
extern "C"
{
#include "OpenCL.cl"
}
MINICL_REGISTER(invsqrt_cl)
MINICL_REGISTER(volume_inside_count)
#endif

// For OpenCL we want the data types available here as well
#ifdef OpenCL
#undef CL_TO_STRING
#define CL_TO_STRING(CL) CL
extern "C"
{
#include "datatypes.h"
}
#endif

int main(){
	cout << "*** OpenCL.cpp ***\n";
	cout << "-> Creating random array of " << Nr << " doubles within range [" << low << "," << high << "].\n";
	unsigned int i,k;
	double* values = new double[Nr];
	__int32_t idum=-8;
	for(i=0; i<Nr; i++) values[i]=low+ran2(idum)*(high-low);
	cout << "<- Done.\n\n";
	
	cout << "-> Calculating inverse square root using reference 1.0/sqrt().\n";
	double* reference = new double[Nr];
	double tstart = clock();
	for(k=0; k<1000; k++){
		reference[k]=0.0; // make sure loop is executed w/ icc
		for(i=0; i<Nr; i++) reference[i]=1.0/sqrt(values[i]);
	}
	double tend = clock();
	cout << "<- Done, took " << (tend-tstart)/CLOCKS_PER_SEC << " ms per " << Nr << " values.\n\n";
	
	double* compare = new double[Nr];
	double maxerror=0.0;
	unsigned int error_max=0;
	cout << "-> Calculating inverse square root using OpenCL.\n";
	cl_context Context=NULL;
	cl_command_queue Queue=NULL;
	cl_device_id deviceID=GetDevice(CL_DEVICE_TYPE_ALL,Context,Queue,true);
	cl_program program = CreateProgram(deviceID,Context,&OpenCL_cl);
/*
 *                                                                 arg
 * Creating kernel: invsqrt_cl(__global double* input,              0
 *                             __global double* output,             1
 *                             const unsigned int count GUID_ARG)   2
 */
	cl_kernel kernel = CreateKernel(program,"invsqrt_cl");
	
	int* Error=NULL;
#ifndef OpenCL
	Error = new int;
#endif
	cl_mem input = clCreateBuffer(Context,CL_MEM_READ_ONLY,sizeof(double)*Nr,NULL,Error);
	cl_mem output = clCreateBuffer(Context,CL_MEM_WRITE_ONLY,sizeof(double)*Nr,NULL,Error);
	if(!input || !output){
		cout << "ERROR: Could not allocate device memory.\n";
		exit(1);
	}
	int Error_Code = clEnqueueWriteBuffer(Queue,input,CL_TRUE,0,sizeof(double)*Nr,values,0,NULL,NULL);
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
	cout << "-> Executing kernel ...";
	cout.flush();
	tstart = clock();
	for(k=0; k<1000; k++){
		Error_Code = clEnqueueNDRangeKernel(Queue,kernel,1,NULL,&global,&local,0,NULL,NULL);
		if(Error_Code!=CL_SUCCESS){
			cout << "ERROR: Failed to execute kernel: " << Error_Code << "\n";
			exit(1);
		}
		clFinish(Queue);
	}
	cout << " Done.\n";
	cout.flush();
	Error_Code = clEnqueueReadBuffer(Queue,output,CL_TRUE,0,sizeof(double)*Nr,compare,0,NULL,NULL);
	if(Error_Code!=CL_SUCCESS){
		cout << "ERROR: Could not read results.\n";
		exit(1);
	}
	tend = clock();
	
	clReleaseKernel(kernel);
	maxerror=0.0;
	error_max=0;
	for(i=0; i<Nr; i++){
		double err=0.0;
		if(reference[i]>EPS) err=fabs(compare[i]-reference[i])/reference[i];
		if(maxerror<err){
			maxerror=fabs(compare[i]-reference[i])/reference[i];
			error_max=i;
		}
	}
	
	cout << "<- Done, took " << (tend-tstart)/CLOCKS_PER_SEC << " ms";
#ifndef IN_WINDOWS
	cout << " (CPU time)";
#endif
	cout << " per " << Nr << " values (maximum relative error: " << maxerror*100 << "% at value " << values[error_max] << ").\n\n";
	clReleaseMemObject(input);
	clReleaseMemObject(output);
	
	delete[] values;
	delete[] reference;
	delete[] compare;
	
	cout << "-> Calculating excluded volume of benzene ring.\n";
	
	unsigned int element_count=12;
	double4* benzene_properties=new double4[element_count];
	benzene_properties[0]=create_double4(-2.3562,1.8209,-1.2822,3.3/2.0);
	benzene_properties[1]=create_double4(-1.1761,1.9604,-2.0142,3.3/2.0);
	benzene_properties[2]=create_double4(-0.5162,0.8292,-2.4966,3.3/2.0);
	benzene_properties[3]=create_double4(-0.7719,2.9478,-2.2081,2.5/2.0);
	benzene_properties[4]=create_double4(-1.0365,-0.4416,-2.2471,3.3/2.0);
	benzene_properties[5]=create_double4(0.4006,0.9376,-3.0654,2.5/2.0);
	benzene_properties[6]=create_double4(-2.2165,-0.5811,-1.5151,3.3/2.0);
	benzene_properties[7]=create_double4(-0.5238,-1.3205,-2.6219,2.5/2.0);
	benzene_properties[8]=create_double4(-2.8764,0.5502,-1.0327,3.3/2.0);
	benzene_properties[9]=create_double4(-2.8688,2.6999,-0.9074,2.5/2.0);
	benzene_properties[10]=create_double4(-3.7933,0.4417,-0.4640,2.5/2.0);
	benzene_properties[11]=create_double4(-2.6208,-1.5685,-1.3212,2.5/2.0);
	
	Frame_Attrib* frame=new Frame_Attrib;
	frame->box_corner=create_double4(-3.614091,-3.672024,-1.764658,0.0);
	frame->box=create_double4(7.228196,7.344105,3.414713,0.0);
	frame->center=create_double4(-1.696325,0.689667,-1.764658,0.0);
	Vec4 aa(0.284026,0.066874,0.956482,2.513526);
	frame->rot=AxisAngle2Quaternion(aa);
	
	Element_Dynamic* benzene_dyn=new Element_Dynamic[element_count];
	Element_Static* benzene_stat=new Element_Static[element_count];
	
	for(unsigned int i=0; i<element_count; i++){
		benzene_dyn[i].center=benzene_properties[i];
		benzene_dyn[i].center.w=0.0;
		benzene_dyn[i].rot=create_double4(0.0);
		benzene_stat[i].saxes=create_double4(benzene_properties[i].w);
	}
	
	tstart = clock();
/*
 *                                                                              arg
 * Creating kernel: volume_inside_count(Frame_Attrib frame,                      0
 *                                      __global double4* random_nrs,            1
 *                                      __global Element_Dynamic* element_dyn,   2
 *                                      __global Element_Static* element_stat,   3
 *                                      __global cl_short* output_count,            4
 *                                      const unsigned int element_count,        5
 *                                      const unsigned int count GUID_ARG)       6
 */
	kernel = CreateKernel(program,"volume_inside_count");
	
	local=0;
	Error_Code = clGetKernelWorkGroupInfo(kernel, deviceID, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local), &local, NULL);
	if(Error_Code!=CL_SUCCESS){
		cout << "ERROR: Could not determine workgroup information.\n";
		exit(1);
	}
	unsigned int NpI=(100000/local)*local;
	cout << "\t-> Using " << NpI << " iteration increment (local workgroup size: " << local << ")\n";
	// create random numbers
	cout << "\t-> Create first set of random numbers. ";
	cout.flush();
	idum=-8;
	double4* random_nrs=new double4[NpI];
	for(unsigned int i=0; i<NpI; i++){
		random_nrs[i].x=ran2(idum);
		random_nrs[i].y=ran2(idum);
		random_nrs[i].z=ran2(idum);
		random_nrs[i].w=0.0;
	}
	cout << "Done.\n";
	
	cl_mem frameattrib = clCreateBuffer(Context,CL_MEM_READ_ONLY,sizeof(Frame_Attrib),NULL,Error);
	cl_mem random = clCreateBuffer(Context,CL_MEM_READ_ONLY,sizeof(double4)*NpI,NULL,Error);
	cl_mem dynamic = clCreateBuffer(Context,CL_MEM_READ_ONLY,sizeof(Element_Dynamic)*element_count,NULL,Error);
	cl_mem stat = clCreateBuffer(Context,CL_MEM_READ_ONLY,sizeof(Element_Static)*element_count,NULL,Error);
	cl_mem countoutput = clCreateBuffer(Context,CL_MEM_WRITE_ONLY,sizeof(cl_short)*NpI,NULL,Error);
	if(!frameattrib || !random || !dynamic || !stat || !countoutput){
		cout << "ERROR: Could not allocate device memory.\n";
		exit(1);
	}
	Error_Code  = clEnqueueWriteBuffer(Queue,frameattrib,CL_TRUE,0,sizeof(Frame_Attrib),frame,0,NULL,NULL);
	Error_Code |= clEnqueueWriteBuffer(Queue,random,CL_TRUE,0,sizeof(double4)*NpI,random_nrs,0,NULL,NULL);
	Error_Code |= clEnqueueWriteBuffer(Queue,dynamic,CL_TRUE,0,sizeof(Element_Dynamic)*element_count,benzene_dyn,0,NULL,NULL);
	Error_Code |= clEnqueueWriteBuffer(Queue,stat,CL_TRUE,0,sizeof(Element_Static)*element_count,benzene_stat,0,NULL,NULL);
	if(Error_Code!=CL_SUCCESS){
		cout << "ERROR: Could not copy values to device.\n";
		exit(1);
	}
	Error_Code = clSetKernelArg(kernel, 0, sizeof(cl_mem), &frameattrib);
	Error_Code|= clSetKernelArg(kernel, 1, sizeof(cl_mem), &random);
	Error_Code|= clSetKernelArg(kernel, 2, sizeof(cl_mem), &dynamic);
	Error_Code|= clSetKernelArg(kernel, 3, sizeof(cl_mem), &stat);
	Error_Code|= clSetKernelArg(kernel, 4, sizeof(cl_mem), &countoutput);
	Error_Code|= clSetKernelArg(kernel, 5, sizeof(unsigned int), &element_count);
	global=NpI;
	count=NpI;
	Error_Code|= clSetKernelArg(kernel, 6, sizeof(unsigned int), &count);
	if(Error_Code!=CL_SUCCESS){
		cout << "ERROR: Failed to set kernel arguments.\n";
		exit(1);
	}
	
	// create output values on our site
	cl_short* outcounts=new cl_short[NpI];
	
	cout << "-> Executing kernel ...";
	cout.flush();
	double inside;
	unsigned int Nin=0;
	unsigned int Niter=0;
	do{
		if(Niter>0){
			Error_Code = clEnqueueWriteBuffer(Queue,random,CL_TRUE,0,sizeof(double4)*NpI,random_nrs,0,NULL,NULL);
			if(Error_Code!=CL_SUCCESS){
				cout << "ERROR: Could not copy values to device.\n";
				exit(1);
			}
		}
		Error_Code = clEnqueueNDRangeKernel(Queue,kernel,1,NULL,&global,&local,0,NULL,NULL);
		if(Error_Code!=CL_SUCCESS){
			cout << "ERROR: Failed to execute kernel: " << Error_Code << "\n";
			exit(1);
		}
		// Use time waiting for kernel to finish usefully
		for(unsigned int i=0; i<NpI; i++){
			random_nrs[i].x=ran2(idum);
			random_nrs[i].y=ran2(idum);
			random_nrs[i].z=ran2(idum);
			random_nrs[i].w=0.0;
		}
		clFinish(Queue);
		Error_Code = clEnqueueReadBuffer(Queue,countoutput,CL_TRUE,0,sizeof(cl_short)*NpI,outcounts,0,NULL,NULL);
		if(Error_Code!=CL_SUCCESS){
			cout << "ERROR: Could not read results.\n";
			exit(1);
		}
		for(unsigned int i=0; i<NpI; i++) Nin+=outcounts[i];
		Niter+=NpI;
		inside=(double)Nin/Niter;
	} while(((1.0-inside)/(Niter*inside)>1E-6) || (Nin<1000*element_count)); // first statement is testing for relative error of 1E-3 (0.1%)
	cout << " Done.\n";
	
	double boundingV=frame->box.x*frame->box.y*frame->box.z;
	double var=boundingV*boundingV*inside*(1-inside)/Niter;
	double origV=inside*boundingV;
	
	cout << "\t-> Excluded volume: " << origV << " +/- " << sqrt(var) << " Angström³\n";
	
	clReleaseKernel(kernel);
	tend = clock();
	cout << "<- Done, took " << (tend-tstart)/CLOCKS_PER_SEC << " s";
#ifndef IN_WINDOWS
	cout << " (CPU time)";
#endif
	cout << " for " << Niter << " iterations (including kernel build, memory transfers, etc.).\n\n";
	
	// MiniCL needs them released before CL memory objects (wierd, I know ...)
	delete frame;
	delete[] benzene_properties;
	delete[] benzene_dyn;
	delete[] benzene_stat;
	delete[] random_nrs;
	delete[] outcounts;
	
	clReleaseMemObject(frameattrib);
	clReleaseMemObject(random);
	clReleaseMemObject(dynamic);
	clReleaseMemObject(stat);
	clReleaseMemObject(countoutput);
	
	clReleaseProgram(program);
	clReleaseCommandQueue(Queue);
	clReleaseContext(Context);
	
	cout << "*** Finished. ***\n\n";
	return 0;
}

