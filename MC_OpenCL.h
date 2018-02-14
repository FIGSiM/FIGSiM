/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * MC_OpenCL.h : include file for using OpenCL (enables use of multicore CPUs and GPUs)
 *
 * Initial version by Andreas Tillack on Aug 13, 2011
 *
 */

#ifndef INCLUDED_MC_OPENCL
#define INCLUDED_MC_OPENCL

#include "Config.h"

#include <iostream>
#include <cstdlib>

#ifdef OpenCL
	#ifdef __APPLE__
		#include "OpenCL/opencl.h"
	#else
		#include "CL/opencl.h"
	#endif
#else
	#include "MiniCL/cl.h"
#endif

#define EPS_CL EPS

inline cl_device_id GetDevice(cl_device_type DeviceType, cl_context &Context, cl_command_queue &Queue, bool double_precision)
{
	std::cout << "-> Initializing OpenCL device (";
	if(double_precision) std::cout << "double "; else std::cout << "single ";
	std::cout << "precision)\n";
	cl_device_id deviceID=NULL;
	
	int Error_Code;
	bool found_device=false;
#ifdef OpenCL
	unsigned int num_platforms=0;
	Error_Code = clGetPlatformIDs(0, NULL, &num_platforms);
	if(Error_Code!=CL_SUCCESS){
		std::cout << "ERROR: Could not determine OpenCL capable platforms.\n";
		exit(1);
	}
	cl_platform_id* platforms = new cl_platform_id[num_platforms];
	Error_Code = clGetPlatformIDs(num_platforms, platforms, NULL);
	if(Error_Code!=CL_SUCCESS){
		std::cout << "ERROR: No OpenCL capable device found.\n";
		exit(1);
	}
	for(unsigned int i=0; i<num_platforms; i++){
		unsigned int num_devices=0;
		Error_Code = clGetDeviceIDs(platforms[i], DeviceType, 0, NULL, &num_devices);
		if(Error_Code!=CL_SUCCESS){
			std::cout << "ERROR: No OpenCL capable device found.\n";
			exit(1);
		}
		cl_device_id* IDs = new cl_device_id[num_devices];
		Error_Code = clGetDeviceIDs(platforms[i], DeviceType, num_devices, IDs, NULL);
		if(Error_Code!=CL_SUCCESS){
			std::cout << "ERROR: Could not access OpenCL capable device.\n";
			exit(1);
		}
		// test which device is capable of double precision
		cl_device_fp_config fp_config=0;
		cl_device_type dev_type=0;
		cl_bool dev_avail=false;
		cl_bool compiler_avail=false;
		for(unsigned int j=0; j<num_devices; j++){
			if(double_precision){
				Error_Code = clGetDeviceInfo(IDs[j],CL_DEVICE_DOUBLE_FP_CONFIG,sizeof(cl_device_fp_config),&fp_config,NULL);
			} else Error_Code = clGetDeviceInfo(IDs[j],CL_DEVICE_SINGLE_FP_CONFIG,sizeof(cl_device_fp_config),&fp_config,NULL);
			Error_Code |= clGetDeviceInfo(IDs[j],CL_DEVICE_TYPE,sizeof(cl_device_type),&dev_type,NULL);
			Error_Code |= clGetDeviceInfo(IDs[j],CL_DEVICE_AVAILABLE,sizeof(cl_bool),&dev_avail,NULL);
			Error_Code |= clGetDeviceInfo(IDs[j],CL_DEVICE_COMPILER_AVAILABLE,sizeof(cl_bool),&compiler_avail,NULL);
			if(Error_Code!=CL_SUCCESS){
				std::cout << "ERROR: Could not get device information.\n";
				exit(2);
			}
			if(fp_config && dev_avail && compiler_avail && (dev_type&DeviceType)){
				deviceID=IDs[j];
				found_device=true;
				break; // use the first available target device with compiler found
			}
		}
		delete[] IDs;
		if(found_device) break;
	}
	delete[] platforms;
#else
	deviceID=0;
	found_device=true;
#endif
	if(!found_device){
		std::cout << "ERROR: No device capable of ";
		if(double_precision) std::cout << "double "; else std::cout << "single ";
		std::cout << "precision floating point operations found.\n";
		exit(3);
	}
	char devicename[256];
	Error_Code = clGetDeviceInfo(deviceID,CL_DEVICE_NAME,256,devicename,NULL);
	if(Error_Code!=CL_SUCCESS){
		std::cout << "ERROR: Could not get device information.\n";
		exit(2);
	}
	std::cout << "\t-> Using device <" << devicename << ">.\n";
	
#ifdef OpenCL
	Context = clCreateContext(0, 1, &deviceID, NULL, NULL, &Error_Code);
#else
	Context = clCreateContextFromType(NULL, DeviceType, NULL, NULL, &Error_Code);
#endif
	if(!Context){
		std::cout << "ERROR: Unable to create compute context on device.\n";
		exit(3);
	}
	Queue = clCreateCommandQueue(Context, deviceID, 0, &Error_Code);
	if(!Queue){
		std::cout << "ERROR: Could not create command queue for device.\n";
		exit(4);
	}
	return deviceID;
}

inline cl_program CreateProgram(cl_device_id &deviceID, cl_context &Context, const char** Code)
{
	int Error_Code=0;
#if DEBUG_LEVEL>1
	std::cout << "-> Compiling program from source code ... ";
	std::cout.flush();
#endif
	cl_program program = clCreateProgramWithSource(Context,1,Code,NULL,&Error_Code);
	if(!program || (Error_Code!=CL_SUCCESS)){
		std::cout << "\nERROR: Could not create compute program.\n";
		exit(1);
	}
	Error_Code = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
#ifdef OpenCL
	char* build_output;
	size_t error_size;
	clGetProgramBuildInfo(program, deviceID, CL_PROGRAM_BUILD_LOG, 0, NULL, &error_size);
	build_output=new char[error_size+1];
	clGetProgramBuildInfo(program, deviceID, CL_PROGRAM_BUILD_LOG, error_size, build_output, NULL);
#if DEBUG_LEVEL>1
	if(error_size>1){
		std::cout << "\n" << build_output;
		if(build_output[error_size]!='\n') std::cout << "\n";
	}
#endif // ifdef OpenCL
#endif
	if(Error_Code!=CL_SUCCESS){
		std::cout << "ERROR: Could not build program (error code: " << Error_Code << ").";
#if DEBUG_LEVEL<=1 && defined(OpenCL)
		std::cout << " Build log:\n" << build_output;
#endif
		std::cout << "\n";
		exit(2);
	}
#ifdef OpenCL
	delete[] build_output;
#endif
#if DEBUG_LEVEL>1
#ifdef OpenCL
	if(error_size>1) std::cout << "<- ";
	std::cout << "Done (";
	if(error_size>1) std::cout << "with warnings"; else std::cout << "no warnings";
	std::cout << ").\n";
#else
	deviceID=0;
	std::cout << "Done.\n";
#endif // ifdef OpenCL
#endif
	return program;
}

inline cl_kernel CreateKernel(cl_program &program, const char* kernel_name)
{
#if DEBUG_LEVEL>1
	std::cout << "\t-> Creating kernel \"" << kernel_name << "\" ... ";
	std::cout.flush();
#endif
	int Error_Code=0;
	cl_kernel kernel = clCreateKernel(program,kernel_name, &Error_Code);
	if(!kernel || (Error_Code!=CL_SUCCESS)){
		std::cout << "\nERROR: Could not create compute kernel <" << kernel_name << "> (error code: " << Error_Code << ").\n";
		exit(1);
	}
#if DEBUG_LEVEL>1
	std::cout << "Done.\n";
	std::cout.flush();
#endif
	return kernel;
}

#endif

