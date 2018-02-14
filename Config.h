/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * AT June 21, 2011
 */
/*!\mainpage Robinson group Monte Carlo simulation code
 * \section OpenCL_sec OpenCL code:
 * OpenCL code is included through CL_code.cl.
 * \section changelog_sec Changelog:
 * \verbinclude CHANGELOG
 */

#ifndef INCLUDED_CONFIG
#define INCLUDED_CONFIG

#define BUILD "AT" // Date of current build. Do not update if only simulation parameters are changed.
#define VERSION "Beta 16.4"

// Special care needs to be taken when compiling under Windows
#if defined(_WIN32) || defined(WIN32) || defined(__CYGWIN__) || defined(__MINGW32__) || defined(__BORLANDC__)
#define IN_WINDOWS
#include <stdint.h>
typedef uint32_t __uint32_t;
typedef uint64_t __uint64_t;
typedef int32_t __int32_t;
typedef int64_t __int64_t;
#endif

// Load code in executable
#ifdef __APPLE__
    #include <mach-o/getsect.h>
    #define EXTERN_OBJLOAD(OBJNAME) extern const unsigned char _ ## OBJNAME __asm("section$start$__DATA$__" #OBJNAME);
    #define EXTERN_OBJDATA(OBJNAME) &_ ## OBJNAME
    #define EXTERN_OBJLENGTH(OBJNAME) (getsectbyname("__DATA","__" #OBJNAME)->size)
#elif (defined IN_WINDOWS)
    #define EXTERN_OBJLOAD(OBJNAME) extern const unsigned char binary_ ## OBJNAME ## _start[];   extern const unsigned char binary_ ## OBJNAME ## _end[];
    #define EXTERN_OBJDATA(OBJNAME) binary_ ## OBJNAME ## _start
    #define EXTERN_OBJLENGTH(OBJNAME) ((binary_ ## OBJNAME ## _end) - (binary_ ## OBJNAME ## _start))
#else
    #define EXTERN_OBJLOAD(OBJNAME) extern const unsigned char _binary_ ## OBJNAME ## _start[]; extern const unsigned char _binary_ ## OBJNAME ## _end[];
    #define EXTERN_OBJDATA(OBJNAME) _binary_ ## OBJNAME ## _start
    #define EXTERN_OBJLENGTH(OBJNAME) ((_binary_ ## OBJNAME ## _end) - (_binary_ ## OBJNAME ## _start))
#endif

#include "tests/use_opencl.h" // either defines OpenCL or not (file is empty if Makefile is not used and/or OpenCL is not installed)

// only uncomment if you know what your doing
// #undef OpenCL // do not compile for OpenCL
// #define OpenCL // do compile for OpenCL (only set manually if <Open>CL/cl.h exists and linking to OpenCL is worked out)

#ifndef OpenCL
#define MINI_CL
#endif

// when USE_OPENCL is defined OpenCL is used for Monte-Carlo calculations (currently unstable, should be left commented out ...)
// #define USE_OPENCL

// when FRIEDMAN_IMAGE is defined the Friedman image charge method (essentially q/R_C) is used for reaction field calculation
// #define FRIEDMAN_IMAGE

/* Level of debug output (cumulative)
 *  0: minimal output and warnings/exit messages
 *  1: informative messages (includes output for default values if certain parameters have not been set)
 *  2: messages for the advanced user
 *  3: algorithm specific output (for example control output, etc.)
 *  4: highest level: function entry/exit notifications and inner loop algorithm output
 */
#define DEBUG_LEVEL 2

#define g0 1.0
// #define SAMPLE_LOD_ULJ
// if defined use weighted average for LOD ellipsoid center calculation
#define LOD_CENTER_WEIGHTED
// if defined calculate epsilon for different g-factors during epsilon matching
// #define EPSILON_G
// if defined add derivatives for epsilon matching
// #define EPSILON_DVLJ
// if defined use lab frame for rotation of elements (note that OpenCL is in element frame by default, for lab frame one needs to change evolve_system.cl manually)
// #define LAB_FRAME_ROTATION
// if defined run a test loop of 0 <= rT <= 10,000 and output data of first ellipsoid
// #define RT_TEST_LOOP

#define theta_oversampling 4 // theta oversampling of theta for x3d texture output (please leave as is)
#define r_oversampling 2 // r oversampling used in calculations

#define theta_res 256
#define phi_res 256
#define r_res 128
#define voxelcount theta_res*phi_res*r_res
#define pi_quarter_fraction (2.0-sqrt(2.0))/4.0 // cos(pi/4) = sqrt(2)/2 => 1-2*u = sqrt(2)/2 => 1-2*i/N = sqrt(2)/2 => i = N*(2-sqrt(2))/4

#define theta_res_half theta_res/2.0
#define phi_res_per_tau phi_res/(2.0*pi)

// comment this line out if code shall use previously used ran2 random number generator
#define USE_CMWC4096

#define R_COMPRESSION 3

// define OUTPUT_EPSILON_OF_RT

#define NEIGHBOR_UPDATE_FREQ 2

// group potentials defined
#define stretch_potential 1
#define bend_potential 2
#define dihedral_potential 3
#define improper_dihedral_potential 4
#define spring_potential 5
#define chainspring_potential 6

// config/trajectory file reading related stuff
#define READBLOCKSIZE 1024
#define MAXCONFSIZE 1024*1024
#define RECURSION_LIMIT 100

#define NUM_V_STORE 4 // number of potential energies stored during each step (currently VmuE, VES, VLJ, and VG)

// define step to exit for debugging purposes (if <0 will not prematurely exit)
#define DEBUG_EXIT_STEP -1

// define size limit for trajectory files (new file will be created after approx. 512 MB)
#define TRAJECTORY_SIZE 512*1024*1024

// use 4 bits for special distances (<bond_distance_factors> for each group)
#define SPECIAL_DISTANCES_BIT_EXPONENT 2 // 2^2 = 4
#define SPECIAL_DISTANCES_BIT_SCREEN ((1<<(1<<SPECIAL_DISTANCES_BIT_EXPONENT))-1)
#define SPECIAL_DISTANCES (SPECIAL_DISTANCES_BIT_SCREEN-1) // first two are used for less than bond_range (0) and larger, but not special (1)

#define SUPERGROUP_ROTATE_STEPS 2000

//Define constants here
#define pi 3.14159265
#define pi_half pi/2.0
#define EPS 1.2e-7
#define NA 6.02214179E23
#define two_to_onesixth 1.12246204831
#define e_in_esu 4.80320427 // 1 e in 10^(-10) stat coulomb (= 10^-(10) esu)
// 1 atm = 101.325 kPa = 1.01325x10^5 J/m^3 ; 1 pErg/Angström^3 = 1E-19 J / (1E-30 m^3) = 1E11 J/m^3 = 1E11 Pa => 1 Pa = 1E-11 pErg/Angström^3 => 1 atm = 1.01325E-6 pErg/Angström^3
#define atm_in_pErg_per_Ang3 1.01325E-6
#define perg_to_kJ_per_mol NA*1E-19/1000
#define perg_to_kcal_per_mol perg_to_kJ_per_mol/4.184
#define kB 1.3806488e-4
// 1 Debye = 10^-21 C*m^2/s / c 
// => perg/Debye = 10^-19 J/(10^-21 C*m^2/s)*2.99792458*10^8*m/s = 2.99792458x10^10 N/C
// 1 MV/m = 10^6 V/m = 10^6 N/C (CV = J => CV/m = N => N/C = V/m)
// => 1 MV/m = 10^6/(2.99792458*10^10) perg/Debye = 3.335646x10^-5 perg/D
#define MV_per_m_to_perg_per_Debye 3.335646E-5
// hbar = 1.05457E-34 Js
#define hbar_perg_ps 1.05457E-34/(1E-19*1E-12)
// 1 amu = 1.66053886E-27 J/(m/s)^2 -> Js^2/m^2
#define amu_to_perg_ps2_per_Ang2 1.66053886E-27/(1E-19*1E-24)*(1E-20)

#endif

