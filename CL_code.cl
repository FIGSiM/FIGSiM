/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
Selectively and automatically include OpenCL code into program for runtime compilation.\n
Code inclusion is governed through compiler directives set by including code part.\n\n
Currently, the following code is included:\n\n
fit2lod: \includelineno CL_epsilon_matching.h
MC_elements.cpp: \includelineno CL_evolve_system.h
*/

#ifdef FIT2LOD
#include "CL_epsilon_matching.h"
#endif
#ifdef MC_ELEMENTS
#include "CL_evolve_system.h"
#endif

/// magic to include code into executables, code is in included CL_*.h files
CL_TO_STRING(
);

