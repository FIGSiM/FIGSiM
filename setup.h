/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * setup.h
 * Created by ljohnson on 5/31/09.
 */

#ifndef INCLUDED_SETUP
#define INCLUDED_SETUP

#ifndef INCLUDED_CONFIG
#include "Config.h"
#endif
#ifndef INCLUDED_MC_ELEMENTS
#include "MC_Elements.h"
#endif

void phys_configuration(Config_Data* configuration); ///< Load physical constants into configuration
void setup_vdw(Config_Data* configuration); ///< Precalculate reusable parameters for LJ interactions and load into configuration
void setup_electrostatics(Config_Data* configuration); ///< Configure electrostatics calculations, including reaction field
void check_box(Config_Data* configuration); ///< Check periodic boundary conditions and calculate box volume
void final_validation(Config_Data* configuration); ///< Check essential simulation parameters to make sure it is runnable (e.g. not zero particles, etc)
void check_oids(MC_Elements &BeOids);
void open_summary (Config_Data* configuration); ///< Open batch summary file and write header data

#endif

