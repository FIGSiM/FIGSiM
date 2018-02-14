/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

CL_TO_STRING(
typedef struct _Frame_Attrib{
	double4 box_corner;
	double4 box;
	double4 center;
	double4 rot;
} Frame_Attrib;

typedef struct _Simulation_Attribs{
	double4 saxes;
	double4 invsaxes2;
	double4 center;
	double rot[9];
	double unit[9];
	double avg_width;
	double inv_avg_area;
	double sphericity;
	
	double weight_factor;
	unsigned int r_overscan;
	unsigned int vdwtype;
	unsigned int LJexp[2];
	double Solvent[3];
	
	unsigned int num_element_types;
	unsigned int touchtrunc;
	unsigned int latticetype;
	unsigned int PBCs[3];
	double boxlength[3];
	double in2;
	double escut2;
	double ljcut2;
	unsigned int cut;
	double rf_correction;
	
	double rp;
	double dphi_3;
	double dcostheta;
	double rmin;
	double rmax;
	double k0;
	double r_r;
	double beta;
	
	double r;
} Simulation_Attribs;

typedef struct _Element_Dynamic{
	double4 center;
	double rot[9];
	double4 dipole;
	
	unsigned int type;
	unsigned int charges_start;
	int group_nr;
} Element_Dynamic;

typedef struct _Element_Static{
	double4 saxes_vdw;
	double4 initial_dipole;
	unsigned int charges_start;
	unsigned int nr_charges;
} Element_Static;

typedef struct _Group_Type{
	unsigned int group_dipole;
	double4 center;
} Group_Type;

typedef struct _One_Move{
	unsigned int counter;
	unsigned int kk;
	double rand;
	double4 deltatrans; // use deltatrans.w to store displacement^2
	double deltarot[10]; // use deltarot[9] to store theta displacement^2
} One_Move;

)

