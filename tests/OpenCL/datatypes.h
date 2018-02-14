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

typedef struct _Epsilon_Attrib{
	double4 saxes;
	double4 center;
	double4 rot;
	unsigned int vdwtype;
	unsigned int LJexp[2];
	double Solvent[3];
	double rmin;
	double rmax;
	double r;
} Epsilon_Attrib;

typedef struct _Element_Dynamic{
	double4 center;
	double4 rot; // use quaternion here
	double4 dipole;
	
	unsigned int type;
	int group_nr;
} Element_Dynamic;

typedef struct _Element_Static{
	double4 saxes;
	double vdw;
	
	unsigned int charges_start;
	unsigned int nr_charges;
} Element_Static;

typedef struct _Group_Type{
	bool group_dipole;
	double4 center;
} Group_Type;

)

