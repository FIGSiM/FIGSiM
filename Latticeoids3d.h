/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

#ifndef INCLUDED_LATTICEOIDS3D
#define INCLUDED_LATTICEOIDS3D

#ifndef INCLUDED_CONFIG
#include "Config.h"
#endif
#ifndef INCLUDED_MC_ELEMENTS
#include "MC_Elements.h"
#endif

class LatticeOids3d
{
public:
	LatticeOids3d(MC_Elements &BeOids); ///< Lattice constructor
	~LatticeOids3d(); ///< destructor
private:
	unsigned int nr_groups; ///< the number of created groups
	Element_Group** created_groups; ///< groups created through constructor
	
	void SimpleCubic(MC_Elements &BeOids); ///< Make a simple cubic lattice
	void GroupPacking(MC_Elements &BeOids); ///< pack groups in lattice
	void SeqLoadGroups(MC_Elements &BeOids); ///< Sequentially load element groups (a.k.a. molecules) into the lattice
	void RandLoadGroups(MC_Elements &BeOids); ///< Randomly load element groups (a.k.a. molecules) into the lattice
	void SeqLoad(MC_Elements &BeOids); ///< Sequentially load elements into the lattice
	void RandLoad(MC_Elements &BeOids); ///< Randomly load elements into the lattice
};

#endif
