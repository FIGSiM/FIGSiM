/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*
testcomplex.cpp
	Test complex vector implementation
	written by Andreas Tillack
*/

#include <iostream>
#include "../../VecMat.h"
using namespace std;

int main(){
	cout << "*** testcomplex.cpp ***\n";

	CVec3 cv(complex<double>(1,2),complex<double>(3,4),complex<double>(5,6));
	Vec3 t;
	t=cv.Re();
	cout << "\tRe(1+2i, 3+4i, 5+6i) = (" << t.V3Str() << ")\n";
	t=cv.Im();
	cout << "\tIm(1+2i, 3+4i, 5+6i) = (" << t.V3Str() << ")\n";
	t=cv.Abs();
	cout << "\tAbs(1+2i, 3+4i, 5+6i) = (" << t.V3Str() << ")\n";

	cout << "*** Finished. ***\n\n";
	return 0;
}
