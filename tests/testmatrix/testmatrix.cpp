/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*
testmatrix.cpp
	Test matrix implementation
	written by Andreas Tillack Jan 22, 2011
*/

#include <iostream>
#include "../../ScalarMat.h"
#include "../../VecMat.h"

using namespace std;

int main(){
	cout << "*** testmatrix.cpp ***\n";
	
	Mat33 A;
	Vec3 v(0,0,0);
	A.mat[0][0]=9;
	A.mat[0][1]=85;
	A.mat[0][2]=5;
	
	A.mat[1][0]=-36;
	A.mat[1][1]=-215;
	A.mat[1][2]=5;
	
	A.mat[2][0]=-27;
	A.mat[2][1]=-210;
	A.mat[2][2]=-6;
	
	cout << "A =\n" << A.M3Str() << "\n\n";
	A.M3GEPP(v);
	cout << "gauss(A) =\n" << A.M3Str() << "\n\nv =\n" << v.V3Str() << "\n\n";
	
	A.mat[0][0]=9;
	A.mat[0][1]=85;
	A.mat[0][2]=0;
	
	A.mat[1][0]=-36;
	A.mat[1][1]=-215;
	A.mat[1][2]=0;
	
	A.mat[2][0]=-27;
	A.mat[2][1]=-210;
	A.mat[2][2]=0;
	
	v=Vec3(5,5,-6);
	cout << "A =\n" << A.M3Str() << "\n\n";
	cout << "v = " << v.V3Str() << "\n\n";
	cout << "Tr(A) = " << A.M3Trace() << "\n";
	cout << "Det(A) = " << A.M3Det() << "\n";
	cout << "A*v=x ; x = " << (A*v).V3Str() << "\n";
	bool multiples[3];
	cout << "A*x=v ; x = " << A.M3LinSolve(v,multiples[0]).V3Str();
	if (multiples[0]!=0){ cout << " (and multiples)"; }
	cout << "\n\n";
	
	CVec3 ev;
/*	Mat33 B; // determine runtime in nanoseconds
	B=A;
	for(int i=0; i<1E9; i++){
		B.mat[(i%9)%3][i%3]+=double(i%11)-5.0;
		ev=B.Eigenvalues();
	}*/
	ev=A.Eigenvalues();
	cout << "eigenvalues(A) = " << ev.CV3Str() << "\n";
/*	cout << "Are eigenvalues correct? ";
	if (((ev.cvec[0]==-2.0) && (ev.cvec[1]==ev.cvec[2])) && (ev.cvec[1]==2.0)){
	    cout << "Yes :-)\n";
	} else{
	    cout << "No :-(\n" ;
//	    return 1;
	}*/
	cout << "Now take eigenvalues and determine eigenvectors\n";
	Vec3 eValues(ev.Re());
	Mat33 eVectors=A.Eigenvectors(eValues,multiples,true);
	cout << "eigenvector to " << real(ev.cvec[0]) << ": " << eVectors.ColumnVec3(0).V3Str();
	if (multiples[0]!=0){ cout << " (and multiples)"; }
	cout << "\n";
	cout << "eigenvector to " << real(ev.cvec[1]) << ": " << eVectors.ColumnVec3(1).V3Str();
	if (multiples[1]!=0){ cout << " (and multiples)"; }
	cout << "\n";
	cout << "eigenvector to " << real(ev.cvec[2]) << ": " << eVectors.ColumnVec3(2).V3Str();
	if (multiples[2]!=0){ cout << " (and multiples)"; }
	cout << "\n\n";
	cout << "Let's try some LU decomposition ...\n A = \n";
	cout << A.M3Str() << "\n\n";
	A.LUDecomposition();
	cout << "L*R =\n" << A.M3Str() << "\n\n";
	cout << "Testing Vec2Rot implementation:\n-> Vec2Rot of random vector and axis, then comparison with rotation matrix times axis vector ...\n";
	unsigned int axis;
	Vec3 av, rotate, rotated;
	Mat33 rot;
	__int32_t idum=-8;
	bool OK=true;
	double tstart = clock();
	for(unsigned int i=0; i<10000000; i++){
		axis=(unsigned int)(3.0*ran2(idum));
		av.vec[0]=0.0; av.vec[1]=0.0; av.vec[2]=0.0;
		av.vec[axis]=1.0;
		rotate.vec[0]=100.0*ran2(idum)-50.0; rotate.vec[1]=100.0*ran2(idum)-50.0; rotate.vec[2]=100.0*ran2(idum)-50.0;
		rot=Vec2Rot(rotate,axis);
		rotated=rot*av;
		if(((fabs(rotated.vec[0]-rotate.vec[0])>EPS) || (fabs(rotated.vec[1]-rotate.vec[1])>EPS)) || (fabs(rotated.vec[2]-rotate.vec[2])>EPS)){
			cout << "ERROR: (" << rotate.V3Str(',') << ")!=(" << rotated.V3Str(',') << ")\n";
			cout << "rotation matrix around axis " << axis << ":\n" << rot.M3Str() << "\n";
			OK=false;
		}
	}
	double tend = clock();
	if(OK) cout << "-> Succesfully finished test.\n"; else cout << "-> Failed test (errors are listed above).\n";
	cout << "<- Done, took " << (tend-tstart)/CLOCKS_PER_SEC*100 << " ns per vector.\n\n";
	cout << "Testing AxisAngle implementation:\n-> Create rotation matrix from random axis vector and angle,\n-> then compare with rotation matrix from recreated Axis and angle ...\n";
	OK=true;
	unsigned int errorcount=0;
	Vec4 axisangle, axisangle2;
	Mat33 rot2;
	tstart = clock();
	for(unsigned int i=0; i<10000000; i++){
		axisangle.vec[0]=100.0*ran2(idum)-50.0; axisangle.vec[1]=100.0*ran2(idum)-50.0; axisangle.vec[2]=100.0*ran2(idum)-50.0;
		axisangle.vec[3]=pi*(2.0*ran2(idum)-1.0); // rotate around axis +/- pi
		rot=AxisAngle2Rot(axisangle);
		axisangle2=Rot2AxisAngle(rot);
		rot2=AxisAngle2Rot(axisangle2);
		bool error=false;
		for(unsigned int n=0; n<3; n++){
			for(unsigned int m=0; m<3; m++){
				if(fabs(rot.mat[m][n]-rot2.mat[m][n])>EPS){
					error=true;
					break;
				}
			}
		}
		if(error){
			cout << "ERROR: rotation matrix from random axis and angle: (" << axisangle.V4Str(',') << ")\nrotation matrix:\n" << rot.M3Str() << "\ndoes not yield equal rotation matrix from recreated axis and angle: (" << axisangle2.V4Str(',') << ")\nrecreated rotation matrix:\n" << rot2.M3Str() << "\n";
			errorcount++;
			OK=false;
		}
	}
	tend = clock();
	if(OK) cout << "-> Succesfully finished test.\n"; else cout << "-> Failed test (" << errorcount << " errors are listed above).\n";
	cout << "<- Done, took " << (tend-tstart)/CLOCKS_PER_SEC*100 << " ns per vector.\n\n";
	cout << "Testing RotAtoB implementation:\n-> Create rotation matrix from two vectors mapping one onto the other ...\n";
	OK=true;
	errorcount=0;
	tstart = clock();
	Vec3 a,b,a2,b2;
	for(unsigned int i=0; i<10000000; i++){
		a.vec[0]=100.0*ran2(idum)-50.0; a.vec[1]=100.0*ran2(idum)-50.0; a.vec[2]=100.0*ran2(idum)-50.0;
		b.vec[0]=100.0*ran2(idum)-50.0; b.vec[1]=100.0*ran2(idum)-50.0; b.vec[2]=100.0*ran2(idum)-50.0;
		rot=RotAtoB(a,b);
		a/=a.V3Norm();
		b/=b.V3Norm();
		a2=rot.M3Transpose()*b;
		b2=rot*a;
		bool error=false;
		for(unsigned int n=0; n<3; n++){
			if((fabs(a2.vec[n]-a.vec[n])>EPS) || (fabs(b2.vec[n]-b.vec[n])>EPS)){
				error=true;
				break;
			}
		}
		if(error){
			cout << "ERROR: mismatch: ";
			if((a2-a).V3Norm()>EPS){
				cout << "a: (" << a.V3Str(',') << ") =/= (" << a2.V3Str(',') << ")";
				if((b2-b).V3Norm()>EPS) cout << ", ";
			}
			if((b2-b).V3Norm()>EPS) cout << "b: (" << b.V3Str(',') << ") =/= (" << b2.V3Str(',') << ")";
			cout << "\n";
			errorcount++;
			OK=false;
		}
	}
	tend = clock();
	if(OK) cout << "-> Succesfully finished test.\n"; else cout << "-> Failed test (" << errorcount << " errors are listed above).\n";
	cout << "<- Done, took " << (tend-tstart)/CLOCKS_PER_SEC*100 << " ns per vector.\n\n";
	cout << "*** Finished. ***\n\n";
	return 0;
}
