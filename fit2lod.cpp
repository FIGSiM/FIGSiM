/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * Robinson group Monte Carlo simulation code setup tool for LOD model
 * Command line utility with one argument, the configuration file
 * Usage example: fit2lod testconfigs/benzene-LOD.conf
 * Created by Andreas Tillack, Jan 2012
 */

#include "MC_OpenCL.h"
#include "MC_Config.h"
#include "ScalarMat.h"
#include "VecMat.h"
#include "setup.h"
#include "x3d_functions.h"

#define EPS_ELEMENT_COUNT 4
#define FIT2LOD // directive to tell compiler which cl code to include (leads to shorter compile time at runtime)

/* Magic to include OpenCL source code
 * from external file into executable
 */
#ifdef OpenCL
#define CL_TO_STRING(CL) #CL
const char* OpenCL_code = "#ifdef KHR_DP_EXTENSION\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n#else\n#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n#endif\n#define EPS_CL 1.2E-7\n#define pi M_PI\n#define GUID_ARG\n"
#include "CL_code.cl"
#else
#define CL_TO_STRING(CL) CL
const char* OpenCL_code = "";
#include "MiniCL/cl_MiniCL_Defs.h"
extern "C"
{
#include "CL_code.cl"
}
MINICL_REGISTER(match_epsilon)
#endif

// For OpenCL we want the data types available here as well
#ifdef OpenCL
#undef CL_TO_STRING
#define CL_TO_STRING(CL) CL
extern "C"
{
#include "CL_datatypes.h"
}
#endif

double EPSofRT(double xval, double* coefficients, unsigned int nr_coefficients)
{
	double fx=0.0;
	double xacc=1.0;
	for(unsigned int j=1; j<nr_coefficients; j++){
		fx+=coefficients[j]*xacc;
		xacc*=xval;
	}
	fx/=qqpwr(coefficients[0]+xval,nr_coefficients-2);
	return fx;
}

double IAadjust(double xval, double* coefficients, unsigned int nr_coefficients)
{
	if(xval>coefficients[0]){
		// a_h+b_h*x+c_h*x^2
		return coefficients[1]+xval*(coefficients[2]+coefficients[3]*xval);
	} else{
		// a_l+b_l*x+c_l*x^2
		// a_l=a_h-x_0^2*(c_h-c_l)
		double x0cd=coefficients[0]*(coefficients[3]-coefficients[4]);
		coefficients[5]=coefficients[1]-coefficients[0]*x0cd;
		// b_l=2.0*x_0*(c_h-c_l)+b_h
		coefficients[6]=2.0*x0cd+coefficients[2];
		return coefficients[5]+xval*(coefficients[6]+coefficients[4]*xval);
	}
}

inline double Chi2(Vec3* datapoints, unsigned int nr_points, double (*fxn)(double,double*,unsigned int),double* coefficients, unsigned int nr_coefficients)
{
	double oxval;
	double xval=datapoints[0].vec[0];
	double fx=fxn(xval,coefficients,nr_coefficients);
	double X2=qqpwr((fx-datapoints[0].vec[1]),2)*(datapoints[0].vec[2]);
	for(unsigned int i=1; i<nr_points; i++){
		// calculate function value
		oxval=xval;
		xval=datapoints[i].vec[0];
		if(fabs(xval-oxval)>EPS){
			fx=fxn(xval,coefficients,nr_coefficients);
		}
		X2+=qqpwr((fx-datapoints[i].vec[1]),2)*(datapoints[i].vec[2]);
	}
	return X2;
}

inline void SortData(Vec3* datapoints, unsigned int nr_points)
{
	Vec3 temp;
	for(unsigned int i=1; i<nr_points-1; i++){
		for(unsigned int j=0; j<nr_points-i; j++){ // last point is already the largest (no need to check again)
			if(datapoints[j].vec[0]>datapoints[j+1].vec[0]){
				temp=datapoints[j+1];
				datapoints[j+1]=datapoints[j];
				datapoints[j]=temp;
			}
		}
	}
}

// datapoints vector is expected to be in the following format:
// vec[0] ... x value
// vec[1] ... y value
// vec[2] ... reciproke variance of y value (1/sigma_i^2)
inline void Coefficient_Fit(Vec3* datapoints, unsigned int nr_points, double (*fxn)(double,double*,unsigned int), double* coefficients, double* coeff_sigmas, unsigned int nr_coefficients)
{
	double X2=Chi2(datapoints,nr_points,fxn,coefficients,nr_coefficients);
	double ChiStart=X2;
	double* coeff_new=new double[nr_coefficients];
	double* best_coeff=new double[nr_coefficients];
	for(unsigned int i=0; i<nr_coefficients; i++) best_coeff[i]=coefficients[i];
	double best_X2=X2;
	double lambda=1.0;
#ifdef USE_CMWC4096
	init_CMWC4096(8);
	zigset();
#else
	__int32_t idum = -8;
#endif
	unsigned int Nit=0;
	unsigned int best_N=0;
	bool iterate=true;
	double tstart = clock();
	unsigned int which;
	do{
#ifdef USE_CMWC4096
		which=(unsigned int)((nr_coefficients-1)*ranQ()+0.5);
		for(unsigned int i=0; i<nr_coefficients; i++) coeff_new[i]=coefficients[i];
		coeff_new[which]=coefficients[which]+(2.0*ranQ()-1.0)*coeff_sigmas[which];
//		for(unsigned int i=0; i<nr_coefficients; i++) coeff_new[i]=coefficients[i]+(2.0*ranQ()-1.0)*coeff_sigmas[i];
#else
		which=(unsigned int)((nr_coefficients-1)*ran2(idum)+0.5);
		for(unsigned int i=0; i<nr_coefficients; i++) coeff_new[i]=coefficients[i];
		coeff_new[which]=coefficients[which]+(2.0*ranQ()-1.0)*coeff_sigmas[which];
//		for(unsigned int i=0; i<nr_coefficients; i++) coeff_new[i]=coefficients[i]+(2.0*ran2(idum)-1.0)*coeff_sigmas[i];
#endif
		double X2new=Chi2(datapoints,nr_points,fxn,coeff_new,nr_coefficients);
#ifdef USE_CMWC4096
		if(exp(X2-X2new)>(1.0-lambda)+lambda*ranQ()){ // likelyhood ratio criterion which goes over to X2new<X2 over the iterations ...
#else
		if(exp(X2-X2new)>(1.0-lambda)+lambda*ran2(idum)){ // likelyhood ratio criterion which goes over to X2new<X2 over the iterations ...
#endif
			X2=X2new;
			for(unsigned int i=0; i<nr_coefficients; i++){
				coefficients[i]=coeff_new[i];
				if(X2<best_X2) best_coeff[i]=coefficients[i];
			}
			if(X2<best_X2){
				iterate=(best_N<1E6);
				best_X2=X2;
				best_N=0;
				lambda=0.001+best_X2/ChiStart;
			}
		}
		Nit++;
		best_N++;
		if(Nit==best_N) iterate&=(best_N<1E6); else iterate&=(best_N<1E7);
	} while(iterate && (Nit<1E8));
	double tend = clock();
	delete[] coeff_new;
	cout << "Fit results: ";
	for(unsigned int i=0; i<nr_coefficients; i++){
		coefficients[i]=best_coeff[i];
		if(i>0) cout << ", ";
		cout << coefficients[i];
	}
	cout << " (Reduced Chi2 = " << best_X2/(nr_points-nr_coefficients-1) << ", " << Nit << " iterations, took " << (tend-tstart)/CLOCKS_PER_SEC << " s)\n";
	delete[] best_coeff;
}

inline unsigned int space_number(string &s)
{
	unsigned int nr=0;
	for(unsigned int i=0; i<s.length(); i++) if(s[i]==' ') nr++;
	return nr+1;
}

inline void DrawCalcSphere(ofstream &vizout, double radius, Vec3 position)
{
	vizout << "\t\t<Transform translation='"<< position.V3Str(' ') << "'>\n";
	vizout << "\t\t\t<Shape>\n";
	vizout << "\t\t\t\t<Appearance>\n";
	vizout << "\t\t\t\t\t<Material diffuseColor='1.000000 1.000000 1.000000' shininess='0.8' transparency='0.8'/>\n";
	vizout << "\t\t\t\t</Appearance>\n";
	vizout << "\t\t\t\t<Sphere radius='" << radius << "'/>\n";
	vizout << "\t\t\t</Shape>\n";
	vizout << "\t\t</Transform>\n";
}

inline void DrawSelfModel(ofstream &vizout, string element_name, Vec3 position, Vec4 rotation)
{
	vizout <<  "\t\t<Transform>\n";
	vizout <<  "\t\t\t<ProtoInstance name='" << element_name << "'>\n";
	vizout <<  "\t\t\t\t<fieldValue name='position' value='" << position.V3Str(' ') << "'/>\n";
	vizout <<  "\t\t\t\t<fieldValue name='transparency' value='0.6'/>\n";
	vizout <<  "\t\t\t\t<fieldValue name='inside_transparency' value='0.6'/>\n";
	vizout <<  "\t\t\t\t<fieldValue name='rotation' value='" << rotation.V4Str(' ') << "'/>\n";
	vizout <<  "\t\t\t</ProtoInstance>\n";
	vizout <<  "\t\t</Transform>\n";
}

inline void DrawSelfModelReduced(ofstream &vizout, Vec3 saxes, Vec3 position, Vec4 rotation, Vec3 color)
{
	vizout <<  "\t\t<Transform translation='" << position.V3Str(' ') << "' rotation='" << rotation.V4Str(' ') << "'>\n";
	vizout <<  "\t\t\t<Transform scale='" << saxes.V3Str(' ') << "'>\n";
	vizout <<  "\t\t\t\t<Shape>\n";
	vizout <<  "\t\t\t\t\t<Appearance>\n";
	vizout <<  "\t\t\t\t\t\t<Material diffuseColor='" << color.V3Str(' ') << "' shininess='0.8' transparency='0.6'/>\n";
	vizout <<  "\t\t\t\t\t</Appearance>\n";
	vizout <<  "\t\t\t\t\t<Sphere radius='1.0'/>\n";
	vizout <<  "\t\t\t\t</Shape>\n";
	vizout <<  "\t\t\t</Transform>\n";
	vizout <<  "\t\t</Transform>\n";
}

inline void DrawData(ofstream &vizout, string datapoints, Vec3 color)
{
	unsigned int nr_points=space_number(datapoints)/3;
	if(nr_points>0){
		vizout << "\t\t<Transform>\n";
		vizout << "\t\t\t<Shape>\n";
		vizout << "\t\t\t\t<Appearance>\n";
		vizout << "\t\t\t\t\t<Material diffuseColor='1 1 1' emissiveColor='" << color.V3Str(' ') << "' shininess='0.8'/>\n";
		vizout << "\t\t\t\t</Appearance>\n";
		vizout << "\t\t\t\t<IndexedLineSet coordIndex='";
		for(unsigned int i=0; i<nr_points; i++){
			if(i>0) vizout << " ";
			vizout << i;
		}
		vizout << "'>\n";
		vizout << "\t\t\t\t\t<Coordinate point='" << datapoints << "'/>\n";
		vizout << "\t\t\t\t</IndexedLineSet>\n";
		vizout << "\t\t\t</Shape>\n";
		vizout << "\t\t</Transform>\n";
	} else cout << "\nWARNING: Could not draw potentials (no data points).\n";
}

inline double VLJ_LOD(Simulation_Attribs* parameters, double rmin, double wA, double wB, double r, bool LJ_adjust_width)
{
	double VLJ;
	double r_over_r=rmin/r;
	if(LJ_adjust_width && (parameters->avg_width>0.0)){
		double dist=r-(rmin-(wA+wB));
		if(dist<0.0) dist=EPS;
		r_over_r=(wA+wB)/dist;
	}
	switch(parameters->vdwtype){
		case 1:
		case 4:
		default:
			r_over_r = qqpwr(r_over_r,parameters->LJexp[0]);
			VLJ=parameters->r*r_over_r*(r_over_r-1.0);
			break;
		case 2:
			// LJ + Bruce correction around r=rmin. Tanh functions have been replaced with a Gaussian for speed
			VLJ=parameters->r*(qqpwr(r_over_r,parameters->LJexp[1])-qqpwr(r_over_r,parameters->LJexp[0])) + parameters->Solvent[1]*exp(parameters->Solvent[0]*(r - parameters->Solvent[2])*(r - parameters->Solvent[2]));
			break;
		case 3:
		case 5:
			VLJ=parameters->r*qqpwr(r_over_r,parameters->LJexp[0]<<1); // bitshift left by one is multiplication by 2
			break;
	}
	return VLJ;
}

inline double VLJatoms(Config_Data* configuration, unsigned int lod_nr, unsigned int lod_level, unsigned int element_nr, Vec3 rvec, double &r_partner)
{
	double VLJatomistic=0.0;
	
	Element_Group* group=configuration->lods[lod_nr]->groups[0];
	Element_Group* lod_group=configuration->lods[lod_nr]->groups[lod_level+1];
	unsigned int element_count=lod_group->Type->LOD->element_groups[lod_group->levelofdetail-1][lod_group->Type->LOD->ellipsoid_counts[lod_group->levelofdetail-1][element_nr]];
	
	for(unsigned int l=0; l<element_count; l++){
		Element* curr=&configuration->group_elements[group->elements[lod_group->Type->LOD->element_groups[lod_group->levelofdetail-1][lod_group->Type->LOD->ellipsoid_counts[lod_group->levelofdetail-1][element_nr]+l+1]-1]];
		Vec3 pv=rvec-curr->center;
		double rc=pv.V3Norm();
		double rminc=EllipsoidRmin(pv,curr->MyType->saxes,curr->rot)+r_partner;
		double disp;
		switch(configuration->vdwtype){
			case 1:
			case 4:
			default:
				disp = qqpwr(rminc/rc,configuration->LJexp[0]);
				VLJatomistic+=sqrt(curr->MyType->Vvdw)*configuration->r*disp*(disp-1.0);
				break;
			case 2:
				disp = rminc/rc;
				// LJ + Bruce correction around r=rmin. Tanh functions have been replaced with a Gaussian for speed
				VLJatomistic+=sqrt(curr->MyType->Vvdw)*configuration->r*(qqpwr(disp,configuration->LJexp[1])-qqpwr(disp,configuration->LJexp[0])) + configuration->Solvent[1]*exp(configuration->Solvent[0]*(rc - configuration->Solvent[2])*(rc - configuration->Solvent[2]));
				break;
			case 3:
			case 5:
				VLJatomistic+=sqrt(curr->MyType->Vvdw)*configuration->r*qqpwr(rminc/rc,configuration->LJexp[0]<<1); // bitshift left by one is multiplication by 2
				break;
		}
	}
	return VLJatomistic;
}

inline double VLJatoms(Config_Data* configuration, unsigned int lod_nr, unsigned int lod_level, unsigned int element_nr, Element* partner_elements, Vec3 rvec)
{
	double VLJatomistic=0.0;
	
	Element_Group* group=configuration->lods[lod_nr]->groups[0];
	Element_Group* lod_group=configuration->lods[lod_nr]->groups[lod_level+1];
	unsigned int element_count=lod_group->Type->LOD->element_groups[lod_group->levelofdetail-1][lod_group->Type->LOD->ellipsoid_counts[lod_group->levelofdetail-1][element_nr]];
	
	for(unsigned int l=0; l<element_count; l++){
		Element* curr=&configuration->group_elements[group->elements[lod_group->Type->LOD->element_groups[lod_group->levelofdetail-1][lod_group->Type->LOD->ellipsoid_counts[lod_group->levelofdetail-1][element_nr]+l+1]-1]];
		for(unsigned int p=0; p<element_count; p++){
			Element* partner=&(partner_elements[p]);
			Vec3 pv=rvec+partner->center-curr->center;
			double rc=pv.V3Norm();
			double rminc=EllipsoidRmin(pv,curr->MyType->saxes,curr->rot)+EllipsoidRmin(pv,partner->MyType->saxes,partner->rot); // fortunately ellipsoids are symmetric around the center, so EllipsoidRmin(pv)=EllipsoidRmin(-pv)
			double disp;
			switch(configuration->vdwtype){
				case 1:
				case 4:
				default:
					disp = qqpwr(rminc/rc,configuration->LJexp[0]);
					VLJatomistic+=sqrt(curr->MyType->Vvdw*partner->MyType->Vvdw)*configuration->r*disp*(disp-1.0);
					break;
				case 2:
					disp = rminc/rc;
					// LJ + Bruce correction around r=rmin. Tanh functions have been replaced with a Gaussian for speed
					VLJatomistic+=sqrt(curr->MyType->Vvdw*partner->MyType->Vvdw)*configuration->r*(qqpwr(disp,configuration->LJexp[1])-qqpwr(disp,configuration->LJexp[0])) + configuration->Solvent[1]*exp(configuration->Solvent[0]*(rc - configuration->Solvent[2])*(rc - configuration->Solvent[2]));
					break;
				case 3:
				case 5:
					VLJatomistic+=sqrt(curr->MyType->Vvdw*partner->MyType->Vvdw)*configuration->r*qqpwr(rminc/rc,configuration->LJexp[0]<<1); // bitshift left by one is multiplication by 2
					break;
			}
		}
	}
	return VLJatomistic;
}

inline void DrawLJPotentials(ofstream &vizout, Config_Data* configuration, unsigned int lod_nr, unsigned int lod_level, unsigned int element_nr, Simulation_Attribs* parameters, double eps_avg, axis_t &axis)
{
	Element_Group* lod_group=configuration->lods[lod_nr]->groups[lod_level+1];
	Element* element=&configuration->group_elements[lod_group->elements[element_nr]];
	double rmin_draw=mind(element->MyType->saxes.vec,3)/2.0;
	
	double VLJamplitude=maxd(element->MyType->saxes.vec,3)*0.5;
	
	// x-Axis
	unsigned int t=theta_res/2; // is cos(pi/2)
	unsigned int p=0;
	double rmin = element->MyType->saxes.vec[0]+parameters->rp;
	
	string coords_VLJatomistic="";
	string coords_VLJ="";
	string coords_VLJ_texeps="";
	
	// Create file
	ofstream data;
	string dataname = lod_group->Type->name+"_"+int2str(element_nr+1)+"_x.dat";
	data.open(dataname.c_str());
	if(data.fail()){
		cout << "Unable to open output file.\n";
		exit(1);
	}
	data << "#r_diff:\t" << element->MyType->r_diff[t*phi_res+p] << "\n";
	data << "#avg_width:\t" << element->MyType->avg_width << "\n";
	data << "#r\tatomistic\tbest-fit epsilon\tepsilon texture\n";
	
	Vec3 e_x, fx;
	e_x=Vec3(1.0,0.0,0.0);
	fx=Vec3(0.0,0.0,1.0);
	double eps_corr=1.0;
	if(configuration->LJ_interaction_area) eps_corr=IA(element->MyType,e_x);
	for(unsigned int k=0; k<=2048; k++){
		double r=(double)k/2048.0*parameters->rmax+rmin_draw;
		if(element->MyType->eps_texture[t*phi_res+p]>EPS){
			// Calculate LOD ellipsoid Lennard-Jones potential
			double VLJ=VLJ_LOD(parameters,rmin,parameters->avg_width,parameters->rp,r,configuration->LJ_adjust_width);
			// Calculate Lennard-Jones potential of fully atomistic model
			double VLJatomistic=VLJatoms(configuration,lod_nr,lod_level,element_nr,element->rot*(e_x*r)+element->center,parameters->rp);
			data << r << "\t" << VLJatomistic << "\t" << VLJ*eps_avg*eps_corr << "\t" << VLJ*eps_avg*element->MyType->eps_texture[t*phi_res+p] << "\n";
			if(((VLJatomistic<2.0) && (r<rmin)) || (VLJatomistic<-0.01)){
				if(coords_VLJatomistic!="") coords_VLJatomistic+=" ";
				coords_VLJatomistic+=(e_x*r+fx*VLJatomistic/eps_avg*VLJamplitude).V3Str(' ');
			}
			if(((VLJ<2.0) && (r<rmin)) || (VLJ<-0.01)){
				if(coords_VLJ!=""){
					coords_VLJ+=" ";
					coords_VLJ_texeps+=" ";
				}
				coords_VLJ+=(e_x*r+fx*eps_corr*VLJ*VLJamplitude).V3Str(' ');
				coords_VLJ_texeps+=(e_x*r+fx*VLJ*element->MyType->eps_texture[t*phi_res+p]*VLJamplitude).V3Str(' ');
				axis.length.vec[0] = r;
			}
		}
	}
	data.close();
	DrawData(vizout,coords_VLJatomistic,Vec3(1.0,0.0,0.0));
	DrawData(vizout,coords_VLJ,Vec3(1.0,1.0,1.0));
	DrawData(vizout,coords_VLJ_texeps,Vec3(0.0,1.0,1.0));
	DrawCalcSphere(vizout,parameters->rp,e_x*axis.length.vec[0]);
	// y-Axis
	t=theta_res/2;
	p=phi_res/4;
	rmin = element->MyType->saxes.vec[1]+parameters->rp;
	
	dataname = lod_group->Type->name+"_"+int2str(element_nr+1)+"_y.dat";
	data.open(dataname.c_str());
	if(data.fail()){
		cout << "Unable to open output file.\n";
		exit(1);
	}
	data << "#r_diff:\t" << element->MyType->r_diff[t*phi_res+p] << "\n";
	data << "#avg_width:\t" << element->MyType->avg_width << "\n";
	data << "#r\tatomistic\tbest-fit epsilon\tepsilon texture\n";
	
	coords_VLJatomistic="";
	coords_VLJ="";
	coords_VLJ_texeps="";
	e_x=Vec3(0.0,1.0,0.0);
	fx=Vec3(0.0,0.0,1.0);
	if(configuration->LJ_interaction_area) eps_corr=IA(element->MyType,e_x);
	for(unsigned int k=0; k<=2048; k++){
		double r=(double)k/2048.0*parameters->rmax+rmin_draw;
		if(element->MyType->eps_texture[t*phi_res+p]>EPS){
			// Calculate LOD ellipsoid Lennard-Jones potential
			double VLJ=VLJ_LOD(parameters,rmin,parameters->avg_width,parameters->rp,r,configuration->LJ_adjust_width);
			// Calculate Lennard-Jones potential of fully atomistic model
			double VLJatomistic=VLJatoms(configuration,lod_nr,lod_level,element_nr,element->rot*(e_x*r)+element->center,parameters->rp);
			data << r << "\t" << VLJatomistic << "\t" << VLJ*eps_avg*eps_corr << "\t" << VLJ*eps_avg*element->MyType->eps_texture[t*phi_res+p] << "\n";
			if(((VLJatomistic<2.0) && (r<rmin)) || (VLJatomistic<-0.01)){
				if(coords_VLJatomistic!="") coords_VLJatomistic+=" ";
				coords_VLJatomistic+=(e_x*r+fx*VLJatomistic/eps_avg*VLJamplitude).V3Str(' ');
			}
			if(((VLJ<2.0) && (r<rmin)) || (VLJ<-0.01)){
				if(coords_VLJ!=""){
					coords_VLJ+=" ";
					coords_VLJ_texeps+=" ";
				}
				coords_VLJ+=(e_x*r+fx*eps_corr*VLJ*VLJamplitude).V3Str(' ');
				coords_VLJ_texeps+=(e_x*r+fx*VLJ*element->MyType->eps_texture[t*phi_res+p]*VLJamplitude).V3Str(' ');
				axis.length.vec[1] = r;
			}
		}
	}
	data.close();
	DrawData(vizout,coords_VLJatomistic,Vec3(1.0,0.0,0.0));
	DrawData(vizout,coords_VLJ,Vec3(1.0,1.0,1.0));
	DrawData(vizout,coords_VLJ_texeps,Vec3(0.0,1.0,1.0));
	DrawCalcSphere(vizout,parameters->rp,e_x*axis.length.vec[1]);
	// z-Axis
	t=theta_res; // special value
	p=0;
	rmin = element->MyType->saxes.vec[2]+parameters->rp;
	
	dataname = lod_group->Type->name+"_"+int2str(element_nr+1)+"_z.dat";
	data.open(dataname.c_str());
	if(data.fail()){
		cout << "Unable to open output file.\n";
		exit(1);
	}
	data << "#r_diff:\t" << element->MyType->r_diff[t*phi_res+p] << "\n";
	data << "#avg_width:\t" << element->MyType->avg_width << "\n";
	data << "#r\tatomistic\tbest-fit epsilon\tepsilon texture\n";
	coords_VLJatomistic="";
	coords_VLJ="";
	coords_VLJ_texeps="";
	e_x=Vec3(0.0,0.0,1.0);
	fx=Vec3(0.0,-1.0,0.0);
	if(configuration->LJ_interaction_area) eps_corr=IA(element->MyType,e_x);
	for(unsigned int k=0; k<=2048; k++){
		double r=(double)k/2048.0*parameters->rmax+rmin_draw;
		if(element->MyType->eps_texture[t*phi_res+p]>EPS){
			// Calculate LOD ellipsoid Lennard-Jones potential
			double VLJ=VLJ_LOD(parameters,rmin,parameters->avg_width,parameters->rp,r,configuration->LJ_adjust_width);
			// Calculate Lennard-Jones potential of fully atomistic model
			double VLJatomistic=VLJatoms(configuration,lod_nr,lod_level,element_nr,element->rot*(e_x*r)+element->center,parameters->rp);
			data << r << "\t" << VLJatomistic << "\t" << VLJ*eps_avg*eps_corr << "\t" << VLJ*eps_avg*element->MyType->eps_texture[t*phi_res+p] << "\n";
			if(((VLJatomistic<2.0) && (r<rmin)) || (VLJatomistic<-0.01)){
				if(coords_VLJatomistic!="") coords_VLJatomistic+=" ";
				coords_VLJatomistic+=(e_x*r+fx*VLJatomistic/eps_avg*VLJamplitude).V3Str(' ');
			}
			if(((VLJ<2.0) && (r<rmin)) || (VLJ<-0.01)){
				if(coords_VLJ!=""){
					coords_VLJ+=" ";
					coords_VLJ_texeps+=" ";
				}
				coords_VLJ+=(e_x*r+fx*eps_corr*VLJ*VLJamplitude).V3Str(' ');
				coords_VLJ_texeps+=(e_x*r+fx*VLJ*element->MyType->eps_texture[t*phi_res+p]*VLJamplitude).V3Str(' ');
				axis.length.vec[2] = r;
			}
		}
	}
	data.close();
	DrawData(vizout,coords_VLJatomistic,Vec3(1.0,0.0,0.0));
	DrawData(vizout,coords_VLJ,Vec3(1.0,1.0,1.0));
	DrawData(vizout,coords_VLJ_texeps,Vec3(0.0,1.0,1.0));
	DrawCalcSphere(vizout,parameters->rp,e_x*axis.length.vec[2]);
}

inline void DrawLJinteraction(ofstream &vizout, Config_Data* configuration, unsigned int lod_nr, unsigned int lod_level, unsigned int element_nr, Simulation_Attribs* parameters, double eps_avg, axis_t &axis)
{
	Element_Group* group=configuration->lods[lod_nr]->groups[0];
	Element_Group* lod_group=configuration->lods[lod_nr]->groups[lod_level+1];
	Element* element=&configuration->group_elements[lod_group->elements[element_nr]];
	unsigned int element_count=lod_group->Type->LOD->element_groups[lod_group->levelofdetail-1][lod_group->Type->LOD->ellipsoid_counts[lod_group->levelofdetail-1][element_nr]];
	double rmin_draw=mind(element->MyType->saxes.vec,3)/2.0;
	
	double VLJamplitude=maxd(element->MyType->saxes.vec,3)*0.5;
	
	// x-Axis (draw x-x and x-y interaction)
	unsigned int t=theta_res/2; // is cos(pi/2)
	unsigned int p=0;
	double rmin = element->MyType->saxes.vec[0];
	
	string coords_VLJ_x="";
	string coords_VLJ_y="";
	string coords_VLJ_z="";
	
	string coords_VLJsingle_x="";
	string coords_VLJsingle_y="";
	string coords_VLJsingle_z="";
	
	string coords_VLJatomistic_x="";
	string coords_VLJatomistic_y="";
	string coords_VLJatomistic_z="";
	
	Vec3 e_x, fx;
	e_x=Vec3(1.0,0.0,0.0);
	fx=Vec3(0.0,0.0,1.0);
	
	Vec4 rot_x(0.0,0.0,1.0,pi);
	Vec4 rot_y(0.0,0.0,1.0,pi+pi/2);
	Vec4 rot_z;
	
	Mat33 rotm_x=AxisAngle2Rot(rot_x);
	Mat33 rotm_y=AxisAngle2Rot(rot_y);
	Mat33 rotm_z;
	
	Element* partner_x = new Element[element_count];
	Element* partner_y = new Element[element_count];
	Element* partner_z = new Element[element_count];
	
	for(unsigned int i=0; i<element_count; i++){
		Element* curr=&configuration->group_elements[group->elements[lod_group->Type->LOD->element_groups[lod_group->levelofdetail-1][lod_group->Type->LOD->ellipsoid_counts[lod_group->levelofdetail-1][element_nr]+i+1]-1]];
		partner_x[i].MyType=curr->MyType;
		partner_y[i].MyType=curr->MyType;
		
		partner_x[i].center=element->rot*(rotm_x*(element->rot.M3Transpose()*(curr->center-element->center)))+element->center;
		partner_y[i].center=element->rot*(rotm_y*(element->rot.M3Transpose()*(curr->center-element->center)))+element->center;
		
		partner_x[i].rot=element->rot*(rotm_x*(element->rot.M3Transpose()*curr->rot));
		partner_y[i].rot=element->rot*(rotm_y*(element->rot.M3Transpose()*curr->rot));
	}

	// Create file
	ofstream data;
	string dataname = lod_group->Type->name+"_"+int2str(element_nr+1)+"_interactions_x.dat";
	data.open(dataname.c_str());
	if(data.fail()){
		cout << "Unable to open output file.\n";
		exit(1);
	}
	data << "#r\tatomistic x-x\tatomistic x-y\tbest-fit epsilon x-x\tbest-fit epsilon x-y\tepsilon texture x-x\tepsilon texture x-y\n";
	
	// idea is that average LJ epsilon depends on the average size of the blob out there (eps_avg) and the interaction area ratio with the average blob area
	Vec3 ia, eps;
	ia.vec[0]=eps_avg*IA(element->MyType,Vec3(1.0,0.0,0.0));
	ia.vec[1]=eps_avg*IA(element->MyType,Vec3(0.0,1.0,0.0));
	ia.vec[2]=eps_avg*IA(element->MyType,Vec3(0.0,0.0,1.0));
	
	// texture intrinsically carry their own information about how much area is seen from the blob out there by how they're set up, base LJ energy again is dependent on average blob size
	// Rules for t and p:
	// t = (1.0-cos(theta))*theta_res/2
	// p = (phi_res+(phi*phi_res/(2*pi)) % phi_res;
	// x face:
	// theta = 90 degree, phi = 0 degree
	// t = (1-cos(pi/2))*theta_res/2 = theta_res/2
	// p = (phi_res + (0*phi_res/(2*pi))) % phi_res = phi_res % phi_res = 0
	// -y face:
	// theta = 90 degree, phi = -90 degree
	// t = (1-cos(pi/2))*theta_res/2 = theta_res/2
	// p = (phi_res + (-pi/2*phi_res/(2*pi))) % phi_res = 3/4*phi_res % phi_res = 3/4*phi_res
	eps.vec[0]=eps_avg*element->MyType->eps_texture[(theta_res/2)*phi_res];
	eps.vec[1]=eps_avg*element->MyType->eps_texture[(theta_res/2)*phi_res+(3*phi_res/4)];
	double source_epsilon=eps_avg*element->MyType->eps_texture[t*phi_res+p];
	
	for(unsigned int k=0; k<=2048; k++){
		double r=(double)k/2048.0*parameters->rmax+rmin_draw;
		
		if(element->MyType->eps_texture[t*phi_res+p]>EPS){
			// Calculate LOD ellipsoid Lennard-Jones potential
			double VLJ_x=VLJ_LOD(parameters,rmin+element->MyType->saxes.vec[0],parameters->avg_width,parameters->avg_width,r,configuration->LJ_adjust_width);
			double VLJ_y=VLJ_LOD(parameters,rmin+element->MyType->saxes.vec[1],parameters->avg_width,parameters->avg_width,r,configuration->LJ_adjust_width);
			if(((VLJ_x<2.0) && (r>0.8*(rmin+element->MyType->saxes.vec[0]))) && (r<2.5*(rmin+element->MyType->saxes.vec[0]))){
				if(coords_VLJ_x!=""){
					coords_VLJ_x+=" ";
					coords_VLJsingle_x+=" ";
				}
				coords_VLJ_x+=(e_x*r+fx*VLJ_x*eps.vec[0]*source_epsilon/(eps_avg*eps_avg)*VLJamplitude).V3Str(' ');
				if(configuration->LJ_interaction_area) coords_VLJsingle_x+=(e_x*r+fx*VLJ_x*ia.vec[0]*ia.vec[0]/(eps_avg*eps_avg)*VLJamplitude).V3Str(' '); else coords_VLJsingle_x+=(e_x*r+fx*VLJ_x*VLJamplitude).V3Str(' ');
				axis.length.vec[0] = r;
			}
			if(((VLJ_y<2.0) && (r>0.8*(rmin+element->MyType->saxes.vec[1]))) && (r<2.5*(rmin+element->MyType->saxes.vec[1]))){
				if(coords_VLJ_y!=""){
					coords_VLJ_y+=" ";
					coords_VLJsingle_y+=" ";
				}
				coords_VLJ_y+=(AxisAngle2Rot(Vec4(e_x.vec[0],e_x.vec[1],e_x.vec[2],pi/2.0))*(e_x*r+fx*VLJ_y*eps.vec[1]*source_epsilon/(eps_avg*eps_avg)*VLJamplitude)).V3Str(' ');
				if(configuration->LJ_interaction_area) coords_VLJsingle_y+=(AxisAngle2Rot(Vec4(e_x.vec[0],e_x.vec[1],e_x.vec[2],pi/2.0))*(e_x*r+fx*VLJ_y*ia.vec[1]*ia.vec[0]/(eps_avg*eps_avg)*VLJamplitude)).V3Str(' '); else coords_VLJsingle_y+=(AxisAngle2Rot(Vec4(e_x.vec[0],e_x.vec[1],e_x.vec[2],pi/2.0))*(e_x*r+fx*VLJ_y*VLJamplitude)).V3Str(' ');
				axis.length.vec[0] = r;
			}
			
			// Calculate Lennard-Jones potential of fully atomistic model
			double VLJatomistic_x=VLJatoms(configuration,lod_nr,lod_level,element_nr,partner_x,element->rot*(e_x*r))/(eps_avg*eps_avg);
			double VLJatomistic_y=VLJatoms(configuration,lod_nr,lod_level,element_nr,partner_y,element->rot*(e_x*r))/(eps_avg*eps_avg);
			if(((VLJatomistic_x<2.0) && (r>0.8*(rmin+element->MyType->saxes.vec[0]))) && (r<2.5*(rmin+element->MyType->saxes.vec[0]))){
				if(coords_VLJatomistic_x!="") coords_VLJatomistic_x+=" ";
				coords_VLJatomistic_x+=(e_x*r+fx*VLJatomistic_x*VLJamplitude).V3Str(' ');
				axis.length.vec[0] = r;
			}
			if(((VLJatomistic_y<2.0) && (r>0.8*(rmin+element->MyType->saxes.vec[1]))) && (r<2.5*(rmin+element->MyType->saxes.vec[1]))){
				if(coords_VLJatomistic_y!="") coords_VLJatomistic_y+=" ";
				coords_VLJatomistic_y+=(AxisAngle2Rot(Vec4(e_x.vec[0],e_x.vec[1],e_x.vec[2],pi/2.0))*(e_x*r+fx*VLJatomistic_y*VLJamplitude)).V3Str(' ');
				axis.length.vec[0] = r;
			}
			data << r << "\t" << eps_avg*eps_avg*VLJatomistic_x << "\t" << eps_avg*eps_avg*VLJatomistic_y << "\t";
			if(configuration->LJ_interaction_area) data << VLJ_x*ia.vec[0]*ia.vec[0] << "\t" << VLJ_y*ia.vec[1]*ia.vec[0] << "\t"; else data << VLJ_x*eps_avg*eps_avg << "\t" << VLJ_y*eps_avg*eps_avg << "\t";
			data << VLJ_x*eps.vec[0]*source_epsilon << "\t" << VLJ_y*eps.vec[1]*source_epsilon << "\n";
		}
	}
	data.close();
	DrawData(vizout,coords_VLJsingle_x,Vec3(1.0,1.0,1.0));
	DrawData(vizout,coords_VLJsingle_y,Vec3(1.0,1.0,1.0));
	
	DrawData(vizout,coords_VLJ_x,Vec3(1.0,1.0,0.0));
	DrawData(vizout,coords_VLJatomistic_x,Vec3(1.0,0.0,0.0));
	DrawData(vizout,coords_VLJ_y,Vec3(0.0,1.0,1.0));
	DrawData(vizout,coords_VLJatomistic_y,Vec3(1.0,0.0,0.0));
	
	double pos=axis.length.vec[0];
	DrawSelfModel(vizout,element->MyType->name+" ("+lod_group->Type->name+")",e_x*pos,rot_x); // x-face
	DrawSelfModelReduced(vizout,element->MyType->saxes,e_x*pos,rot_y,Vec3(0.0,1.0,1.0)); // y-face (aquamarin)
	
	// y-Axis (draw y-y and y-z interaction)
	t=theta_res/2; // is cos(pi/2)
	p=phi_res/4; // also cos(pi/2)
	rmin = element->MyType->saxes.vec[1];
	
	dataname = lod_group->Type->name+"_"+int2str(element_nr+1)+"_interactions_y.dat";
	data.open(dataname.c_str());
	if(data.fail()){
		cout << "Unable to open output file.\n";
		exit(1);
	}
	data << "#r\tatomistic y-y\tatomistic y-z\tbest-fit epsilon y-y\tbest-fit epsilon y-z\tepsilon texture y-y\tepsilon texture y-z\n";
	
	coords_VLJ_y="";
	coords_VLJ_z="";
	
	coords_VLJsingle_y="";
	coords_VLJsingle_z="";
	
	coords_VLJatomistic_y="";
	coords_VLJatomistic_z="";
	
	e_x=Vec3(0.0,1.0,0.0);
	fx=Vec3(0.0,0.0,1.0);
	
	rot_y=Vec4(0.0,0.0,1.0,pi);
	rot_z=Vec4(1.0,0.0,0.0,pi/2);
	
	rotm_y=AxisAngle2Rot(rot_y);
	rotm_z=AxisAngle2Rot(rot_z);
	
	for(unsigned int i=0; i<element_count; i++){
		Element* curr=&configuration->group_elements[group->elements[lod_group->Type->LOD->element_groups[lod_group->levelofdetail-1][lod_group->Type->LOD->ellipsoid_counts[lod_group->levelofdetail-1][element_nr]+i+1]-1]];
		partner_y[i].MyType=curr->MyType;
		partner_z[i].MyType=curr->MyType;
		
		partner_y[i].center=element->rot*(rotm_y*(element->rot.M3Transpose()*(curr->center-element->center)))+element->center;
		partner_z[i].center=element->rot*(rotm_z*(element->rot.M3Transpose()*(curr->center-element->center)))+element->center;
		
		partner_y[i].rot=element->rot*(rotm_y*(element->rot.M3Transpose()*curr->rot));
		partner_z[i].rot=element->rot*(rotm_z*(element->rot.M3Transpose()*curr->rot));
	}
	
	// Rules for t and p:
	// t = (1.0-cos(theta))*theta_res/2
	// p = (phi_res+(phi*phi_res/(2*pi)) % phi_res;
	// y face:
	// theta = 90 degree, phi = pi/2 degree
	// t = (1-cos(pi/2))*theta_res/2 = theta_res/2
	// p = (phi_res + (pi/2*phi_res/(2*pi))) % phi_res = 5/4*phi_res % phi_res = phi_res/4
	// -z face:
	// theta = 180 degree, phi = 0 degree
	// t = (1-cos(pi))*theta_res/2 = theta_res
	// p = (phi_res + (0*phi_res/(2*pi))) % phi_res = phi_res % phi_res = 0
	eps.vec[1]=eps_avg*element->MyType->eps_texture[(theta_res/2)*phi_res+(phi_res/4)];
	// poles are special values
	eps.vec[2]=eps_avg*element->MyType->eps_texture[theta_res*phi_res+1];
	source_epsilon=eps_avg*element->MyType->eps_texture[t*phi_res+p];
	
	for(unsigned int k=0; k<=2048; k++){
		double r=(double)k/2048.0*parameters->rmax+rmin_draw;
		
		if(element->MyType->eps_texture[t*phi_res+p]>EPS){
			// Calculate LOD ellipsoid Lennard-Jones potential
			double VLJ_y=VLJ_LOD(parameters,rmin+element->MyType->saxes.vec[1],parameters->avg_width,parameters->avg_width,r,configuration->LJ_adjust_width);
			double VLJ_z=VLJ_LOD(parameters,rmin+element->MyType->saxes.vec[2],parameters->avg_width,parameters->avg_width,r,configuration->LJ_adjust_width);
			if(((VLJ_y<2.0) && (r>0.8*(rmin+element->MyType->saxes.vec[1]))) && (r<2.5*(rmin+element->MyType->saxes.vec[1]))){
				if(coords_VLJ_y!=""){
					coords_VLJ_y+=" ";
					coords_VLJsingle_y+=" ";
				}
				coords_VLJ_y+=(e_x*r+fx*VLJ_y*eps.vec[1]*source_epsilon/(eps_avg*eps_avg)*VLJamplitude).V3Str(' ');
				if(configuration->LJ_interaction_area) coords_VLJsingle_y+=(e_x*r+fx*VLJ_y*ia.vec[1]*ia.vec[1]/(eps_avg*eps_avg)*VLJamplitude).V3Str(' '); else coords_VLJsingle_y+=(e_x*r+fx*VLJ_y*VLJamplitude).V3Str(' ');
				axis.length.vec[1] = r;
			}
			if(((VLJ_z<2.0) && (r>0.8*(rmin+element->MyType->saxes.vec[2]))) && (r<2.5*(rmin+element->MyType->saxes.vec[2]))){
				if(coords_VLJ_z!=""){
					coords_VLJ_z+=" ";
					coords_VLJsingle_z+=" ";
				}
				coords_VLJ_z+=(AxisAngle2Rot(Vec4(e_x.vec[0],e_x.vec[1],e_x.vec[2],pi/2.0))*(e_x*r+fx*VLJ_z*eps.vec[2]*source_epsilon/(eps_avg*eps_avg)*VLJamplitude)).V3Str(' ');
				if(configuration->LJ_interaction_area) coords_VLJsingle_z+=(AxisAngle2Rot(Vec4(e_x.vec[0],e_x.vec[1],e_x.vec[2],pi/2.0))*(e_x*r+fx*VLJ_z*ia.vec[2]*ia.vec[1]/(eps_avg*eps_avg)*VLJamplitude)).V3Str(' '); else coords_VLJsingle_z+=(AxisAngle2Rot(Vec4(e_x.vec[0],e_x.vec[1],e_x.vec[2],pi/2.0))*(e_x*r+fx*VLJ_z*VLJamplitude)).V3Str(' ');
				axis.length.vec[1] = r;
			}
			
			// Calculate Lennard-Jones potential of fully atomistic model
			double VLJatomistic_y=VLJatoms(configuration,lod_nr,lod_level,element_nr,partner_y,element->rot*(e_x*r))/(eps_avg*eps_avg);
			double VLJatomistic_z=VLJatoms(configuration,lod_nr,lod_level,element_nr,partner_z,element->rot*(e_x*r))/(eps_avg*eps_avg);
			if(((VLJatomistic_y<2.0) && (r>0.8*(rmin+element->MyType->saxes.vec[1]))) && (r<2.5*(rmin+element->MyType->saxes.vec[1]))){
				if(coords_VLJatomistic_y!="") coords_VLJatomistic_y+=" ";
				coords_VLJatomistic_y+=(e_x*r+fx*VLJatomistic_y*VLJamplitude).V3Str(' ');
				axis.length.vec[1] = r;
			}
			if(((VLJatomistic_z<2.0) && (r>0.8*(rmin+element->MyType->saxes.vec[2]))) && (r<2.5*(rmin+element->MyType->saxes.vec[2]))){
				if(coords_VLJatomistic_z!="") coords_VLJatomistic_z+=" ";
				coords_VLJatomistic_z+=(AxisAngle2Rot(Vec4(e_x.vec[0],e_x.vec[1],e_x.vec[2],pi/2.0))*(e_x*r+fx*VLJatomistic_z*VLJamplitude)).V3Str(' ');
				axis.length.vec[1] = r;
			}
			data << r << "\t" << eps_avg*eps_avg*VLJatomistic_y << "\t" << eps_avg*eps_avg*VLJatomistic_z << "\t";
			if(configuration->LJ_interaction_area) data << VLJ_y*ia.vec[1]*ia.vec[1] << "\t" << VLJ_z*ia.vec[2]*ia.vec[1] << "\t"; else data << VLJ_y*eps_avg*eps_avg << "\t" << VLJ_z*eps_avg*eps_avg << "\t";
			data << VLJ_y*eps.vec[1]*source_epsilon << "\t" << VLJ_z*eps.vec[2]*source_epsilon << "\n";
		}
	}
	data.close();
	DrawData(vizout,coords_VLJsingle_y,Vec3(1.0,1.0,1.0));
	DrawData(vizout,coords_VLJsingle_z,Vec3(1.0,1.0,1.0));
	
	DrawData(vizout,coords_VLJ_y,Vec3(0.0,1.0,1.0));
	DrawData(vizout,coords_VLJatomistic_y,Vec3(1.0,0.0,0.0));
	DrawData(vizout,coords_VLJ_z,Vec3(0.0,1.0,0.0));
	DrawData(vizout,coords_VLJatomistic_z,Vec3(1.0,0.0,0.0));
	
	pos=axis.length.vec[1];
	DrawSelfModel(vizout,element->MyType->name+" ("+lod_group->Type->name+")",e_x*pos,rot_y); // y-face
	DrawSelfModelReduced(vizout,element->MyType->saxes,e_x*pos,rot_z,Vec3(0.0,1.0,0.0)); // z-face
	
	// z-Axis (draw z-z and x-z interaction)
	t=theta_res; // is special value
	p=0; // also cos(0)
	rmin = element->MyType->saxes.vec[2];
	
	dataname = lod_group->Type->name+"_"+int2str(element_nr+1)+"_interactions_z.dat";
	data.open(dataname.c_str());
	if(data.fail()){
		cout << "Unable to open output file.\n";
		exit(1);
	}
	data << "#r\tatomistic z-z\tatomistic x-z\tbest-fit epsilon z-z\tbest-fit epsilon x-z\tepsilon texture z-z\tepsilon texture x-z\n";
	
	coords_VLJ_x="";
	coords_VLJ_z="";
	
	coords_VLJsingle_x="";
	coords_VLJsingle_z="";
	
	coords_VLJatomistic_x="";
	coords_VLJatomistic_z="";
	
	e_x=Vec3(0.0,0.0,1.0);
	fx=Vec3(0.0,-1.0,0.0);
	
	rot_x=Vec4(0.0,1.0,0.0,pi/2);
	rot_z=Vec4(0.0,1.0,0.0,pi);
	
	rotm_x=AxisAngle2Rot(rot_x);
	rotm_z=AxisAngle2Rot(rot_z);
	
	for(unsigned int i=0; i<element_count; i++){
		Element* curr=&configuration->group_elements[group->elements[lod_group->Type->LOD->element_groups[lod_group->levelofdetail-1][lod_group->Type->LOD->ellipsoid_counts[lod_group->levelofdetail-1][element_nr]+i+1]-1]];
		partner_x[i].MyType=curr->MyType;
		partner_z[i].MyType=curr->MyType;
		
		partner_x[i].center=element->rot*(rotm_x*(element->rot.M3Transpose()*(curr->center-element->center)))+element->center;
		partner_z[i].center=element->rot*(rotm_z*(element->rot.M3Transpose()*(curr->center-element->center)))+element->center;
		
		partner_x[i].rot=element->rot*(rotm_x*(element->rot.M3Transpose()*curr->rot));
		partner_z[i].rot=element->rot*(rotm_z*(element->rot.M3Transpose()*curr->rot));
	}
	
	// Rules for t and p:
	// t = (1.0-cos(theta))*theta_res/2
	// p = (phi_res+(phi*phi_res/(2*pi)) % phi_res;
	// x face:
	// theta = 90 degree, phi = 0 degree
	// t = (1-cos(pi/2))*theta_res/2 = theta_res/2
	// p = (phi_res + (0*phi_res/(2*pi))) % phi_res = phi_res % phi_res = 0
	// z face:
	// theta = 0 degree, phi = 0 degree
	// t = (1-cos(0))*theta_res/2 = 0
	// p = (phi_res + (0*phi_res/(2*pi))) % phi_res = phi_res % phi_res = 0
	eps.vec[0]=eps_avg*element->MyType->eps_texture[(theta_res/2)*phi_res];
	eps.vec[2]=eps_avg*element->MyType->eps_texture[theta_res*phi_res];
	source_epsilon=eps_avg*element->MyType->eps_texture[t*phi_res+p];
	
	for(unsigned int k=0; k<=2048; k++){
		double r=(double)k/2048.0*parameters->rmax+rmin_draw;
		
		if(element->MyType->eps_texture[t*phi_res+p]>EPS){
			// Calculate LOD ellipsoid Lennard-Jones potential
			double VLJ_x=VLJ_LOD(parameters,rmin+element->MyType->saxes.vec[0],parameters->avg_width,parameters->avg_width,r,configuration->LJ_adjust_width);
			double VLJ_z=VLJ_LOD(parameters,rmin+element->MyType->saxes.vec[2],parameters->avg_width,parameters->avg_width,r,configuration->LJ_adjust_width);
			if(((VLJ_z<2.0) && (r>0.8*(rmin+element->MyType->saxes.vec[2]))) && (r<2.5*(rmin+element->MyType->saxes.vec[2]))){
				if(coords_VLJ_z!=""){
					coords_VLJ_z+=" ";
					coords_VLJsingle_z+=" ";
				}
				coords_VLJ_z+=(e_x*r+fx*VLJ_z*eps.vec[2]*source_epsilon/(eps_avg*eps_avg)*VLJamplitude).V3Str(' ');
				if(configuration->LJ_interaction_area) coords_VLJsingle_z+=(e_x*r+fx*VLJ_z*ia.vec[2]*ia.vec[2]/(eps_avg*eps_avg)*VLJamplitude).V3Str(' '); else coords_VLJsingle_z+=(e_x*r+fx*VLJ_z*VLJamplitude).V3Str(' ');
				axis.length.vec[2] = r;
			}
			if(((VLJ_x<2.0) && (r>0.8*(rmin+element->MyType->saxes.vec[0]))) && (r<2.5*(rmin+element->MyType->saxes.vec[0]))){
				if(coords_VLJ_x!=""){
					coords_VLJ_x+=" ";
					coords_VLJsingle_x+=" ";
				}
				coords_VLJ_x+=(AxisAngle2Rot(Vec4(e_x.vec[0],e_x.vec[1],e_x.vec[2],pi/2.0))*(e_x*r+fx*VLJ_x*eps.vec[0]*source_epsilon/(eps_avg*eps_avg)*VLJamplitude)).V3Str(' ');
				if(configuration->LJ_interaction_area) coords_VLJsingle_x+=(AxisAngle2Rot(Vec4(e_x.vec[0],e_x.vec[1],e_x.vec[2],pi/2.0))*(e_x*r+fx*VLJ_x*ia.vec[0]*ia.vec[2]/(eps_avg*eps_avg)*VLJamplitude)).V3Str(' '); else coords_VLJsingle_x+=(AxisAngle2Rot(Vec4(e_x.vec[0],e_x.vec[1],e_x.vec[2],pi/2.0))*(e_x*r+fx*VLJ_x*VLJamplitude)).V3Str(' ');
				axis.length.vec[2] = r;
			}
			
			// Calculate Lennard-Jones potential of fully atomistic model
			double VLJatomistic_x=VLJatoms(configuration,lod_nr,lod_level,element_nr,partner_x,element->rot*(e_x*r))/(eps_avg*eps_avg);
			double VLJatomistic_z=VLJatoms(configuration,lod_nr,lod_level,element_nr,partner_z,element->rot*(e_x*r))/(eps_avg*eps_avg);
			if(((VLJatomistic_z<2.0) && (r>0.8*(rmin+element->MyType->saxes.vec[2]))) && (r<2.5*(rmin+element->MyType->saxes.vec[2]))){
				if(coords_VLJatomistic_z!="") coords_VLJatomistic_z+=" ";
				coords_VLJatomistic_z+=(e_x*r+fx*VLJatomistic_z*VLJamplitude).V3Str(' ');
				axis.length.vec[2] = r;
			}
			if(((VLJatomistic_x<2.0) && (r>0.8*(rmin+element->MyType->saxes.vec[0]))) && (r<2.5*(rmin+element->MyType->saxes.vec[0]))){
				if(coords_VLJatomistic_x!="") coords_VLJatomistic_x+=" ";
				coords_VLJatomistic_x+=(AxisAngle2Rot(Vec4(e_x.vec[0],e_x.vec[1],e_x.vec[2],pi/2.0))*(e_x*r+fx*VLJatomistic_x*VLJamplitude)).V3Str(' ');
				axis.length.vec[2] = r;
			}
			data << r << "\t" << eps_avg*eps_avg*VLJatomistic_z << "\t" << eps_avg*eps_avg*VLJatomistic_x << "\t";
			if(configuration->LJ_interaction_area) data << VLJ_z*ia.vec[2]*ia.vec[2] << "\t" << VLJ_x*ia.vec[0]*ia.vec[2] << "\t"; else data << VLJ_z*eps_avg*eps_avg << "\t" << VLJ_x*eps_avg*eps_avg << "\t";
			data << VLJ_z*eps.vec[2]*source_epsilon << "\t" << VLJ_x*eps.vec[0]*source_epsilon << "\n";
		}
	}
	data.close();
	DrawData(vizout,coords_VLJsingle_z,Vec3(1.0,1.0,1.0));
	DrawData(vizout,coords_VLJsingle_x,Vec3(1.0,1.0,1.0));
	
	DrawData(vizout,coords_VLJ_z,Vec3(0.0,1.0,0.0));
	DrawData(vizout,coords_VLJ_x,Vec3(0.0,1.0,1.0));
	DrawData(vizout,coords_VLJatomistic_z,Vec3(1.0,0.0,0.0));
	DrawData(vizout,coords_VLJatomistic_x,Vec3(1.0,0.0,0.0));
	
	pos=axis.length.vec[2];
	DrawSelfModel(vizout,element->MyType->name+" ("+lod_group->Type->name+")",e_x*pos,rot_z); // z-face
	DrawSelfModelReduced(vizout,element->MyType->saxes,e_x*pos,rot_x,Vec3(1.0,1.0,0.0)); // x-face
}

void epsilon_matching(MC_Config* config, Config_Data* configuration, unsigned int lod_nr, unsigned int lod_level, unsigned int element_nr, bool output_x3d, string &texture_string, string &partner_string)
{
	configuration->lods[lod_nr]->visual_original[lod_level]=true;
	Element_Group* group=configuration->lods[lod_nr]->groups[0];
	Element_Group* lod_group=configuration->lods[lod_nr]->groups[lod_level+1];
	Element* element=&configuration->group_elements[lod_group->elements[element_nr]];
	
	// Set up Lennard-Jones related parameters from configuration
	Simulation_Attribs* parameters=new Simulation_Attribs;
	parameters->saxes=create_double4(element->MyType->saxes);
	parameters->invsaxes2=create_double4(element->MyType->invsaxes2);
	parameters->center=create_double4(element->center);
	for(unsigned int j=0; j<9; j++){
		parameters->rot[j]=element->rot.mat[j/3][j%3];
		parameters->unit[j]=0.0;
	}
	parameters->avg_width=-1.0;
	parameters->inv_avg_area=-1.0;
	if(configuration->LJ_adjust_width) parameters->avg_width=element->MyType->avg_width;
	if(configuration->LJ_interaction_area) parameters->inv_avg_area=pi*element->MyType->inv_avg_area;
	parameters->sphericity=element->MyType->sphericity;
	parameters->unit[0]=1.0; parameters->unit[4]=1.0; parameters->unit[8]=1.0;
	parameters->vdwtype=configuration->vdwtype;
	parameters->LJexp[0]=configuration->LJexp[0]; parameters->LJexp[1]=configuration->LJexp[1];
	parameters->Solvent[0]=configuration->Solvent[0]; parameters->Solvent[1]=configuration->Solvent[1]; parameters->Solvent[2]=configuration->Solvent[2];
	
	unsigned int count_t=theta_res;
	unsigned int count_p=phi_res;
	unsigned int count_r=r_res;
	
	parameters->r_overscan=r_oversampling;
	parameters->r=configuration->r;
	
	// Set up current LOD ellipsoid's fully atomistic model
	unsigned int element_count=lod_group->Type->LOD->element_groups[lod_group->levelofdetail-1][lod_group->Type->LOD->ellipsoid_counts[lod_group->levelofdetail-1][element_nr]];
	Element_Dynamic* element_dyn=new Element_Dynamic[element_count];
	Element_Static* element_stat=new Element_Static[element_count];
	
	for(unsigned int i=0; i<element_count; i++){
		Element* curr=&configuration->group_elements[group->elements[lod_group->Type->LOD->element_groups[lod_group->levelofdetail-1][lod_group->Type->LOD->ellipsoid_counts[lod_group->levelofdetail-1][element_nr]+i+1]-1]];
		
		element_dyn[i].center=create_double4(curr->center);
		for(unsigned int j=0; j<9; j++) element_dyn[i].rot[j]=curr->rot.mat[j/3][j%3];
		
		element_stat[i].saxes_vdw=create_double4(curr->MyType->saxes);
		element_stat[i].saxes_vdw.w=sqrt(curr->MyType->Vvdw); // only use sqrt in calculation
	}
	parameters->weight_factor=2.0*sqrt(element->MyType->Vvdw)*configuration->beta; // = sum sqrt(eps) * 1/kT
//	if(!configuration->LJ_interaction_area) parameters->weight_factor=0.0;
	
	size_t global=count_t*count_p+2;
	cl_double* epsilons=new cl_double[(voxelcount+2*r_res)*EPS_ELEMENT_COUNT];
	for(unsigned int i=0; i<(voxelcount+2*r_res)*EPS_ELEMENT_COUNT; i++) epsilons[i]=0.0;
	
	// voxel partitioned linear in phi, cos(theta), and r^6
	parameters->dphi_3=1.0/(3.0*count_p);
	parameters->dcostheta=1.0/count_t; // individual theta=acos(1.0-t*dcostheta)
	
	// Set up OpenCL
	cl_context Context=NULL;
	cl_command_queue Queue=NULL;
	cl_device_type dev=CL_DEVICE_TYPE_ALL;
	if(configuration->use_gpu) dev=CL_DEVICE_TYPE_GPU;
	cl_device_id deviceID=GetDevice(dev,Context,Queue,true);
	cl_program program = CreateProgram(deviceID,Context,&OpenCL_code);
/*
 *                                                                              arg
 * Creating kernel: match_epsilon(__global Simulation_Attribs* attrib,           0
 *                                __global Element_Dynamic* element_dyn,         1
 *                                __global Element_Static* element_stat,         2
 *                                __global double* r_diff,                       3
 *                                __global double* output,                       4
 *                                const unsigned int k,                          5
 *                                const unsigned int element_count,              6
 *                                const unsigned int count_t,                    7
 *                                const unsigned int count_p,                    8
 *                                const unsigned int count_r GUID_ARG)           9
 */
	cl_kernel kernel = CreateKernel(program,"match_epsilon");
	
	int* Error=NULL;
#ifndef OpenCL
	Error = new int;
#endif
	cl_mem eps_attrib = clCreateBuffer(Context,CL_MEM_READ_ONLY,sizeof(Simulation_Attribs),NULL,Error);
	cl_mem dynamic = clCreateBuffer(Context,CL_MEM_READ_ONLY,sizeof(Element_Dynamic)*element_count,NULL,Error);
	cl_mem stat = clCreateBuffer(Context,CL_MEM_READ_ONLY,sizeof(Element_Static)*element_count,NULL,Error);
	cl_mem rd = clCreateBuffer(Context,CL_MEM_READ_ONLY,sizeof(cl_double)*(count_t*count_p+2),NULL,Error);
	cl_mem output = clCreateBuffer(Context,CL_MEM_READ_WRITE,sizeof(cl_double)*(voxelcount+2*r_res)*EPS_ELEMENT_COUNT,NULL,Error);
	if(!eps_attrib || !output || !dynamic || !stat){
		cout << "ERROR: Could not allocate device memory.\n";
		exit(1);
	}
	
	int Error_Code = clEnqueueWriteBuffer(Queue,dynamic,CL_TRUE,0,sizeof(Element_Dynamic)*element_count,element_dyn,0,NULL,NULL);
	Error_Code |= clEnqueueWriteBuffer(Queue,stat,CL_TRUE,0,sizeof(Element_Static)*element_count,element_stat,0,NULL,NULL);
	Error_Code |= clEnqueueWriteBuffer(Queue,rd,CL_TRUE,0,sizeof(cl_double)*(count_t*count_p+2),element->MyType->r_diff,0,NULL,NULL);
	if(Error_Code!=CL_SUCCESS){
		cout << "ERROR: Could not copy values to device.\n";
		exit(1);
	}
	
	Error_Code = clSetKernelArg(kernel, 0, sizeof(cl_mem), &eps_attrib);
	Error_Code|= clSetKernelArg(kernel, 1, sizeof(cl_mem), &dynamic);
	Error_Code|= clSetKernelArg(kernel, 2, sizeof(cl_mem), &stat);
	Error_Code|= clSetKernelArg(kernel, 3, sizeof(cl_mem), &rd);
	Error_Code|= clSetKernelArg(kernel, 4, sizeof(cl_mem), &output);
	Error_Code|= clSetKernelArg(kernel, 6, sizeof(unsigned int), &element_count);
	Error_Code|= clSetKernelArg(kernel, 7, sizeof(unsigned int), &count_t);
	Error_Code|= clSetKernelArg(kernel, 8, sizeof(unsigned int), &count_p);
	Error_Code|= clSetKernelArg(kernel, 9, sizeof(unsigned int), &count_r);
	if(Error_Code!=CL_SUCCESS){
		cout << "ERROR: Failed to set kernel arguments.\n";
		exit(1);
	}
	
	size_t local=0;
#ifdef CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE
	Error_Code = clGetKernelWorkGroupInfo(kernel, deviceID, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(local), &local, NULL);
#else
	Error_Code = clGetKernelWorkGroupInfo(kernel, deviceID, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local), &local, NULL);
#endif
	if(Error_Code!=CL_SUCCESS){
		cout << "ERROR: Could not determine workgroup information.\n";
		exit(1);
	}
	if(local>global) global=local;
	cout << "-> Preferred work group size is " << local << "\n";
	double newglobal=(global/local+(global%local>0))*local; // make global worksize a multiple of local one
	if(global!=newglobal){
		cout << "\t-> Adjusted number of work units from " << global << " to " << newglobal << "\n";
		global=newglobal;
	}
	
	double wsum, eps_avg, eps_old, vareps;
	eps_avg=2.0*sqrt(element->MyType->Vvdw);
	
	cout << "-> Running kernel: 0%";
	cout.flush();
	double tstart = clock();
	parameters->rp=100.0/0.95;
	unsigned int NrT=0;
	do{
		if(parameters->rp>1.0) parameters->rp-=parameters->rp/20.0; else parameters->rp-=0.05;
		if(parameters->rp<0.0) parameters->rp=0.0;
		NrT++;
	} while(parameters->rp>0.0);
	
	Vec3* eps_rT = new Vec3[NrT];
	double* tex_avg_rT = new double[NrT];
	
	parameters->rp=100.0/0.95;
	unsigned int percentage=0;
	unsigned int Nit=0;
	unsigned int Nrp=0;
	unsigned int r_idx=0;
	bool drawingrun=false;
	bool endrun=false;
	bool first=true;
	do{
		if(!drawingrun){
			if(parameters->rp>1.0) parameters->rp-=parameters->rp/20.0; else parameters->rp-=0.05;
			if(parameters->rp<0.0) parameters->rp=0.0;
		} else{
			drawingrun=false;
			endrun=true;
		}
		first=true;
		unsigned int Ncr=0;
		double texture_average;
		do{
			Ncr++;
			parameters->rmax=(parameters->rp+maxd(element->MyType->saxes.vec,3))*pow(0.01,-1.0/parameters->LJexp[0]); // maximum distance before (rmin/r)^LJexponent is less than LJ cutoff (0.01)
			parameters->rmin=mind(element->MyType->saxes.vec,3)/2.0+parameters->rp;
			parameters->r_r=pow(parameters->rmin/parameters->rmax,1.0/(double)R_COMPRESSION);
			parameters->k0=(unsigned int)(parameters->r_r/(1.0-parameters->r_r)*count_r);
			for(unsigned int i=0; i<(voxelcount+2*r_res)*EPS_ELEMENT_COUNT; i++) epsilons[i]=0.0;
			Error_Code = clEnqueueWriteBuffer(Queue,eps_attrib,CL_TRUE,0,sizeof(Simulation_Attribs),parameters,0,NULL,NULL);
			Error_Code |= clEnqueueWriteBuffer(Queue,output,CL_TRUE,0,sizeof(cl_double)*(voxelcount+2*r_res)*EPS_ELEMENT_COUNT,epsilons,0,NULL,NULL);
			if(Error_Code!=CL_SUCCESS){
				cout << "ERROR: Failed to set kernel arguments.\n";
				exit(1);
			}
			r_idx=0;
			do{
				Error_Code = clSetKernelArg(kernel, 5, sizeof(unsigned int), &r_idx);
				if(Error_Code!=CL_SUCCESS){
					cout << "ERROR: Failed to set kernel argument #3.\n";
					exit(1);
				}
				Error_Code = clEnqueueNDRangeKernel(Queue,kernel,1,NULL,&global,&local,0,NULL,NULL);
				if(Error_Code!=CL_SUCCESS){
					cout << "ERROR: Failed to execute kernel: " << Error_Code << "\n";
					exit(1);
				}
				clFinish(Queue);
				
				if(first){
					Nit++;
					double percent=100.0*double(Nit)/double(NrT*count_r*r_oversampling);
					if(percent>=percentage+5){
						percentage=(unsigned int)floor(percent/5.0)*5;
						if(percentage%20==0) cout << percentage << "%"; else cout << ".";
						cout.flush();
					}
				}
				r_idx++;
			} while(r_idx<count_r*r_oversampling);
			Error_Code = clEnqueueReadBuffer(Queue,output,CL_TRUE,0,sizeof(cl_double)*(voxelcount+2*r_res)*EPS_ELEMENT_COUNT,epsilons,0,NULL,NULL);
			if(Error_Code!=CL_SUCCESS){
				cout << "ERROR: Could not read results.\n";
				exit(1);
			}
			eps_old=eps_avg;
			// go over all datapoints to get average epsilon
			wsum=0.0;
			eps_avg=0.0;
			vareps=0.0;
			texture_average=0.0;
			for(unsigned int t=0; t<count_t; t++){
				for(unsigned int p=0; p<count_p; p++){
					double tex_wsum=0.0;
					double tex_eps_avg=0.0;
					double tex_vareps=0.0;
					for(int k=count_r-1; k>=0; k--){
						unsigned int idx=((t*count_p+p)*count_r+(unsigned int)k)*EPS_ELEMENT_COUNT;
						
						double c_eps=epsilons[idx]; // sum of w*eps
						double c_eps2=epsilons[idx+1]; // sum of w*eps*eps
						double ws=epsilons[idx+2]; // sum of w
						double curr_var=(c_eps2-c_eps*c_eps/ws)/ws;
						if(fabs(curr_var)<EPS) curr_var=0.0;
						if(ws>EPS){ // sanity checks
							double tex_curr_avg=tex_eps_avg+c_eps;
							double tex_curr_var=tex_vareps+c_eps2;
							tex_curr_avg/=(tex_wsum+ws);
							tex_curr_var/=(tex_wsum+ws);
							tex_curr_var-=tex_curr_avg*tex_curr_avg;
							if(4.0*tex_curr_var<tex_curr_avg*tex_curr_avg){
								wsum+=ws;
								tex_wsum+=ws;
								tex_eps_avg+=c_eps;
								tex_vareps+=c_eps2;
								eps_avg+=c_eps;
								vareps+=c_eps2;
							} else break; // abort if per direction std.dev. would be more than half the average value
						}
					}
					if(wsum>EPS) texture_average+=tex_eps_avg/tex_wsum;
				}
			}
			eps_avg/=wsum;
			texture_average/=count_t*count_p*eps_avg;
			vareps/=wsum;
			vareps-=eps_avg*eps_avg;
			parameters->weight_factor=eps_avg*configuration->beta; // = sum sqrt(eps) * 1/kT
/*			if(!configuration->LJ_interaction_area){
				parameters->weight_factor=0.0;
				eps_old=eps_avg; // no Boltzmann-weighting is going on
			}*/
			first=false;
		} while((fabs(eps_avg-eps_old)/eps_avg>0.05) && (Ncr<100));
		if(!endrun){
			eps_rT[Nrp].vec[0]=parameters->rp;
			eps_rT[Nrp].vec[1]=eps_avg;
			eps_rT[Nrp].vec[2]=1.0/vareps;
			tex_avg_rT[Nrp]=texture_average;
			Nrp++;
		}
		
		if(parameters->rp==0.0){
			drawingrun=true;
			parameters->rp=configuration->rT;
			parameters->rp=element->MyType->rT;
		}
	} while((parameters->rp>0.0) && !endrun);
	double tend = clock();
	cout << "\n\t-> " << Nit*global << " iterations, took " << (tend-tstart)/CLOCKS_PER_SEC << " s\n";
	
	if(element->MyType->eps_texture) delete[] element->MyType->eps_texture; // get rid of old texture if it exists
	element->MyType->eps_texture=new double[count_t*count_p+4];
	double tex_min=1.0;
	double tex_max=1.0;
	
	Vec3 direction;
	double phi;
	double texture_average=0.0;
	Vec3* texture_data=NULL;
	if(configuration->LJ_interaction_area) texture_data=new Vec3[count_t*count_p+2];
	for(unsigned int t=0; t<count_t; t++){
		for(unsigned int p=0; p<count_p; p++){
			phi=2.0*pi*(double)(p+0.5)/(double)phi_res;
			double cost=1.0-2.0*(double)(t+0.5)/(double)theta_res;
			double sint=sqrt(1.0-cost*cost);
			direction.vec[0]=sint*cos(phi);
			direction.vec[1]=sint*sin(phi);
			direction.vec[2]=cost;
			wsum=0.0;
			double tex_eps_avg=0.0;
			double tex_vareps=0.0;
			bool add=true;
			for(int k=count_r-1; k>=0; k--){
				unsigned int idx=((t*count_p+p)*count_r+(unsigned int)k)*EPS_ELEMENT_COUNT;
				
				double c_eps=epsilons[idx]; // sum of w*eps
				double c_eps2=epsilons[idx+1]; // sum of w*eps*eps
				double ws=epsilons[idx+2]; // sum of w
				double curr_var=(c_eps2-c_eps*c_eps/ws)/ws;
				if(fabs(curr_var)<EPS) curr_var=0.0;
				
				if((ws>EPS) && add){ // sanity checks
					epsilons[idx]/=ws;
					epsilons[idx+1]=sqrt(curr_var);
					double tex_curr_avg=tex_eps_avg+c_eps;
					double tex_curr_var=tex_vareps+c_eps2;
					tex_curr_avg/=(wsum+ws);
					tex_curr_var/=(wsum+ws);
					tex_curr_var-=tex_curr_avg*tex_curr_avg;
					if(4.0*tex_curr_var<tex_curr_avg*tex_curr_avg){ // do not accept inward values if std.dev. is more than half the average
						wsum+=ws;
						tex_eps_avg+=c_eps;
						tex_vareps+=c_eps2;
					} else add=false;
				} else{
					epsilons[idx]=0.0;
					epsilons[idx+1]=0.0;
				}
			}
			if(wsum>EPS){
				tex_eps_avg/=wsum;
				tex_vareps/=wsum;
				tex_vareps-=tex_eps_avg*tex_eps_avg;
				element->MyType->eps_texture[t*count_p+p]=tex_eps_avg/eps_avg;
				if(tex_eps_avg>EPS){ // find maximum and minimum texture value
					if(tex_eps_avg/eps_avg>tex_max) tex_max=tex_eps_avg/eps_avg;
					if(tex_eps_avg/eps_avg<tex_min) tex_min=tex_eps_avg/eps_avg;
				}
			} else element->MyType->eps_texture[t*count_p+p]=0.0;
			texture_average+=element->MyType->eps_texture[t*count_p+p];
			if(configuration->LJ_interaction_area){
				texture_data[t*count_p+p].vec[0]=IA(element->MyType,direction);
				texture_data[t*count_p+p].vec[1]=element->MyType->eps_texture[t*count_p+p];
//				if(wsum>EPS) texture_data[t*count_p+p].vec[2]=1.0/tex_vareps;
				texture_data[t*count_p+p].vec[2]=1.0;
			}
		}
	}
	// take care of poles
	for(unsigned int i=0; i<2; i++){
		wsum=0.0;
		double tex_eps_avg=0.0;
		double tex_vareps=0.0;
		bool add=true;
		for(int k=count_r-1; k>=0; k--){
			unsigned int idx=((count_t*count_p+i)*count_r+(unsigned int)k)*EPS_ELEMENT_COUNT;
			
			double c_eps=epsilons[idx]; // sum of w*eps
			double c_eps2=epsilons[idx+1]; // sum of w*eps*eps
			double ws=epsilons[idx+2]; // sum of w
			double curr_var=(c_eps2-c_eps*c_eps/ws)/ws;
			if(fabs(curr_var)<EPS) curr_var=0.0;
			
			if((ws>EPS) && add){ // sanity checks
				epsilons[idx]/=ws;
				epsilons[idx+1]=sqrt(curr_var);
				double tex_curr_avg=tex_eps_avg+c_eps;
				double tex_curr_var=tex_vareps+c_eps2;
				tex_curr_avg/=(wsum+ws);
				tex_curr_var/=(wsum+ws);
				tex_curr_var-=tex_curr_avg*tex_curr_avg;
				if(4.0*tex_curr_var<tex_curr_avg*tex_curr_avg){ // do not accept inward values if std.dev. is more than half the average
					wsum+=ws;
					tex_eps_avg+=c_eps;
					tex_vareps+=c_eps2;
				} else add=false;
			} else{
				epsilons[idx]=0.0;
				epsilons[idx+1]=0.0;
			}
		}
		if(wsum>EPS){
			tex_eps_avg/=wsum;
			tex_vareps/=wsum;
			tex_vareps-=tex_eps_avg*tex_eps_avg;
			element->MyType->eps_texture[count_t*count_p+i]=tex_eps_avg/eps_avg;
			if(tex_eps_avg>EPS){ // find maximum and minimum texture value
				if(tex_eps_avg/eps_avg>tex_max) tex_max=tex_eps_avg/eps_avg;
				if(tex_eps_avg/eps_avg<tex_min) tex_min=tex_eps_avg/eps_avg;
			}
		} else element->MyType->eps_texture[count_t*count_p+i]=0.0;
		texture_average+=element->MyType->eps_texture[count_t*count_p+i];
		if(configuration->LJ_interaction_area){
			Vec3 direction(0.0,0.0,1.0);
			texture_data[count_t*count_p+i].vec[0]=IA(element->MyType,direction);
			texture_data[count_t*count_p+i].vec[1]=element->MyType->eps_texture[count_t*count_p+i];
//			if(wsum>EPS) texture_data[count_t*count_p+i].vec[2]=1.0/tex_vareps;
			texture_data[count_t*count_p+i].vec[2]=1.0;
		}
	}
	texture_average/=count_t*count_p+2;
	element->MyType->eps_texture[count_t*count_p+2]=tex_min;
	element->MyType->eps_texture[count_t*count_p+3]=tex_max;
	
	cout << "\t-> For a test sphere with radius equal to average interaction width (" << parameters->rp << " Angström)\n\t   the average Lennard-Jones epsilon is: " << eps_avg*eps_avg << " +/- " << sqrt(qqpwr(2.0*eps_avg,2)*vareps) << " (between " << tex_min*tex_min*eps_avg*eps_avg << " and " << tex_max*tex_max*eps_avg*eps_avg << ")\n";
	if(element->MyType->nr_Vvdw_coefficients>0){
		double* coeff_sigmas=new double[8];
		for(unsigned int i=0; i<8; i++) coeff_sigmas[i]=0.1*element->MyType->Vvdw_coefficients[i];
		coeff_sigmas[0]=0.0;
		coeff_sigmas[7]=0.0;
		cout << "\t-> Fitting coefficients to average sqrt(epsilon) with respect to test sphere size.\n\t\t-> ";
		Coefficient_Fit(eps_rT,Nrp,&EPSofRT,element->MyType->Vvdw_coefficients,coeff_sigmas,8);
		for(unsigned int i=0; i<8; i++){
			if(i>0) partner_string+=",";
			partner_string+=double2str(element->MyType->Vvdw_coefficients[i]);
		}
		delete[] coeff_sigmas;
	}
	// Create epsilon data file
	ofstream data;
	string dataname = lod_group->Type->name+"_"+int2str(element_nr+1)+".epsilon.dat";
	data.open(dataname.c_str());
	if(data.fail()){
		cout << "Unable to open output file.\n";
		exit(1);
	}
	if(element->MyType->nr_Vvdw_coefficients>0){
		data << "#coefficients:";
		for(unsigned int i=0; i<8; i++) data << "\t" << element->MyType->Vvdw_coefficients[i];
		data << "\n";
	}
	data << "# sigma_{test}\t<sqrt(eps)>\tstd.dev.\ttexture average\n";
	for(unsigned int i=0; i<NrT; i++) data << eps_rT[i].vec[0] << "\t" << eps_rT[i].vec[1] << "\t" << sqrt(1.0/eps_rT[i].vec[2]) << "\t" << tex_avg_rT[i] << "\n"; // eps_rT contains variance ...
	data.close();
	if(configuration->LJ_interaction_area && configuration->LJ_interaction_area_fit){
		SortData(texture_data,count_t*count_p+2);
		element->MyType->IA_coefficients=new double[7];
		element->MyType->IA_coefficients[0]=1.0;
		element->MyType->IA_coefficients[1]=0.0;
		element->MyType->IA_coefficients[2]=1.0;
		element->MyType->IA_coefficients[3]=0.0;
		element->MyType->IA_coefficients[4]=0.0;
		double* coeff_sigmas=new double[7];
		coeff_sigmas[0]=0.01;
		coeff_sigmas[1]=0.05;
		coeff_sigmas[2]=0.05;
		coeff_sigmas[3]=0.01;
		coeff_sigmas[4]=0.01;
		coeff_sigmas[5]=0.0;
		coeff_sigmas[6]=0.0;
		cout << "\t-> Fitting IA coefficients to average texture values (may take a while).\n\t\t-> ";
		Coefficient_Fit(texture_data,count_t*count_p+2,&IAadjust,element->MyType->IA_coefficients,coeff_sigmas,7);
		delete[] texture_data;
		// Create output data file
		dataname = lod_group->Type->name+"_"+int2str(element_nr+1)+".IA_vs_texture.dat";
		data.open(dataname.c_str());
		if(data.fail()){
			cout << "Unable to open output file.\n";
			exit(1);
		}
		data << "#coefficients:";
		for(unsigned int i=0; i<7; i++) data << "\t" << element->MyType->IA_coefficients[i];
		data << "\n";
		data << "#cos(theta)\tphi\tIA\ttexture\n";
		double IAadj_average=0.0;
		for(unsigned int ct=0; ct<count_t; ct++){
			double cost=1.0-2.0*(ct+0.5)/count_t;
			double sint=sqrt(1.0-cost*cost);
			for(unsigned int j=0; j<count_p; j++){
				double phi=2.0*pi*(j+0.5)/count_p;
				Vec3 direction(sint*cos(phi),sint*sin(phi),cost);
				data << cost << "\t" << phi << "\t" << Ellipsoid_Cross_Section(element->MyType->invsaxes2,direction)*element->MyType->inv_avg_area << "\t" << double2str(element->MyType->eps_texture[ct*count_p+j]) << "\n";
				element->MyType->eps_texture[ct*count_p+j]=IA(element->MyType,direction);
				IAadj_average+=element->MyType->eps_texture[ct*count_p+j];
			}
		}
		element->MyType->eps_texture[count_t*count_p]=IA(element->MyType,Vec3(0.0,0.0,1.0));
		element->MyType->eps_texture[count_t*count_p+1]=IA(element->MyType,Vec3(0.0,0.0,1.0));
		IAadj_average+=2*element->MyType->eps_texture[count_t*count_p];
		IAadj_average/=count_t*count_p+2;
		cout << IAadj_average << "\n";
		data.close();
	}
	
	if(output_x3d){
		cout << "-> Generating Lennard-Jones potential X3D file ... ";
		// Create file
		ofstream vizout;
		string vizname = lod_group->Type->name+"_"+int2str(element_nr+1)+".x3d";
		vizout.open(vizname.c_str());
		if (vizout.fail()){
			cout << "Unable to open output file.\n";
			exit(1);
		}
		// Create model file
		ofstream model;
		string modelname = lod_group->Type->name+"_"+int2str(element_nr+1)+".model.x3d";
		model.open(modelname.c_str());
		if(model.fail()){
			cout << "Unable to open output file.\n";
			exit(1);
		}
		
		// Configure viewpoint -- want to have left-handed coordinate system (x to the right, y to the back, z up)
		viewpoint_t viewpoint;
		viewpoint.Center = Vec3(0,0,0);
		double lookdownangle=-pi/7.854;
		double viewradius=parameters->rmax/3.0;
		Mat33 viewrot=AxisAngle2Rot(Vec4(0,0,1,pi*0.25))*AxisAngle2Rot(Vec4(1,0,0,pi*0.5+lookdownangle));
		viewpoint.Orientation = Rot2AxisAngle(viewrot);
		Vec3 viewpos(0.0,0.0,viewradius);
		viewpoint.Position = viewrot*viewpos;
		viewpoint.Description = "Default viewpoint";
		
		// Output file header info
		output_header(vizout, vizname, viewpoint);
		output_header(model, modelname, viewpoint);
		
		// Load and define element types
		vizout << "\t\t<!-- Define element types -->\n";
		model << "\t\t<!-- Define element types -->\n";
		
		// define LOD ellipsoid archetype (need not forget transparency)
		define_element(vizout,element->MyType,lod_group,false,0,configuration);
		define_element(model,element->MyType,lod_group,false,0,configuration);
		vizout << "\t\t<!-- done -->\n";
		model << "\t\t<!-- done -->\n";
		
		vizout <<  "\t\t<Transform>\n";
		vizout <<  "\t\t\t<ProtoInstance name='" << element->MyType->name << " (" << lod_group->Type->name << ")' DEF='" << lod_group->Type->name << "_" << element_nr+1 << "'/>\n";
		vizout <<  "\t\t</Transform>\n";
		model <<  "\t\t<Transform>\n";
		model <<  "\t\t\t<ProtoInstance name='" << element->MyType->name << " (" << lod_group->Type->name << ")' DEF='" << lod_group->Type->name << "_" << element_nr+1 << "'/>\n";
		model <<  "\t\t</Transform>\n";
		model << "\t</Scene>\n";
		model << "</X3D>\n";
		model.close();
		
		//Configure footer - replace this with zvis read and defaults
		axis_t axis;
		axis.DrawAxes = true; //draw axis markers
		
		if(element->MyType->nr_Vvdw_coefficients>0) eps_avg=EPSofRT(element->MyType->rT,element->MyType->Vvdw_coefficients,element->MyType->nr_Vvdw_coefficients); // *sqrt(element->MyType->sphericity);
		DrawLJPotentials(vizout,configuration,lod_nr,lod_level,element_nr,parameters,eps_avg,axis);
		
		axis.ColorX = Vec3(1,0,0); //color of x axis marker
		axis.ColorY = Vec3(0.1,0.5,1); //color of y axis marker
		axis.ColorZ = Vec3(0,1,0); //you've probably already guessed what this is (green? -- AT)
		
		// Write footer
		output_footer(vizout,axis);
		vizout.close();
		
		cout << "Done.\n";
		cout << "-> Generating Lennard-Jones interaction with self X3D file ... ";
		// Create file
		vizname = lod_group->Type->name+"_"+int2str(element_nr+1)+"_interaction.x3d";
		vizout.open(vizname.c_str());
		if (vizout.fail()){
			cout << "Unable to open output file.\n";
			exit(1);
		}
		
		// Configure viewpoint -- want to have left-handed coordinate system (x to the right, y to the back, z up)
		viewpoint.Center = Vec3(0,0,0);
		lookdownangle=-pi/7.854;
		viewradius=parameters->rmax/2.5;
		viewrot=AxisAngle2Rot(Vec4(0,0,1,pi*0.25))*AxisAngle2Rot(Vec4(1,0,0,pi*0.5+lookdownangle));
		viewpoint.Orientation = Rot2AxisAngle(viewrot);
		viewpos=Vec3(0.0,0.0,viewradius);
		viewpoint.Position = viewrot*viewpos;
		viewpoint.Description = "Default viewpoint";
		// Output file header info
		output_header(vizout, vizname, viewpoint);
		// Load and define element types
		vizout << "\t\t<!-- Define element types -->\n";
		
		// define LOD ellipsoid archetype (need not forget transparency)
		define_element(vizout,element->MyType,lod_group,false,0,configuration,true);
		vizout << "\t\t<!-- done -->\n";
		
		vizout <<  "\t\t<Transform>\n";
		vizout <<  "\t\t\t<ProtoInstance name='" << element->MyType->name << " (" << lod_group->Type->name << ")' DEF='" << lod_group->Type->name << "_" << element_nr+1 << "'/>\n";
		vizout <<  "\t\t</Transform>\n";
		
		DrawLJinteraction(vizout,configuration,lod_nr,lod_level,element_nr,parameters,eps_avg,axis);
		
		// Write footer
		output_footer(vizout,axis);
		vizout.close();
		
		cout << "Done.\n";
	}
	
	element->MyType->Vvdw=eps_avg*eps_avg; // we're outputting epsilon - not its square root
	// output texture string
	for(unsigned int i=0; i<count_t*count_p+4; i++){
		if(i>0) texture_string+=",";
		texture_string+=double2str(element->MyType->eps_texture[i]);
	}
	
	// clean up
	clReleaseKernel(kernel);
	
#ifndef OpenCL
	delete Error;
#endif
	delete parameters;
	delete[] element_dyn;
	delete[] element_stat;
	delete[] epsilons;
	delete[] eps_rT;
	
	clReleaseMemObject(dynamic);
	clReleaseMemObject(stat);
	clReleaseMemObject(output);
	
	clReleaseProgram(program);
	clReleaseCommandQueue(Queue);
	clReleaseContext(Context);
}

int main(int argc, char* argv[])
{
	cout << "\nRobinson Group Level-of-Detail Fitting Tool\n";
	cout << "Compiled " << __DATE__ << " (Build " << BUILD << ", " << VERSION << ")\n";
	
	MC_Config config;
	Config_Data* configuration = &config.parameters;
	configuration->fit2lod=true;
	
	string conffile="";
	// Check for command line parameters
	if(argc>1){ // yes, there are some -- only parameter accepted is configuration filename
		conffile=argv[1]; // easy enough
	} else{
		cout << "Syntax: " << argv[0] << " <configuration file> <optional: group name>\n\n";
		exit(1);
	}
	
	// load configuration from file
#ifndef USE_CMWC4096
	configuration->idum = new __int32_t; // so we don't get segfaults ...
	*configuration->idum = -8;
#endif
	config.GetFromFile(conffile.c_str());
	phys_configuration(configuration);
	
	cout << "\n";
	if(configuration->num_levelofdetail==0){
		cout << "ERROR: Specified configuation file does not contain any LOD groups.\n";
		exit(1);
	} else{
		cout << "-> Found " << configuration->num_levelofdetail << " Level-of-Detail group";
		if(configuration->num_levelofdetail>1) cout << "s";
		cout << ".\n";
	}
	
	int group_nr=-1;
	if(argc>2){ // second parameter: particular group
		for(unsigned int i=0; i<configuration->num_levelofdetail; i++){
			if(compare_strings(argv[2],configuration->lods[i]->groups[0]->Type->name.c_str())){
				group_nr=(int)i;
				break;
			}
		}
		if(group_nr<0){
			cout << "ERROR: Group <" << argv[2] << "> is nowhere to be found ...\n";
			exit(2);
		} else{
			cout << "\t-> Group <" << configuration->lods[(unsigned int)group_nr]->groups[0]->Type->name << "> found\n";
		}
	}
	cout << "\n";
	
	double cr = (double)(configuration->LJexp[1])/(double)(configuration->LJexp[0]);
	configuration->r=pow(cr,cr/(cr-1.0))/(cr-1.0);
	
	string storage_filename=configuration->fileout+".conf.dat";
	fstream output;
	output.open(storage_filename.c_str(), fstream::out | fstream::trunc);
	if(output.fail()){
		cout << "Unable to open configuration output file.\n";
		exit(1);
	}
	output << "// This is an automatically created file.\n// Please do not edit this file manually.\n// -- Use \"fit2lod <configuration file>\"\n\n";
	
	unsigned int i_s=0;
	if(group_nr>=0) i_s=(unsigned int)group_nr;
	for(unsigned int i=i_s; (i<configuration->num_levelofdetail) && ((group_nr<0) || ((group_nr>=0) && (i==(unsigned int)group_nr))); i++){
		Element_Group* group=configuration->lods[i]->groups[0];
		bool output_detail=true;
			
		for(unsigned int j=0; j<configuration->lods[i]->levels; j++){
			Element_Group* newgroup=configuration->lods[i]->groups[j+1];
			if(newgroup->number>0){
				if(output_detail){
					output << "[Detail: " << group->Type->name << "]\n";
					output_detail=false; // header needs to be output only once
				}
				cout << "Setting up <" << newgroup->Type->name << ">:\n";
				string model_filename=newgroup->Type->name+".model.conf";
				fstream model;
				model.open(model_filename.c_str(), fstream::out | fstream::trunc);
				if(model.fail()){
					cout << "Unable to open model output file.\n";
					exit(1);
				}
				model << "// This file was automatically created by fit2lod.\n";
				if(newgroup->Type->LOD->match_epsilon){
					cout << "-> Adjusting Lennard-Jones epsilons to match fully atomistic potential energies.\n";
					cout.flush();
				}
				string element_list="";
				string connectivity="";
				string element_positions="";
				string element_rotations="";
				string connection_points="";
				// Output what we already know
				string volume_string="lod_volumes."+int2str(j+1)+" = {";
				string epsilon_string="lod_epsilons."+int2str(j+1)+" = {";
				string IA_string="lod_IA."+int2str(j+1)+" = {";
				bool IA_adj_used=false;
				string texture_string="lod_textures."+int2str(j+1)+" = {";
				string partner_string="lod_epsilon_of_rT."+int2str(j+1)+" = {";
				for(unsigned int k=0; k<newgroup->nr_elements; k++){
					Element* element=&configuration->group_elements[newgroup->elements[k]];
					cout << "\n-> Calculating properties for element <" << element->MyType->name << ">\n";
					model << "\n[Element: " << element->MyType->name << "]\n";
					element_positions+="position."+int2str(k+1)+" = ("+element->center.V3Str(',')+")\n";
					element_rotations+="rotation."+int2str(k+1)+" = ("+Rot2AxisAngle(element->rot).V4Str(',')+")\n";
					bool comma=false;
					for(unsigned int l=0; l<element->nr_interactions; l++){
						unsigned int partner=element->interactions[l].partner;
						if(partner>k){
							if(comma) connectivity+=","; else comma=true;
							connectivity+=int2str(partner+1);
						}
						if(element->interactions[l].location) connection_points+="connection_point."+int2str(k+1)+"-"+int2str(partner+1)+" = ("+element->interactions[l].initial_location->V3Str(',')+")\n";
					}
					if(k+1<newgroup->nr_elements) connectivity+="|";
					if(k>0) element_list+=", ";
					element_list+=element->MyType->name;
					if(k>0){
						volume_string+=", ";
						epsilon_string+=", ";
						IA_string+="|";
						texture_string+="|";
						partner_string+="|";
					}
					if(k<newgroup->Type->LOD->nr_ellipsoids[j]){
						volume_string+=double2str(element->MyType->saxes.vec[0]*element->MyType->saxes.vec[1]*element->MyType->saxes.vec[2]*4.0/3.0*pi);
						if(newgroup->Type->LOD->match_epsilon) epsilon_matching(&config,configuration,i,j,k,true,texture_string,partner_string); // epsilon matching happens here
						epsilon_string+=double2str(element->MyType->Vvdw);
						if(element->MyType->IA_coefficients){
							IA_adj_used=true;
							for(unsigned int l=0; l<7; l++){
								if(l>0) IA_string+=",";
								IA_string+=double2str(element->MyType->IA_coefficients[l]);
							}
						}
					}
					model << "r0 = (" << element->MyType->saxes.vec[0] << ", " << element->MyType->saxes.vec[1] << ", " << element->MyType->saxes.vec[2] << ")\n";
					model << "vdw = " << element->MyType->Vvdw << "\n";
					if(element->MyType->avg_width>EPS) model << "avg_width = " << element->MyType->avg_width << "\n";
					if(element->MyType->mass>EPS) model << "mass = " << element->MyType->mass << "\n";
					if(element->MyType->hasmu) model << "dipole = (" << element->MyType->initial_dipole.V3Str(',') << ")\n";
					if(element->MyType->nr_charges>0){
						model << "charges = {";
						string locations="";
						for(unsigned int l=0; l<element->MyType->nr_charges; l++){
							if(l>0) model << ",";
							model << element->MyType->q[l]/e_in_esu;
							if(element->MyType->q_pos[l]!=Vec3(0.0)) locations+="charge_location."+int2str(l+1)+" = ("+element->MyType->q_pos[l].V3Str(',')+")\n";
						}
						model << "}\n";
						if(locations!="") model << locations;
					}
				}
				volume_string+="}\n";
				epsilon_string+="}\n";
				texture_string+="}\n";
				partner_string+="}\n";
				
				output << volume_string;
				output << epsilon_string;
				if(IA_adj_used) output << IA_string << "}\n";
				output << partner_string;
				output << texture_string;
				
				model << "\n[Group: " << newgroup->Type->name << "]\n";
				model << "elements = " << element_list << "\n";
				if(connectivity!="") model << "connectivity = {" << connectivity << "}\n";
				model << "\n" << element_positions << element_rotations;
				if(connection_points!="") model << connection_points;
				model.close();
				cout << "<- Done.\n\n";
				output << "\n";
			}
		}
	}
	output << "\n";
	output.close();
	
	return 0;
}

