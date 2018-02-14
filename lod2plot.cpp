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

inline double VLJ_LOD(Config_Data* configuration, double rmin, double wA, double wB, double r)
{
	double VLJ;
	double r_over_r=rmin/r;
	if(configuration->LJ_adjust_width){
		double dist=r-(rmin-(wA+wB));
		if(dist<0.0) dist=EPS;
		r_over_r=(wA+wB)/dist;
	}
	switch(configuration->vdwtype){
		case 1:
		case 4:
		default:
			r_over_r = qqpwr(r_over_r,configuration->LJexp[0]);
			VLJ=configuration->r*r_over_r*(r_over_r-1.0);
			break;
		case 2:
			// LJ + Bruce correction around r=rmin. Tanh functions have been replaced with a Gaussian for speed
			VLJ=configuration->r*(qqpwr(r_over_r,configuration->LJexp[1])-qqpwr(r_over_r,configuration->LJexp[0])) + configuration->Solvent[1]*exp(configuration->Solvent[0]*(r - configuration->Solvent[2])*(r - configuration->Solvent[2]));
			break;
		case 3:
		case 5:
			VLJ=configuration->r*qqpwr(r_over_r,configuration->LJexp[0]<<1); // bitshift left by one is multiplication by 2
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

inline void DrawLJPotentials(ofstream &vizout, Config_Data* configuration, unsigned int lod_nr, unsigned int lod_level, unsigned int element_nr, double eps_avg, axis_t &axis)
{
/*	Element_Group* lod_group=configuration->lods[lod_nr]->groups[lod_level+1];
	Element* element=&configuration->group_elements[lod_group->elements[element_nr]];
	double rmin_draw=mind(element->MyType->saxes.vec,3)/2.0;
	
	double VLJamplitude=maxd(element->MyType->saxes.vec,3)*0.5;
	
	// x-Axis
	unsigned int t=theta_res/2; // is cos(pi/2)
	unsigned int p=0;
	double rmin = element->MyType->saxes.vec[0]+element->MyType->rT;
	
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
	DrawCalcSphere(vizout,parameters->rp,e_x*axis.length.vec[2]);*/
}

inline void DrawLJinteraction(ofstream &vizout, Config_Data* configuration, unsigned int lod_nr, unsigned int lod_level, unsigned int element_nr, double eps_avg, axis_t &axis)
{
/*	Element_Group* group=configuration->lods[lod_nr]->groups[0];
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
	
	cout << ia.vec[0]*ia.vec[0] << ", " << ia.vec[2]*ia.vec[2] << "\n";
	
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
	DrawSelfModelReduced(vizout,element->MyType->saxes,e_x*pos,rot_x,Vec3(1.0,1.0,0.0)); // x-face*/
}


int main(int argc, char* argv[])
{
	cout << "\nRobinson Group Level-of-Detail Data Plotting Tool\n";
	cout << "Compiled " << __DATE__ << " (Build " << BUILD << ", " << VERSION << ")\n";
	
	MC_Config config;
	Config_Data* configuration = &config.parameters;
	configuration->fit2lod=true;
	
	string conffile="";
	// Check for command line parameters
	if(argc>2){ // yes, there are some
		conffile=argv[1]; // easy enough
	} else{
		cout << "Syntax: " << argv[0] << " <configuration file> <LOD group name>\n\n";
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
	
	Element_Group* lod_group=configuration->lods[(unsigned int)group_nr]->groups[1];
	
	double cr = (double)(configuration->LJexp[1])/(double)(configuration->LJexp[0]);
	configuration->r=pow(cr,cr/(cr-1.0))/(cr-1.0);
	
	cout << "-> Generating LOD model X3D file ... ";
	// Create model file
	ofstream model;
	string modelname = lod_group->Type->name+".model.x3d";
	model.open(modelname.c_str());
	if(model.fail()){
		cout << "Unable to open output file.\n";
		exit(1);
	}
	
	// Configure viewpoint -- want to have left-handed coordinate system (x to the right, y to the back, z up)
	viewpoint_t viewpoint;
	viewpoint.Center = Vec3(0,0,0);
	double lookdownangle=-pi/7.854;
	double viewradius=configuration->group_radii[lod_group->type];
	Mat33 viewrot=AxisAngle2Rot(Vec4(0,0,1,pi*0.25))*AxisAngle2Rot(Vec4(1,0,0,pi*0.5+lookdownangle));
	viewpoint.Orientation = Rot2AxisAngle(viewrot);
	Vec3 viewpos(0.0,0.0,viewradius);
	viewpoint.Position = viewrot*viewpos;
	viewpoint.Description = "Default viewpoint";
	
	// Output file header info
	output_header(model, modelname, viewpoint);
	for(unsigned int element_nr=0; element_nr<lod_group->nr_elements; element_nr++){
		Element* element=&configuration->group_elements[lod_group->elements[element_nr]];
		double transparency=element->MyType->transparency;
		if(lod_group->Type->transparency>=0.0) transparency=lod_group->Type->transparency;
		if(lod_group->Type->LOD && (lod_group->levelofdetail==0)){
			if((lod_group->Type->LOD->nr_components>0) && (lod_group->Type->LOD->element_in_component[element_nr]>=0)){
				if(lod_group->Type->LOD->component_transparency[lod_group->Type->LOD->element_in_component[element_nr]]>=0.0){
					transparency=lod_group->Type->LOD->component_transparency[lod_group->Type->LOD->element_in_component[element_nr]];
				}
			}
		}
		define_element(model,element->MyType,lod_group,false,0,configuration,true,false,element->center,Rot2AxisAngle(element->rot),transparency);
	}
	model << "\t</Scene>\n";
	model << "</X3D>\n";
	model.close();
	
	cout << "<- Done.\n\n";
	return 0;
}

