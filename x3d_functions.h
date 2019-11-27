/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

// Viewpoint info for header
typedef struct _viewpoint_t {
	Vec3 Position;
	Vec3 Center;
	Vec4 Orientation;
	string Description;
} viewpoint_t;

// Axis info for footer
typedef struct _axis_t {
	bool DrawAxes;
	Vec3 ColorX;
	Vec3 ColorZ;
	Vec3 ColorY;
	Vec3 length;
} axis_t;

/// Output header info to x3d file
inline void output_header(ofstream &vizout, string &vizname, viewpoint_t &viewpoint, bool animate, unsigned int steps)
{
	vizout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	vizout << "<!DOCTYPE X3D PUBLIC \"ISO//Web3D//DTD X3D 3.1//EN\" \"http://www.web3d.org/specifications/x3d-3.1.dtd\">\n";
	vizout << "<X3D profile='Interchange' version='3.1' xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' xsd:noNamespaceSchemaLocation=' http://www.web3d.org/specifications/x3d-3.1.xsd '>\n";
	vizout << "\t<head>\n";
	vizout << "\t\t<meta name='title' content='" << vizname << "'/>\n";
	vizout << "\t\t<meta name='description' content='MC simulation visualization'/>\n";
	vizout << "\t\t<meta name='generator' content='Robinson group'/>\n";
	vizout << "\t\t<meta name='generator' content='Department of Chemistry, University of Washington'/>\n";
	vizout << "\t</head>\n";
	vizout << "\t<Scene>\n";
	vizout << "\t\t<Viewpoint position='" << viewpoint.Position.V3Str(' ') << "' centerOfRotation='" << viewpoint.Center.V3Str(' ') << "' ";
	vizout << "orientation='" << viewpoint.Orientation.V4Str(' ') << "' " << "fieldOfView='0.7854' ";
	vizout << "description='" << viewpoint.Description << "'/>\n";
	vizout << "\t\t<NavigationInfo headlight='true' blendingSort='3D'/>\n";
	vizout << "\t\t<!-- Define line -->\n";
	vizout << "\t\t<ProtoDeclare name='line'>\n";
	vizout << "\t\t\t<ProtoInterface>\n";
	vizout << "\t\t\t\t<field accessType='inputOutput' name='AB' type='MFVec3f' value='0 0 0 1 0 0'/>\n";
	vizout << "\t\t\t\t<field accessType='inputOutput' name='position' type='SFVec3f' value='0 0 0'/>\n";
	vizout << "\t\t\t\t<field accessType='inputOutput' name='rotation' type='SFRotation' value='0 0 0 0'/>\n";
	vizout << "\t\t\t\t<field accessType='inputOutput' name='color' type='SFColor' value='1 1 1'/>\n";
	vizout << "\t\t\t\t<field accessType='inputOutput' name='transparency' type='SFFloat' value='0.2'/>\n";
	vizout << "\t\t\t</ProtoInterface>\n";
	vizout << "\t\t\t<ProtoBody>\n";
	vizout << "\t\t\t\t<Group>\n";
	vizout << "\t\t\t\t\t<Transform>\n";
	vizout << "\t\t\t\t\t\t<IS>\n";
	vizout << "\t\t\t\t\t\t\t<connect nodeField='translation' protoField='position'/>\n";
	vizout << "\t\t\t\t\t\t\t<connect nodeField='rotation' protoField='rotation'/>\n";
	vizout << "\t\t\t\t\t\t</IS>\n";
	vizout << "\t\t\t\t\t\t<Shape>\n";
	vizout << "\t\t\t\t\t\t\t<Appearance>\n";
	vizout << "\t\t\t\t\t\t\t\t<Material shininess='0.8' specularColor='0.8 0.8 0.8'>\n";
	vizout << "\t\t\t\t\t\t\t\t\t<IS>\n";
	vizout << "\t\t\t\t\t\t\t\t\t\t<connect nodeField='diffuseColor' protoField='color'/>\n";
	vizout << "\t\t\t\t\t\t\t\t\t\t<connect nodeField='emissiveColor' protoField='color'/>\n";
	vizout << "\t\t\t\t\t\t\t\t\t\t<connect nodeField='transparency' protoField='transparency'/>\n";
	vizout << "\t\t\t\t\t\t\t\t\t</IS>\n";
	vizout << "\t\t\t\t\t\t\t\t</Material>\n";
	vizout << "\t\t\t\t\t\t\t</Appearance>\n";
	vizout << "\t\t\t\t\t\t\t<IndexedLineSet coordIndex='0 1'>\n";
	vizout << "\t\t\t\t\t\t\t\t<Coordinate>\n";
	vizout << "\t\t\t\t\t\t\t\t\t<IS>\n";
	vizout << "\t\t\t\t\t\t\t\t\t\t<connect nodeField='point' protoField='AB'/>\n";
	vizout << "\t\t\t\t\t\t\t\t\t</IS>\n";
	vizout << "\t\t\t\t\t\t\t\t</Coordinate>\n";
	vizout << "\t\t\t\t\t\t\t</IndexedLineSet>\n";
	vizout << "\t\t\t\t\t\t</Shape>\n";
	vizout << "\t\t\t\t\t</Transform>\n";
	vizout << "\t\t\t\t</Group>\n";
	vizout << "\t\t\t</ProtoBody>\n";
	vizout << "\t\t</ProtoDeclare>\n";
	vizout << "\t\t<!-- Define 3D vector -->\n";
	vizout << "\t\t<ProtoDeclare name='Vector3D'>\n";
	vizout << "\t\t\t<ProtoInterface>\n";
	vizout << "\t\t\t\t<field accessType='inputOutput' name='position' type='SFVec3f' value='0 0 0'/>\n";
	vizout << "\t\t\t\t<field accessType='inputOutput' name='length' type='SFFloat' value='1.0'/>\n";
	vizout << "\t\t\t\t<field accessType='inputOutput' name='tip_vector' type='SFVec3f' value='0 0.5 0'/>\n";
	vizout << "\t\t\t\t<field accessType='inputOutput' name='radius' type='SFFloat' value='0.05'/>\n";
	vizout << "\t\t\t\t<field accessType='inputOutput' name='tip_width' type='SFFloat' value='0.2'/>\n";
	vizout << "\t\t\t\t<field accessType='inputOutput' name='tip_height' type='SFFloat' value='0.4'/>\n";
	vizout << "\t\t\t\t<field accessType='inputOutput' name='rotation_from_y' type='SFRotation' value='0 0 0 0'/>\n";
	vizout << "\t\t\t\t<field accessType='inputOutput' name='color' type='SFColor' value='1 1 1'/>\n";
	vizout << "\t\t\t\t<field accessType='inputOutput' name='transparency' type='SFFloat' value='0'/>\n";
	if(animate){
		vizout << "\t\t\t\t<field accessType='inputOutput' name='key' type='MFFloat'/>\n";
		vizout << "\t\t\t\t<field accessType='inputOutput' name='position_values' type='MFVec3f'/>\n";
		vizout << "\t\t\t\t<field accessType='inputOutput' name='rotation_values' type='MFRotation'/>\n";
		vizout << "\t\t\t\t<field accessType='inputOutput' name='length_values' type='MFVec3f'/>\n";
	}
	vizout << "\t\t\t</ProtoInterface>\n";
	vizout << "\t\t\t<ProtoBody>\n";
	vizout << "\t\t\t\t<Group>\n";
	if(animate) vizout << "\t\t\t\t\t<Transform DEF='overall_transformation'>\n"; else vizout << "\t\t\t\t\t<Transform>\n";
	vizout << "\t\t\t\t\t\t<IS>\n";
	vizout << "\t\t\t\t\t\t\t<connect nodeField='translation' protoField='position'/>\n";
	vizout << "\t\t\t\t\t\t\t<connect nodeField='rotation' protoField='rotation_from_y'/>\n";
	vizout << "\t\t\t\t\t\t</IS>\n";
	vizout << "\t\t\t\t\t\t<Shape>\n";
	vizout << "\t\t\t\t\t\t\t<Appearance>\n";
	vizout << "\t\t\t\t\t\t\t\t<Material shininess='0.8' specularColor='0.8 0.8 0.8'>\n";
	vizout << "\t\t\t\t\t\t\t\t\t<IS>\n";
	vizout << "\t\t\t\t\t\t\t\t\t\t<connect nodeField='diffuseColor' protoField='color'/>\n";
	vizout << "\t\t\t\t\t\t\t\t\t\t<connect nodeField='specularColor' protoField='color'/>\n";
	vizout << "\t\t\t\t\t\t\t\t\t\t<connect nodeField='transparency' protoField='transparency'/>\n";
	vizout << "\t\t\t\t\t\t\t\t\t</IS>\n";
	vizout << "\t\t\t\t\t\t\t\t</Material>\n";
	vizout << "\t\t\t\t\t\t\t</Appearance>\n";
	vizout << "\t\t\t\t\t\t\t<Cylinder>\n";
	vizout << "\t\t\t\t\t\t\t\t<IS>\n";
	vizout << "\t\t\t\t\t\t\t\t\t<connect nodeField='radius' protoField='radius'/>\n";
	vizout << "\t\t\t\t\t\t\t\t\t<connect nodeField='height' protoField='length'/>\n";
	vizout << "\t\t\t\t\t\t\t\t</IS>\n";
	vizout << "\t\t\t\t\t\t\t</Cylinder>\n";
	vizout << "\t\t\t\t\t\t</Shape>\n";
	vizout << "\t\t\t\t\t\t<Transform>\n";
	vizout << "\t\t\t\t\t\t\t<IS>\n";
	vizout << "\t\t\t\t\t\t\t\t<connect nodeField='translation' protoField='tip_vector'/>\n";
	vizout << "\t\t\t\t\t\t\t</IS>\n";
	vizout << "\t\t\t\t\t\t\t<Shape>\n";
	vizout << "\t\t\t\t\t\t\t\t<Appearance>\n";
	vizout << "\t\t\t\t\t\t\t\t\t<Material shininess='0.8' specularColor='0.8 0.8 0.8'>\n";
	vizout << "\t\t\t\t\t\t\t\t\t\t<IS>\n";
	vizout << "\t\t\t\t\t\t\t\t\t\t\t<connect nodeField='diffuseColor' protoField='color'/>\n";
	vizout << "\t\t\t\t\t\t\t\t\t\t\t<connect nodeField='specularColor' protoField='color'/>\n";
	vizout << "\t\t\t\t\t\t\t\t\t\t\t<connect nodeField='transparency' protoField='transparency'/>\n";
	vizout << "\t\t\t\t\t\t\t\t\t\t</IS>\n";
	vizout << "\t\t\t\t\t\t\t\t\t</Material>\n";
	vizout << "\t\t\t\t\t\t\t\t</Appearance>\n";
	vizout << "\t\t\t\t\t\t\t\t<Cone>\n";
	vizout << "\t\t\t\t\t\t\t\t\t<IS>\n";
	vizout << "\t\t\t\t\t\t\t\t\t\t<connect nodeField='bottomRadius' protoField='tip_width'/>\n";
	vizout << "\t\t\t\t\t\t\t\t\t\t<connect nodeField='height' protoField='tip_height'/>\n";
	vizout << "\t\t\t\t\t\t\t\t\t</IS>\n";
	vizout << "\t\t\t\t\t\t\t\t</Cone>\n";
	vizout << "\t\t\t\t\t\t\t</Shape>\n";
	vizout << "\t\t\t\t\t\t</Transform>\n";
	vizout << "\t\t\t\t\t</Transform>\n";
	vizout << "\t\t\t\t</Group>\n";
	if(animate){
		vizout << "\t\t\t\t<TimeSensor DEF='SimGrStep' cycleInterval='" << steps << "' loop='true'/>\n";
		vizout << "\t\t\t\t<PositionInterpolator DEF='positions'>\n";
		vizout << "\t\t\t\t\t<IS>\n";
		vizout << "\t\t\t\t\t\t<connect nodeField='key' protoField='key'/>\n";
		vizout << "\t\t\t\t\t\t<connect nodeField='keyValue' protoField='position_values'/>\n";
		vizout << "\t\t\t\t\t</IS>\n";
		vizout << "\t\t\t\t</PositionInterpolator>\n";
		vizout << "\t\t\t\t<OrientationInterpolator DEF='rotations'>\n";
		vizout << "\t\t\t\t\t<IS>\n";
		vizout << "\t\t\t\t\t\t<connect nodeField='key' protoField='key'/>\n";
		vizout << "\t\t\t\t\t\t<connect nodeField='keyValue' protoField='rotation_values'/>\n";
		vizout << "\t\t\t\t\t</IS>\n";
		vizout << "\t\t\t\t</OrientationInterpolator>\n";
		vizout << "\t\t\t\t<PositionInterpolator DEF='lengths'>\n";
		vizout << "\t\t\t\t\t<IS>\n";
		vizout << "\t\t\t\t\t\t<connect nodeField='key' protoField='key'/>\n";
		vizout << "\t\t\t\t\t\t<connect nodeField='keyValue' protoField='length_values'/>\n";
		vizout << "\t\t\t\t\t</IS>\n";
		vizout << "\t\t\t\t</PositionInterpolator>\n";
		vizout << "\t\t\t\t<ROUTE fromNode='SimGrStep' fromField='fraction_changed' toNode='positions' toField='set_fraction'/>\n";
		vizout << "\t\t\t\t<ROUTE fromNode='SimGrStep' fromField='fraction_changed' toNode='rotations' toField='set_fraction'/>\n";
		vizout << "\t\t\t\t<ROUTE fromNode='SimGrStep' fromField='fraction_changed' toNode='lengths' toField='set_fraction'/>\n";
		vizout << "\t\t\t\t<ROUTE fromNode='positions' fromField='value_changed' toNode='overall_transformation' toField='translation'/>\n";
		vizout << "\t\t\t\t<ROUTE fromNode='rotations' fromField='value_changed' toNode='overall_transformation' toField='rotation'/>\n";
		vizout << "\t\t\t\t<ROUTE fromNode='lengths' fromField='value_changed' toNode='overall_transformation' toField='scale'/>\n";
	}
	vizout << "\t\t\t</ProtoBody>\n";
	vizout << "\t\t</ProtoDeclare>\n";
}

inline void output_header(ofstream &vizout, string &vizname, viewpoint_t &viewpoint){
	output_header(vizout,vizname,viewpoint,false,0);
}

inline string Byte2Hex(unsigned short b)
{
	string result="";
	char hex[16]={'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};
	
	result+=hex[(b>>4)&0xF];
	result+=hex[b&0xF];
	
	return result;
}

#define atan_step_max 4.0
#define inv_atan_max 1.0/atan(atan_step_max)

#define min_color {255,0,255}
#define avg_color {0,0,255}
#define max_color {0,255,0}

inline string ColorMapping(double value, double min, double max, double avg)
{
	string result="0x";
	
	unsigned short ac[3]=avg_color;
	
	unsigned short r,g,b;
	if((value<min) || (value>max)){
		r=127;
		g=127;
		b=127;
	} else{
		if(value<avg){
			unsigned short mc[3]=min_color;
			double frac=atan(atan_step_max*(avg-value)/(avg-min))*inv_atan_max;
			if((fabs(avg-min)<EPS) || (frac>1.0)) frac=1.0;
			
			r=(unsigned short)((1.0-frac)*ac[0]+frac*mc[0]);
			g=(unsigned short)((1.0-frac)*ac[1]+frac*mc[1]);
			b=(unsigned short)((1.0-frac)*ac[2]+frac*mc[2]);
		} else{
			unsigned short mc[3]=max_color;
			double frac=atan(atan_step_max*(value-avg)/(max-avg))*inv_atan_max;
			if((fabs(max-avg)<EPS) || (frac>1.0)) frac=1.0;
			
			r=(unsigned short)((1.0-frac)*ac[0]+frac*mc[0]);
			g=(unsigned short)((1.0-frac)*ac[1]+frac*mc[1]);
			b=(unsigned short)((1.0-frac)*ac[2]+frac*mc[2]);
		}
	}
	
	result+=Byte2Hex(r)+Byte2Hex(g)+Byte2Hex(b);
	
	return result;
}

inline Vec3 ColorMappingVec(double value, double min, double max, double avg)
{
	unsigned short ac[3]=avg_color;
	
	double r,g,b;
	if((value<min) || (value>max)){
		r=0.5;
		g=0.5;
		b=0.5;
	} else{
		if(value<avg){
			unsigned short mc[3]=min_color;
			double frac=atan(atan_step_max*(avg-value)/(avg-min))*inv_atan_max;
			if((fabs(avg-min)<EPS) || (frac>1.0)) frac=1.0;
			
			r=((1.0-frac)*ac[0]+frac*mc[0])/255.0;
			g=((1.0-frac)*ac[1]+frac*mc[1])/255.0;
			b=((1.0-frac)*ac[2]+frac*mc[2])/255.0;
		} else{
			unsigned short mc[3]=max_color;
			double frac=atan(atan_step_max*(value-avg)/(max-avg))*inv_atan_max;
			if((fabs(max-avg)<EPS) || (frac>1.0)) frac=1.0;
			
			r=((1.0-frac)*ac[0]+frac*mc[0])/255.0;
			g=((1.0-frac)*ac[1]+frac*mc[1])/255.0;
			b=((1.0-frac)*ac[2]+frac*mc[2])/255.0;
		}
	}
	
	return Vec3(r,g,b);
}

inline Vec3 ColorMappingVecES(double value, double min, double max, double avg)
{
	unsigned short ac[3]={0,0,0};
	
	double r,g,b,q;
	if(fabs(min)>fabs(max)){
		q=min;
	} else q=max;
	if((value<-q) || (value>q)){
		r=0.5;
		g=0.5;
		b=0.5;
	} else{
		if(value<avg){
			unsigned short mc[3]={0,0,255};
			double frac=atan(atan_step_max*(avg-value)/(q-avg))*inv_atan_max;
			if((fabs(avg-q)<EPS) || (frac>1.0)) frac=1.0;
			
			r=((1.0-frac)*ac[0]+frac*mc[0])/255.0;
			g=((1.0-frac)*ac[1]+frac*mc[1])/255.0;
			b=((1.0-frac)*ac[2]+frac*mc[2])/255.0;
		} else{
			unsigned short mc[3]={255,0,0};
			double frac=atan(atan_step_max*(value-avg)/(q-avg))*inv_atan_max;
			if((fabs(q-avg)<EPS) || (frac>1.0)) frac=1.0;
			
			r=((1.0-frac)*ac[0]+frac*mc[0])/255.0;
			g=((1.0-frac)*ac[1]+frac*mc[1])/255.0;
			b=((1.0-frac)*ac[2]+frac*mc[2])/255.0;
		}
	}
	
	return Vec3(r,g,b);
}

/// Write element definitions to x3d file
inline void define_element(ofstream &vizout, Element_Type* type, Element_Group* group, bool animate, unsigned int steps, Config_Data* configuration, bool output_directly, bool inside_transparency, Vec3 position, Vec4 rotation, double transparency)
{
	string pretabs="\t\t";
	if(output_directly) pretabs="";
	string pixeltexture="";
	if(type->eps_texture){
		pixeltexture+=int2str(phi_res)+" "+int2str(theta_res*theta_oversampling)+" 3";
		double dt=pi/(theta_res*theta_oversampling);
		for(unsigned int j=0; j<theta_res*theta_oversampling; j++){
			for(unsigned int i=0; i<phi_res; i++){
				if(i==0) pixeltexture+="\n"; else pixeltexture+=" ";
				unsigned int p=(i+(phi_res>>2))%phi_res; // phi shift by 90 degrees
				unsigned int t=(unsigned int)(theta_res*(1.0-(1.0-cos(j*dt))/2.0)); // theta reversed
				pixeltexture+=ColorMapping(type->eps_texture[t*phi_res+p],type->eps_texture[phi_res*theta_res+2],type->eps_texture[phi_res*theta_res+3],1.0);
			}
		}
	}
	string groupname="";
	if(group){
		groupname=" ("+group->Type->name+")";
	}
	if(!output_directly){
		vizout << pretabs << "<ProtoDeclare name='" << type->name+groupname << "'>\n";
		vizout << pretabs << "\t<ProtoInterface>\n";
		vizout << pretabs << "\t\t<field accessType='inputOutput' name='position' type='SFVec3f' value='0 0 0'/>\n";
		vizout << pretabs << "\t\t<field accessType='inputOutput' name='rotation' type='SFRotation' value='0 0 0 0'/>\n";
		vizout << pretabs << "\t\t<field accessType='inputOutput' name='transparency' type='SFFloat' value='0.2'/>\n";
		if(inside_transparency) vizout << pretabs << "\t\t<field accessType='inputOutput' name='inside_transparency' type='SFFloat' value='0.0'/>\n";
		if(animate){
			vizout << pretabs << "\t\t<field accessType='inputOutput' name='key' type='MFFloat'/>\n";
			vizout << pretabs << "\t\t<field accessType='inputOutput' name='position_values' type='MFVec3f'/>\n";
			vizout << pretabs << "\t\t<field accessType='inputOutput' name='rotation_values' type='MFRotation'/>\n";
		}
		vizout << pretabs << "\t</ProtoInterface>\n";
		vizout << pretabs << "\t<ProtoBody>\n";
	}
	vizout << pretabs << "\t\t<Group>\n";
	bool show_original=false;
	if(!group || (group && ((group->levelofdetail==0) || ((group->levelofdetail>0) && (!group->Type->LOD->visual_original[group->levelofdetail-1]))))){
		if(output_directly){
			vizout << pretabs << "\t\t\t<Transform DEF='overall_transformation' scale='" << type->saxes.V3Str(' ') << "' translation='" << position.V3Str(' ') << "' rotation='" << rotation.V4Str(' ') << "'>\n";
		} else{
			vizout << pretabs << "\t\t\t<Transform DEF='overall_transformation' scale='" << type->saxes.V3Str(' ') << "'>\n";
			vizout << pretabs << "\t\t\t\t<IS>\n";
			vizout << pretabs << "\t\t\t\t\t<connect nodeField='translation' protoField='position'/>\n";
			vizout << pretabs << "\t\t\t\t\t<connect nodeField='rotation' protoField='rotation'/>\n";
			vizout << pretabs << "\t\t\t\t</IS>\n";
		}
		if(type->texture_is_dipole && (!type->eps_texture && type->texture!="")){
			Mat33 rotation=AxisAngle2Rot(Vec4(1.0,0.0,0.0,pi/2.0)); // base rotation
			Vec3 dipole=type->initial_dipole; // get overall dipole moment of element (including charges)
			for(unsigned i=0; i<type->nr_charges; i++) dipole+=type->q_pos[i]*type->q[i];
			Mat33 dipole_z_rotation=Vec2Rot(dipole,2); // rotates z axis onto dipole moment
			rotation=dipole_z_rotation*rotation;
			vizout << pretabs << "\t\t\t\t<Transform rotation='" << Rot2AxisAngle(rotation).V4Str(' ') << "'>\n";
		} else{
			vizout << pretabs << "\t\t\t\t<Transform rotation='1 0 0 " << pi/2.0 << "'>\n";
		}
		vizout << pretabs << "\t\t\t\t\t<Shape>\n";
		vizout << pretabs << "\t\t\t\t\t\t<Appearance>\n";
		if(output_directly){
			if(type->eps_texture || type->texture!="") vizout << pretabs << "\t\t\t\t\t\t\t<Material transparency='" << transparency << "'>\n"; else vizout << pretabs << "\t\t\t\t\t\t\t<Material diffuseColor='" << type->color.V3Str(' ') << "' shininess='0.8' specularColor='0.8 0.8 0.8' transparency='" << transparency << "'>\n";
		} else{
			if(type->eps_texture || type->texture!="") vizout << pretabs << "\t\t\t\t\t\t\t<Material>\n"; else vizout << pretabs << "\t\t\t\t\t\t\t<Material diffuseColor='" << type->color.V3Str(' ') << "' shininess='0.8' specularColor='0.8 0.8 0.8'>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t\t<IS>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t\t\t<connect nodeField='transparency' protoField='transparency'/>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t\t</IS>\n";
		}
		vizout << pretabs << "\t\t\t\t\t\t\t</Material>\n";
		
		if(type->eps_texture){
			vizout << pretabs << "\t\t\t\t\t\t\t<PixelTexture image='" << pixeltexture << "'>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t<TextureProperties magnificationFilter='NEAREST_PIXEL'/>\n</PixelTexture>\n";
		} else{
			if(type->texture!="") vizout << pretabs << "\t\t\t\t\t\t\t<ImageTexture url='\"" << type->texture << "\"'/>\n";
		}
		vizout << pretabs << "\t\t\t\t\t\t</Appearance>\n";
		if(!group)
			vizout << pretabs << "\t\t\t\t\t\t<Sphere radius='" << type->scale_visual_vdw << "'/>\n";
		else
			vizout << pretabs << "\t\t\t\t\t\t<Sphere radius='" << group->Type->scale_visual_vdw << "'/>\n";
		vizout << pretabs << "\t\t\t\t\t</Shape>\n";
		vizout << pretabs << "\t\t\t\t</Transform>\n";
		bool show_labels=((type->label_elements) && (type->transparency>0.0));
		if(group){
			show_labels=((group->Type->label_elements) && ((group->Type->transparency>0.0) || ((group->Type->transparency<0.0) && show_labels)));
		}
		if(show_labels){
			vizout << pretabs << "\t\t\t\t<Transform rotation='1 0 0 " << pi/2.0 << "'>\n";
			vizout << pretabs << "\t\t\t\t\t<Billboard axisOfRotation='0 0 0'>\n";
			vizout << pretabs << "\t\t\t\t\t\t<Shape>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t<Appearance>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t\t<Material diffuseColor='0.9 0.9 0.9' emissiveColor='0.9 0.9 0.9'>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t\t\t<IS>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t<connect nodeField='transparency' protoField='transparency'/>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t\t\t</IS>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t\t</Material>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t</Appearance>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t<Text string='\"" << type->name << "\"'>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t\t<FontStyle justify='\"MIDDLE\" \"MIDDLE\"' size='0.8'/>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t</Text>\n";
			vizout << pretabs << "\t\t\t\t\t\t</Shape>\n";
			vizout << pretabs << "\t\t\t\t\t</Billboard>\n";
			vizout << pretabs << "\t\t\t\t</Transform>\n";
		}
		vizout << pretabs << "\t\t\t</Transform>\n";
	} else{
		if(group && ((group->levelofdetail>0) && (group->Type->LOD->visual_original[group->levelofdetail-1]))){
			Element_Group* original_group=group->Type->LOD->groups[0];
			unsigned int element_nr=0;
			bool found=false;
			while(element_nr<group->nr_elements){
				if(type->thistype==configuration->group_elements[group->elements[element_nr]].mytype){
					found=true;
					break;
				}
				element_nr++;
			}
			if(found){ // "magic" starts to happen here ...
				Element* thiselement=&configuration->group_elements[group->elements[element_nr]];
				Mat33 antirot=thiselement->rot.M3Transpose(); // for rotation matrix R^(-1) = R^T
				Vec3 axis(1.0,0.0,0.0);
				double theta=pi/2.0;
				show_original=true;
				if(output_directly){
					vizout << pretabs << "\t\t\t<Transform DEF='inside_transformation' translation='" << position.V3Str(' ') << "' rotation='" << rotation.V4Str(' ') << "'>\n";
				} else{
					vizout << pretabs << "\t\t\t<Transform DEF='inside_transformation'>\n";
					vizout << pretabs << "\t\t\t\t<IS>\n";
					vizout << pretabs << "\t\t\t\t\t<connect nodeField='translation' protoField='position'/>\n";
					vizout << pretabs << "\t\t\t\t\t<connect nodeField='rotation' protoField='rotation'/>\n";
					vizout << pretabs << "\t\t\t\t</IS>\n";
				}
				vizout << pretabs << "\t\t\t\t<Group>\n";
				Vec3 ydir(0.0,1.0,0.0);
				for(unsigned int i=0; (int)i<group->Type->LOD->element_groups[group->levelofdetail-1][group->Type->LOD->ellipsoid_counts[group->levelofdetail-1][element_nr]]; i++){
					Element* element=&configuration->group_elements[original_group->elements[group->Type->LOD->element_groups[group->levelofdetail-1][group->Type->LOD->ellipsoid_counts[group->levelofdetail-1][element_nr]+i+1]-1]];
					Vec3 translation=antirot*(element->center-thiselement->center);
					Mat33 rotation=element->rot*antirot;
					vizout << pretabs << "\t\t\t\t\t<Transform scale='" << element->MyType->saxes.V3Str(' ') << "' translation='" << translation.V3Str(' ') << "'>\n";
					vizout << pretabs << "\t\t\t\t\t\t<Transform rotation='" << Rot2AxisAngle(rotation*AxisAngle2Rot(axis,theta)).V4Str(' ') << "'>\n";
					vizout << pretabs << "\t\t\t\t\t\t\t<Shape>\n";
					vizout << pretabs << "\t\t\t\t\t\t\t\t<Appearance>\n";
					vizout << pretabs << "\t\t\t\t\t\t\t\t\t<Material diffuseColor='";
					if(group->Type->LOD->internal_charge_colors) vizout << ColorMappingVecES(element->MyType->charge,group->Type->LOD->groups[0]->Type->mincharge,group->Type->LOD->groups[0]->Type->maxcharge,group->Type->LOD->groups[0]->Type->charge).V3Str(' '); else vizout << element->MyType->color.V3Str(' ');
					vizout << "' shininess='0.8' specularColor='0.8 0.8 0.8'";
					if(inside_transparency){
						vizout << ">\n";
						vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t<IS>\n";
						vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t<connect nodeField='transparency' protoField='inside_transparency'/>\n";
						vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t</IS>\n";
						vizout << pretabs << "\t\t\t\t\t\t\t\t\t</Material>\n";
					} else{
						vizout << " transparency='" << element->MyType->transparency << "'/>\n";
					}
					if(element->MyType->texture!="") vizout << pretabs << "\t\t\t\t\t\t\t\t\t<ImageTexture url='\"" << element->MyType->texture << "\"'/>\n";
					vizout << pretabs << "\t\t\t\t\t\t\t\t</Appearance>\n";
					vizout << pretabs << "\t\t\t\t\t\t\t\t<Sphere radius='" << group->Type->LOD->scale_original_vdw << "'/>\n";
					vizout << pretabs << "\t\t\t\t\t\t\t</Shape>\n";
					vizout << pretabs << "\t\t\t\t\t\t</Transform>\n";
					vizout << pretabs << "\t\t\t\t\t</Transform>\n";
					for(unsigned int j=0; (j<element->nr_interactions && group->Type->LOD->show_internal_bonds); j++){ // take care of bonds
						Interaction* bond=&element->interactions[j];
						if(group->Type->LOD->element_groups[group->levelofdetail-1][group->Type->LOD->ellipsoid_counts[group->levelofdetail-1][element_nr]+i+1]-1<(int)bond->partner){
							Vec3 direction=antirot*(configuration->group_elements[original_group->elements[bond->partner]].center-thiselement->center)-translation;
							Vec3 middle=translation+direction*0.5;
							double bond_length=direction.V3Norm();
							direction/=bond_length;
							Vec4 rot=Rot2AxisAngle(RotAtoB(ydir,direction));
							if(bond->bond_order<2.0){ // single bond
								vizout << pretabs << "\t\t\t\t\t<Transform translation='" << middle.V3Str(' ') << "'>\n";
								vizout << pretabs << "\t\t\t\t\t\t<Transform rotation='" << rot.V4Str(' ') << "'>\n";
								vizout << pretabs << "\t\t\t\t\t\t\t<Shape>\n";
								vizout << pretabs << "\t\t\t\t\t\t\t\t<Appearance>\n";
								vizout << pretabs << "\t\t\t\t\t\t\t\t\t<Material diffuseColor='0.5 0.5 0.5' shininess='0.8' specularColor='0.8 0.8 0.8'";
								if(inside_transparency){
									vizout << ">\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t<IS>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t<connect nodeField='transparency' protoField='inside_transparency'/>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t</IS>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t</Material>\n";
								} else{
									vizout << " transparency='" << element->MyType->transparency << "'/>\n";
								}
								vizout << pretabs << "\t\t\t\t\t\t\t\t</Appearance>\n";
								vizout << pretabs << "\t\t\t\t\t\t\t\t<Cylinder radius='0.05' top='false' bottom='false' height='" << bond_length << "'/>\n";
								vizout << pretabs << "\t\t\t\t\t\t\t</Shape>\n";
								vizout << pretabs << "\t\t\t\t\t\t</Transform>\n";
								vizout << pretabs << "\t\t\t\t\t</Transform>\n";
							} else{
								if((int)qround(bond->bond_order)==2){ // double bond
									vizout << pretabs << "\t\t\t\t\t<Transform translation='" << middle.V3Str(' ') << "'>\n";
									vizout << pretabs << "\t\t\t\t\t\t<Transform rotation='" << rot.V4Str(' ') << "'>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t<Transform translation='0.1 0 0'>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t<Shape>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t<Appearance>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t<Material diffuseColor='0.5 0.5 0.5' shininess='0.8' specularColor='0.8 0.8 0.8'";
									if(inside_transparency){
										vizout << ">\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t<IS>\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t\t<connect nodeField='transparency' protoField='inside_transparency'/>\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t</IS>\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t</Material>\n";
									} else{
										vizout << " transparency='" << element->MyType->transparency << "'/>\n";
									}
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t</Appearance>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t<Cylinder radius='0.05' top='false' bottom='false' height='" << bond_length << "'/>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t</Shape>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t</Transform>\n";
									vizout << pretabs << "\t\t\t\t\t\t</Transform>\n";
									vizout << pretabs << "\t\t\t\t\t</Transform>\n";
									vizout << pretabs << "\t\t\t\t\t<Transform translation='" << middle.V3Str(' ') << "'>\n";
									vizout << pretabs << "\t\t\t\t\t\t<Transform rotation='" << rot.V4Str(' ') << "'>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t<Transform translation='-0.1 0 0'>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t<Shape>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t<Appearance>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t<Material diffuseColor='0.5 0.5 0.5' shininess='0.8' specularColor='0.8 0.8 0.8'";
									if(inside_transparency){
										vizout << ">\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t<IS>\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t\t<connect nodeField='transparency' protoField='inside_transparency'/>\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t</IS>\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t</Material>\n";
									} else{
										vizout << " transparency='" << element->MyType->transparency << "'/>\n";
									}
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t</Appearance>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t<Cylinder radius='0.05' top='false' bottom='false' height='" << bond_length << "'/>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t</Shape>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t</Transform>\n";
									vizout << pretabs << "\t\t\t\t\t\t</Transform>\n";
									vizout << pretabs << "\t\t\t\t\t</Transform>\n";
								} else{ // triple bond
									vizout << pretabs << "\t\t\t\t\t<Transform translation='" << middle.V3Str(' ') << "'>\n";
									vizout << pretabs << "\t\t\t\t\t\t<Transform rotation='" << rot.V4Str(' ') << "'>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t<Transform translation='0 0 0.12'>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t<Shape>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t<Appearance>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t<Material diffuseColor='0.5 0.5 0.5' shininess='0.8' specularColor='0.8 0.8 0.8'";
									if(inside_transparency){
										vizout << ">\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t<IS>\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t\t<connect nodeField='transparency' protoField='inside_transparency'/>\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t</IS>\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t</Material>\n";
									} else{
										vizout << " transparency='" << element->MyType->transparency << "'/>\n";
									}
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t</Appearance>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t<Cylinder radius='0.05' top='false' bottom='false' height='" << bond_length << "'/>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t</Shape>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t</Transform>\n";
									vizout << pretabs << "\t\t\t\t\t\t</Transform>\n";
									vizout << pretabs << "\t\t\t\t\t</Transform>\n";
									vizout << pretabs << "\t\t\t\t\t<Transform translation='" << middle.V3Str(' ') << "'>\n";
									vizout << pretabs << "\t\t\t\t\t\t<Transform rotation='" << rot.V4Str(' ') << "'>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t<Transform translation='0.1039231 0 -0.06'>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t<Shape>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t<Appearance>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t<Material diffuseColor='0.5 0.5 0.5' shininess='0.8' specularColor='0.8 0.8 0.8'";
									if(inside_transparency){
										vizout << ">\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t<IS>\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t\t<connect nodeField='transparency' protoField='inside_transparency'/>\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t</IS>\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t</Material>\n";
									} else{
										vizout << " transparency='" << element->MyType->transparency << "'/>\n";
									}
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t</Appearance>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t<Cylinder radius='0.05' top='false' bottom='false' height='" << bond_length << "'/>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t</Shape>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t</Transform>\n";
									vizout << pretabs << "\t\t\t\t\t\t</Transform>\n";
									vizout << pretabs << "\t\t\t\t\t</Transform>\n";
									vizout << pretabs << "\t\t\t\t\t<Transform translation='" << middle.V3Str(' ') << "'>\n";
									vizout << pretabs << "\t\t\t\t\t\t<Transform rotation='" << rot.V4Str(' ') << "'>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t<Transform translation='-0.1039231 0 -0.06'>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t<Shape>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t<Appearance>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t<Material diffuseColor='0.5 0.5 0.5' shininess='0.8' specularColor='0.8 0.8 0.8'";
									if(inside_transparency){
										vizout << ">\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t<IS>\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t\t<connect nodeField='transparency' protoField='inside_transparency'/>\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t\t</IS>\n";
										vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t</Material>\n";
									} else{
										vizout << " transparency='" << element->MyType->transparency << "'/>\n";
									}
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t</Appearance>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t\t<Cylinder radius='0.05' top='false' bottom='false' height='" << bond_length << "'/>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t\t</Shape>\n";
									vizout << pretabs << "\t\t\t\t\t\t\t</Transform>\n";
									vizout << pretabs << "\t\t\t\t\t\t</Transform>\n";
									vizout << pretabs << "\t\t\t\t\t</Transform>\n";
								}
							}
						}
					}
				}
				vizout << pretabs << "\t\t\t\t</Group>\n";
				vizout << pretabs << "\t\t\t</Transform>\n";
			}
		}
		if(output_directly){
			vizout << pretabs << "\t\t\t<Transform DEF='overall_transformation' translation='" << position.V3Str(' ') << "' rotation='" << rotation.V4Str(' ') << "'>\n";
		} else{
			vizout << pretabs << "\t\t\t<Transform DEF='overall_transformation'>\n";
			vizout << pretabs << "\t\t\t\t<IS>\n";
			vizout << pretabs << "\t\t\t\t\t<connect nodeField='translation' protoField='position'/>\n";
			vizout << pretabs << "\t\t\t\t\t<connect nodeField='rotation' protoField='rotation'/>\n";
			vizout << pretabs << "\t\t\t\t</IS>\n";
		}
		vizout << pretabs << "\t\t\t\t<Transform scale='" << type->saxes.V3Str(' ') << "'>\n";
		if(type->texture_is_dipole && (!type->eps_texture && type->texture!="")){
			Mat33 rotation=AxisAngle2Rot(Vec4(1.0,0.0,0.0,pi/2.0)); // base rotation
			Vec3 dipole=type->initial_dipole;
			Mat33 dipole_z_rotation=Vec2Rot(dipole,2); // rotates z axis onto dipole moment
			for(unsigned i=0; i<type->nr_charges; i++) dipole+=type->q_pos[i]*type->q[i];
			rotation=dipole_z_rotation*rotation;
			vizout << pretabs << "\t\t\t\t\t<Transform rotation='" << Rot2AxisAngle(rotation).V4Str(' ') << "'>\n";
		} else{
			vizout << pretabs << "\t\t\t\t\t<Transform rotation='1 0 0 " << pi/2.0 << "'>\n";
		}
		vizout << pretabs << "\t\t\t\t\t\t<Shape>\n";
		vizout << pretabs << "\t\t\t\t\t\t\t<Appearance>\n";
		string ts=" transparency='"+double2str(transparency)+"'";
		if(type->eps_texture || type->texture!=""){
			vizout << pretabs << "\t\t\t\t\t\t\t\t<Material" << ts << ">\n";
		} else{
			if(group && ((group->levelofdetail>0) && (group->Type->LOD->visual_original[group->levelofdetail-1]))) vizout << pretabs << "\t\t\t\t\t\t\t<Material" << ts << " diffuseColor='" << type->color.V3Str(' ') << "'>\n"; else vizout << pretabs << "\t\t\t\t\t\t\t<Material" << ts << " diffuseColor='" << type->color.V3Str(' ') << "' shininess='0.8' specularColor='0.8 0.8 0.8'>\n";
		}
		if(!output_directly){
			vizout << pretabs << "\t\t\t\t\t\t\t\t\t<IS>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t\t\t\t<connect nodeField='transparency' protoField='transparency'/>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t\t\t</IS>\n";
		}
		vizout << pretabs << "\t\t\t\t\t\t\t\t</Material>\n";
		
		if(type->eps_texture){
			vizout << pretabs << "\t\t\t\t\t\t\t\t<PixelTexture image='" << pixeltexture << "'>\n";
			vizout << pretabs << "\t\t\t\t\t\t\t\t<TextureProperties magnificationFilter='NEAREST_PIXEL'/>\n</PixelTexture>\n";
		} else{
			if(type->texture!="") vizout << pretabs << "\t\t\t\t\t\t\t\t<ImageTexture url='\"" << type->texture << "\"'/>\n";
		}
		if(!output_directly && (group && (group->Type->transparency<0.5))) vizout << pretabs << "\t\t\t\t\t\t\t\t<BlendMode srcFactor='src_alpha' destFactor='one_minus_src_color'/>\n"; else vizout << pretabs << "\t\t\t\t\t\t\t\t<BlendMode srcFactor='src_alpha' destFactor='one_minus_src_alpha'/>\n";
		vizout << pretabs << "\t\t\t\t\t\t\t</Appearance>\n";
		vizout << pretabs << "\t\t\t\t\t\t\t<Sphere radius='" << group->Type->scale_visual_vdw << "'/>\n";
		vizout << pretabs << "\t\t\t\t\t\t</Shape>\n";
		vizout << pretabs << "\t\t\t\t\t</Transform>\n";
		vizout << pretabs << "\t\t\t\t</Transform>\n";
		vizout << pretabs << "\t\t\t</Transform>\n";
	}
	vizout << pretabs << "\t\t</Group>\n";
	if(animate){
		vizout << pretabs << "\t\t<TimeSensor DEF='SimGrStep' cycleInterval='" << steps << "' loop='true'/>\n";
		vizout << pretabs << "\t\t<PositionInterpolator DEF='positions'>\n";
		vizout << pretabs << "\t\t\t<IS>\n";
		vizout << pretabs << "\t\t\t\t<connect nodeField='key' protoField='key'/>\n";
		vizout << pretabs << "\t\t\t\t<connect nodeField='keyValue' protoField='position_values'/>\n";
		vizout << pretabs << "\t\t\t</IS>\n";
		vizout << pretabs << "\t\t</PositionInterpolator>\n";
		vizout << pretabs << "\t\t<OrientationInterpolator DEF='rotations'>\n";
		vizout << pretabs << "\t\t\t<IS>\n";
		vizout << pretabs << "\t\t\t\t<connect nodeField='key' protoField='key'/>\n";
		vizout << pretabs << "\t\t\t\t<connect nodeField='keyValue' protoField='rotation_values'/>\n";
		vizout << pretabs << "\t\t\t</IS>\n";
		vizout << pretabs << "\t\t</OrientationInterpolator>\n";
		vizout << pretabs << "\t\t<ROUTE fromNode='SimGrStep' fromField='fraction_changed' toNode='positions' toField='set_fraction'/>\n";
		vizout << pretabs << "\t\t<ROUTE fromNode='SimGrStep' fromField='fraction_changed' toNode='rotations' toField='set_fraction'/>\n";
		vizout << pretabs << "\t\t<ROUTE fromNode='positions' fromField='value_changed' toNode='overall_transformation' toField='translation'/>\n";
		vizout << pretabs << "\t\t<ROUTE fromNode='rotations' fromField='value_changed' toNode='overall_transformation' toField='rotation'/>\n";
		if(show_original){
			vizout << pretabs << "\t\t<ROUTE fromNode='positions' fromField='value_changed' toNode='inside_transformation' toField='translation'/>\n";
			vizout << pretabs << "\t\t<ROUTE fromNode='rotations' fromField='value_changed' toNode='inside_transformation' toField='rotation'/>\n";
		}
	}
	if(!output_directly){
		vizout << pretabs << "\t</ProtoBody>\n";
		vizout << pretabs << "</ProtoDeclare>\n";
	}
}

/// Wrapper for above
inline void define_element(ofstream &vizout, Element_Type* type, Element_Group* group, bool animate, unsigned int steps, Config_Data* configuration, bool inside_transparency)
{
	define_element(vizout,type,group,animate,steps,configuration,false,inside_transparency,Vec3(0.0),Vec4(0.0,0.0,0.0,0.0),0.0);
}

/// Wrapper for above
inline void define_element(ofstream &vizout, Element_Type* type, Element_Group* group, bool animate, unsigned int steps, Config_Data* configuration)
{
	define_element(vizout,type,group,animate,steps,configuration,false,false,Vec3(0.0),Vec4(0.0,0.0,0.0,0.0),0.0);
}


/// Output visualization file footer, including information about axes
inline void output_footer(ofstream &vizout, axis_t &axis)
{
	//Draw axes
	if (axis.DrawAxes == 1) {
		vizout << "\t\t<!-- Draw coordinate system -->\n";
		vizout << "\t\t<ProtoInstance name='line'>\n";
		vizout << "\t\t\t<fieldValue name='AB' value='0 0 0 " << axis.length.vec[0] << " 0 0'/>\n";
		vizout << "\t\t\t<fieldValue name='color' value='" << axis.ColorX.V3Str(' ') << "'/>\n";
		vizout << "\t\t</ProtoInstance>\n";
		vizout << "\t\t<ProtoInstance name='Vector3D'>\n";
		vizout << "\t\t\t<fieldValue name='position' value='" << axis.length.vec[0]/2.0 << " 0 0'/>\n";
		vizout << "\t\t\t<fieldValue name='tip_vector' value='0 " << axis.length.vec[0]/2.0 << " 0'/>\n";
		vizout << "\t\t\t<fieldValue name='length' value='" << axis.length.vec[0] << "'/>\n";
		vizout << "\t\t\t<fieldValue name='rotation_from_y' value='0 0 1 " << -pi/2.0 << "'/>\n";
		vizout << "\t\t\t<fieldValue name='color' value='" << axis.ColorX.V3Str(' ') << "'/>\n";
		vizout << "\t\t</ProtoInstance>\n";
		vizout << "\t\t<Transform rotation='1 0 0 " << pi/2.0 << "' translation='" << axis.length.vec[0]+0.2 << " 0 0'>\n";
		vizout << "\t\t\t<Billboard axisOfRotation='0 0 0'>\n";
		vizout << "\t\t\t\t<Shape>\n";
		vizout << "\t\t\t\t\t<Appearance>\n";
		vizout << "\t\t\t\t\t\t<Material diffuseColor='0.9 0.9 0.9' emissiveColor='0.9 0.9 0.9'/>\n";
		vizout << "\t\t\t\t\t</Appearance>\n";
		vizout << "\t\t\t\t\t<Text string='\"x\"'>\n";
		vizout << "\t\t\t\t\t\t<FontStyle justify='\"MIDDLE\" \"MIDDLE\" \"MIDDLE\"' size='0.8'/>\n";
		vizout << "\t\t\t\t\t</Text>\n";
		vizout << "\t\t\t\t</Shape>\n";
		vizout << "\t\t\t</Billboard>\n";
		vizout << "\t\t</Transform>\n";
		vizout << "\t\t<ProtoInstance name='line'>\n";
		vizout << "\t\t\t<fieldValue name='AB' value='0 0 0 0 " << axis.length.vec[1] << " 0'/>\n";
		vizout << "\t\t\t<fieldValue name='color' value='" << axis.ColorY.V3Str(' ') << "'/>\n";
		vizout << "\t\t</ProtoInstance>\n";
		vizout << "\t\t<ProtoInstance name='Vector3D'>\n";
		vizout << "\t\t\t<fieldValue name='position' value='0 " << axis.length.vec[1]/2.0 << " 0'/>\n";
		vizout << "\t\t\t<fieldValue name='tip_vector' value='0 " << axis.length.vec[1]/2.0 << " 0'/>\n";
		vizout << "\t\t\t<fieldValue name='length' value='" << axis.length.vec[1] << "'/>\n";
		vizout << "\t\t\t<fieldValue name='color' value='" << axis.ColorY.V3Str(' ') << "'/>\n";
		vizout << "\t\t</ProtoInstance>\n";
		vizout << "\t\t<Transform rotation='1 0 0 " << pi/2.0 << "' translation='0 " << axis.length.vec[1]+0.2 << " 0'>\n";
		vizout << "\t\t\t<Billboard axisOfRotation='0 0 0'>\n";
		vizout << "\t\t\t\t<Shape>\n";
		vizout << "\t\t\t\t\t<Appearance>\n";
		vizout << "\t\t\t\t\t\t<Material diffuseColor='0.9 0.9 0.9' emissiveColor='0.9 0.9 0.9'/>\n";
		vizout << "\t\t\t\t\t</Appearance>\n";
		vizout << "\t\t\t\t\t<Text string='\"y\"'>\n";
		vizout << "\t\t\t\t\t\t<FontStyle justify='\"MIDDLE\" \"MIDDLE\" \"MIDDLE\"' size='0.8'/>\n";
		vizout << "\t\t\t\t\t</Text>\n";
		vizout << "\t\t\t\t</Shape>\n";
		vizout << "\t\t\t</Billboard>\n";
		vizout << "\t\t</Transform>\n";
		vizout << "\t\t<ProtoInstance name='line'>\n";
		vizout << "\t\t\t<fieldValue name='AB' value='0 0 0 0 0 " << axis.length.vec[2] << "'/>\n";
		vizout << "\t\t\t<fieldValue name='color' value='" << axis.ColorZ.V3Str(' ') << "'/>\n";
		vizout << "\t\t</ProtoInstance>\n";
		vizout << "\t\t<ProtoInstance name='Vector3D'>\n";
		vizout << "\t\t\t<fieldValue name='position' value='0 0 " << axis.length.vec[2]/2.0 << "'/>\n";
		vizout << "\t\t\t<fieldValue name='tip_vector' value='0 " << axis.length.vec[2]/2.0 << " 0'/>\n";
		vizout << "\t\t\t<fieldValue name='length' value='" << axis.length.vec[2] << "'/>\n";
		vizout << "\t\t\t<fieldValue name='rotation_from_y' value='1 0 0 " << pi/2.0 << "'/>\n";
		vizout << "\t\t\t<fieldValue name='color' value='" << axis.ColorZ.V3Str(' ') << "'/>\n";
		vizout << "\t\t</ProtoInstance>\n";
		vizout << "\t\t<Transform rotation='1 0 0 " << pi/2.0 << "' translation='0 0 " << axis.length.vec[2]+0.2 << "'>\n";
		vizout << "\t\t\t<Billboard axisOfRotation='0 0 0'>\n";
		vizout << "\t\t\t\t<Shape>\n";
		vizout << "\t\t\t\t\t<Appearance>\n";
		vizout << "\t\t\t\t\t\t<Material diffuseColor='0.9 0.9 0.9' emissiveColor='0.9 0.9 0.9'/>\n";
		vizout << "\t\t\t\t\t</Appearance>\n";
		vizout << "\t\t\t\t\t<Text string='\"z\"'>\n";
		vizout << "\t\t\t\t\t\t<FontStyle justify='\"MIDDLE\" \"MIDDLE\" \"MIDDLE\"' size='0.8'/>\n";
		vizout << "\t\t\t\t\t</Text>\n";
		vizout << "\t\t\t\t</Shape>\n";
		vizout << "\t\t\t</Billboard>\n";
		vizout << "\t\t</Transform>\n";
		vizout << "\t\t<!-- done -->\n";
	}
	vizout << "\t</Scene>\n";
	vizout << "</X3D>\n";
}

