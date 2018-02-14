/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

/*!\file
 * AT May 12, 2011
*/

#ifndef INCLUDED_CONFIGREADER
#define INCLUDED_CONFIGREADER

#include <string.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>

#ifndef INCLUDED_CONFIG
#include "Config.h"
#endif
#ifndef INCLUDED_SCALARMAT
#include "ScalarMat.h"
#endif

template <class T> bool from_string(T &t, std::string &s);
template <class T> bool from_string(T &t, std::string &s){
	std::stringstream stream(s);
	return (stream >> t).fail();
}

template <class T> bool from_string(T &t, const char* s);
template <class T> bool from_string(T &t, const char* s){
	std::stringstream stream(s);
	return (stream >> t).fail();
}

template <class T> string to_string(T t);
template <class T> string to_string(T t){
	std::stringstream stream;
	stream << t;
	return stream.str();
}

typedef struct _readfile{
	ifstream file;
	string filename;
	string directory;
} readfile;

class ConfigReader
{
public:
	char* constant_block; ///< block to look for constants if item could not be found in current section
	
	ConfigReader();
	virtual ~ConfigReader();
	
	int GetNrSections(readfile* file, string name, unsigned int pos){ return GetNrSections(file,name,pos,true); };
	int GetNrSections(readfile* file, string name, unsigned int pos, bool combine_sections);
	
	char* GetSection(readfile* file, string name, unsigned int &pos){ string subname; return GetSection(file,name,pos,subname); };
	char* GetSection(readfile* file, string name, unsigned int &pos, string &subname){ unsigned int include_position=0; return GetSection(file,name,pos,include_position,subname,true); };
	char* GetSection(readfile* file, string name, unsigned int &pos, unsigned int &include_position, string &subname, bool use_includes){ return GetSection(file,name,pos,include_position,subname,use_includes,true); };
	char* GetSection(readfile* file, string name, unsigned int &pos, unsigned int &include_position, string &subname, bool use_includes, bool combine_sections){ return GetSection(file,name,pos,include_position,subname,use_includes,combine_sections,NULL); };
	char* GetSection(readfile* file, string name, unsigned int &pos, unsigned int &include_position, string &subname, bool use_includes, bool combine_sections, string* inc_files);
	
	char* GetExclusiveSection(readfile* file, string name, unsigned int &pos, string &subname){ unsigned int include_position=0; return GetSection(file,name,pos,include_position,subname,true,false,NULL); };
	
	char* GetItemValue(char* conf, string item){ return GetItemValue(conf,item,false); };
	char* GetItemValue(char* conf, string item, bool item_is_string){ return GetItemValue(conf,item,item_is_string,true); };
	char* GetItemValue(char* conf, string item, bool item_is_string, bool get_last){ unsigned int pos=0; return GetItemValue(conf,item,pos,item_is_string,get_last); };
	char* GetItemValue(char* conf, string item, unsigned int &pos, bool item_is_string, bool get_last);
	
	double* get_flex_tupel(const char* item, char* conf){ return get_flex_tupel(item,conf,NULL,NULL,0); };
	double* get_flex_tupel(const char* item, char* conf, char** item_names, unsigned int* element_count, unsigned int nr_elements);
	double* get_MxN_tupel(const char* item, char* conf, int m, int n, double* default_value){ return get_MxN_tupel(item,conf,m,n,default_value,true); };
	double* get_MxN_tupel(const char* item, char* conf, int m, int n, double* default_value, bool warn);
	
	//! set parameters from text configuration file, if item is not found use default_value
	void SetParam(string &s, const char* item, char* conf, string default_value);
	void SetParam(bool &b, const char* item, char* conf, bool default_value);
	void SetParam(unsigned int &i, const char* item, char* conf, int default_value);
	void SetParam(int &i, const char* item, char* conf, int default_value);
	void SetParam(double &d, const char* item, char* conf, double default_value){ return SetParam(d,item,conf,default_value,false); };
	void SetParam(double &d, const char* item, char* conf, double default_value, bool quiet);
	void SetParam(double* v, int len, const char* item, char* conf, double* default_value);
	void SetParam(int* v, int len, const char* item, char* conf, int* default_value);
	void SetParam(bool* v, int len, const char* item, char* conf, bool* default_value);
	void SetParam(double v[3], const char* item, char* conf, Vec3 default_value);
	void SetParam(int v[3], const char* item, char* conf, Vec3 default_value);
	void SetParam(bool v[3], const char* item, char* conf, Vec3 default_value);
	void SetParam(double A[3][3], const char* item, char* conf, Mat33 default_value);
	void SetParam(Mat33 &A, const char* item, char* conf, Mat33 default_value);
	
	//! versions below bail out if item is not specified
	void SetParam(string &s, const char* item, char* conf);
	void SetParam(bool &b, const char* item, char* conf);
	void SetParam(unsigned int &i, const char* item, char* conf);
	void SetParam(int &i, const char* item, char* conf);
	void SetParam(double &d, const char* item, char* conf){ SetParam(d,item,conf,false); };
	bool SetParam(double &d, const char* item, char* conf, bool quiet);
	void SetParam(double* v, int len, const char* item, char* conf);
	void SetParam(int* v, int len, const char* item, char* conf);
	void SetParam(bool* v, int len, const char* item, char* conf);
	void SetParam(double v[3], const char* item, char* conf);
	void SetParam(int v[3], const char* item, char* conf);
	void SetParam(bool v[3], const char* item, char* conf);
protected: //! All of the stuff below is private because it should only be called by class methods (but not from the outside)
	int recursion_count; ///< just to be save limit recursive lookup to a depth of 100 (should be unrestrictive enough ...)
	
	char* GetPreamble(readfile* file);
	virtual char* IncludeFile(readfile* file, char* include_type, unsigned int &pos, unsigned int &include_position, string name, string &subname, bool combine_sections, string* inc_files);
	
	char* ItemPos(char* conf, string item, unsigned int &pos, bool get_last);
	char* GetItemLine(char* conf, string item, unsigned int &pos);
	
	char** string2UPN(const char* value, int &items); ///< takes a string for a calculation (e.g. (1+2)/3) and converts it into a string array (of items-length) in reverse polish notation
	double UPN2double(char** upn, int items, char* conf, bool &success); ///< reverse polish notation string array to double
	bool do_number_magic(double &d, const char* item, const char* value, char* conf); ///< basically takes string value and returns corresponding double
	bool do_int_number_magic(int &i, const char* item, const char* value, char* conf); ///< overloaded version of above
};

/// function to compare if strings are equal (uppercase)
inline bool compare_strings(const char* s1, const char* s2)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	unsigned int len=strlen(s1);
	if(len==strlen(s2)){
		for(unsigned int i=0; i<len; i++){
			if(toupper(s1[i])!=toupper(s2[i])){
#if DEBUG_LEVEL>3
				cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
				return false;
			}
		}
#if DEBUG_LEVEL>3
		cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
		return true;
	}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return false;
}

inline char* find_string(char* content, string search)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	char* pos=content;
	if(pos){
		int len=search.length();
		if(len>0){
			bool found=false;
			while((pos) && ((*pos!='\0') && (!found))){
				pos=strchr(pos,search[0]);
				found=true;
				if(pos){
					for(int i=0; i<len; i++){
						if(*(pos+i)=='\0'){
							found=false;
							break;
						}
						if(toupper(*(pos+i))!=toupper(search[i])){
							pos+=i;
							found=false;
							break;
						}
					}
				}
			}
		}
	}
	if(pos && (*pos=='\0')) pos=NULL;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return pos;
}

inline char* stringNcopy(const char* source, int num) /// I grew tired of typing this ;-)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	char* s=NULL;
	if((source!=NULL) && (num>0)){
		s = new char[num+1];
		strncpy(s,source,num);
		s[num]='\0';
	}
	
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return s;
}

inline string WindowsFilter(string filename)
{
#ifdef IN_WINDOWS
	string result="";
	
	for(unsigned int i=0; i<filename.length(); i++){
		if(filename[i]!='/') result+=filename[i]; else result+="\\";
	}
	
	return result;
#else
	return filename;
#endif
}

inline string GetDirectory(string filename)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	string fn=filename;
#ifdef IN_WINDOWS
	fn=WindowsFilter(filename);
#endif
	string directory="";
	unsigned int i=fn.length()-1;
	char delimiter='/';
// are we compiling for Windows?
#ifdef IN_WINDOWS
	delimiter='\\';
#endif
	while((i>0) && (fn[i]!=delimiter)) i--;
	for(unsigned int j=0; j<i; j++) directory+=fn[j];
	if(directory!="") directory+=delimiter;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return directory;
}

inline void GetDirectory(readfile &file)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
#ifdef IN_WINDOWS
	file.directory=WindowsFilter(file.directory);
	file.filename=WindowsFilter(file.filename);
#endif
	string filename="";
	string directory="";
	unsigned int i=file.filename.length()-1;
// are we compiling for Windows?
#ifdef IN_WINDOWS
	char delimiter='\\';
#else
	char delimiter='/';
#endif
	while((i>0) && (file.filename[i]!=delimiter)) i--;
	for(unsigned int j=0; j<i; j++) directory+=file.filename[j];
	for(unsigned int j=i; j<file.filename.length(); j++) if(file.filename[j]!=delimiter) filename+=file.filename[j];
	
	bool absolute=false;
	if(directory!=""){
#ifdef IN_WINDOWS
		// as usual, Windows a bit long and winding ...
		absolute=(((directory[0]==delimiter) && (directory[1]==delimiter)) || ((directory[1]==':') && (directory[2]==delimiter)));
#else
		absolute=(directory[0]==delimiter);
#endif
		directory+=delimiter;
	}
	
	file.filename=filename;
	if(absolute) file.directory=directory; else file.directory+=directory;
	
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

#endif

