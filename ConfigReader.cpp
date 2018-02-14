/****************************************************/
/* This file is distributed under the               */
/* University of Illinois/NCSA Open Source License. */
/* See LICENSE file in top directory for details.   */
/*                                                  */
/* Copyright (c) 2016 FIGSiM developers             */
/****************************************************/

#include <list>
#include "ConfigReader.h"
#define IN_CONFIGREADER
#include "MC_Elements.h"

/*!
Nothing to see here but a default constructor and destructor. Move along.
Not anymore, look, there's something going on :-) -- AT, Feb 8, 2011
*/
ConfigReader::ConfigReader(){
	recursion_count=0;
	constant_block=NULL;
}

ConfigReader::~ConfigReader()
{
}

char* ConfigReader::ItemPos(char* conf, string item, unsigned int &pos, bool get_last)
{
	char* itemline=GetItemLine(conf,item,pos);

	if(get_last){
		char* newline;
		// try to find more (and use last one only)
		newline=itemline;
		while(newline!=NULL){
			newline=GetItemLine(conf,item,pos);
			if(newline!=NULL){
				delete[] itemline; // delete previous line
				itemline=newline;
			}
		}
	}
	return itemline;
}

/// returns value of item (the string between "=" and the end of the line (or comment start)
/// if item_is_string is true then no whitespace is removed in the entry (means strings don't need to be defined within quotes "")
char* ConfigReader::GetItemValue(char* conf, string item, unsigned int &pos, bool item_is_string, bool get_last)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	char *line, *newline, *temp, *temp2;
	string temp_string="";
	char* itemline=ItemPos(conf,item,pos,get_last);
	line=itemline;

	// line now contains the item line (or NULL)
	if(line!=NULL){ // only return value between "=" and end of line or comments
		newline=strchr(line,'=');
		if (newline==NULL){ // no "=" after item is fine by me as well
			newline=itemline+item.length();
		}
		newline++; // skip to end of time or over "=" sign
		// erase leading whitespaces (I'll count " and ' as well ...)
		while((((*(newline)==' ') || (*(newline)=='\t')) || (*(newline)=='\"')) || (*(newline)=='\'')){
			newline++;
			if((*(newline-1)=='\"') || (*(newline-1)=='\'')){
				item_is_string=true;
				break; // don't touch strings
			}
		}
		temp=strchr(newline,'#');
		temp2=strstr(newline,"//");
		if((temp!=NULL) && (temp2!=NULL)){ // found both comments after value
			if(temp2<temp) temp=temp2; // cut after the closest one
		} else{
			if(temp2!=NULL){ // only "//" found
				temp=temp2;
			} else if(temp==NULL) temp=strchr(newline,'\0'); // no comments, value till end of line
			// else temp stays at "#" comment sign
		}
		// erase trailing whitespaces
		while((((*(temp-1)==' ') || (*(temp-1)=='\t')) || (*(temp-1)=='\"')) || (*(temp-1)=='\'')){
			temp--;
			if((*(temp)=='\"') || (*(temp)=='\'')) break; // don't touch strings (item_is_string should be true already, if not fun occurs down the road ...)
		}
		if(temp>newline){
			// finally get rid of whitespaces in between (if item is not a string)
			for(temp2=newline; temp2<temp; temp2++)
				// '\r' takes care of carriage return found on Windows
				if((*(temp2)!='\r') && (((*(temp2)!=' ') && (*(temp2)!='\t')) || (item_is_string))) temp_string+=*(temp2);
			line = new char[temp_string.length()+1];
			strcpy(line,temp_string.c_str());
		} else line = NULL;
		delete[] itemline;
	}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return line;
}

char* ConfigReader::GetItemLine(char* conf, string item, unsigned int &pos) /// returns line containing item (and position of linestart)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	char* line=NULL;
	if(conf+pos<conf+strlen(conf)){
		char* oldline=conf+pos;
		line = strchr(oldline,'\n');
		char* temp;
		bool found, whitespace;
		int i, j;
		int itemlength=item.length();
		while(line!=NULL){
			if(line-oldline>itemlength+1){ // in order to find item, line needs to be at least as unsigned int as item+1
				found=false;
				for(i=0;i<line-oldline-itemlength; i++){ // find item in line (all uppercase) ...
					whitespace=false;
					for(j=0;j<itemlength; j++){ // ... by going through entire line, until
						if(toupper(item[j])!=toupper(*(oldline+i+j))){
							if((*(oldline+i+j)==' ') || (*(oldline+i+j)=='\t')){
								whitespace=true;
							}
							break;
						}
					}
					if(j==itemlength){ // either item has been found,
						if(((*(oldline+i+j)==' ') || (*(oldline+i+j)=='\t')) || (*(oldline+i+j)=='=')) found=true;
						break;
					}
					if(!whitespace) break;
					if((*(oldline+i)=='#') || (*(oldline+i)=='=')) break; // a comment sign, or "=" sign is encountered (after "=" no items are defined),
					if((*(oldline+i)=='/') && (*(oldline+i+1)=='/')) break; // or a comment starts
				}
				if (found==true){ // item has been found
					temp = new char[line-oldline+1];
					strncpy(temp,oldline,line-oldline);
					temp[line-oldline]='\0';
					pos=line-conf; // set new position to line end
					line=temp;
					break; // done, found item line
				}
			}
			oldline=line+1;
			line=strchr(oldline,'\n');
		}
	}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return line;
}

char* ConfigReader::GetPreamble(readfile* file)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	file->file.seekg(0); // go to beginning of file
	char* conf=new char[READBLOCKSIZE+1];
	file->file.read(conf,READBLOCKSIZE); // read in READBLOCKSIZE byte from file (or less if EOF is reached)
	file->file.clear();
	unsigned int readcount=file->file.gcount();
	conf[readcount]='\0';
	
	char* line = NULL;
	if(readcount>0){
		line = strchr(conf,'[');
		if(line!=conf) while((line) && (*(line-1)!='\n')) line = strchr(conf,'[');
		if(line){
			readcount=(line-conf);
			if(readcount>0){
				line=new char[readcount+1]; line[readcount]='\0';
				file->file.seekg(0);
				file->file.read(line,readcount);
				file->file.clear();
			} else line=NULL;
		} else{
			cout << "Configuration file \"" << file->filename << "\" preamble, i.e. for \"&include\" statements, is too big (>" << READBLOCKSIZE/1024 << " kB)\n";
			exit(1);
		}
	} else cout << "WARNING: Configuration file \"" << file->filename << "\" is of zero size.\n";
	delete[] conf;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return line;
}

char* ConfigReader::IncludeFile(readfile* file, char* include_type, unsigned int &pos, unsigned int &include_position, string name, string &subname, bool combine_sections, string* inc_files)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	char* include=NULL;
	if(compare_strings(include_type,"CONF")){
#if DEBUG_LEVEL>2
		cout << "Including file: " << file->directory+file->filename << "\n";
#endif
		recursion_count++;
		if(recursion_count>RECURSION_LIMIT){
			cout << "\nToo deep includes - Round and round go the circular includes?\n";
			exit(4);
		}
		file->file.open((file->directory+file->filename).c_str(),ios::binary);
		if(file->file.fail()){
			cout << "Could not open included file " << file->filename << ".\n";
			exit(1);
		}
		// Get file size
		file->file.seekg(0,ifstream::end);
		unsigned int filesize=file->file.tellg();
		include=GetSection(file,name,pos,include_position,subname,true,combine_sections,inc_files);
		
		if(combine_sections || !include) include_position+=filesize;
		
		file->file.close();
		recursion_count--;
	} else{
		cout << "WARNING: Can not include configuration file \"" << file->filename << "\" as type <" << include_type << ">.\n";
	}
#if DEBUG_LEVEL>3
		cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return include;
}

/// returns content of section matching uppercase name (file parsing version)
char* ConfigReader::GetSection(readfile* file, string name, unsigned int &pos, unsigned int &include_position, string &subname, bool use_includes, bool combine_sections, string* inc_files)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	string sn_save=subname;
	// Get includes at beginning of file first
	char* included_content=NULL;
	char* conf = NULL;
	if(use_includes) conf=GetPreamble(file);
	unsigned int include_first_section_pos=0;
	if(conf){
		char* include_file;
		unsigned int inc_pos=0;
		while((include_file=GetItemValue(conf,"&include",inc_pos,true,false))){
			// found &include line, now determine include type
			char* hastype=strstr(include_file,"as");
			while((hastype) && (!(((*(hastype-1)==' ') || (*(hastype-1)=='\t')) && ((*(hastype+2)==' ') || (*(hastype+2)=='\t'))))){
				char* temp=strstr(hastype,"as");
				if(temp!=hastype) hastype=temp; else hastype=NULL;
			}
			char* include_filename=include_file;
			char* include_type=new char[5];
			memcpy(include_type,"conf\0",5);
			if(hastype){
				char* fn=hastype-1;
				while((*fn==' ') || (*fn=='\t')){
					fn--;
					if(*fn=='\"'){
						fn--;
						break;
					}
				}
				fn++;
				char* it=hastype+3;
				while((*it==' ') || (*it=='\t')){
					it++;
					if(*it=='\"'){
						it++;
						break;
					}
				}
				include_filename=new char[fn-include_file+1]; include_filename[fn-include_file]='\0';
				memcpy(include_filename,include_file,fn-include_file);
				delete[] include_type;
				include_type=new char[strlen(include_file)-(it-include_file)+1]; include_type[strlen(include_file)-(it-include_file)]='\0';
				memcpy(include_type,it,strlen(include_file)-(it-include_file));
				delete[] include_file;
			}
			readfile inc_file;
			inc_file.filename=include_filename;
			inc_file.directory=file->directory;
			GetDirectory(inc_file);
			if(inc_files){
				if(*inc_files!="") *inc_files+="\n";
				if(inc_file.directory=="") *inc_files+="./";
				*inc_files+=inc_file.directory+"\t"+inc_file.filename;
			}
			char* includefile_content=IncludeFile(&inc_file,include_type,pos,include_position,name,subname,combine_sections,inc_files); // subname set first stays throughout
			if(includefile_content){
				if(include_first_section_pos==0) include_first_section_pos=pos;
				if(included_content){ // append included configuration files sequentially in order of appearance
					char* temp=new char[strlen(included_content)+strlen(includefile_content)+1]; temp[strlen(included_content)+strlen(includefile_content)]='\0';
					memcpy(temp,included_content,strlen(included_content));
					memcpy(temp+strlen(included_content),includefile_content,strlen(includefile_content));
					delete[] includefile_content;
					delete[] included_content;
					included_content=temp;
				} else{
					included_content=includefile_content;
					if(!combine_sections) return included_content;
				}
			}
			delete[] include_filename;
			delete[] include_type;
		}
		delete[] conf;
	}
	char* line=NULL;
	
	string subname_search="";
	if(subname!="") subname_search=subname;
	
	// Get file size
	file->file.seekg(0,ifstream::end);
	unsigned int filesize=file->file.tellg();
	
	int fileposition=pos-include_position;
	if(fileposition<0) fileposition=0;
	if(fileposition<int(filesize)){ // with includes one needs to check if position is in the right file
		bool filestart=false;
		if(fileposition==0) filestart=true;
		file->file.seekg(fileposition); // go to position pos in file
		conf=new char[READBLOCKSIZE+1];
		file->file.read(conf,READBLOCKSIZE); // read in READBLOCKSIZE byte from file (or less if EOF is reached)
		file->file.clear();
		unsigned int readcount=file->file.gcount();
		conf[readcount]='\0';
		subname="";
		if(readcount>0){
			line = strchr(conf,'[');
			while(!line && (readcount==READBLOCKSIZE)){ // read file until '[' found or EOF reached
				filestart=false; // reading is definitely not from the filestart anymore
				fileposition=file->file.tellg(); // set fileposition to current block's start
				file->file.read(conf,READBLOCKSIZE); // read in READBLOCKSIZE byte from file (or less if EOF is reached)
				file->file.clear();
				readcount=file->file.gcount();
				conf[readcount]='\0';
				line = strchr(conf,'[');
			}
			if(line-conf>1){ // now that '[' is found, make sure we read in enough text after it
				fileposition+=(line-conf)-1; // new position is two characters before '['
				file->file.seekg(fileposition); // go to position two characters before '[' in file
				file->file.read(conf,READBLOCKSIZE); // read in READBLOCKSIZE byte from file (or less if EOF is reached)
				file->file.clear();
				readcount=file->file.gcount();
				conf[readcount]='\0';
				line = strchr(conf,'[');
			}
			char *temp, *temp2;
			bool found;
			unsigned int i;
			while(line!=NULL){ // char* are NULL-terminated - that's why ;-)
				if(((line==conf) && (filestart)) || ((*line=='[') && (*(line-1)=='\n'))){
					temp=strchr(line+1,']');
					temp2=strchr(line+1,'\n');
					if(temp<temp2){ // found section (remember: "FILESTART or \n[section name]\n)"
						// Is this the right one?
						found=true; // As in live, one starts enthusiastically ...
						line++;
						while((*line==' ') || (*line=='\t')) line++; // take care of whitespace before section name
						for(i=0; i<name.length(); i++){
							if(toupper(name[i])!=toupper(*line)){
								found=false; // ... and is either disappointed ...
								break;
							}
							line++;
						}
						if(found){ // ... or not :-)
							while((*line==' ') || (*line=='\t')) line++; // take care of whitespace *after* section name
							if(*line==':'){ // subname found
								// take care of whitespace after ":"
								line++;
								while((*line==' ') || (*line=='\t')) line++;
								// take care of whitespace between subname and "]"
								temp--;
								while((*temp==' ') || (*temp=='\t')) temp--;
								// subname is between line and temp
								while(line<=temp){
									subname+=*line;
									line++;
								}
							}
							if(subname_search!=""){
								if(!compare_strings(subname.c_str(),subname_search.c_str())){
									subname="";
									found=false;
								}
							}
						}
						if(found){
							// Section content starts right after "[Section name]" (temp2),
							// this includes '\n' -- makes things easier down the road
							
							// Where does the section end? -- At the next "\n[" ...
							unsigned int posinc=(temp2-conf);
							pos=fileposition+posinc;
							temp=strchr(temp2,'['); // look in current buffer if next section starts
							bool entertest=false;
							if(temp && (temp-conf>0)) entertest=(*(temp-1)=='\n');
							while((!temp && (readcount==READBLOCKSIZE)) || (temp && (!entertest))){ // read file until '[' found or EOF reached
								if(!temp){
									filestart=false; // reading is definitely not from the filestart anymore
									fileposition+=posinc;
									file->file.seekg(fileposition);
									if(conf[READBLOCKSIZE-1]=='\n') entertest=true; else entertest=false;
									file->file.read(conf,READBLOCKSIZE); // read in READBLOCKSIZE byte from file (or less if EOF is reached)
									file->file.clear();
									readcount=file->file.gcount();
									conf[readcount]='\0';
									posinc=readcount;
									temp = strchr(conf,'[');
									entertest&=(temp==conf);
									if(temp-conf>0) entertest=(*(temp-1)=='\n');
								}
								if(temp && !entertest){
									temp = strchr(temp+1,'[');
									if(temp) entertest=(*(temp-1)=='\n');
								}
							}
							if(!temp) temp=strchr(temp2,'\0'); // if next section could not be found, output ends at string's end
							readcount=fileposition+(temp-conf)-pos;
							unsigned int included=0;
							// find identical sections within same file and append at end of section (done here to get correct size for output)
							unsigned int moretemp=pos;
							string subtemp=subname; // needed -- otherwise the next line upon not finding another similar section will overwrite subname
							char* more=NULL;
							if(combine_sections) more=GetSection(file,name,moretemp,included,subtemp,false); // combine this files sections, no loop needed - is recursive ...
							moretemp=0;
							if(more){
								moretemp=strlen(more);
								if(!compare_strings(subtemp.c_str(),subname.c_str())){ // only if subnames are matching do we have a match
									delete[] more;
									more=NULL;
								}
							}
							included=0;
							if(included_content) included=strlen(included_content);
							line = new char[readcount+included+moretemp+1]; line[readcount+included+moretemp]='\0';
							file->file.seekg(pos);
							pos+=include_position;
							// append included section to beginning of section
							file->file.read(line+included,readcount);
							file->file.clear();
							if(included_content){ // includes at beginning
								memcpy(line,included_content,included);
								delete[] included_content;
							}
							if(more){ // later sections in same file at end
								memcpy(line+included+readcount,more,moretemp);
								delete[] more;
							}
							break;
						}
					}
					line=temp2;
					if(*line=='[') line--; // just to be save
				} else{
					line=strchr(line+1,'[');
					bool entertest=false;
					while(!line && (readcount==READBLOCKSIZE)){ // read file until '[' found or EOF reached
						filestart=false; // reading is definitely not from the filestart anymore
						fileposition=file->file.tellg(); // set fileposition to current block's start
						if(conf[READBLOCKSIZE-1]=='\n') entertest=true; else entertest=false;
						file->file.read(conf,READBLOCKSIZE); // read in READBLOCKSIZE byte from file (or less if EOF is reached)
						file->file.clear();
						readcount=file->file.gcount();
						conf[readcount]='\0';
						line = strchr(conf,'[');
						entertest&=(line==conf);
					}
					if((line-conf>1) || entertest){ // now that '[' is found, make sure we read in enough text after it
						fileposition+=(line-conf)-1; // new position is two characters before '['
						file->file.seekg(fileposition); // go to position two characters before '[' in file
						file->file.read(conf,READBLOCKSIZE); // read in READBLOCKSIZE byte from file (or less if EOF is reached)
						file->file.clear();
						readcount=file->file.gcount();
						conf[readcount]='\0';
						line = strchr(conf,'[');
					}
				}
			}
		}
		delete[] conf;
	}
	if(include_first_section_pos>0) pos=include_first_section_pos;
	if((included_content) && (!line)){
		subname=subname_search;
		line=included_content;
	} else if(!line) subname=sn_save;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return line;
}

int ConfigReader::GetNrSections(readfile* file, string name, unsigned int pos, bool combine_sections) /// returns number of sections matching uppercase name (file parsing version)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	int nrsections=0;
	char* conf;
	string subname="";
	list<string> names;
	list<string>::iterator i;
	unsigned int position=pos;
	unsigned int inc_pos=0;
	while((conf=GetSection(file,name,position,inc_pos,subname,true,combine_sections))){
		inc_pos=0;
		bool found=false;
		if(combine_sections){
			for(i=names.begin(); i!=names.end(); ++i){
				if(compare_strings(i->c_str(),subname.c_str())){
					found=true;
					break;
				}
			}
		}
		if(!found){
			names.push_back(subname);
			nrsections++;
		}
		subname="";
		delete[] conf;
	}
	names.clear();
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return nrsections;
}

bool ConfigReader::do_int_number_magic(int &i, const char* item, const char* value, char* conf)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double d;
	bool good;
	good=do_number_magic(d,item,value,conf);
	i=(int)qround(d);
	if(fabs(i-d)>EPS){
		cout << "WARNING: <" << item << "> should be an integer.\n";
	}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return good;
}

char** ConfigReader::string2UPN(const char* value, int &items)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
// maximum number of items (including operators) in reverse polish notation (umgekehrte polnische notation - UPN) is length of string
	items=0;
	const char math_key[]="+-*/(eE\0"; // ")" not in list on purpose (each "(" should have corresponding ")" which will be searched for later)
	const char numbers[]="0123456789";
	char** storage = new char*[strlen(value)];
	char** outsourced;
	int i,j;
	char *temp, *temp2;
	char operation='\0';
	char operation2='\0';

	char* old_pos=(char*)value;
	char* pos=strpbrk(old_pos,math_key);

	while(pos!=NULL){
		switch(*pos){
			case '+':
			case '-':
				if(pos-old_pos>0){ // may be zero after ()
					storage[items]=stringNcopy(old_pos,pos-old_pos);
					items++;
				}
				if(operation!='\0'){
					storage[items]=new char[2];
					storage[items][0]=operation;
					storage[items][1]='\0';
					items++;
				}
				if(operation2!='\0'){
					storage[items]=new char[2];
					storage[items][0]=operation2;
					storage[items][1]='\0';
					items++;
					operation2='\0';
				}
				operation=*pos;
				old_pos=pos+1;
				break;
			case '*':
			case '/':
				if(pos-old_pos>0){ // may be zero after ()
					storage[items]=stringNcopy(old_pos,pos-old_pos);
					items++;
				}
				if(operation!='\0'){
					if((operation=='*') || (operation=='/')){
						storage[items]=new char[2];
						storage[items][0]=operation;
						storage[items][1]='\0';
						items++;
					} else operation2=operation;
				}
				operation=*pos;
				old_pos=pos+1;
				break;
			case 'e':
			case 'E':
				if((strchr(numbers,*(pos-1)) && strchr(numbers,*(pos+1))) || (strchr(numbers,*(pos-1)) && (((*(pos+1)=='+')) || (*(pos+1)=='-')))){
					if(pos-old_pos>0){ // may be zero after ()
						storage[items]=stringNcopy(old_pos,pos-old_pos);
						items++;
					}
					if(operation!='\0'){
						if((operation=='e') || (operation=='E')){
							storage[items]=new char[2];
							storage[items][0]=operation;
							storage[items][1]='\0';
							items++;
						} else operation2=operation;
					}
					operation='e';
					if(*(pos+1)=='-') operation='i';
					if(((*(pos+1)=='+')) || (*(pos+1)=='-')) pos++;
					old_pos=pos+1;
				}
				break;
			case '(': // "outsource" to self
				// need to find corresponding ")" first ...
				i=0;
				temp=strchr(pos+1,')');
				if(temp==NULL){
					cout << "Could not evaluate \"" << value << "\", missing \")\".\n";
					exit(4);
				}
				temp2=strchr(pos+1,'(');
				while((temp2<temp) && (temp2!=NULL)){ // if there are any '('s are there any more? if yes, how many?
					i++;
					temp2=strchr(temp2+1,'(');
				}
				if (i>0){ // i more nested levels found
					temp2=strchr(temp+1,')');
					while((temp2!=NULL) && (i>0)){
						i--;
						if(i==0){ // found corresponding ")"
							temp=temp2;
							break;
						}
						temp2=strchr(temp2+1,')');
					}
					if(i>0){
						cout << "Could not evaluate \"" << value << "\", missing \")\".\n";
						exit(4);
					}
				}
				char* outsource_temp=stringNcopy(pos+1,temp-pos-1);
				outsourced=string2UPN(outsource_temp,j);
				delete[] outsource_temp;
				for(i=0; i<j; i++){
					storage[items]=outsourced[i];
					delete[] outsourced[i];
					items++;
				}
				delete[] outsourced;
				pos=temp;
				old_pos=pos+1;
				break;
		}
		pos=strpbrk(pos+1,math_key);
	}

	// finish up and write last element and operations (if existing, respectively)
	if(strlen(old_pos)>0){
		storage[items]=stringNcopy(old_pos,strlen(old_pos));
		items++;
	}
	if(operation!='\0'){
		storage[items]=new char[2];
		storage[items][0]=operation;
		storage[items][1]='\0';
		items++;
	}
	if(operation2!='\0'){
		storage[items]=new char[2];
		storage[items][0]=operation2;
		storage[items][1]='\0';
		items++;
	}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return storage;
}

double ConfigReader::UPN2double(char** upn, int items, char* conf, bool &success)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	success=true;
	recursion_count++;
	if(recursion_count>RECURSION_LIMIT){
		cout << "\nToo many recursive value lookups - who wants to ride the infinity loop?\n";
		exit(4);
	}
	double d=0.0;
	int i,j;
	if(items==0){ // An error that should not happen ...
		cout << "Cannot calculate from an empty string. I could make up values but that may not be desired ...\n";
		exit(5);
	}
	if(items==1){ // variable name lookup case
		if((toupper(upn[0][0])=='P') && (toupper(upn[0][1])=='I')){
			d = pi;
		} else{
			Vec3 v;
			char c[2];
			c[1]='\0';
			if(upn[0][strlen(upn[0])-1]==']'){ // last char is ']'
				if(strrchr(upn[0],'[')-upn[0]==(int)strlen(upn[0])-3){
					if(upn[0][strlen(upn[0])-4]==']'){
						if(upn[0][strlen(upn[0])-6]=='['){
							char* temp=stringNcopy(upn[0],strlen(upn[0])-6);
							double* tupel=get_MxN_tupel(temp,conf,3,3,NULL);
							delete[] temp;
							c[0]=upn[0][strlen(upn[0])-5];
							i=atoi(c);
							c[0]=upn[0][strlen(upn[0])-2];
							j=atoi(c);
							d=tupel[i*3+j];
							delete[] tupel;
						}
					} else{
						char* temp=stringNcopy(upn[0],strlen(upn[0])-3);
						double* tupel=get_MxN_tupel(temp,conf,3,1,NULL);
						delete[] temp;
						c[0]=upn[0][strlen(upn[0])-2];
						i=atoi(c);
						d=tupel[i];
						delete[] tupel;
					}
				} else{
					cout << "Vector/Matrix index cannot be more than one character.\n";
					exit(6);
				}
			} else{
				char* line;
				if((line=GetItemValue(conf,upn[0]))){
					delete[] line;
					success=SetParam(d,upn[0],conf,true);
				} else success=SetParam(d,upn[0],constant_block,true);
			}
		}
	} else{
		double* storage = new double[items]; // storage for values
		int pos=0;
		for(i=0; i<items; i++){
			if(upn[i][0]=='+'){
				if (pos<1){
					cout << "Error evaluating calculation.\n"; // This error should not happen, but you never know ...
					success=false;
					break;
				}
				if (pos>1) storage[pos-2]+=storage[pos-1]; // else storage[0] stays positive
				pos--;
			} else{
				if(upn[i][0]=='-'){
					if (pos<1){
						cout << "Error evaluating calculation.\n";
						success=false;
						break;
					}
					if (pos>1) storage[pos-2]-=storage[pos-1]; else storage[pos-1]*=-1.0;
					pos--;
				} else{
					if(upn[i][0]=='*'){
						if (pos<2){
							cout << "Error evaluating calculation.\n";
							success=false;
							break;
						}
						storage[pos-2]*=storage[pos-1];
						pos--;
					} else{
						if(upn[i][0]=='/'){
							if (pos<2){
								cout << "Error evaluating calculation.\n";
								success=false;
								break;
							}
							storage[pos-2]/=storage[pos-1];
							pos--;
						} else{
							if((upn[i][0]=='e') || (upn[i][0]=='E')){ // exponential
								if (pos<2){
									cout << "Error evaluating calculation.\n";
									success=false;
									break;
								}
								storage[pos-2]*=pow(10.0,storage[pos-1]);
								pos--;
							} else{
								if((upn[i][0]=='i') || (upn[i][0]=='I')){ // inverse exponential ...
									if (pos<2){
										cout << "Error evaluating calculation.\n";
										success=false;
										break;
									}
									storage[pos-2]/=pow(10.0,storage[pos-1]);
									pos--;
								} else{
									success=do_number_magic(storage[pos],"__evaluation__",upn[i],conf);
									pos++;
								}
							}
						}
					}
				}
			}
		}
		d=storage[0];
		delete[] storage;
	}
	recursion_count--;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return d;
}

bool ConfigReader::do_number_magic(double &d, const char* item, const char* value, char* conf)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	bool good=true;
	char fp_number_key[]="1234567890.";

	int value_length=strlen(value);
	int i=strspn(value,fp_number_key);
	if(i==value_length){
		good=!from_string(d,value);
	} else{
		int upn_items;
		char** upn=string2UPN(value,upn_items);
#if DEBUG_LEVEL>2
		if(upn_items>1){
			if(!compare_strings(item,"__evaluation__")){
				for(i=0; i<recursion_count; i++) cout << "\t";
				cout << "<" << item << "> : ";
			}
			cout << "UPN = |";
			for(i=0; i<upn_items; i++) cout << upn[i] << "|";
			cout << "\n";
			for(i=0; i<recursion_count; i++) cout << "\t";
		}
#endif
		// Now all we need to do is actually calculate the value ...
		d=UPN2double(upn,upn_items,conf,good);
#if DEBUG_LEVEL>2
		if(upn_items>1){
			for(i=0; i<recursion_count; i++) cout << "\t";
			cout << "= " << d << "\n";
		}
#endif
		// clean up
		for(int ii=0; ii<upn_items; ii++) delete[] upn[ii];
		delete[] upn;
	}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return good;
}

void ConfigReader::SetParam(string &s, const char* item, char* conf,string default_value)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	char* value_string=GetItemValue(conf,item,true);
	if (value_string!=NULL){
		s=value_string;
		delete[] value_string;
	} else{
#if DEBUG_LEVEL>0
		if(default_value!="") cout << "Parameter <" << item << "> empty, using default (" << default_value << ")\n";
#endif
		s=default_value;
	}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(string &s, const char* item, char* conf)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	char* value_string=GetItemValue(conf,item,true);
	if (value_string!=NULL){
		s=value_string;
		delete[] value_string;
	} else{
		cout << "Parameter <" << item << "> not set.\n";
		exit(3);
	}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(bool &b, const char* item, char* conf,bool default_value)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	int i;
	SetParam(i,item,conf,int(default_value));
	if(i>0) b=true; else b=false;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(bool &b, const char* item, char* conf)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	int i;
	SetParam(i,item,conf);
	if(i>0) b=true; else b=false;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(unsigned int &i, const char* item, char* conf, int default_value)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	int si;
	SetParam(si,item,conf,default_value);
	if(si<0){
		cout << "Expecting non-negative number for <" << item << ">.\n";
		exit(3);
	} else i=si;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(int &i, const char* item, char* conf, int default_value)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	char* value_string=GetItemValue(conf,item);
	if (value_string!=NULL){
		if(!do_int_number_magic(i,item,value_string,conf)){
			cout << "Cannot read parameter <" << item << ">. \"" << value_string << "\" cannot be evaluated.\n";
			exit(3);
		}
		delete[] value_string;
	} else{
#if DEBUG_LEVEL>0
		cout << "Parameter <" << item << "> not set, using default (" << default_value << ")\n";
#endif
		i=default_value;
	}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(unsigned int &i, const char* item, char* conf)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	int si;
	SetParam(si,item,conf);
	if(si<0){
		cout << "Expecting non-negative number for <" << item << ">.\n";
		exit(3);
	} else i=si;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(int &i, const char* item, char* conf)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	char* value_string=GetItemValue(conf,item);
	if (value_string!=NULL){
		if(!do_int_number_magic(i,item,value_string,conf)){
			cout << "Cannot read parameter <" << item << ">. \"" << value_string << "\" cannot be evaluated.\n";
			exit(3);
		}
		delete[] value_string;
	} else{
		cout << "Parameter <" << item << "> not set.\n";
		exit(3);
	}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(double &d, const char* item, char* conf, double default_value, bool quiet)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	char* value_string=GetItemValue(conf,item);
	if (value_string!=NULL){
		if(!do_number_magic(d,item,value_string,conf)){
			cout << "Cannot read parameter <" << item << ">. \"" << value_string << "\" cannot be evaluated.\n";
			exit(3);
		}
		delete[] value_string;
	} else{
#if DEBUG_LEVEL>0
		if(!quiet) cout << "Parameter <" << item << "> not set, using default (" << default_value << ")\n";
#endif
		d=default_value;
	}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

bool ConfigReader::SetParam(double &d, const char* item, char* conf, bool quiet)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	char* value_string=GetItemValue(conf,item);
	if(value_string){
		if(!do_number_magic(d,item,value_string,conf)){
			cout << "Cannot read parameter <" << item << ">. \"" << value_string << "\" cannot be evaluated.\n";
			exit(3);
		}
		delete[] value_string;
		return true;
	} else{
		if(!quiet){
			cout << "Parameter <" << item << "> not set.\n";
			exit(3);
		}
		return false;
	}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

inline bool could_it_be_tupel(char* tupel_string, int minimum_length) ///  needs properly whitespace-removed string
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	bool now_could_it=true;
	int last = strlen(tupel_string)-1;
	if(last<minimum_length-1) now_could_it=false;
		 // let's not be picky what constitutes brackets
		else if(!((((tupel_string[0]=='(') && (tupel_string[last]==')')) || ((tupel_string[0]=='{') && (tupel_string[last]=='}'))) || (((tupel_string[0]=='[') && (tupel_string[last]==']')) || ((tupel_string[0]=='<') && (tupel_string[last]==']'))))) now_could_it=false;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return now_could_it;
}

double* ConfigReader::get_flex_tupel(const char* item, char* conf, char** item_names, unsigned int* element_count, unsigned int nr_elements)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	char* values=GetItemValue(conf,item);
	char* value_string=values;
	double* result_tupel = NULL;
	char *pos;
	unsigned int i,j;
	bool element_tupels=(item_names && element_count && (nr_elements>0));
	if(value_string){ // Now that we got ourselves a nice value string ...
		// guess how many entities we'll need is string length+2
		unsigned int howmany=strlen(value_string)+2;
		result_tupel=new double[howmany];
		if(strpbrk(value_string,"({[<")) value_string=strpbrk(value_string,"({[<")+1;
		if(strpbrk(value_string,",;|)}]>")) pos=strpbrk(value_string,",;|)}]>"); else pos=value_string+strlen(value_string);
		i=0; // count position
		result_tupel[i]=0;
		bool range=false;
		j=1; // running nr (first element 0 is first count)
		while(pos!=NULL){
			if(pos>value_string){ // make sure empty entries are represented as such
				char* temp=stringNcopy(value_string,pos-value_string);
				bool part_of_list=false;
				if(!do_number_magic(result_tupel[j],item,temp,conf)){
					if(!element_tupels){
						if((compare_strings(temp,"...") || compare_strings(temp,"..")) && ((j>0) && (!range || ((j>0) && (fabs(qround(result_tupel[j-1])-result_tupel[j-1])<EPS))))){
							part_of_list=true;
						} else{
							if((compare_strings(temp,"NaN") || compare_strings(temp,"-NaN")) || (compare_strings(temp,"inf") || compare_strings(temp,"-inf"))){
								if(compare_strings(temp,"NaN") || compare_strings(temp,"-NaN")) result_tupel[j]=sin(1.0/0.0); else result_tupel[j]=1.0/0.0;
							} else{
								cout << "ERROR: <" << item << "> is not a valid flex tupel (e.g. [1,...,a|1,...,b|1,...,c]):\n";
								value_string=strstr(values,temp);
								if(!value_string) value_string=values;
								cout << item << " = " << values << "\n";
								for(unsigned int i=0; i<(value_string-values)+strlen(item)+3; i++) cout << " ";
								cout << "^\n";
								exit(2);
							}
						}
					} else{
						if(item_names && (nr_elements>0)){
							char* name=temp;
							char* number=NULL;
							int nr=-1;
							// find last occurence of ":" in (group) element name (is ignored if only looking for element names)
							char* temp2=strrchr(name,':');
							if(temp2 && element_count){
								name=stringNcopy(temp,temp2-temp);
								number=stringNcopy(temp2+1,strlen(temp)-(temp2-temp)-1);
							}
							if(element_count && number){
								if(!do_int_number_magic(nr,item,number,conf)){
									cout << "ERROR: Could not determine element number.\nSyntax for special elements: <item name>:<element number>\n.";
									exit(3);
								}
								if(nr<1){
									cout << "ERROR: Element number needs to be integer >= 1.\n";
									exit(3);
								}
							}
							unsigned int k=0;
							bool found=false;
							for(k=0; k<nr_elements; k++){
								if(compare_strings(name,item_names[k])){
									found=true;
									break;
								}
								if(number){
									string test2=name;
									test2+=":";
									test2+=number;
									if(compare_strings(test2.c_str(),item_names[k])){
										delete[] name;
										delete[] number;
										name=stringNcopy(test2.c_str(),test2.length());
										number=NULL;
										found=true;
										break;
									}
								}
							}
							if(found){
								double element_nr_fraction=0.0;
								if(element_count && number){
									if(nr>(int)element_count[k]){
										cout << "ERROR: For <" << name << ">, only up to " << element_count[k] << " items are specifiable.\n";
										exit(3);
									}
									// element_nr_fraction is in [1/element count,element count/(element count+1.0)]
									element_nr_fraction=(double)(nr)/(double)(element_count[k]+1.0);
								}
								if(temp2 && number){
									delete[] name;
									delete[] number;
								}
								// negative numbers indicate successful item lookups
								result_tupel[j]=-1.0*(double)(k+element_nr_fraction);
							} else{
								cout << "ERROR: For <" << item << ">, could not attribute <" << temp << "> to anything.\n";
								exit(3);
							}
						} else{
							cout << "ERROR: No elements specified.\n";
							exit(2);
						}
					}
				} else{
					if(element_tupels && ((result_tupel[j]<1.0) || (fabs(qround(result_tupel[j])-result_tupel[j])>EPS))){
						cout << "ERROR: <" << item << "> element tupels need to be integers >0.\n";
						exit(4);
					}
				}
				delete[] temp;
				if(part_of_list){
					range=true;
				} else{
					if(range){
						if(fabs(qround(result_tupel[j])-result_tupel[j])>EPS){
							cout << "<" << item << "> is not a valid flex tupel (e.g. [1,...,a|1,...,b|1,...,c]), \"...\" needs to be followed by an integer.\n";
							exit(2);
						}
						unsigned int needed=abs((int)result_tupel[j]-(int)result_tupel[j-1]);
						double* new_values=new double[howmany+needed];
						for(unsigned int k=0; k<j; k++) new_values[k]=result_tupel[k];
						for(unsigned int k=1; k<needed; k++) new_values[j+k-1]=result_tupel[j-1]+(double)k*((int)needed/((int)result_tupel[j]-(int)result_tupel[j-1]));
						new_values[i]+=needed;
						new_values[j+needed-1]=result_tupel[j];
						j+=needed;
						delete[] result_tupel;
						result_tupel=new_values;
						howmany+=needed;
						range=false;
					} else{
						result_tupel[i]+=1;
						j++;
					}
				}
			}
			if(*pos=='|'){ // "|" is only divider for ranges
				if(range){
					cout << "<" << item << "> is not a valid flex tupel (e.g. [1,...,a|1,...,b|1,...,c]), \"...\" needs to be followed by an integer.\n";
					exit(2);
				}
				i=j;
				result_tupel[i]=0;
				j++;
			}
			if(pos+1<value_string+strlen(value_string)){
				value_string=pos+1;
				if(strpbrk(value_string,",;|)}]>")) pos=strpbrk(value_string,",;|)}]>"); else pos=value_string+strlen(value_string);
			} else pos=NULL;
		}
		result_tupel[j]=-1;
		delete[] values;
	} else{
#if DEBUG_LEVEL>2
		cout << "Parameter <" << item << "> not set.\n"; // don't bail out here, take care of that later ...
#endif
	}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return result_tupel;
}

double* ConfigReader::get_MxN_tupel(const char* item, char* conf, int m, int n, double* default_value, bool warn)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double* result_tupel = new double[n*m];
	char* values=GetItemValue(conf,item);
	char* value_string=values;
	char* temp;
	char *pos;
	int i,j;
	if(value_string!=NULL){ // Now that we got ourselves a nice value string ...
// shortest MxN tupel is (1,...,m;1,...,m;...n times), which has 2 brackets, m-1 commas, n-1 semicolons, and m*n at least single digit numbers
		if(could_it_be_tupel(value_string,m+n+m*n)){ // 2-1-1 = 0 ...
			value_string=strpbrk(value_string,"({[<")+1;
			pos=strpbrk(value_string,",;|");
			if(!pos && (m==1) && (n==1)){ // only one value
				temp=stringNcopy(value_string,strchr(value_string,'\0')-value_string-1);
				do_number_magic(result_tupel[0],item,temp,conf);
				delete[] temp;
			} else{
				i=0; // down
				j=0; // across
				while(pos){
					if((*pos==',') || ((n==1) && ((*pos==';') || (*pos=='|')))){
						temp=stringNcopy(value_string,pos-value_string);
						do_number_magic(result_tupel[i*m+j],item,temp,conf);
						delete[] temp;
						j++;
						if(j==m-1){ // end of m-vector (can end in ";" or closing bracket -- because we already checked there must be at least a bracket still left)
							if(i<n-1){
								temp=stringNcopy(pos+1,strpbrk(pos+1,";|")-pos-1);
								do_number_magic(result_tupel[i*m+j],item,temp,conf);
								delete[] temp;
							} else{
								temp=stringNcopy(pos+1,strchr(pos,'\0')-pos-2);
								do_number_magic(result_tupel[i*m+j],item,temp,conf);
								delete[] temp;
							}
						}
					} else if((*pos==';') || (*pos=='|')){
							j=0;
							i++;
						}
					value_string=pos+1;
					pos=strpbrk(pos+1,",;|");
				}
				if((i!=n-1) || (j!=m-1)){
					cout << i << ", " << j << "\n";
					cout << "<" << item << "> is not a valid " << m << "*" << n << " tupel, found " << i*m+j+1 << " elements.\n";
					exit(2);
				}
			}
		} else{
			cout << "<" << item << "> = " << value_string << " is not a " << m << "*" << n << " tupel (e.g. [1,...,m;1,...,m;...n times]).\n";
			exit(2);
		}
		delete[] values;
	} else{
		if(default_value!=NULL){
#if DEBUG_LEVEL>0
			cout << "Parameter <" << item << "> not set, using default ";
			if(n*m<=9){
				cout << "[";
				for(i=0; i<n; i++){
					for(j=0; j<m; j++){
						cout << default_value[i*m+j];
						if(j+1<m) cout << ", ";
					}
					if(i+1<n) cout << "; ";
				}
				cout << "].\n";
			} else cout << m << "*" << n << " tupel (read documentation).\n";
#endif
			memcpy(result_tupel,default_value,sizeof(double)*n*m);
		} else{
			if(warn) cout << "Parameter <" << item << "> not set.\n";
			delete[] result_tupel;
			return NULL;
		}
	}
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	return result_tupel;
}

void ConfigReader::SetParam(double* v, int len, const char* item, char* conf, double* default_value)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double* result_vec=get_MxN_tupel(item,conf,len,1,default_value);
	for(int i=0; i<len; i++) v[i]=result_vec[i];
	delete[] result_vec;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(int* v, int len, const char* item, char* conf, int* default_value)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double* def_val=new double[len];
	for(int i=0; i<len; i++) def_val[i]=default_value[i];
	double* result_vec=get_MxN_tupel(item,conf,len,1,def_val);
	for(int i=0; i<len; i++){
		v[i]=(int)qround(result_vec[i]);
		if(fabs(v[i]-result_vec[i])>EPS){
			cout << "WARNING: <" << item << "> should be an integer vector.\n";
		}
	}
	delete[] result_vec;
	delete[] def_val;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(bool* v, int len, const char* item, char* conf, bool* default_value)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double* def_val=new double[len];
	for(int i=0; i<len; i++) def_val[i]=default_value[i];
	double* result_vec=get_MxN_tupel(item,conf,len,1,def_val);
	for(int i=0; i<len; i++){
		v[i]=!((int)qround(result_vec[i])==0);
		if(fabs(v[i]-result_vec[i])>EPS){
			cout << "WARNING: <" << item << "> should be a boolean vector.\n";
		}
	}
	delete[] result_vec;
	delete[] def_val;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(double* v, int len, const char* item, char* conf)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double* result_vec=get_MxN_tupel(item,conf,len,1,NULL);
	for(int i=0; i<len; i++) v[i]=result_vec[i];
	delete[] result_vec;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(int* v, int len, const char* item, char* conf)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double* result_vec=get_MxN_tupel(item,conf,len,1,NULL);
	for(int i=0; i<len; i++){
		v[i]=(int)qround(result_vec[i]);
		if(fabs(v[i]-result_vec[i])>EPS){
			cout << "WARNING: <" << item << "> should be an integer vector.\n";
		}
	}
	delete[] result_vec;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(bool* v, int len, const char* item, char* conf)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double* result_vec=get_MxN_tupel(item,conf,len,1,NULL);
	for(int i=0; i<len; i++){
		v[i]=!((int)qround(result_vec[i])==0);
		if(fabs(v[i]-result_vec[i])>EPS){
			cout << "WARNING: <" << item << "> should be a boolean vector.\n";
		}
	}
	delete[] result_vec;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(double v[3], const char* item, char* conf, Vec3 default_value)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double* result_vec=get_MxN_tupel(item,conf,3,1,default_value.vec);
	for(int i=0; i<3; i++) v[i]=result_vec[i];
	delete[] result_vec;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(int v[3], const char* item, char* conf, Vec3 default_value)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double* result_vec=get_MxN_tupel(item,conf,3,1,default_value.vec);
	for(int i=0; i<3; i++){
		v[i]=(int)qround(result_vec[i]);
		if(fabs(v[i]-result_vec[i])>EPS){
			cout << "WARNING: <" << item << "> should be an integer vector.\n";
		}
	}
	delete[] result_vec;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(bool v[3], const char* item, char* conf, Vec3 default_value)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double* result_vec=get_MxN_tupel(item,conf,3,1,&default_value.vec[0]);
	for(int i=0; i<3; i++){
		v[i]=!((int)qround(result_vec[i])==0);
		if(fabs(v[i]-result_vec[i])>EPS){
			cout << "WARNING: <" << item << "> should be a boolean vector.\n";
		}
	}
	delete[] result_vec;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(double v[3], const char* item, char* conf)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double* result_vec=get_MxN_tupel(item,conf,3,1,NULL);
	for(int i=0; i<3; i++) v[i]=result_vec[i];
	delete[] result_vec;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(int v[3], const char* item, char* conf)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double* result_vec=get_MxN_tupel(item,conf,3,1,NULL);
	for(int i=0; i<3; i++){
		v[i]=(int)qround(result_vec[i]);
		if(fabs(v[i]-result_vec[i])>EPS){
			cout << "WARNING: <" << item << "> should be an integer vector.\n";
		}
	}
	delete[] result_vec;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(bool v[3], const char* item, char* conf)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double* result_vec=get_MxN_tupel(item,conf,3,1,NULL);
	for(int i=0; i<3; i++){
		v[i]=!((int)qround(result_vec[i])==0);
		if(fabs(v[i]-result_vec[i])>EPS){
			cout << "WARNING: <" << item << "> should be a boolean vector.\n";
		}
	}
	delete[] result_vec;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(double A[3][3], const char* item, char* conf, Mat33 default_value)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double* result_vec=get_MxN_tupel(item,conf,3,3,&default_value.mat[0][0]);
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++) A[i][j]=result_vec[i*3+j];
	delete[] result_vec;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

void ConfigReader::SetParam(Mat33 &A, const char* item, char* conf, Mat33 default_value)
{
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
	double* result_vec=get_MxN_tupel(item,conf,3,3,&default_value.mat[0][0]);
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++) A.mat[i][j]=result_vec[i*3+j];
	delete[] result_vec;
#if DEBUG_LEVEL>3
	cout << "*** " << __FUNCTION__ << " *** (" << __FILE__ << ": " << __LINE__ << ")\n";
#endif
}

