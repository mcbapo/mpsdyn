#include "Properties.h"
#include <fstream>

char Properties::keysep='=';

// auxiliary function to strip whitespace from values and names
const string strip(const string str);

Properties::Properties(const char* filename):table_(){
  loadFromFile(filename);
}

Properties::~Properties(){
  // Empty the table
  table_.clear();
}

const string Properties::getProperty(const string& name) const{
  map<string,string>::const_iterator iter=table_.find(name);
  if(iter!=table_.end()){
    return iter->second;
  }
  return "";
}

int Properties::getIntProperty(const string& name) const{
  string val=getProperty(name);
  if(val.empty())
    return -1; // WARNING: Should be exception
  else
    return atoi(val.c_str());
}

double Properties::getDoubleProperty(const string& name) const{
  string val=getProperty(name);
  if(val.empty())
    return -1; // WARNING: Should be exception
  else
    return atof(val.c_str());
}

#include <sstream>
vector<double> Properties::getVectorDoubleProperty(const string& name) const{
  string val=getProperty(name);
  vector<double> ret;
  if(val.empty())
    return ret; // WARNING: Should be exception
  else{
    istringstream input(val,ios_base::in);
    double aux;
    int cnt=0;
    while(!input.eof()){
      input>>aux;
      ret.push_back(aux);
      cnt++;
    }
    return ret;
  }
}

vector<int> Properties::getVectorIntProperty(const string& name) const{
  string val=getProperty(name);
  vector<int> ret;
  if(val.empty())
    return ret; // WARNING: Should be exception
  else{
    istringstream input(val,ios_base::in);
    int aux;
    int cnt=0;
    while(!input.eof()){
      input>>aux;
      ret.push_back(aux);
      cnt++;
    }
    return ret;
  }
}

void Properties::addProperty(const string& name,const string& value){
  // Erase the property, if it existed
  table_.erase(name);
  cout<<"Inserting("<<name<<"-"<<value<<")"<<endl;
  table_.insert(make_pair(name,value));
}
  

#define MAXLINE 240

void Properties::loadFromFile(const char* filename){
  // Parse the file
  ifstream input(filename);
  if(!input.is_open()){
    cout<<"Error: could not open file "<<filename<<" to read properties"
	<<endl;
    exit(1);
  }
  char aux[MAXLINE];
  while(true){
    input.getline(aux,MAXLINE);
    if(input.good()){
      string tmp(aux); // handle it with more practical string
      if(tmp.length()!=0){ // ignore empty line
	if(tmp[0]!='#'){ // Ignore comment lines
	  parseString(tmp);
	}
      }
    }
    else{
      cout<<"No more lines in "<<filename<<endl;
      input.close();
      return;
    }
  }


}

void Properties::loadProperties(int nr,const char* args[]){
  for(int k=0;k<nr;k++){
    cout<<"Parsing arg "<<k<<": "<<args[k]<<endl;
    string tmp(args[k]);
    if(tmp[0]=='-'){ // remove leading dash
      tmp.erase(0,1);
      parseString(tmp);
    }
    else{
      cout<<"Argument "<<args[k]<<"not a valid property: ignoring"<<endl;
    }
  }

}

void Properties::parseString(string tmp){
  // look for separator
  int pos=tmp.find(keysep);
  //cout<<"In "<<tmp<<", = identified at pos "<<pos<<endl;
  if(pos>=0&&pos<tmp.length()){
    string name=tmp.substr(0,pos);
    string value=tmp.substr(pos+1,tmp.length()-pos-1);
    // strip white space at the beginning and the end
    addProperty(strip(name),strip(value));
    //cout<<"Identified key="<<name<<" and value="<<value<<endl;
  }
  else{
    cout<<"Separator "<<keysep<<" not found in line "<<tmp
	<<". Ignoring line."<<endl;
  }
}

const string strip(const string str){
  const char* whites=" \t";
  int pos1=str.find_first_not_of(whites);
  int pos2=str.find_last_not_of(whites);
  return str.substr(pos1,pos2-pos1+1);
}
