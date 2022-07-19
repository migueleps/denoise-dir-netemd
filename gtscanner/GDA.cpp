#include "GDA.h"

#include "Error.h"

#include <fstream>
#include <sstream>
#include <string>

#include <cstring>

using namespace std;

int constCharToint(const char* value);
vector<string> &splits(const string &s, char delim, vector<string> &elems);
vector<string> splits(const string &s, char delim);

void GDA::newGDD(const char* file){
  GDD gdd1 = readOrbitFile(file);
  gdds.push_back(gdd1);
}

GDD GDA::readOrbitFile(const char* file){
  GDD this_gdd;	
  
  std::ifstream infile(file);
  
  vector<string> tokens;
  vector<string> freqs;
  vector<string> freq;
	
  std::string line;
  while (std::getline(infile, line)){
	std::map<int, float> orbit_vector;
	unsigned int k;
	unsigned int d;
	
	unsigned int N;
	float T = 0;
	
	tokens = splits(line.c_str(), '-');
	if(tokens.size() == 2){
	  
	  freqs = splits(tokens[1].c_str(), ',');
	  for(int i = 0; i < freqs.size(); i++){
	    freq = splits(freqs[i].c_str(), ':'); 
	    
	    k = constCharToint(freq[0].c_str());
	    d = constCharToint(freq[1].c_str());
	    
	    orbit_vector[k] = float(d)/float(k);
	    T += orbit_vector[k];
	  }
	  
	  for(map<int, float>::iterator it = orbit_vector.begin(); it != orbit_vector.end(); ++it){
	    orbit_vector[it->first] = it->second/T;
	  }
    } else{
	  orbit_vector[0] = -1;
	}
	this_gdd.push_back(orbit_vector);
  }
  
  return this_gdd;
}

float GDA::computeGDA(int gdd_i1, int gdd_i2){
  std::vector<GDD> gdds = GDA::getGDDs();	
  int max_k;	
	
  GDD gdd1 = gdds[gdd_i1];
  GDD gdd2 = gdds[gdd_i2];

  GDD joint_gdd;

  vector<float> A;
  float this_gda = 0;
  
  if(gdd1.size() != gdd2.size())
	Error::msg("Files do not have the same number of orbits!");
  
  for(int o = 0; o < gdd1.size(); o++){
	joint_gdd.push_back(gdd1[o]);
	
	float dif    = 0;
	float this_D = 0;
	
	int max_k = std::max((gdd1[o].rbegin())->first, (gdd2[o].rbegin())->first);
    if(max_k == 0) continue;
	
	for(map<int, float>::iterator it = gdd1[o].begin(); it != gdd1[o].end(); ++it){
	  //if(it->second != -1){
	    joint_gdd[o][it->first] = pow(it->second - gdd2[o][it->first], 2);
	  
	    dif += joint_gdd[o][it->first];
	    //printf("%d: [%.4f,%.4f]%d: %.5f == %.5f\n", o, gdd1[o][it->first], gdd2[o][it->first], it->first, pow(gdd1[o][it->first] - gdd2[o][it->first], 2), dif);
      //}
	}
	for(map<int, float>::iterator it = gdd2[o].begin(); it != gdd2[o].end(); ++it){
	  if(gdd1[o][it->first] == 0 /*&& it->second != -1*/){
		joint_gdd[o][it->first] = pow(gdd1[o][it->first] - it->second, 2);
	  
	    dif += joint_gdd[o][it->first];
	    //printf("%d: [%.4f,%.4f]%d: %.5f == %.5f\n", o, gdd2[o][it->first], gdd1[o][it->first], it->first, pow(gdd2[o][it->first] - it->second, 2), dif);
      }
	}
	
	//printf("%d--> dif: %.5f\n", o, dif);
	this_D = (1.0/sqrt(2.0)) * pow(dif, 0.5);
	A.push_back(1-this_D);
	
	//printf("<Agreement Orbit %d: %.5f>\n", o, 1-this_D);
  }
  
  for(int o = 0; o < A.size(); o++){
    this_gda += A[o];
  }
  
  return this_gda/A.size();
}

/***********************************/

vector<string> &splits(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> splits(const string &s, char delim) {
    vector<string> elems;
    splits(s, delim, elems);
    return elems;
}

int constCharToint(const char* value){
	stringstream strValue;
	strValue << value;

	unsigned int intValue;
	strValue >> intValue;
	
	return intValue;
}
