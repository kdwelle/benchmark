
#include <string>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstring>
#include <limits>
#include "simbox.hpp"

using namespace std;


Simbox::Simbox(int sideLength, int SLZ, string filename): sideLength(sideLength), SLZ(SLZ)
{
  gam = 1;
  readInput(filename);
  initialize();
}

int Simbox::readInput(string filename){
  ifstream fin;
  fin.open(filename); // open a file
  if (!fin.good()){
    return 1; // exit if file not found
  }
  //read first line (how many objects)
  char buf[100];
  fin.getline(buf, 100);
  numObjs = atoi(strtok(buf, " "));
  // read rest of the file

  vector<float> itemp;
  itemp.resize(3);
  position.resize(numObjs,itemp);
  charge.resize(numObjs);
  int objIndex = 0;

  while (!fin.eof())
  {
    fin.getline(buf, 100); // read an entire line into memory
    // parse the line into blank-delimited tokens
    int n = 0; // a for-loop index
    int maxTokens = 5;
    const char* token[maxTokens] = {}; // array to store memory addresses of the tokens in buf

    // parse the line
    token[0] = strtok(buf, " "); // first token
    if (token[0]){ //if there are any tokens
      for (n = 1; n < maxTokens; n++)
      {
        token[n] = strtok(NULL, " "); // subsequent tokens
        if (!token[n]) break; // no more tokens
      }

    position[objIndex][0] = atof(token[1]);  // x-position
    position[objIndex][1] = atof(token[2]);  // y-position
    position[objIndex][2] = atof(token[3]);  // z-position
    charge[objIndex] = atof(token[4]);       //charge

    objIndex++;
    }
  }

  for(int i=0; i<numObjs; ++i){
    cout << position[i][0] << " " << position[i][1] << " " << position[i][2] << endl;
    cout << "charge is " << charge[i] << endl;
  }
  cout << numObjs << endl;
  if (numObjs == objIndex){
    return 0; //whee everythig's good
  }
  else{
    return 2; //mismatch between declared number of objects and actual number
  }
}

void Simbox::initialize(){
  inf = std::numeric_limits<double>::infinity();

  numPairs = 0;
  vector<int> pairindTemp;
  vector<int> pairind2Temp(numObjs);
  pairind2.reserve(numObjs);
  for(int i = 0; i < numObjs; i++){ // Build list of ion pairs
    pairind2.push_back(pairind2Temp);
    for(int j = i; j < numObjs; j++){
      pairindTemp.clear();
      pairindTemp.push_back(i);
      pairindTemp.push_back(j);
      pairind.push_back(pairindTemp);
      pairind2[i][j] = numPairs;
      numPairs++;
    }
  }

  for(int i = 0 ; i < numObjs-1 ; i++){ // mirror pairind2
    for(int j = i+1 ; j < numObjs ; j++){
      pairind2[j][i] = pairind2[i][j];
    }
  }

  position.reserve(numObjs);
  vector<double> posTemp(3);

  drpair.reserve(numPairs); //initialize drpairs
  for (int i=0; i<numPairs; ++i){
    drpair.push_back(posTemp);
  }

  get_drpair0();

}

void Simbox::get_drpair0(){
  vector<double> dr(3);
  int ind1,ind2;
  for(int i = 0 ; i < numPairs ; i++){
    ind1=pairind[i][0];
    ind2=pairind[i][1];
    if (charge[ind1] && charge[ind2]){ //if they both have charge
      for(int k = 0 ; k < 3 ; k++){
        dr[k]=position[ind1][k]-position[ind2][k];
        if(k<2){
        	if(dr[k] > sideLength/2.) dr[k]-=sideLength;
        	if(dr[k] < -sideLength/2.) dr[k]+=sideLength;
        }else{ //nearest neighbor convention for z-direction
          if(dr[k] > SLZ) dr[k]-=SLZ*2;
          if(dr[k] < -SLZ) dr[k]+=SLZ*2;
        }
      }
    }else{ //otherwise they are "infinitely far apart"
      for(int k = 0 ; k < 3 ; k++){
        dr[k] = inf;
      }
    }
    drpair[i] = dr;
    cout << dr[0] << " " << dr[1] << " " << dr[2] << " " << endl;
  }
}

void Simbox::get_drpair1(int ind1){

  vector<double> dr(3);
  int pairIndex;
  for(int ind2 = 0 ; ind2 < numObjs ; ind2++){
    pairIndex=pairind2[ind1][ind2];
    if (charge[ind1] && charge[ind2]){ //if they both have charge
      for(int k = 0 ; k < 3 ; k++){
      	dr[k]=position[ind1][k]-position[ind2][k];
      	if(k<2){ //nearest neighbor convention
      	  if(dr[k] > sideLength/2.) dr[k]-=sideLength;
      	  if(dr[k] < -sideLength/2.) dr[k]+=sideLength;
      	}else{ //nearest neighbor convention for z-direction
          if(dr[k] > SLZ) dr[k]-=SLZ*2;
          if(dr[k] < -SLZ) dr[k]+=SLZ*2;
        }
      }
    }else{
      for(int k = 0 ; k < 3 ; k++){
        dr[k] = inf;
      }
    }
    drpair[pairIndex]=dr;
  }
}

void Simbox::translate_particle(int index, float dx, float dy, float dz){
  // translate particle index by dx,dy,dz
  position[index][0] += dx;
  position[index][1] += dy;
  position[index][2] += dz;

  get_drpair1(index);
}

void Simbox::set_position(int index, float x, float y, float z){
  // move particle index to x,y,z
  position[index][0] = x;
  position[index][1] = y;
  position[index][2] = z;

  get_drpair1(index);
}

void Simbox::change_charge(int index, float dq){
  charge[index] += dq;
}

void Simbox::set_charge(int index, float newq){
  charge[index] = newq;
}
