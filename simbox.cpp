
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstring>
#include "simbox.hpp"

using namespace std;


Simbox::Simbox(int sideLength, int SLZ, string filename): sideLength(sideLength), SLZ(SLZ)
{
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
  int numObjs = atoi(strtok(buf, " "));
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
// need to set up drpair

}
