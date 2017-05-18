#include "simbox.cpp"
#include "noPeriodicity.cpp"
#include "sweepLtR.cpp"


int main(int argc, char *argv[]){

//************************* Set Parameters **************************************
// Default parameters
  int sideLength = 10; // dimension in x and y dimension
  int SLZ = 10;        // distance between electrodes
  int imageChargesExist = 1;

  // float zend = 3.0*SLZ/2.0;


  switch(argc){ //fallthrough means all previous args also get assigned
    case 4: imageChargesExist = atoi( argv[3] );
    case 3: SLZ = atoi( argv[2] );
    case 2: sideLength = atoi( argv[1] );
  }

  cout << "image charges: " << imageChargesExist << endl;

  // float zbegin = (imageChargesExist) ? SLZ/2.0 : 0;

  // Make a lattice object:
  Simbox sample(sideLength, SLZ, imageChargesExist, "input.in");
  cout << get_energy(sample) << endl;

  // for (int i=0; i<99; ++i){
  //   double z = zbegin+0.1+(i*1.0/SLZ);
  //   cout << z << " ";
  //   sample.set_position(0,0,0,z);
  //   cout << get_energy(sample) << endl;
  // }


}
