#include "simbox.cpp"
#include "noPeriodicity.cpp"
#include "sweepLtR.cpp"


int main(int argc, char *argv[]){

//************************* Set Parameters **************************************
// Default parameters
  int sideLength = 10; // dimension in x and y dimension
  int SLZ = 10;        // distance between electrodes
  bool imageChargesExist = true;


  switch(argc){ //fallthrough means all previous args also get assigned
    case 4: imageChargesExist = atoi( argv[3] );
    case 3: SLZ = atoi( argv[2] );
    case 2: sideLength = atoi( argv[1] );
  }


  // Make a lattice object:
  Simbox sample(sideLength, SLZ, imageChargesExist, "input.in");
  cout << get_energy(sample) << endl;

  for (int i=0; i<100; ++i){
    double z = i*1.0/SLZ;
    cout << z << " ";
    sample.set_position(0,0,0,z);
    cout << get_energy(sample) << endl;
  }


}
