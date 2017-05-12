#include "simbox.cpp"
#include "noPeriodicity.cpp"
#include "sweepLtR.cpp"


int main(int argc, char *argv[]){

//************************* Set Parameters **************************************
// Default parameters
  int sideLength = 10; // dimension in x and y dimension
  int SLZ = 10;        // distance between electrodes


  switch(argc){ //fallthrough means all previous args also get assigned
    case 3: SLZ = atoi( argv[2] );
    case 2: sideLength = atoi( argv[1] );
  }


  // Make a lattice object:
  Simbox sample(sideLength, SLZ, "input.in");
  cout << get_energy(sample) << endl;

  


}
