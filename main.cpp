#include "simbox.cpp"
#include "noPeriodicity.cpp"
int main(int argc, char *argv[]){
//delta E simulation
//************************* Set Parameters **************************************
// Default parameters
  int SideLength = 10; // number of lattice sites in x and y dimension
  int SLZ = 10;        // number of lattice sites between electrodes -- has to be even


  switch(argc){ //fallthrough means all previous args also get assigned
    case 3: SLZ = atoi( argv[2] );
    case 2: SideLength = atoi( argv[1] );
  }


  // Make a lattice object:
  Simbox sample(SideLength, SLZ, "input.in");

}
