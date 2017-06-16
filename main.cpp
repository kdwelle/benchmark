#include "periodicImages.cpp"
#include "simbox.cpp"
// #include "noPeriodicity.cpp"
// #include "twoDSlab.cpp"
#include "ewald.cpp"


int main(int argc, char *argv[]){

//************************* Set Parameters **************************************
// Default parameters
  int sideLength = 10; // dimension in x and y dimension
  int SLZ = 10;        // distance between electrodes
  int imageCharges = 0;

  // float zend = 3.0*SLZ/2.0;


  switch(argc){ //fallthrough means all previous args also get assigned
    case 4: imageCharges = atoi( argv[3] );
    case 3: SLZ = atoi( argv[2] );
    case 2: sideLength = atoi( argv[1] );
  }

  cout << "image charges: " << imageCharges << endl;



  // Make a lattice object:
  Simbox sample(sideLength, SLZ, imageCharges, "input.in");
  // and an image vector object:
  PeriodicImages imageItem = get_images(sideLength,SLZ);
  cout << get_energy(sample,imageItem) << endl;

  float zbegin = 0;
  sample.output_analysis();
  for (int i=0; i<97; ++i){
    double z = zbegin+0.2+(i*1.0/SLZ);
    cout << z << " ";
    sample.set_position(0,0,0,z);
    cout << get_energy(sample,imageItem) << endl;
  }
  sample.output_analysis();

}
