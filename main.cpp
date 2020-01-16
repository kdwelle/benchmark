#include "periodicImages.cpp"
#include "simbox.cpp"
// #include "noPeriodicity.cpp"
// #include "twoDSlab.cpp"
#include "ewald.cpp"
#include <stdio.h>
#include <string>


int main(int argc, char *argv[]){

//************************* Set Parameters **************************************
// Default parameters
  double SLX = 68.8558; 
  double SLY = 31.0000;
  double SLZ = 31.0000; //image plane dimension if image charges
  int imageCharges = 0;

  // set some binning/sampling parameters
  float zbegin = -34;
  int nbins = 100;
  float dz = (SLX)/nbins;
  float x1begin = -0.21;
  float x2begin = -0.21;
  int nbinsx = 2;
  float dx1 = (SLY)/nbinsx;
  float dx2 = (SLZ)/nbinsx;
  int nTimesteps = 2;


  switch(argc){ //fallthrough means all previous args also get assigned
    case 5: imageCharges = atoi( argv[4] );
    case 4: SLZ = atoi( argv[3] );
    case 3: SLY = atoi( argv[2] );
    case 2: SLX = atoi( argv[1] );
  }

  // cout << "image charges: " << imageCharges << endl;


  for (int i=0;i<=nTimesteps;++i){
    // Make a lattice object:
    char number[3];
    sprintf(number,"%.3d",i);
    std::string number2 = number;
    std::string name = "input/DumpOut2-2-";
    name += number2;

    Simbox sample(SLX,SLY,SLZ, imageCharges, name);
    // cout << name << " loaded" << endl;
    // and an image vector object:
    PeriodicImages imageItem = get_images(SLX,SLY,SLZ);
    for (int i=0; i<nbins; ++i){
      double z = zbegin+(i*dz);
      cout << z << " ";
      for (int j=1; j<=nbinsx; ++j){
        for (int k=1; k<=nbinsx; ++k){
          double x1 = x1begin + (j*dx1);
          double x2 = x2begin + (k*dx2);
          cout << get_potential(sample,imageItem, z,x1,x2) << " ";
        }
      }
      cout << endl;
    }
  }

}
