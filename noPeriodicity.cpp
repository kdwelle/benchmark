
#include <math.h>

using namespace std;

double get_energy(const Simbox& config){
  double eng=0.0;
  double ftemp;
  float q1,q2;
  int x1,y1,z1;
  double rEnergy = 0.0;
  double r_ij;

  for (int i=0; i<config.numPairs; ++i){
    q1=config.charge[config.pairind[i][0]]; //charge of the ions
    q2=config.charge[config.pairind[i][1]];
    if(q1 && q2){
      x1 = config.drpair[i][0]; //non-periodic distance
      y1 = config.drpair[i][1];
      z1 = config.drpair[i][2];
      r_ij = sqrt(x1*x1+y1*y1+z1*z1);
      if(r_ij > 0){
        ftemp=q1*q2/r_ij;
        rEnergy += ftemp;
      }
    }
  }
  eng = rEnergy/config.gam;
  return eng;
}
