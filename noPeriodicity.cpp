
#include <math.h>

using namespace std;

double get_energy(const Simbox& config){
  double eng=0.0;
  double ftemp;
  
  for(int i=0; i<config.numReal; ++i){ //sum over madeleung potentials for real ions
    int ionIndex = config.realCharges[i];
    ftemp = get_mad_potential(config,ionIndex)*config.charge[ionIndex];
    eng += ftemp;
  }
  return eng/2;
}

double get_mad_potential(const Simbox& config, int ionIndex){
  double pot=0.0;
  double ftemp;
  float q1,q2,q;
  double x1,y1,z1;
  int pairIndex;

  double rEnergy = 0.0; //real-space energy
  double r_ij;

  for (int i=0; i<config.numObjs; ++i){ //sum over all ions except ionIndex (not including own periodic images)
    if(i != ionIndex){ //madeleung potential doesn't include own interaction
      pairIndex = config.pairind2[ionIndex][i];
      q1=config.charge[config.pairind[pairIndex][0]]; //0 index is smaller index than 1 index
      q2=config.charge[config.pairind[pairIndex][1]]; //charge of the ions
      if(i < ionIndex){
        q=q1;
      }else if (i > ionIndex){
        q=q2;
      }
      if(q){ //only need if interacting ion has charge
        x1 = config.drpair[pairIndex][0]; //non-periodic distance
        y1 = config.drpair[pairIndex][1];
        z1 = config.drpair[pairIndex][2];
        r_ij = sqrt(x1*x1+y1*y1+z1*z1);
        if(r_ij > 0){
          ftemp=q/r_ij;
          rEnergy += ftemp;
        }
      }
    }
  }
  pot = rEnergy/config.gam;
  return pot;
}
