#include <math.h>
#include <complex>

using namespace std;
double alpha = 6.67; //constant for 2-D slab summation


double get_energy(const Simbox& config, const PeriodicImages& imageItem){
  double eng=0.0;
  double ftemp;

  for(int i=0; i<config.numReal; ++i){ //sum over madeleung potentials for real ions
    int ionIndex = config.realCharges[i];
    ftemp = get_mad_potential(config,ionIndex,imageItem)*config.charge[ionIndex];
    eng += ftemp;
  }
  return eng/2;
}


double get_mad_potential(const Simbox& config, int ionIndex, const PeriodicImages& imageItem){
  //gets the madeleung potential of an ion
  double alpha = 3.2/config.rsCuttoff;
  double pot=0.0;
  double ftemp;
  float q1,q2,q;
  double dx,dy,dz;
  double x1,y1,z1;
  int pairIndex;

  //Ewald Summation
  double rEnergy = 0.0; //real-space energy
  double r_ij;
  double siEnergy = 0.0; //self-interaction energy

  // REAL SPACE Part
  for (int i=0; i<config.numObjs; ++i){ //sum over all ions except ionIndex (not including own periodic images)
    pairIndex = config.pairind2[ionIndex][i];
    q1=config.charge[config.pairind[pairIndex][0]]; //0 index is smaller index than 1 index
    q2=config.charge[config.pairind[pairIndex][1]]; //charge of the ions
    if(i < ionIndex){
      q=q1;
    }else{
      q=q2;
    }
    if(q){ //only need if interacting ion has charge
      x1 = config.drpair[pairIndex][0]; //non-periodic distance
      y1 = config.drpair[pairIndex][1];
      z1 = config.drpair[pairIndex][2];
      for (int img = 0; img < imageItem.rCount; img++){ //loop over all periodic images
        dx = x1 + imageItem.rspace[img][0];
        dy = y1 + imageItem.rspace[img][1];
        dz = z1 + imageItem.rspace[img][2];
        r_ij = sqrt(dx*dx+dy*dy+dz*dz);
        if(r_ij > 0 && r_ij < config.rsCuttoff){
          ftemp=q*erfc(sqrt(alpha)*r_ij)/r_ij;
          rEnergy += ftemp;
        }
      }
    }
  }

  // Interactions with self in periodic images
  q1 = config.charge[ionIndex];
  if(q1){ //only need to run if q1 has charge
    siEnergy = -sqrt(alpha/M_PI)*q1*2;  //self-interaction term here for efficiency
    // End self-interaction term

    for (int img = 0; img < imageItem.rCount; img++){ //loop over all real-space periodic images
      dx = imageItem.rspace[img][0];
      dy = imageItem.rspace[img][1];
      dz = imageItem.rspace[img][2];
      if (dx != 0 || dy != 0 || dz != 0){ //don't want dist = 0 interaction
        r_ij = sqrt(dx*dx+dy*dy+dz*dz);
        ftemp=q1*erfc(sqrt(alpha)*r_ij)/r_ij;
        rEnergy += ftemp;
      }
    }
  }

  //FOURIER SPACE Part
  complex<double> rho;
  double ftEnergy = 0.0;
  double k1,k2,k3,k_2;
  complex<double> dot;
  for (int img = 0; img < imageItem.ftCount; img++){ //loop over all periodic images
    k1=imageItem.fspace[img][0];
    k2=imageItem.fspace[img][1];
    k3=imageItem.fspace[img][2];
    k_2 = k1*k1+k2*k2+k3*k3; //k^2
    if (k_2 > .000001){
      rho = complex<double>(0.0,0.0);
      for (int i=0; i<config.numObjs; ++i){
        q2 = config.charge[i];
        if (q2){ //don't need if no charge
          pairIndex = config.pairind2[ionIndex][i];
          dot = complex<double>(0.0,k1*config.drpair[pairIndex][0]+k2*config.drpair[pairIndex][1]+k3*config.drpair[pairIndex][2]); //r_ij dot ik
          rho = exp(dot);
          ftemp = 4*M_PI*q2/k_2*real(rho)*exp(-k_2/(4*alpha));
          ftEnergy += ftemp;
        }
      }
    }
  }
  ftEnergy = ftEnergy/(SideLength*SideLength*SideLength); //factor of 1/V

  pot = (ftEnergy+siEnergy+rEnergy)/gam;
  return pot; //*q1*0.5;
}

PeriodicImages get_images(double sideLength, double SLZ){
  int dimen = 3;  // this is a 3D ewald code so 3 dimensions
  int n_max = 1;  // only need one real image
  int m_max = 5; // seems like 3D slab converges around 5 k-space images
  PeriodicImages imageItem(dimen,n_max,m_max,sideLength,SLZ);
  return imageItem;
}
