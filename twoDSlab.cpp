#include <math.h>
#include <complex>

using namespace std;
double functionF(double,double,double);
double alpha = 6.67; //constant for 2-D slab summation



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
  //gets the madeleung potential of an ion
  double pot=0.0;
  double ftemp;
  float q,q1,q2;
  int dx,dy;
  int x1,y1,z1,pairIndex;

  //Ewald Summation
  double rEnergy = 0.0;
  double r_ij;
  double siEnergy = 0.0;
  double cosTerm;
  double functionF;
  double functiong;
  double ftEnergy = 0.0;
  double h1,h2,h_2,h;
  complex<double> dot;

  for (int i=0; i<config.numObjs; ++i){
    pairIndex = config.pairind2[ionIndex][i];
    q1=config.charge[config.pairind[pairIndex][0]]; //0 index is smaller index than 1 index
    q2=config.charge[config.pairind[pairIndex][1]]; //charge of the ions
    if(i < ionIndex){
      q=q1;
    }else if (i > ionIndex){
      q=q2;
    }
    if(q){
      x1 = config.drpair[i][0]; //non-periodic distance
      y1 = config.drpair[i][1];
      z1 = config.drpair[i][2];
      for (int img = 0; img < config.rCount; img++){ //loop over all periodic images
        dx = x1 + config.rspace[img][0];
        dy = y1 + config.rspace[img][1];
        r_ij = sqrt(dx*dx+dy*dy+z1*z1);
        if(r_ij < config.rsCuttoff && r_ij > 0){
          ftemp=q*erfc(alpha*r_ij)/r_ij;
          rEnergy += ftemp;
        }
      }
      for(int img = 0; img<config.ftCount; ++img){ //FOURIER SPACE Part
        h1=config.fspace[img][0];
        h2=config.fspace[img][1];
        h_2 = h1*h1+h2*h2;
        h=sqrt(h_2);
        if (h_2 > .00001){
          cosTerm = cos(h1*x1+h2*y1);
          functionF = exp(h*z1)*erfc(alpha*z1+(h/(2*alpha))) + exp(-1*h*z1)*erfc(-1*alpha*z1+(h/(2*alpha)));
          functiong = z1*erfc(alpha*z1)+ exp(-1*(z1*alpha*z1*alpha))/(alpha*M_PI);
          ftEnergy += q*M_PI/(config.sideLength*config.sideLength)*cosTerm*functionF-functiong;
        }
      }
    }
  }

  for (int i=0; i<config.numObjs; ++i){
    q1 = config.charge[i];
    if(q1){
      //SELF-INTERACTION
      ftemp = -1*alpha/sqrt(M_PI)*q1*q1;
      siEnergy += ftemp;
    }
  }

  pot = (ftEnergy+siEnergy+rEnergy)/config.gam;
  return pot; //*q1*0.5;
}


// double functionF(double u, double z, double alpha){
//   return (exp(2*u*z)*erfc(alpha*u+(z/alpha)) + exp(-2*u*z)*erfc(alpha*u+(z/alpha)))/(2*u);
//
// }
