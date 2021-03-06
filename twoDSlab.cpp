#include <math.h>

using namespace std;
double alpha = 6.67; //constant for 2-D slab summation
// double alpha = 1000;


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
  double pot=0.0;
  double ftemp;
  float q,q1,q2;
  double dx,dy;
  double x1,y1,z1;
  int pairIndex;

  //Ewald Summation
  double rEnergy = 0.0;  // real-space energy
  double r_ij;
  double siEnergy = 0.0;  //self-interaction energy
  double cosTerm;
  double functionF;
  double functiong;
  double ftEnergy = 0.0;
  double h1,h2,h_2,h;

  for (int i=0; i<config.numObjs; ++i){
    pairIndex = config.pairind2[ionIndex][i];
    q1=config.charge[config.pairind[pairIndex][0]]; //0 index is smaller index than 1 index
    q2=config.charge[config.pairind[pairIndex][1]]; //charge of the ions
    if(i < ionIndex){
      q=q1;
    }else{
      q=q2;
    }
    if(q){
      x1 = config.drpair[pairIndex][0]; //non-periodic distance
      y1 = config.drpair[pairIndex][1];
      z1 = config.drpair[pairIndex][2];
      for (int img = 0; img < imageItem.rCount; img++){ //loop over all periodic images
        dx = x1 + imageItem.rspace[img][0];
        dy = y1 + imageItem.rspace[img][1];
        r_ij = sqrt(dx*dx+dy*dy+z1*z1);
        if(r_ij > 0 && r_ij < config.rsCuttoff){
          ftemp=q*erfc(r_ij/alpha)/r_ij;
          rEnergy += ftemp;
        }
      }
      for(int img = 0; img<imageItem.ftCount; ++img){ //FOURIER SPACE Part
        h1=imageItem.fspace[img][0];
        h2=imageItem.fspace[img][1];
        h_2 = h1*h1+h2*h2;
        h=sqrt(h_2);
        if (h_2 > .00001){
          cosTerm = cos((2*M_PI*h1*x1)+(2*M_PI*h2*y1))/(2*h);
          functionF = exp(2*h*z1)*erfc(alpha*h+(z1/alpha)) + exp(-2*h*z1)*erfc(alpha*h-(z1/alpha));
          ftemp = q*cosTerm*functionF/(config.sideLength*config.sideLength);
          ftEnergy += ftemp;
        }
      }
      functiong = z1*erf(z1/alpha) + alpha*exp(-1*(z1*z1)/(alpha*alpha))/sqrt(M_PI);
      ftEnergy += -2*M_PI*q*functiong/(config.sideLength*config.sideLength);
    }
  }

  q1=config.charge[ionIndex];
  siEnergy = -2*q1/(alpha*sqrt(M_PI));

  pot = (ftEnergy+siEnergy+rEnergy)/config.gam;
  // cout << endl << "rEnergy: " << rEnergy << " siEnergy: " << siEnergy << " ftEnergy: " << ftEnergy << endl;
  // cout <<  ionIndex << " " << pot << endl;
  return pot; //*q1*0.5;
}

PeriodicImages get_images(double sideLength, double SLZ){
  int dimen = 2;  // this is a 2D slab code so 2 dimensions
  int n_max = 1;  // only need one real image
  int m_max = 20; // seems like 2D slab converges around 20 k-space images
  PeriodicImages imageItem(dimen,n_max,m_max,sideLength,SLZ);
  return imageItem;
}
