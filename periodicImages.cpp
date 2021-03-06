#include <vector>
#include <math.h>

#include "periodicImages.hpp"


using namespace std;

PeriodicImages::PeriodicImages(int dimen, int n_max=0, int m_max=0, double SLX=0, double SLY=0, double SLZ=0): n_max(n_max), m_max(m_max), SLX(SLX), SLY(SLY), SLZ(SLZ) {
  if (dimen == 2){
    initialize2D();
  }else if (dimen == 3){
    initialize3D();
  }else{ //placeholder
    rCount=0;
    ftCount=0;
  }
}

void PeriodicImages::initialize2D(){
  //initialize the real and fourier space periodic image vectors
  ftCount = (m_max*2+1)*(m_max*2+1); //total number of fourier space periodic images
  fspace.resize(ftCount);  //fourier space
  int count = 0;
  double kx = 2*M_PI/SLY; //lattice vector dimensions
  for (int i = -m_max; i <= m_max; ++i){
    for (int j = -m_max; j <= m_max; ++j){
      fspace[count].resize(2); // {0,0}
      fspace[count][0]=i*kx;
      fspace[count][1]=j*kx;
      count++;
    }
  }
  rCount = (n_max*2+1)*(n_max*2+1); // total number of real-space periodic images
  rspace.resize(rCount); //real space
  count = 0;
  for (int i = -n_max; i <= n_max; ++i){
    for (int j = -n_max; j <= n_max; ++j){
      rspace[count].resize(2); // {0,0}
      rspace[count][0]=i*SLY;
      rspace[count][1]=j*SLY;
      count++;
    }
  }
}

void PeriodicImages::initialize3D(){
  //initialize the real and fourier space periodic image vectors
  ftCount = (m_max*2+1)*(m_max*2+1)*(m_max*2+1); //total number of fourier space periodic images
  fspace.resize(ftCount);  //fourier space
  int count = 0;
  double kx = 2*M_PI/SLX; //lattice vector dimensions
  double ky = 2*M_PI/SLY;
  double kz = 2*M_PI/SLZ;
  for (int i = -m_max; i <= m_max; ++i){
    for (int j = -m_max; j <= m_max; ++j){
      for (int k = -m_max; k <= m_max; ++k){
        fspace[count].resize(3); // {0,0,0}
        fspace[count][0]=i*kx;
        fspace[count][1]=j*ky;
        fspace[count][2]=k*kz;
        count++;
      }
    }
  }
  rCount = (n_max*2+1)*(n_max*2+1)*(n_max*2+1); // total number of real-space periodic images
  rspace.resize(rCount); //real space
  count = 0;
  for (int i = -n_max; i <= n_max; ++i){
    for (int j = -n_max; j <= n_max; ++j){
      for (int k = -n_max; k <= n_max; ++k){
        rspace[count].resize(3); // {0,0,0}
        rspace[count][0]=i*SLX;
        rspace[count][1]=j*SLY;
        rspace[count][2]=k*SLZ;
        count++;
      }
    }
  }
}
