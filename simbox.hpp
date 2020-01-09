#ifndef SIMBOX_HPP
#define SIMBOX_HPP
#include <string>

class Simbox{
/* Class for the overhead for any generic simulation. Stores all the
  distances, ions, etc.

   Naming conventions:
   functions_have_underscores()
   variablesUseCamelCase
*/
  public:
    Simbox(double, double, double, int, std::string);
    void translate_particle(int, float, float, float); //move particle by x,y,z
    void set_position(int, float, float, float);
    void change_charge(int,float);
    void set_charge(int,float);
    void output_analysis();

    friend double get_energy(const Simbox&, const PeriodicImages&);
    friend double get_mad_potential(const Simbox&, int, const PeriodicImages&);

    int   readInput(std::string);
    void  initialize();
    void  get_drpair0();
    void  get_drpair1(int);
    float get_image(float,int);
    float get_image_charge(float,int);

    double SLX;
    double SLY;
    double SLZ;
    int    numImageReflections;  // how many image charges to include (0 is none, 1 is as many image charges as real charges)
    float  zbegin;
    float  zend;
    int    numObjs;              // total number of objects in the simulation (including image charges)
    int    numReal;              // total number of "real" (i.e. non-image) charge in the simulation
    int    numPairs;
    float  gam;

    double inf;                          // infinity
    int                rsCuttoff;        // real space cutoff radius


    std::vector<int>                   realCharges; // a list of all the non-image charges in the syste, by index
    std::vector< std::vector<float> >  position;    // positoin[ionIndex][x/y/z] = position in x/y/z
    std::vector<float>                 charge;
    std::vector< std::vector<double> > drpair;     // drpair[ionPairIndex][x/y/z] = 3-D distance between ions in ion pair
    std::vector< std::vector<int> >    pairind;    // pairind[ionPairIndex][0/1] = index of 1st/2nd ion in pair
    std::vector< std::vector<int> >    pairind2;   // pairind2[ionIndex 1][ionIndex 2] = ionPairIndex
};

#endif
