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
    Simbox(int, int, std::string);
    friend double get_energy(const Simbox&);
    friend double get_mad_potential(const Simbox&, int);

  private:
    int  readInput(std::string);
    void initialize();
    void get_drpair0();
    void get_drpair1(int);

    int    sideLength;
    int    SLZ;
    int    numObjs;
    int    numPairs;
    float  gam;

    double inf;              // infinity

    std::vector< std::vector<float> >  position;
    std::vector<float>                 charge;
    std::vector< std::vector<double> > drpair;     // drpair[ionPairIndex][x/y/z] = 3-D distance between ions in ion pair
    std::vector< std::vector<int> >    pairind;    // pairind[ionPairIndex][0/1] = index of 1st/2nd ion in pair
    std::vector< std::vector<int> >    pairind2;   // pairind2[ionIndex 1][ionIndex 2] = ionPairIndex

};

#endif
