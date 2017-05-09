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



  private:
    int   readInput(std::string);

    int sideLength;
    int SLZ;

    std::vector< std::vector<float> > position;
    std::vector<float>             charge;





};

#endif
