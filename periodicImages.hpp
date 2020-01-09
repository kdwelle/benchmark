#ifndef IMAGES_HPP
#define IMAGES_HPP

class PeriodicImages{
/* Class for a container to store the vectors of real and reciprical space
images used to calculate various flavors of Ewald summations
*/

  public:
    PeriodicImages(int,int,int,double,double,double);
    void initialize2D();
    void initialize3D();

    int   rCount;           // total number of real periodic images
    int   ftCount;          // total number of fourier periodic images

    std::vector< std::vector<float> >  fspace;           //reciprical (fourier) space repeating images
    std::vector< std::vector<int> >    rspace;           //real space repeating images


  private:
    int n_max;           // cells in real space
    int m_max;           // cells in fourier space
    double SLX;
    double SLY;
    double SLZ;          // dimension normal to elecrodes; =0 for 2D-slab in which case leave default
};

#endif
