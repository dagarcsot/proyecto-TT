
#ifndef PROYECTO_AZELPA_H
#define PROYECTO_AZELPA_H


#include "Matrix.h"

class AzElPa {

private:
    double azimuth;
    double El;
    Matrix *dAds;
    Matrix *dEds;
public:
    AzElPa(Matrix s);
    double getAzimuth();
    double getElevation();
    Matrix *getPartialsAzimuth();
    Matrix *getPartialsElevation();

};


#endif //PROYECTO_AZELPA_H
