
#ifndef PROYECTO_GLOBALES_H
#define PROYECTO_GLOBALES_H

#include "Matrix.h"

#include "Matrix.h"

struct AuxParam {
    double Mjd_TT;
    double Mjd_UTC;
    int n;
    int m;
    bool sun;
    bool moon;
    bool  planets;
};

extern Matrix Cnm;
extern Matrix Snm;
extern Matrix eopdata;
extern AuxParam auxParam;


#endif //PROYECTO_GLOBALES_H
