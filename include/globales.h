
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

extern  double Cnm[181][181];
extern double Snm[181][181];
extern Matrix eopdata;
extern AuxParam auxParam;
extern Matrix PC;


#endif //PROYECTO_GLOBALES_H
