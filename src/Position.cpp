#include "../include/Position.h"
#include "../include/SAT_Const.h"
#include <cmath>
/*
 *  Purpose:
 *      Position vector (r [m]) from geodetic coordinates (Longitude [rad],
 *      latitude [rad], altitude [m])
 *
 */


Matrix Position(double lon, double lat, double h) {
    double R_equ = R_Earth;
    double f = f_Earth;

    double e2 = f * (2.0 - f);         // Square of eccentricity3
    double CosLat = cos(lat);       // (Co)sine of geodetic latitude
    double SinLat = sin(lat);

    // Position vector
    double N = R_equ / sqrt(1.0 - e2 * SinLat * SinLat);

    double v[] = {(N + h) * CosLat * cos(lon), (N + h) * CosLat * sin(lon), ((1.0 - e2) * N + h) * SinLat};
    Matrix r = Matrix(1, 3, v, 3);

    return r;
}
