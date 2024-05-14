
#include <iostream>
#include <cmath>
#include "../include/Geodetic.h"
#include "../include/SAT_Const.h"

/*
 * Purpose:
 *       geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
 *       from given position vector (r [m])
 */


void Geodetic(Matrix &r, double &lon, double &lat, double &h) {

    double R_equ = R_Earth;
    double f = f_Earth;

    double epsRequ = eps * R_equ;       // Convergence criterion
    double e2 = f * (2.0 - f);          // Square of eccentricity

    double X = r(1, 1);             // Cartesian cordenates
    double Y = r(1, 2);
    double Z = r(1, 3);
    double rho2 = X * X + Y * Y;              // Square of distance from z-axis

    // Check validity of input data
    if (r.norma() == 0.0) {
        std::cout << "invalid input in Geodetic constructor\n";
        lon = 0.0;
        lat = 0.0;
        h = -R_Earth;

    } else {
        //Iteration
        double dZ = e2 * Z;

        double ZdZ, Nh, SinPhi, N, dZ_new;

        while (true){
            ZdZ = Z + dZ;
            Nh = std::sqrt(rho2 + ZdZ * ZdZ);
            SinPhi = ZdZ / Nh;                          // Sine of geodetic latitude
            N = R_equ / std::sqrt(1.0 - e2 * SinPhi * SinPhi);
            dZ_new = N * e2 * SinPhi;
            if(abs(dZ - dZ_new) < epsRequ){
                break;
            }
            dZ = dZ_new;
        }

        // Longitude, latitude , altitude
        lon = atan2(Y, X);
        lat = atan2(ZdZ, sqrt(rho2));
        h = Nh - N;
    }

}