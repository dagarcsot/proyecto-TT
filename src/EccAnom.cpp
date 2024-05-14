
#include <cmath>
#include <iostream>
#include "../include/EccAnom.h"
#include "../include/SAT_Const.h"

/*
 * Purpose:
 *      Computes the eccentric anomaly for elliptic orbits
 *  Inputs:
 *      M         Mean anomaly in [rad]
 *      e         Eccentricity of the orbit [0,1]
 *  Output:
 *      Eccentric anomaly in [rad]
 *
 */

double EccAnom(double M, double e) {
    int maxit = 15;
    int i = 1;

    //Starting value
    M = fmod(M, pi2);

    double E;
    if (e < 0.8) {
        E = M;
    } else {
        E = pi;
    }

    double f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));

    //Iteration

    while (abs(f) > 100 * eps) {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i++;
        if (i == maxit) {
            std::cout << " convergence problems in EccAnom";
            exit(-1);
        }
    }
    return E;

}