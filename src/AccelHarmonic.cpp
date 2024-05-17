

#include "../include/AccelHarmonic.h"
#include "../include/Legendre.h"
#include "../include/globales.h"
#include <cmath>

/*
 * Purpose:
 *     Computes the acceleration due to the harmonic gravity field of the
 *      central body
 *
 *  Inputs:
 *      r       Satellite position vector in the inertial system
 *      E       Transformation matrix to body-fixed system
 *      n_max   Maximum degree
 *      m_max   Maximum order (m_max<=n_max; m_max=0 for zonals, only)
 *
 *  Output:
 *      a       Acceleration (a=d^2r/dt^2)
 */
Matrix AccelHarmonic(Matrix &r, Matrix &E, double n_max, double m_max) {

    double r_ref = 6378.1363e3;   // Earth's radius [m]; GGM03S
    double gm = 398600.4415e9; // [m^3/s^2]; GGM03S

    //Body-fixed-position
    Matrix r_bf = E * r;

    //Auxiliary quantities
    double d = r_bf.norma();
    double latgc = asin(r_bf(1, 3));
    double lon = atan2(r_bf(1, 2), r_bf(1, 1));

    Matrix pnm(n_max + 1, m_max + 1);
    Matrix dpnm(n_max + 1, m_max + 1);
    Legendre(n_max, m_max, latgc, pnm, dpnm);

    double dUdr = 0.0;
    double dUdlatgc = 0.0;
    double dUdlon = 0.0;
    double q3 = 0.0;
    double q2 = q3;
    double q1 = q2;

    double b1, b2, b3;

    for (int n = 0; n <= n_max; n++) {
        b1 = (-gm / (d * d)) * pow(r_ref / d, n) * (n + 1);
        b2 = (gm / d) * pow(r_ref / d, n);
        b3 = (gm / d) * pow(r_ref / d, n);
/*
        for (int m = 0; m <= m_max; m++) {
            q1 = q1 + pnm(n + 1, m + 1) * (Cnm[n][m] * cos(m * lon) + Snm[n][m]* sin(m * lon));
            q2 = q2 + dpnm(n + 1, m + 1) * (Cnm[n][m] * cos(m * lon) +Snm[n][m] * sin(m * lon));
            q3 = q3 + m * pnm(n + 1, m + 1) * (Snm[n][m] * cos(m * lon) - Cnm[n][m] * sin(m * lon));
        }
        */
        dUdr = dUdr + q1 * b1;
        dUdlatgc = dUdlatgc + q2 * b2;
        dUdlon = dUdlon + q3 * b3;
        q3 = 0.0;
        q2 = q3;
        q1 = q2;
    }
    // Body-fixed acceleration
    double r2xy = pow(r_bf(1, 1), 2) * pow(r_bf(1, 2), 2);

    double ax = (1 / d * dUdr - r_bf(1, 3) / (d * d * sqrt(r2xy)) * dUdlatgc) * r_bf(1, 1) -
                (1 / r2xy * dUdlon) * r_bf(1, 2);
    double ay = (1 / d * dUdr - r_bf(1, 3) / (d * d * sqrt(r2xy)) * dUdlatgc) * r_bf(1, 2) +
                (1 / r2xy * dUdlon) * r_bf(1, 1);
    double az = 1 / d * dUdr * r_bf(1, 3) + sqrt(r2xy) / (d * d) * dUdlatgc;

    double v[] = {ax, ay, az};
    Matrix a_bf(3, 1, v, 3);

    // Inertial acceleration
    Matrix a =  E.transpuesta() * a_bf;
    return a;

}