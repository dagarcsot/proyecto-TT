#include "../include/AzElPa.h"
#include "../include/SAT_Const.h"
#include <cmath>
#include <iostream>

/*
 * Purpose:
 *      Computes azimuth, elevation and partials from local tangent coordinates.
 * Input:
 *      s      Topocentric local tangent coordinates (East-North-Zenith frame)
 * Outputs:
 *      A      Azimuth [rad]
 *      E      Elevation [rad]
 *      dAds   Partials of azimuth w.r.t. s
 *      dEds   Partials of elevation w.r.t. s
 */
void AzElPa(Matrix &s, double &Az, double &El, Matrix &dAds, Matrix &dEds) {
    double rho = sqrt(s(1, 1) * s(1, 1) + s(1, 2) * s(1, 2));

    // Angles
    Az = atan2(s(1, 1), s(1, 2));
    if (Az < 0.0) {
        Az += pi2;
    }
    El = atan(s(1, 3) / rho);
    // Partials
    double v[] = {s(1, 2) / (rho * rho), -s(1, 1) / (rho * rho), 0.0};
    dAds = Matrix(1, 3, v, 3);

    double dotS = s.dot(s);
    double w[] = {(-s(1, 1) * s(1, 3) / rho) / dotS, (-s(1, 2) * s(1, 3) / rho) / dotS, rho / dotS};
    dEds = Matrix(1, 3, w, 3);
}
