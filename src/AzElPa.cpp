//
// Created by dagarcsot on 03/05/2024.
//


/*
*  Purpose:
  Computes azimuth, elevation and partials from local tangent coordinates

 Input:
   s      Topocentric local tangent coordinates (East-North-Zenith frame)

 Outputs:
   A      Azimuth [rad]
   E      Elevation [rad]
   dAds   Partials of azimuth w.r.t. s
   dEds   Partials of elevation w.r.t. s



 */
#include "../include/AzElPa.h"
#include <math.h>


double AzElPa(Matrix s){


    double pi2 = 2.0 * M_PI;
    double rho = sqrt(s(0, 0) * s(0, 0) + s(0, 1) * s(0, 1));
    double Az = atan2(s(0, 0), s(0, 1));
    if(Az<0.0){
        Az += pi2;
    }
    double El = atan(s(0,3)/rho);
    Matrix *dAds = new Matrix(1, 3);
    (*dAds)(0, 0) = s(1, 0)/(rho*rho);
    (*dAds)(0, 1) = -s(0, 0)/(rho*rho);
    (*dAds)(0, 2) = 0.0;

    Matrix *dEds = new Matrix(1, 3);
    (*dEds)(0, 0) = s(0, 0)*s(2,0)/rho;
    (*dEds)(0, 1) = -s(0, 0)/(rho*rho);
    (*dEds)(0, 2) = 0.0;





}