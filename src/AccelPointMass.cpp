//
// Created by dagarcsot on 19/04/2024.
//
#include "../include/AccelPointMass.h"
#include <vector>
/*
AccelPointMass: Computes the perturbational acceleration due to a point mass

 Inputs:
   r           Satellite position vector
   s           Point mass position vector
   GM          Gravitational coefficient of point mass

 Output:
   a    		Acceleration (a=d^2r/dt^2)

 */

double AccelPointMass(const std::vector<double>&  r, const std::vector<double>&  s, double GM)){
    std::vector<double> d = r - s;
    return -GM * (d/(norm(d)^3) + s/(norm(s)^3) );
}