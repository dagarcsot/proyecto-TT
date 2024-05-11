

#include "../include/AccelPointMass.h"
#include <cmath>


/*
 * AccelPointMass: Computes the perturbational acceleration due to a point mass
 * Inputs:
 *     r           Satellite position vector
 *     s           Point mass position vector
 *     GM          Gravitational coefficient of point mass
 * Output:
 *     a    		Acceleration (a=d^2r/dt^2)
 */

Matrix AccelPointMass(Matrix  &r, Matrix  &s, double GM){

    //Relative position vector of satellite w.r.t. point mass
    Matrix d = r - s;

    double norm3 = pow(d.norma(),3);

    //Acceleration
    return (d * (1/norm3) + s * (1/norm3)) * -GM;

}