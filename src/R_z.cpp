//
// Created by dagarcsot on 12/04/2024.
//
#include <math.h>
#include "../include/Matrix.h"
#include "../include/R_z.h"

/*
      input:
        angle       - angle of rotation [rad]

      output:
        rotmat      - vector result
*/

Matrix R_z(double angle)
{
    double C, S;
    C = cos(angle);
    S = sin(angle);

    Matrix rotmat(3,3);

    rotmat(1,1) =      C;  rotmat(1,2) =   S;  rotmat(1,3) = 0.0;
    rotmat(2,1) = -1.0*S;  rotmat(2,2) =   C;  rotmat(2,3) = 0.0;
    rotmat(3,1) =    0.0;  rotmat(3,2) = 0.0;  rotmat(3,3) = 1.0;

    return rotmat;
}

