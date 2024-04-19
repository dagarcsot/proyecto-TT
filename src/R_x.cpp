//
// Created by dagarcsot on 19/04/2024.
//
#include <math.h>
#include "../include/Matrix.h"
#include "../include/R_x.h"

/*
      input:
        angle       - angle of rotation [rad]

      output:
        rotmat      - vector result
*/

Matrix R_x(double angle)
{
    double C, S;
    C = cos(angle);
    S = sin(angle);

    Matrix rotmat(3,3);

    rotmat(1,1) = 1.0;  rotmat(1,2) =      0.0;  rotmat(1,3) = 0.0;
    rotmat(2,1) = 0.0;  rotmat(2,2) =        C;  rotmat(2,3) = S;
    rotmat(3,1) = 0.0;  rotmat(3,2) = -1.0 * S;  rotmat(3,3) = C;

    return rotmat;
}

