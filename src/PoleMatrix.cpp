



#include "../include/PoleMatrix.h"
#include "../include/Matrix.h"
#include "../include/R_y.h"
#include "../include/R_x.h"



/*
 * Purpose:
 * Transformation from pseudo Earth-fixed to Earth-fixed coordinates
 * for a given date
 *
 * Input:
 *      Pole coordinate(xp,yp)
 *
 * Output:
 *      PoleMate Pole matrix
 *
 */

Matrix PoleMatrix(double xp, double yp) {
    Matrix PoleMat = R_y(-xp) * R_x(-yp);
    return PoleMat;
}

