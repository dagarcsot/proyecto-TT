

#include "../include/NutMatrix.h"
#include "../include/Matrix.h"
#include "../include/NutAngles.h"
#include "../include/R_x.h"
#include "../include/R_z.h"
#include "../include/MeanObliquity.h"
/*
 * Purpose:
 *      Transformation from mean to true equator and equinox
 * Input:
 *       Mjd_TT    Modified Julian Date (Terrestrial Time)
 * Output:
 *       NutMat    Nutation matrix
 *
 */

Matrix NutMatrix(double Mjd_TT) {

    double eps = MeanObliquity(Mjd_TT);


    double dpsi, deps;
    NutAngles(Mjd_TT, dpsi, deps);

    Matrix NutMat = R_x(-eps - deps) * R_z(-dpsi) * R_x(eps);

    return NutMat;
}
