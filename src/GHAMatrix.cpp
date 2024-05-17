

#include "../include/GHAMatrix.h"
#include "../include/gast.h"
#include "../include/R_z.h"


/*
 * Purpose:
 *      Transformation from true equator and equinox to Earth equator and
 *      Greenwich meridian system
 *
 * Input:
 *      Mjd_UT1   Modified Julian Date UT1
 *
 * Output:
 *      GHAmat    Greenwich Hour Angle matrix
 *
 */
Matrix GHAMatrix(double Mjd_UT1){
    return R_z(gast(Mjd_UT1));
}