
#include "../include/gast.h"
#include "../include/SAT_Const.h"
#include "../include/gmst.h"
#include "../include/EqnEquinox.h"
#include <cmath>

/*
 * Purpose:
 *      Greenwich Apparent Sidereal Time
 *
 * Input:
 *      Mjd_UT1   Modified Julian Date UT1
 *
 * Output:
 *      gstime    GAST in [rad]
 *
 */

double gast(double Mjd_UT1){
    return fmod((gmst(Mjd_UT1)) +  EqnEquinox(Mjd_UT1) , pi2);
}