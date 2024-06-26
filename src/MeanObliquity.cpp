
#include "../include/MeanObliquity.h"
#include "../include/SAT_Const.h"

/*
 * Purpose:
 *      Computes the mean obliquity of the ecliptic
 * Input:
 *      Mjd_TT    Modified Julian Date (Terrestrial Time)
 * Output:
 *      MOblq     Mean obliquity of the ecliptic [rad]
 *
 */

double MeanObliquity(double Mjd_TT) {
    double T = (Mjd_TT - MJD_J2000) / 36525.0;
    double MOblq = Rad * (84381.448 / 3600.0 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0);
    return MOblq;
}