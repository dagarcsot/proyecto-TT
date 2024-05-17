

#include "../include/EqnEquinox.h"
#include "../include/NutAngles.h"
#include "../include/MeanObliquity.h"
#include <cmath>
/*
 * Purpose:
 *      Computation of the equation of the equinoxes
 *
 * Input:
 *       Mjd_TT    Modified Julian Date (Terrestrial Time)
 *
 * Output:
 *      EqE      Equation of the equinoxes
 *
 * Notes:
 *      The equation of the equinoxes dpsi*cos(eps) is the right ascension of
 *      the mean equinox referred to the true equator and equinox and is equal
 *       to the difference between apparent and mean sidereal time.
 */

double EqnEquinox(double Mjd_TT){

    //Nutation in longitude and obliquity
    double dpsi,deps;
    NutAngles (Mjd_TT,dpsi,deps);

    //Equation of the equinoxes
    double EqE = dpsi * cos (MeanObliquity(Mjd_TT));
    return EqE;

}