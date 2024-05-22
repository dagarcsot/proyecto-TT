

#include <iostream>
#include <fstream>
#include <cstring>
#include "include/Matrix.h"
#include "include/IERS.h"
#include "include/globales.h"
#include "include/SAT_Const.h"
#include "include/Position.h"
#include "include/Mjday.h"

using namespace std;
/*
 * Initial Orbit Determination using Gauss and Extended Kalman Filter methods
 *
 * References:
 *      O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
 *      Applications", Springer Verlag, Heidelberg, 2000.
 *
 *      D. Vallado, "Fundamentals of Astrodynamics and Applications",
 *      4th Edition, 2013.
 *
 *
 *  Last modified:   2020/03/16   Meysam Mahooti
 *
 */

int main() {



    FILE *fid = fopen ("../data/GGM03S.txt","r");

    int n1, n2;
    double d1,d2,d3,d4;

    for (int n = 0; n <= 180; ++n) {
        for (int m = 0; m <= n; ++m) {

            if (fscanf(fid, "%d %d %lf %lf %lf %lf", &n1, &n2, &d1, &d2, &d3, &d4) != 6) {
                cerr << "Error reading data from GGM03S.txt" << endl;
                fclose(fid);
                return 1;
            }
            Cnm[n][m] = d1;
            Snm[n][m] = d2;
        }
    }
    fclose(fid);

    // Model parameters


    //read Earth orientation parameters
    fid = fopen("../data/eop19620101.txt","r");

    /*
     *------------------------------------------------------------------------------------------------------
     *|  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
     *|(0h UTC)           "         "          s          s          "        "          "         "     s
     *------------------------------------------------------------------------------------------------------
     */

    //eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
    fclose(fid);

    int nobs = 46;
    Matrix obs(nobs,4);

    // read observations
    fid = fopen("../data/GEOS3.txt","r");



    char tline[100];



    int Y, M, D, h, m;
    double s, az, el, Dist;

    for (int i = 1; i <= nobs; ++i) {


        if (fscanf(fid, "%d/%d/%d %d:%d:%lf %lf %lf %lf", &Y, &M, &D, &h, &m, &s, &az, &el, &Dist) != 9) {
            cerr << "Error reading data from GGM03S.txt" << endl;
            fclose(fid);
            return 1;
        }

        obs(i,1) = Mjday(Y,M,D,h,m,s);
        obs(i,2) = Rad*az;
        obs(i,3) = Rad*el;
        obs(i,4) = 1e3*Dist;

    }

    fclose(fid);



    double sigma_range = 92.5;          // [m]
    double sigma_az = 0.0224*Rad; // [rad]
    double sigma_el = 0.0139*Rad; // [rad]

    // Kaena Point station
    double lat = Rad*21.5748;     // [rad]
    double lon = Rad*(-158.2706); // [rad]
    double alt = 300.20;                // [m]

    Matrix Rs = Position(lon, lat, alt).transpuesta();


    double Mjd1 = obs(1,1);
    double Mjd2 = obs(9,1);
    double Mjd3 = obs(18,1);






    return 0;
}
