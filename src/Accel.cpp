
#include "../include/Accel.h"
#include "../include/globales.h"
#include "../include/IERS.h"
#include "../include/timediff.h"
#include "../include/SAT_Const.h"
#include "../include/NutMatrix.h"
#include "../include/Mjday_TDB.h"
#include "../include/AccelHarmonic.h"
#include "../include/AccelPointMass.h"
#include "../include/PoleMatrix.h"
#include "../include/PrecMatrix.h"
#include "../include/GHAMatrix.h"

/*
 * Purpose:
 *     Computes the acceleration of an Earth orbiting satellite due to
 *     - the Earth's harmonic gravity field,
 *     - the gravitational perturbations of the Sun and Moon
 *     - the solar radiation pressure and
 *     - the atmospheric drag
 *
 * Inputs:
 *     Mjd_TT      Terrestrial Time (Modified Julian Date)
 *     Y           Satellite state vector in the ICRF/EME2000 system
 *
 * Output:
 *     dY		    Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
 */

Matrix Accel(double x, Matrix &Y){


    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;

    IERS(eopdata,auxParam.Mjd_TT +x/86400,'l',x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);


    double UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
    timediff(UT1_UTC,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);

    double Mjd_UT1 = auxParam.Mjd_UTC + x/86400 + UT1_UTC/86400;
    double Mjd_TT = auxParam.Mjd_UTC + x/86400 + TT_UTC/86400;

    Matrix P = PrecMatrix(MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    double MJD_TDB = Mjday_TDB(Mjd_TT);
    Matrix r_Mercury(1, 3), r_Venus(1, 3), r_Earth(1, 3), r_Mars(1, 3), r_Jupiter(1, 3),
            r_Saturn(1, 3), r_Uranus(1, 3), r_Neptune(1, 3), r_Pluto(1, 3), r_Moon(1, 3),
            r_Sun(1, 3);
    JPL_Eph_DE430(MJD_TDB,r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun);

    //Acceleration due to harmonic gravity field
    Matrix Yaux(1,3);
    Yaux(1,1) = Y(1,1);
    Yaux(1,2) = Y(2,1);
    Yaux(1,3) = Y(3,1);

    Matrix a = AccelHarmonic(Yaux,E,auxParam.n,auxParam.m);

    // Luni-solar perturbartions

    if(auxParam.sun){
        a = a + AccelPointMass(Yaux,r_Sun,GM_Sun);
    }
    if(auxParam.moon){
        a = a + AccelPointMass(Yaux,r_Moon,GM_Moon);
    }
    if(auxParam.planets){
        a = a + AccelPointMass(Yaux,r_Mercury,GM_Mercury);
        a = a + AccelPointMass(Yaux,r_Venus,GM_Venus);
        a = a + AccelPointMass(Yaux,r_Mars,GM_Mars);
        a = a + AccelPointMass(Yaux,r_Jupiter,GM_Jupiter);
        a = a + AccelPointMass(Yaux,r_Saturn,GM_Saturn);
        a = a + AccelPointMass(Yaux,r_Uranus,GM_Uranus);
        a = a + AccelPointMass(Yaux,r_Neptune,GM_Neptune);
        a = a + AccelPointMass(Yaux,r_Pluto,GM_Pluto);
    }

    //Planetary perturbations



    Matrix dY(2,3);
    dY(1,1) = Y(4,1);
    dY(1,2) = Y(5,1);
    dY(1,3) = Y(6,1);

    dY(2,1) = a(1,1);
    dY(2,2) = a(1,2);
    dY(2,3) = a(1,3);

    return dY;










}