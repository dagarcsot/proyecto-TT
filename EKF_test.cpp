
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include "./include/Matrix.h"
#include "include/SAT_Const.h"
#include "include/R_x.h"
#include "include/R_z.h"
#include "include/R_y.h"
#include "include/NutMatrix.h"
#include "include/Position.h"
#include "include/IERS.h"
#include "include/TimeUpdate.h"
#include "include/unit.h"
#include "include/timediff.h"
#include "include/sign_.h"
#include "include/Frac.h"
#include "include/AccelPointMass.h"
#include "include/AccelHarmonic.h"
#include "include/Legendre.h"
#include "include/MeasUpdate.h"

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

using namespace std;


int proMat_01() {
    double v1[] = {1.0, 2.0, 3.0, 4.0};
    double v2[] = {1.0, 0.0, 0.0, 1.0};
    Matrix m1(2, 2, v1, 4);
    Matrix m2(2, 2, v2, 4);
    Matrix sol(2, 2);

    sol = m1 * m2;

    _assert(sol(1, 1) == m1(1, 1) && sol(1, 2) == m1(1, 2) && sol(2, 1) == m1(2, 1) && sol(2, 2) == m1(2, 2));

    return 0;
}

int R_x_01() {
    Matrix m1 = R_x(cos(pi / 3));

    double v[] = {1.0, 0, 0, 0, 0.877582561890373, 0.479425538604203, 0, -0.479425538604203, 0.877582561890373};
    Matrix m2 = Matrix(3, 3, v, 9);

    _assert(m1.equals(m2, 10));

    return 0;
}

int R_y_01() {
    Matrix m1 = R_y(cos(pi / 3));

    double v[] = {0.877582561890373, 0, -0.479425538604203, 0, 1.0, 0, 0.479425538604203, 0, 0.877582561890373};
    Matrix m2 = Matrix(3, 3, v, 9);

    _assert(m1.equals(m2, 10));

    return 0;
}

int R_z_01() {
    Matrix m1 = R_z(cos(pi / 3));

    double v[] = {0.877582561890373, 0.479425538604203, 0, -0.479425538604203, 0.877582561890373, 0, 0, 0, 1};
    Matrix m2 = Matrix(3, 3, v, 9);

    _assert(m1.equals(m2, 10));

    return 0;
}

int Position_01() {
    double lat = Rad * 21.5748;
    double lon = Rad * (-158.2706);
    double alt = 300.20;

    double v[] = {-5.512567840036068e+06, -2.196994446669333e+06, 2.330804966146887e+06};
    Matrix m1 = Matrix(1, 3, v, 3);
    Matrix m2 = Position(lon, lat, alt);


    _assert(m1.equals(m2, 9));

    return 0;
}

int timediff_01(){


    double UT1_UTC = 1.0;
    double TAI_UTC = 2.0;
    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;

    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

    _assert(UT1_TAI == -1.0);
    _assert(UTC_GPS == 17.0);
    _assert(UT1_GPS == 18.0);
    _assert(TT_UTC == 34.184);
    _assert(GPS_UTC == -17.0);

    return 0;


}


int TimeUpdate_01() {

    double v[] = {1.0, 2.0, 3.0,
                  4.0, 5.0, 6.0,
                  7.0, 8.0, 9.0};
    Matrix P = Matrix(3, 3, v, 9);

    double w[] = {0.5, 1.0, 1.5,
                  2.0, 2.5, 3.0,
                  3.5, 4.0, 4.5};
    Matrix Phi = Matrix(3, 3, w, 9);

    P = TimeUpdate(P, Phi);

    double u[] = {57.0, 138.0, 219.0,
                  129.0, 311.25, 493.5,
                  201.0, 484.5, 768.0};
    Matrix Res(3, 3, u, 9);


    _assert(P.equals(Res,10));

    return 0;
}


int Unit_01() {
    double v[] = {7.0, 2.0, 5.0};
    Matrix M = Matrix(1, 3, v, 3);

    M = unit(M);

    double w[] = { 0.792593923901217, 0.226455406828919 , 0.566138517072298};
    Matrix Res = Matrix(1, 3, w, 3);


    _assert(M.equals(Res, 10));

    return 0;
}

int sign_01() {



    _assert(sign_(-2,1) ==2);
    _assert(sign_(-2,-1) ==-2);
    _assert(sign_(2,-1) ==-2);
    _assert(sign_(2,1) ==2);

    return 0;
}


int Frac_01() {

    _assert(Frac(7.25) == 0.25);

    return 0;
}

int AccelPointMass_01() {


    double v[] = {7.0, 2.0, 5.0};
    Matrix r = Matrix(1, 3, v, 3);

    double w[] = {5.5, 13.0, 1.5};
    Matrix s = Matrix(1, 3, w, 3);

    Matrix apm = AccelPointMass(r,s,GM_Earth);


    double u[] = {-1145527739233.2365, 968224951909.6899, -1093531800186.5676};
    Matrix res(1,3,u,3);


    _assert(res.equals(apm,8));

    return 0;
}

int Legendre_01(){
   int n = 2;
   int m = 2;
   double phi = 0.15;
   Matrix pnm(n+1,n+1);
    Matrix dpnm(n+1,n+1);

    Legendre(n,m,phi,pnm,dpnm);



    return 0;
}


int JPL_Eph_DE430_01(){
    double Mjd_TDB;
    /*JPL_Eph_DE430(Mjd_TDB, Matrix &r_Mercury, Matrix &r_Venus, Matrix &r_Earth, Matrix &r_Mars, Matrix &r_Jupiter,
            Matrix &r_Saturn, Matrix &r_Uranus, Matrix &r_Neptune, Matrix &r_Pluto, Matrix &r_Moon, Matrix &r_Sun)
    */
    return 0;
}


int MeasUpdate_01(){

    int n = 6;
    double v[] = {5738566.5776918,3123975.34092958, 3727114.48156063,
                  5199.63329181126,-2474.43881044665,-7195.16752553894};
    Matrix x(n, 1,v,6);

    double w[] ={ 9.59123748602943e-08,2.16050345227544e-07,-3.27382770920699e-07};

    double z = 1.0559084894933;
    double g = 1.05892995381513;
    double s = 0.00039095375244673;
    Matrix G(1, n,w,6);

    double u[] = {
            101453348.207917, 120429.109355752, 148186.145010685, 39372.9209771494, 3284.21674871861, 4014.15727499737,
            120429.109355752, 101309543.076737, 84141.6477319108, 3284.34773346912, 35369.9224485894, 2255.66799443683,
            148186.145010685, 84141.6477319108, 101344434.103716, 4014.41933457261, 2255.72532464628, 36274.7873567659,
            39372.9209771494, 3284.34773346912, 4014.41933457261, 1001.21615369033, 1.32096249010467, 1.60455480925268,
            3284.21674871861, 35369.9224485894, 2255.72532464628, 1.32096249010467, 999.576829597177, 0.892927374761907,
            4014.15727499737, 2255.66799443683, 36274.7873567659, 1.60455480925268, 0.892927374761907, 999.924178045209};

    Matrix P(n, n,u, 36);
    for (int i = 1; i <= n; ++i) {
        P(i, i) = 1.0;
    }
    Matrix K(n, 1);

    MeasUpdate(x, z, g, s, G, P, n, K);
/*
    K.print();
    x.print();
    P.print();

 */
    return 0;
}


int IERS_01() {
    Matrix eop(13, 2);
    eop(1, 1) = 2022;
    eop(2, 1) = 1;
    eop(3, 1) = 1;
    eop(4, 1) = 59544;
    eop(5, 1) = 0.196;
    eop(6, 1) = 0.232;
    eop(7, 1) = -0.228;
    eop(8, 1) = 0.1;
    eop(9, 1) = 0.2;
    eop(10, 1) = 0.3;
    eop(11, 1) = 0.4;
    eop(12, 1) = 0.5;
    eop(13, 1) = 37;

    eop(1, 2) = 2022;
    eop(2, 2) = 1;
    eop(3, 2) = 2;
    eop(4, 2) = 59545;
    eop(5, 2) = 0.197;
    eop(6, 2) = 0.233;
    eop(7, 2) = -0.227;
    eop(8, 2) = 0.1;
    eop(9, 2) = 0.2;
    eop(10, 2) = 0.3;
    eop(11, 2) = 0.4;
    eop(12, 2) = 0.5;
    eop(13, 2) = 37;


    double Mjd_UTC = 59544.5;
    char interp = 'l';

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;

    IERS(eop, Mjd_UTC, interp, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);

    double tolerancia = 1e-10;

    _assert(fabs(x_pole - 9.526589e-07) < tolerancia);
    _assert(fabs(y_pole - 1.127192e-06) < tolerancia);
    _assert(fabs(UT1_UTC - (-0.227500)) < tolerancia);
    _assert(fabs(LOD - 0.100000) < tolerancia);
    _assert(fabs(dpsi - 9.696274e-07) < tolerancia);
    _assert(fabs(deps - 1.454441e-06) < tolerancia);
    _assert(fabs(dx_pole - 1.939255e-06) < tolerancia);
    _assert(fabs(dy_pole - 2.424068e-06) < tolerancia);
    _assert(fabs(TAI_UTC - 37.000000) < tolerancia);

    interp = 'n';

    IERS(eop, Mjd_UTC, interp, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);

    _assert(fabs(x_pole - 9.502348e-07) < tolerancia);
    _assert(fabs(y_pole - 1.124768e-06) < tolerancia);
    _assert(fabs(UT1_UTC - (-0.22800)) < tolerancia);
    _assert(fabs(LOD - 0.100000) < tolerancia);
    _assert(fabs(dpsi - 9.696274e-07) < tolerancia);
    _assert(fabs(deps - 1.454441e-06) < tolerancia);
    _assert(fabs(dx_pole - 1.939255e-06) < tolerancia);
    _assert(fabs(dy_pole - 2.424068e-06) < tolerancia);
    _assert(fabs(TAI_UTC - 37.000000) < tolerancia);


    return 0;
}


int all_tests() {
    _verify(proMat_01);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(Position_01);
    _verify(TimeUpdate_01);
    _verify(timediff_01);
    _verify(Unit_01);
    _verify(sign_01);
    _verify(Frac_01);
    _verify(AccelPointMass_01);
    _verify(Legendre_01);
    _verify(MeasUpdate_01);
    _verify(IERS_01);


    return 0;
}


int main() {

    int result = all_tests();
    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
