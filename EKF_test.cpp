
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

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

using namespace std;


int proMat_01()
{
    double v1[] = {1.0, 2.0, 3.0, 4.0};
    double v2[] = {1.0, 0.0, 0.0, 1.0};
    Matrix m1(2, 2, v1, 4);
    Matrix m2(2, 2, v2, 4);
    Matrix sol(2, 2);

    sol = m1 * m2;

    _assert(sol(1,1) == m1(1,1) && sol(1,2) == m1(1,2) && sol(2,1) == m1(2,1) && sol(2,2) == m1(2,2));

    return 0;
}

int R_x_01() {
    Matrix m1 = R_x(cos(pi/3));

    double v[] = {1.0, 0, 0,0, 0.877582561890373,  0.479425538604203,0, -0.479425538604203 ,  0.877582561890373};
    Matrix m2 = Matrix(3, 3, v, 9);

    _assert(m1.equals(m2, 10));

    return 0;
}

int R_y_01() {
    Matrix m1 = R_y(cos(pi/3));

    double v[] = { 0.877582561890373, 0, -0.479425538604203,0, 1.0,  0,0.479425538604203, 0 ,  0.877582561890373};
    Matrix m2 = Matrix(3, 3, v, 9);

    _assert(m1.equals(m2, 10));

    return 0;
}

int R_z_01() {
    Matrix m1 = R_z(cos(pi/3));

    double v[] = {0.877582561890373, 0.479425538604203 , 0,-0.479425538604203, 0.877582561890373,  0,0, 0 ,  1};
    Matrix m2 = Matrix(3, 3, v, 9);

    _assert(m1.equals(m2, 10));

    return 0;
}

int Position_01() {
    double lat = Rad*21.5748;
    double lon = Rad*(-158.2706);
    double alt = 300.20;

    double v[] = {-5.512567840036068e+06 , -2.196994446669333e+06,2.330804966146887e+06};
    Matrix m1= Matrix(1, 3, v, 3);
    Matrix m2 = Position(lon,lat,alt);


    _assert(m1.equals(m2, 9));

    return 0;
}


int IERS_01() {
    // Datos de entrada
    Matrix eop(13, 2);
    eop(1, 1) = 2022; eop(2, 1) = 1; eop(3, 1) = 1; eop(4, 1) = 59544; eop(5, 1) = 0.196;
    eop(6, 1) = 0.232; eop(7, 1) = -0.228; eop(8, 1) = 0.1; eop(9, 1) = 0.2; eop(10, 1) = 0.3;
    eop(11, 1) = 0.4; eop(12, 1) = 0.5; eop(13, 1) = 37;

    eop(1, 2) = 2022; eop(2, 2) = 1; eop(3, 2) = 2; eop(4, 2) = 59545; eop(5, 2) = 0.197;
    eop(6, 2) = 0.233; eop(7, 2) = -0.227; eop(8, 2) = 0.1; eop(9, 2) = 0.2; eop(10, 2) = 0.3;
    eop(11, 2) = 0.4; eop(12, 2) = 0.5; eop(13, 2) = 37;


    double Mjd_UTC = 59544.5;
    char interp = 'l';

    // Variables de salida
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;

    // Llamada a la función
    IERS(eop, Mjd_UTC, interp, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);

    // Impresión de resultados
    std::cout << "x_pole: " << x_pole << std::endl;
    std::cout << "y_pole: " << y_pole << std::endl;
    std::cout << "UT1_UTC: " << UT1_UTC << std::endl;
    std::cout << "LOD: " << LOD << std::endl;
    std::cout << "dpsi: " << dpsi << std::endl;
    std::cout << "deps: " << deps << std::endl;
    std::cout << "dx_pole: " << dx_pole << std::endl;
    std::cout << "dy_pole: " << dy_pole << std::endl;
    std::cout << "TAI_UTC: " << TAI_UTC << std::endl<<endl;

    double tolerancia = 1e-10;

    _assert(fabs(x_pole - 9.526589e-07)<tolerancia);
    _assert(fabs(y_pole - 1.127192e-06)<tolerancia);
    _assert(fabs(UT1_UTC -(-0.227500))<tolerancia);
    _assert(fabs(LOD -0.100000)<tolerancia);
    _assert(fabs(dpsi - 9.696274e-07)<tolerancia);
    _assert(fabs(deps -1.454441e-06)<tolerancia);
    _assert(fabs(dx_pole- 1.939255e-06)<tolerancia);
    _assert(fabs(dy_pole - 2.424068e-06)<tolerancia);
    _assert(fabs(TAI_UTC - 37.000000)<tolerancia);

    interp = 'n';


    IERS(eop, Mjd_UTC, interp, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);


    std::cout << "x_pole: " << x_pole << std::endl;
    std::cout << "y_pole: " << y_pole << std::endl;
    std::cout << "UT1_UTC: " << UT1_UTC << std::endl;
    std::cout << "LOD: " << LOD << std::endl;
    std::cout << "dpsi: " << dpsi << std::endl;
    std::cout << "deps: " << deps << std::endl;
    std::cout << "dx_pole: " << dx_pole << std::endl;
    std::cout << "dy_pole: " << dy_pole << std::endl;
    std::cout << "TAI_UTC: " << TAI_UTC << std::endl;

    _assert(fabs(x_pole - 9.502348e-07)<tolerancia);
    _assert(fabs(y_pole -1.124768e-06)<tolerancia);
    _assert(fabs(UT1_UTC -(-0.22800))<tolerancia);
    _assert(fabs(LOD -0.100000)<tolerancia);
    _assert(fabs(dpsi - 9.696274e-07)<tolerancia);
    _assert(fabs(deps -1.454441e-06)<tolerancia);
    _assert(fabs(dx_pole- 1.939255e-06)<tolerancia);
    _assert(fabs(dy_pole - 2.424068e-06)<tolerancia);
    _assert(fabs(TAI_UTC - 37.000000)<tolerancia);





    return 0;


}






int all_tests()
{
    _verify(proMat_01);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(Position_01);
    _verify(IERS_01);

    return 0;
}


int main()
{

    int result = all_tests();
    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
