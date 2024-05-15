//
// Created by dagarcsot on 27/03/2024.
//
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





int all_tests()
{
    _verify(proMat_01);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(Position_01);

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
