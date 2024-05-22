using namespace std;
#include <iostream>
#include "include/Matrix.h"
#include "include/AzElPa.h"
#include "include/AccelPointMass.h"
#include "include/Frac.h"
#include "include/MeanObliquity.h"
#include "include/Mjday.h"
#include "include/PoleMatrix.h"
#include "include/Position.h"
#include "include/sign_.h"
#include "include/unit.h"
#include "include/NutMatrix.h"
#include "include/Cheb3D.h"
#include "include/EccAnom.h"
#include <cmath>
#include <iomanip>
#include "include/Geodetic.h"
#include "include/Mjday_TDB.h"
#include "include/IERS.h"
#include "include/Legendre.h"
#include "include/R_x.h"
#include "include/TimeUpdate.h"
#include "include/PrecMatrix.h"


int main() {
    double v[] = {7.0, 2.0, 5.0,6.0};
    Matrix r = Matrix(2, 2, v, 4);
    r.print();
    Matrix s(1,1);
    s(1,1) = 4;
    s.inversa().print();





    return 0;

}