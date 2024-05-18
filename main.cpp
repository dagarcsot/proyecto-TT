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
    // Ejemplo de uso de la función TimeUpdate
    Matrix eop(13, 2);
    eop(1, 1) = 2022; eop(2, 1) = 5; eop(3, 1) = 1; eop(4, 1) = 59544; eop(5, 1) = 0.196;
    eop(6, 1) = 0.232; eop(7, 1) = -0.228; eop(8, 1) = 0.1; eop(9, 1) = 0.2; eop(10, 1) = 0.3;
    eop(11, 1) = 0.4; eop(12, 1) = 0.5; eop(13, 1) = 37;

    eop(1, 2) = 2022; eop(2, 2) = 1; eop(3, 2) = 2; eop(4, 2) = 59545; eop(5, 2) = 0.197;
    eop(6, 2) = 0.233; eop(7, 2) = -0.227; eop(8, 2) = 0.1; eop(9, 2) = 0.2; eop(10, 2) = 0.3;
    eop(11, 2) = 0.4; eop(12, 2) = 0.5; eop(13, 2) = 37;

    eop.print();

    Matrix a = eop.getPrimeraFil(1,3);
    a.print();
    return 0;

}