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
#include "include/sign.h"
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
   double P = 53005.0;
   double Phi = 53015.0;



    // Llamar a la función TimeUpdate
    Matrix updated_P = PrecMatrix(P, Phi);

    // Mostrar la matriz P actualizada
    updated_P.print();

    return 0;

}