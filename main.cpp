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
#include <cmath>
#include <iomanip>

using namespace std;

int main() {

    // Ejemplo de uso
    int t = 8;
    int N = 3;
    double Ta = 7.0;
    double Tb = 11.0;

    double cx_data[] = {1.0, 2.0, 3.0};
    double cy_data[] = {4.0, 5.0, 6.0};
    double cz_data[] = {7.0, 8.0, 9.0};

    Matrix Cx(1, 3, cx_data, 3);
    Matrix Cy(1, 3, cy_data, 3);
    Matrix Cz(1, 3, cz_data, 3);

    Matrix m = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);
    m.print();

    return 0;
}
