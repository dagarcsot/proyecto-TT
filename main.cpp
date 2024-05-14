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


int main() {
    // Definir los parámetros de la aproximación de Chebyshev
    double t = 0.5; // Valor de tiempo
    double N = 3; // Número de coeficientes
    double Ta = 0.0; // Inicio del intervalo
    double Tb = 1.0; // Fin del intervalo
    double Cxd[] = {1.0, 2.0, 3.0, 4.0}; // Coeficientes para la coordenada x
    double Cyd[] = {5.0, 6.0, 7.0, 8.0}; // Coeficientes para la coordenada y
    double Czd[] = {9.0, 10.0, 11.0, 12.0}; // Coeficientes para la coordenada z

    Matrix Cx(1,4,Cxd,4);
    Matrix Cy(1,4,Cyd,4);
    Matrix Cz(1,4,Czd,4);


    // Calcular la aproximación de Chebyshev
    Matrix result = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);

    // Imprimir el resultado
    std::cout << "Aproximación de Chebyshev en t = " << t << ": "<<endl;
    result.print();

    return 0;

}