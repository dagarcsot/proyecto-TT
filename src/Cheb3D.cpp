
#include <cstdlib>
#include <iostream>
#include "../include/Cheb3D.h"

Matrix Cheb3D(int t, int N, double Ta, double Tb, Matrix &Cx, Matrix &Cy, Matrix &Cz) {

    //Check validity
    if (t < Ta || Tb < t) {
        std::cout << "'ERROR: Time out of range in Cheb3D. \n";
        exit(-1);
    }

    //Clenshaw algorithm
    double tau = (2 * t - Ta - Tb) / (Tb - Ta);

    Matrix f1(1,3);
    Matrix f2(1,3);
    Matrix old_f1(1,3);
    Matrix aux(1,3);


    for(int i = N; i >= 2; i--) {

        aux(1,1) = Cx(1,i);
        aux(1,2) = Cy(1,i);
        aux(1,3) = Cz(1,i);
        old_f1 = f1;
        f1 = f1-f2+aux*2*tau;
        f2 = old_f1;

    }


    aux(1,1) = Cx(1,1);
    aux(1,2) = Cy(1,1);
    aux(1,3) = Cz(1,1);



    Matrix m = f1 * tau - f2 +  aux;
    return m;

}