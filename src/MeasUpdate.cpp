
#include <iostream>
#include "../include/MeasUpdate.h"

void MeasUpdate(Matrix &x,double z,double g, double s , Matrix &G, Matrix &P, int n, Matrix &K ) {

      double Inv_W = s * s;    // Inverse weight (measurement covariance)

    double aux = ( G * P * G.transpuesta())(1,1);
    double v[] ={aux +Inv_W};
    Matrix Aux(1,1,v,1);

    //Kalman gain
    Aux.print();
    Aux.inversa();
    //K = P * G.transpuesta() * (Aux).inversa();
/*
    // State update
    x = x + K*(z-g);
    //
    // Covariance update
    P = (Matrix::identidad(n)-K*G)*P;

*/
}