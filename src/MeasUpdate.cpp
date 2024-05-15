
#include "../include/MeasUpdate.h"

void MeasUpdate(Matrix &z,Matrix &s){
    int m = z.getNumCol();
    Matrix Inv_W(m,m);

    for(int i =1; i<=m;i++){
        Inv_W(i,i) = s(1,i)*s(1,i);    // Inverse weight (measurement covariance)
    }

    //Kalman gain
    //K = P*G'*inv(Inv_W+G*P*G');
    //
    // State update
    //x = x + K*(z-g);
    //
    // Covariance update
    //P = (eye(n)-K*G)*P;

}