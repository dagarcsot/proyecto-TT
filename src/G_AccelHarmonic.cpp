

#include "../include/G_AccelHarmonic.h"
/*
 * Purpose:
 *      Computes the gradient of the Earth's harmonic gravity field
 * Inputs:
 *      r           Satellite position vector in the true-of-date system
 *      U           Transformation matrix to body-fixed system
 *      n           Gravity model degree
 *      m 			Gravity model order
 * Output:
 *      G    		Gradient (G=da/dr) in the true-of-date system
 */


Matrix G_AccelHarmonic(){
    double d = 1.0;
    Matrix G (3,3);
    Matrix dr(3,1);

    // Gradient
    for (int i=1;i<=3;i++){
        // Set offset in i-th component of the position vector
        dr(:) = 0.0;
        dr(i) = d;
        // Acceleration difference
        da = AccelHarmonic ( r+dr/2,U, n_max, m_max ) -AccelHarmonic ( r-dr/2,U, n_max, m_max );
        // Derivative with respect to i-th axis
        G(:,i) = da/d;
        end
    }

}