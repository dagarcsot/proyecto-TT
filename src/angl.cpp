
#include "../include/angl.h"
#include <cmath>


/*
 * inputs:
 *      vec1       -vector 1
 *      vec2       -vector 2
 *
 * output:
 *      theta      -angle between the two vectors  -pi to pi
 */
double angl(Matrix &vec1, Matrix &vec2) {

    double small = 0.00000001;
    double undefined = 999999.1;

    double magv1 = vec1.norma();
    double magv2 = vec2.norma();

    double theta, temp;
    if (magv1 * magv2 > pow(small,2)) {
        temp = vec1.dot(vec2) / (magv1 * magv2);
        if (std::abs(temp) > 1.0) {
            temp = temp > 0 ? 1.0 : -1.0;
        }
        theta = acos(temp);
    } else {
        theta = undefined;
    }

    return theta;

}

