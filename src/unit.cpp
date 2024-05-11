
#include "../include/unit.h"

/*
 *  this function calculates a unit vector given the original vector. if a
 *  zero vector is input, the vector is set to zero.
 *
 *   input:
 *          m1         - vector (use of matrix instead)
 *
 *   output:
 *          m2         - unit vector
 *
 */

Matrix unit(Matrix &vec) {
    double small = 0.000001;
    double magv = vec.norma();
    Matrix outvec(1, 3);

    if (magv > small) {
        for (int i = 1; i <= 3; i++) {
            outvec(1, i) = vec(1, i) / magv;
        }
    } else {
        for (int i = 1; i <= 3; i++) {
            outvec(1, i) = 0.0;
        }
    }
    return outvec;
}

