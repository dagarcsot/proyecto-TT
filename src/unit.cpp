//
// Created by dagarcsot on 24/04/2024.
//

#include "../include/unit.h"

Matrix unit(Matrix m1) {
    double small = 0.000001;
    double magv = m1.norma();
    Matrix *m2 = new Matrix(1, 3);

    if (magv > small) {
        for (int i = 0; i < 3; i++) {
            (*m2)(0, i) = m1(0, i) / magv;
        }
    } else {
        for (int i = 0; i < 3; i++) {
            (*m2)(0, i) = 0.0;
        }
    }
    return *m2;

}
