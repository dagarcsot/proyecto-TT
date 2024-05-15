
#include "../include/TimeUpdate.h"

Matrix TimeUpdate(Matrix& P, Matrix& Phi) {
    Matrix aux =   Phi * P * Phi.transpuesta();

    return aux;
}
