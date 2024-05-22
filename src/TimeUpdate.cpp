
#include "../include/TimeUpdate.h"

Matrix TimeUpdate(Matrix& P, Matrix& Phi) {

    return Phi * P * Phi.transpuesta();
}
