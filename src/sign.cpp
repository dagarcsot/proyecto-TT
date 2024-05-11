
#include "../include/sign.h"

//sign: returns absolute value of a with sign of b

double sign(double a, double b){
    return (b >= 0) ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
