
#include <cmath>
#include "../include/sign_.h"

//sign: returns absolute value of a with sign of b

double sign_(double a, double b){
    double result;
    if (b>=0.0){
        result = abs(a);
    }

    else{
        result = - abs(a);
    }
    return result;
}
