//
// Created by dagarcsot on 24/04/2024.
//

/*
 * Fractional part of a number (y=x-[x])
 */

double Frac(double x){
    int e = x/1;
    return x - e;
}