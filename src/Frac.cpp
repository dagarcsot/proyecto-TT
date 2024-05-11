

/*
 * Fractional part of a number (y=x-[x])
 */

double Frac(double x){
    int e = (int)x;
    return x - e;
}