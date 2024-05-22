

#include "../include/JPL_Eph_DE430.h"
#include "../include/Cheb3D.h"

/*
 * JPL_Eph_DE430: Computes the sun, moon, and nine major planets' equatorial
 *                position using JPL Ephemerides
 * Inputs:
 *       Mjd_TDB         Modified julian date of TDB
 * Output:
 *       r_Earth(solar system barycenter (SSB)),r_Mars,r_Mercury,r_Venus,
 *       r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,
 *        r_Sun(geocentric equatorial position ([m]) referred to the
 *        International Celestial Reference Frame (ICRF))
 * Notes:
 *       Light-time is already taken into account
 */

void
JPL_Eph_DE430(double Mjd_TDB, Matrix &r_Mercury, Matrix &r_Venus, Matrix &r_Earth, Matrix &r_Mars, Matrix &r_Jupiter,
              Matrix &r_Saturn, Matrix &r_Uranus, Matrix &r_Neptune, Matrix &r_Pluto, Matrix &r_Moon, Matrix &r_Sun) {

//global
     extern Matrix PC;


    double JD = Mjd_TDB + 2400000.5;

    int i;
    for (i = 0; i < PC.getNumFil(); ++i) {
        if (PC(i + 1, 1) <= JD && JD <= PC(i + 1, 2)) {
            break;
        }
    }
    Matrix PCtemp = PC.getFil(i);

    double t1 = PCtemp(1, 1) - 2400000.5; // MJD at start of interval
    double dt = Mjd_TDB - t1;


    //------------------EARTH----------------------

    int temp[] = {231, 244, 257, 270};


    Matrix Cx_Earth = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
    Matrix Cy_Earth = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
    Matrix Cz_Earth = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);

    for (int i = 0; i < 4; i++) {
        temp[i] += 39;
    }

    Matrix Cx = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
    Matrix Cy = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
    Matrix Cz = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);

    Cx_Earth = Cx_Earth.concatenar(Cx);
    Cy_Earth = Cy_Earth.concatenar(Cy);
    Cz_Earth = Cz_Earth.concatenar(Cz);

    int j;
    double Mjd0;
    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    } else if (16 < dt && dt <= 32) {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }


    Matrix CxAux = Cx_Earth.getPrimeraFil(13 * j + 1, 13 * j + 13);
    Matrix CyAux = Cy_Earth.getPrimeraFil(13 * j + 1, 13 * j + 13);
    Matrix CzAux = Cz_Earth.getPrimeraFil(13 * j + 1, 13 * j + 13);
    r_Earth = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 16, CxAux, CyAux, CzAux).transpuesta() * 1e3;


    //------------------MOON----------------------


    temp[0] = 441;
    temp[1] = 454;
    temp[2] = 467;
    temp[3] = 480;

    Matrix Cx_Moon = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
    Matrix Cy_Moon = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
    Matrix Cz_Moon = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);

    for (int i = 1; i <= 7; i++) {
        for (int j = 0; j < 4; j++) {
            temp[j] += 39;
        }

        Cx = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
        Cy = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
        Cz = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);

        Cx_Moon = Cx_Moon.concatenar(Cx);
        Cy_Moon = Cy_Moon.concatenar(Cy);
        Cz_Moon = Cz_Moon.concatenar(Cz);
    }

    if (0 <= dt && dt <= 4) {
        j = 0;
        Mjd0 = t1;
    } else if (dt > 4 && dt <= 8) {
        j = 1;
        Mjd0 = t1 + 4 * j;
    } else if (dt > 8 && dt <= 12) {
        j = 2;
        Mjd0 = t1 + 4 * j;
    } else if (dt > 12 && dt <= 16) {
        j = 3;
        Mjd0 = t1 + 4 * j;
    } else if (dt > 16 && dt <= 20) {
        j = 4;
        Mjd0 = t1 + 4 * j;
    } else if (dt > 20 && dt <= 24) {
        j = 5;
        Mjd0 = t1 + 4 * j;
    } else if (dt > 24 && dt <= 28) {
        j = 6;
        Mjd0 = t1 + 4 * j;
    } else if (dt > 28 && dt <= 32) {
        j = 7;
        Mjd0 = t1 + 4 * j;
    }

    CxAux = Cx_Moon.getPrimeraFil(13 * j + 1, 13 * j + 13);
    CyAux = Cy_Moon.getPrimeraFil(13 * j + 1, 13 * j + 13);
    CzAux = Cz_Moon.getPrimeraFil(13 * j + 1, 13 * j + 13);
    r_Moon = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 4, CxAux, CyAux, CzAux).transpuesta() * 1e3;


    //------------------SUN----------------------

    temp[0] = 753;
    temp[1] = 764;
    temp[2] = 775;
    temp[3] = 786;


    Matrix Cx_Sun = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
    Matrix Cy_Sun = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
    Matrix Cz_Sun = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);
    for (int i = 0; i < 4; j++) {
        temp[i] += 33;
    }
    Cx = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
    Cy = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
    Cz = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);

    Cx_Sun = Cx_Sun.concatenar(Cx);
    Cy_Sun = Cy_Sun.concatenar(Cy);
    Cz_Sun = Cz_Sun.concatenar(Cz);
    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    } else if (16 < dt && dt <= 32) {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }


    CxAux = Cx_Moon.getPrimeraFil(11 * j + 1, 11 * j + 11);
    CyAux = Cy_Moon.getPrimeraFil(11 * j + 1, 11 * j + 11);
    CzAux = Cz_Moon.getPrimeraFil(11 * j + 1, 11 * j + 11);
    r_Sun = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 16, CxAux, CyAux, CzAux).transpuesta() * 1e3;



    //------------------MERCURY----------------------


    temp[0] = 3;
    temp[1] = 17;
    temp[2] = 31;
    temp[3] = 45;

    Matrix Cx_Mercury = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
    Matrix Cy_Mercury = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
    Matrix Cz_Mercury = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);
    for (int i = 1; i <= 3; i++) {
        for (int i = 0; i < 4; j++) {
            temp[i] += 42;
        }
        Cx = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
        Cy = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
        Cz = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);

        Cx_Mercury = Cx_Mercury.concatenar(Cx);
        Cy_Mercury = Cy_Mercury.concatenar(Cy);
        Cz_Mercury = Cz_Mercury.concatenar(Cz);

    }

    if (0 <= dt && dt <= 8) {
        j = 0;
        Mjd0 = t1;
    } else if (8 < dt && dt <= 16) {
        j = 1;
        Mjd0 = t1 + 8 * j;
    } else if (16 < dt && dt <= 24) {
        j = 2;
        Mjd0 = t1 + 8 * j;
    } else if (24 < dt && dt <= 32) {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }


    CxAux = Cx_Mercury.getPrimeraFil(14 * j + 1, 14 * j + 14);
    CyAux = Cy_Mercury.getPrimeraFil(14 * j + 1, 14 * j + 14);
    CzAux = Cz_Mercury.getPrimeraFil(14 * j + 1, 14 * j + 14);
    r_Mercury = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0 + 8, CxAux, CyAux, CzAux).transpuesta() * 1e3;


    //------------------VENUS----------------------


    temp[0] = 171;
    temp[1] = 181;
    temp[2] = 191;
    temp[3] = 201;

    Matrix Cx_Venus = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
    Matrix Cy_Venus = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
    Matrix Cz_Venus = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);
    for (int i = 1; i <= 2; i++) {
        for (int j = 0; j < 4; j++) {
            temp[j] += 33;
        }

        Cx = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
        Cy = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
        Cz = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);

        Cx_Venus = Cx_Venus.concatenar(Cx);
        Cy_Venus = Cy_Venus.concatenar(Cy);
        Cz_Venus = Cz_Venus.concatenar(Cz);
    }

    if (0 <= dt && dt <= 8) {
        j = 0;
        Mjd0 = t1;
    } else if (8 < dt && dt <= 16) {
        j = 1;
        Mjd0 = t1 + 8 * j;
    } else if (16 < dt && dt <= 24) {
        j = 2;
        Mjd0 = t1 + 8 * j;
    } else if (24 < dt && dt <= 32) {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }
    CxAux = Cx_Venus.getPrimeraFil(8 * j + 1, 8 * j + 8);
    CyAux = Cy_Venus.getPrimeraFil(8 * j + 1, 8 * j + 8);
    CzAux = Cz_Venus.getPrimeraFil(8 * j + 1, 8 * j + 8);
    r_Venus = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0 + 32, CxAux, CyAux, CzAux).transpuesta() * 1e3;


    //------------------MARS----------------------

    temp[0] = 309;
    temp[1] = 320;
    temp[2] = 331;
    temp[3] = 342;

    Matrix Cx_Mars = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
    Matrix Cy_Mars  = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
    Matrix Cz_Mars  = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);

    j = 0;
    Mjd0 = t1;

    CxAux = Cx_Mars.getPrimeraFil(11 * j + 1, 11* j + 11);
    CyAux = Cy_Mars.getPrimeraFil(11 * j + 1, 11* j + 11);
    CzAux = Cz_Mars.getPrimeraFil(11 * j + 1, 11 * j + 11);
    r_Mars = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 32, CxAux,CyAux, CzAux).transpuesta() *1e3;


    //------------------JUPITER----------------------

    temp[0] = 342;
    temp[1] = 350;
    temp[2] = 358;
    temp[3] = 366;

    Matrix Cx_Jupiter = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
    Matrix Cy_Jupiter = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
    Matrix Cz_Jupiter = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);

    j = 0;
    Mjd0 = t1;

    CxAux = Cx_Jupiter.getPrimeraFil(8 * j + 1, 8 * j + 8);
    CyAux = Cy_Jupiter.getPrimeraFil(8 * j + 1, 8 * j + 8);
    CzAux = Cz_Jupiter.getPrimeraFil(8 * j + 1, 8 * j + 8);
    r_Jupiter = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0 + 32, CxAux, CyAux, CzAux).transpuesta() * 1e3;


    //------------------SATURN----------------------

    temp[0] = 366;
    temp[1] = 373;
    temp[2] = 380;
    temp[3] = 387;

    Matrix Cx_Saturn = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
    Matrix Cy_Saturn = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
    Matrix Cz_Saturn = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);

    j = 0;
    Mjd0 = t1;

    CxAux = Cx_Saturn.getPrimeraFil(7 * j + 1, 7 * j + 7);
    CyAux = Cy_Saturn.getPrimeraFil(7 * j + 1, 7 * j + 7);
    CzAux = Cz_Saturn.getPrimeraFil(7 * j + 1, 7 * j + 7);
    r_Saturn = Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0 + 32, CxAux, CyAux, CzAux).transpuesta() * 1e3;


    //------------------URANUS----------------------


    temp[0] = 387;
    temp[1] = 393;
    temp[2] = 399;
    temp[3] = 405;


    Matrix Cx_Uranus = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
    Matrix Cy_Uranus = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
    Matrix Cz_Uranus = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);

    j = 0;
    Mjd0 = t1;

    CxAux = Cx_Uranus.getPrimeraFil(6 * j + 1, 6 * j + 6);
    CyAux = Cy_Uranus.getPrimeraFil(6 * j + 1, 6 * j + 6);
    CzAux = Cz_Uranus.getPrimeraFil(6 * j + 1, 6 * j + 6);

    r_Uranus = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, CxAux,CyAux, CzAux).transpuesta() * 1e3;


    //------------------NEPTUNE----------------------


    temp[0] = 405;
    temp[1] = 411;
    temp[2] = 417;
    temp[3] = 423;

    Matrix Cx_Neptune = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
    Matrix Cy_Neptune = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
    Matrix Cz_Neptune = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);

    j = 0;
    Mjd0 = t1;

    CxAux = Cx_Neptune.getPrimeraFil(6 * j + 1, 6 * j + 6);
    CyAux = Cy_Neptune.getPrimeraFil(6 * j + 1, 6 * j + 6);
    CzAux = Cz_Neptune.getPrimeraFil(6 * j + 1, 6 * j + 6);
    r_Neptune =  Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, CxAux,CyAux, CzAux).transpuesta() * 1e3;


    //------------------PLUTO----------------------


    temp[0] = 423;
    temp[1] = 429;
    temp[2] = 435;
    temp[3] = 441;

    Matrix Cx_Pluto = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
    Matrix Cy_Pluto = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
    Matrix Cz_Pluto = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);

    j = 0;
    Mjd0 = t1;

    CxAux = Cx_Neptune.getPrimeraFil(6 * j + 1, 6 * j + 6);
    CyAux = Cy_Neptune.getPrimeraFil(6 * j + 1, 6 * j + 6);
    CzAux = Cz_Neptune.getPrimeraFil(6 * j + 1, 6 * j + 6);
    r_Pluto = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, CxAux,CyAux, CzAux).transpuesta() * 1e3;

/*
    //------------------NUTATIONS----------------------

    temp[0] = 819;
    temp[1] = 829;
    temp[2] = 839;

    Matrix Cx_Nutations = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
    Matrix Cy_Nutations = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);

    for (int i = 1; i <= 3; i++) {
        for (int j = 0; j < 3; j++) {
            temp[i] += 20;
        }
        Cx = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
        Cy = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);

        Cx_Nutations = Cx_Nutations.concatenar(Cx);
        Cy_Nutations = Cy_Nutations.concatenar(Cy);
    }

    if (0 <= dt && dt <= 8) {
        j = 0;
        Mjd0 = t1;
    } else if (8 < dt && dt <= 16) {
        j = 1;
        Mjd0 = t1 + 8 * j;
    } else if (16 < dt && dt <= 24) {
        j = 2;
        Mjd0 = t1 + 8 * j;
    } else if (24 < dt && dt <= 32) {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }

    CxAux = Cx_Nutations.getPrimeraFil(10 * j + 1, 10 * j + 10);
    CyAux = Cy_Nutations.getPrimeraFil(10 * j + 1, 10 * j + 10);
    CzAux = Matrix(10,1);
    Matrix Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 8, CxAux,CyAux,CzAux).transpuesta();



    temp[0] = 899;
    temp[1] = 909;
    temp[2] = 919;
    temp[2] = 929;

    Matrix Cx_Librations = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
    Matrix Cy_Librations = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
    Matrix Cz_Librations = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);

    for (i = 1; i <= 3; i++) {

        for (int j = 0; j < 4; j++) {
            temp[i] += 30;
        }
        Cx = PCtemp.getPrimeraFil(temp[0], temp[1] - 1);
        Cy = PCtemp.getPrimeraFil(temp[1], temp[2] - 1);
        Cy = PCtemp.getPrimeraFil(temp[2], temp[3] - 1);

        Cx_Librations = Cx_Librations.concatenar(Cx);
        Cy_Librations = Cy_Librations.concatenar(Cy);
        Cz_Librations = Cz_Librations.concatenar(Cz);
    }
    if (0 <= dt && dt <= 8) {
        j = 0;
        Mjd0 = t1;
    } else if (8 < dt && dt <= 16) {
        j = 1;
        Mjd0 = t1 + 8 * j;
    } else if (16 < dt && dt <= 24) {
        j = 2;
        Mjd0 = t1 + 8 * j;
    } else if (24 < dt && dt <= 32) {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }

    CxAux = Cx_Librations .getPrimeraFil(10 * j + 1, 10 * j + 10);
    CyAux = Cy_Librations .getPrimeraFil(10 * j + 1, 10 * j + 10);
    CyAux = Cz_Librations .getPrimeraFil(10 * j + 1, 10 * j + 10);
    Matrix Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 8, CxAux,CyAux, CzAux).transpuesta();


 */

    double EMRAT = 81.30056907419062; // DE430
    double EMRAT1 = 1 / (1 + EMRAT);
    r_Earth = r_Earth -  r_Moon * EMRAT1 ;
    r_Mercury = r_Mercury - r_Earth;
    r_Venus = r_Venus - r_Earth;
    r_Mars = r_Mars - r_Earth;
    r_Jupiter = r_Jupiter - r_Earth;
    r_Saturn = r_Saturn - r_Earth;
    r_Uranus = r_Uranus - r_Earth;
    r_Neptune = r_Neptune - r_Earth;
    r_Pluto = r_Pluto - r_Earth;
    r_Sun = r_Sun - r_Earth;

}
