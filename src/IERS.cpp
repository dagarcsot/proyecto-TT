
#include "../include/IERS.h"
#include "../include/Frac.h"
#include "../include/SAT_Const.h"
#include <cmath>
#include <iostream>

/*
 * IERS: Management of IERS time and polar motion data
 */

void IERS(Matrix &eop, double Mjd_UTC, char interp, double &x_pole, double &y_pole, double &UT1_UTC,
          double &LOD, double &dpsi, double &deps, double &dx_pole, double &dy_pole, double &TAI_UTC) {


    if (interp != 'n' && interp != 'l') {
        interp = 'n';
    }
    int i;
    double mjd;
    if (interp == 'l') {
        // linear interpolation
        mjd = floor(Mjd_UTC);

        int index = -1;
        double aux;
        for (int i = 0; i < eop.getNumCol(); ++i) {
            aux = floor(eop(4, i));
            if (mjd == aux) {
                index = i;
                break;
            }
        }

        if (index == -1) {
            std::cout<<"Mjd_UTC was not found on eop";
            exit(-1);
        }

        Matrix preeop = eop.getCol(index);
        Matrix nexteop = eop.getCol(index + 1);
        preeop.print();
        nexteop.print();

        double mfme = 1440 * Frac(Mjd_UTC);
        double fixf = mfme / 1440;
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole = preeop(1, 5) + (nexteop(1, 5) - preeop(1, 5)) * fixf;
        y_pole = preeop(1, 6) + (nexteop(1, 6) - preeop(1, 6)) * fixf;
        UT1_UTC = preeop(1, 7) + (nexteop(1, 7) - preeop(1, 7)) * fixf;
        LOD = preeop(1, 8) + (nexteop(1, 8) - preeop(1, 8)) * fixf;
        dpsi = preeop(1, 9) + (nexteop(1, 9) - preeop(1, 9)) * fixf;
        deps = preeop(1, 10) + (nexteop(1, 10) - preeop(1, 10)) * fixf;
        dx_pole = preeop(1, 11) + (nexteop(1, 11) - preeop(1, 11)) * fixf;
        dy_pole = preeop(1, 12) + (nexteop(1, 12) - preeop(1, 12)) * fixf;
        TAI_UTC = preeop(1, 13);

        x_pole /= Arcs; // Pole coordinate [rad]
        y_pole /= Arcs; // Pole coordinate [rad]
        dpsi /= Arcs;
        deps /= Arcs;
        dx_pole /= Arcs; // Pole coordinate [rad]
        dy_pole /= Arcs; // Pole coordinate [rad]

    } else {
        mjd = floor(Mjd_UTC);

        int index = -1;
        for (int i = 0; i < eop.getNumCol(); ++i) {
            if (mjd == eop(4, i)) {
                index = i;
                break;
            }
        }

        if (index == -1) {
            std::cout << "Mjd_UTC was not found on eop";
            exit(-1);
        }

        Matrix preeop = eop.getCol(index);
        Matrix nexteop = eop.getCol(index);

        //Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])

        x_pole = eop(1, 5) / Arcs; // Pole coordinate [rad]
        y_pole = eop(1, 6) / Arcs; // Pole coordinate [rad]
        UT1_UTC = eop(1, 7);      // UT1-UTC time difference [s]
        LOD = eop(1, 8);          //  Length of day [s]
        dpsi = eop(1, 9) / Arcs;
        deps = eop(1, 10) / Arcs;
        dx_pole = eop(1, 11) / Arcs; // Pole coordinate [rad]
        dy_pole = eop(1, 12) / Arcs; // Pole coordinate [rad]
        TAI_UTC = eop(1, 13);       // TAI-UTC time difference [s]
    }

}

