import math

from com.nuaa.gps.nrc.common_func.geo_ecef import geo_ccef
from com.nuaa.gps.nrc.common_func.trace_velo import trace_velo

import numpy as np
from math import cos
from math import sin
from math import pi


def user_positions(t, T, atti, atti_rate, veloB, acceB, posi):
    Re = 6378137.0;
    f = 1 / 298.257223563;

    long = posi(1, 1) * math.pi / 180.0;
    lati = posi(2, 1) * math.pi / 180.0;
    heig = posi(3, 1);

    Rm = Re * (1 - 2 * f + 3 * f * math.sin(lati) * math.sin(lati));
    Rn = Re * (1 + f * math.sin(lati) * math.sin(lati));

    t, atti, atti_rate, veloB, acceB = trace_velo(t, T, atti, atti_rate, veloB, acceB)
    roll = atti(1, 1) * math.pi / 180.0;
    pitch = atti(2, 1) * math.pi / 180.0;
    head = atti(3, 1) * math.pi / 180.0;

    Cbn = np.matrix([[cos(roll) * cos(head) + sin(roll) * sin(pitch) * sin(head),
                      -cos(roll) * sin(head) + sin(roll) * sin(pitch) * cos(head), -sin(roll) * cos(pitch)],
                     [cos(pitch) * sin(head), cos(pitch) * cos(head), sin(pitch)],
                     [sin(roll) * cos(head) - cos(roll) * sin(pitch) * sin(head),
                      -sin(roll) * sin(head) - cos(roll) * sin(pitch) * cos(head), cos(roll) * cos(pitch)]])

    veloN = Cbn.transpose() * veloB;
    acceN = Cbn.transpose() * acceB;
    Ve = veloN(1, 1);
    Vn = veloN(2, 1);
    Vu = veloN(3, 1);

    heig = heig + T * Vu;
    lati = lati + T * (Vn / (Rm + heig));
    long = long + T * (Ve / ((Rn + heig) * cos(lati)));

    posi[1, 1] = long * 180.0 / pi;
    posi[2, 1] = lati * 180.0 / pi;
    posi[3, 1] = heig;
    posiE = geo_ccef(posi);

    Cne = np.matrix([[-sin(long), cos(long), 0],
                     [-sin(lati) * cos(long), -sin(lati) * sin(long), cos(lati)],
                     [cos(lati) * cos(long), cos(lati) * sin(long), sin(lati)]])

    veloE = Cne.transpose() * veloN;

    Xr = posiE;
    Vr = veloE;

    return Xr, Vr
