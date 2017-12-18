import numpy as np


def trace_velo(t, T, atti, atti_rate, veloB, acceB):
    if t > 0:
        acceB[2, 1] = 0;
        veloB[2, 1] = 100;
    veloB[2, 1] = veloB[2, 1] + acceB[2, 1] * T;
    atti_rate = np.zeros[3, 1];

    atti[1, 1] = atti[1, 1] + atti_rate[1, 1] * T;
    atti[2, 1] = atti[2, 1] + atti_rate[2, 1] * T;
    atti[3, 1] = atti[3, 1] + atti_rate[3, 1] * T;
    return t, T, atti, atti_rate, veloB, acceB
