import math
import numpy as np


def sate_posivelo(t, orbit):
    GM = 3.986005e14;  # 地球引力常数
    wie = 7.2921151467e-5;  # 地球自转角速度
    Crs = orbit[1];
    delt_n = orbit[2];
    M0 = orbit[3];
    Cuc = orbit[4];
    e = orbit[5];
    Cus = orbit[6];
    a = np.square(orbit[7]);
    toe = orbit[8];

    Cic = orbit[9];
    omega0 = orbit[10];
    Cis = orbit[11];
    i0 = orbit[12];

    Crc = orbit[13];
    w0 = orbit[14];
    delt_omega = orbit[15];
    d_i = orbit[16];
    af0 = orbit[17];
    af1 = orbit[18];
    af2 = orbit[19];
    toc = orbit[20];

    n = math.sqrt(GM / (a * a * a)) + delt_n;
    delt_t = t - toe;
    M = M0 + n * delt_t;

    E0 = M;
    E = E0;
    delt_M = 1;
    while ~(abs(delt_M) < (1e-12)):
        M1 = E0 - e * math.sin(E0);
        delt_M = M - M1;
        delt_E = delt_M / (1 - e * math.cos(E0));
        E = E0 + delt_E;
        E0 = E;

    r = a * (1 - e * math.cos(E));
    f = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(E / 2));
    fai0 = w0 + f;
    delt_u = Cus * math.sin(2 * fai0) + Cuc * math.cos(2 * fai0);
    delt_r = Crs * math.sin(2 * fai0) + Crc * math.cos(2 * fai0);
    delt_i = Cis * math.sin(2 * fai0) + Cic * math.cos(2 * fai0);
    fai = fai0 + delt_u;
    r = r + delt_r;
    i = i0 + delt_i + d_i * delt_t;
    na = omega0 + (delt_omega - wie) * delt_t - wie * toe;
    X0 = np.dot(r, np.array([math.cos(fai), math.sin(fai), 0]))

    C = np.matrix([[math.cos(na), -math.sin(na) * math.cos(i), math.sin(na) * math.sin(i)],
                   [math.sin(na), math.cos(na) * math.cos(i), -math.cos(na) * math.sin(i)],
                   [0, math.sin(i), math.cos(i)]])

    satp = np.dot(C, X0);
    d_E = n / (1 - e * math.cos(E));
    d_fai0 = d_E * math.sqrt((1 + e) / (1 - e)) * (math.cos(f / 2)) ^ 2 / (math.cos(E / 2)) ^ 2;
    d_u = 2 * d_fai0 * (Cus * math.cos(2 * fai0) - Cuc * math.sin(2 * fai0));
    d_r = 2 * d_fai0 * (Crs * math.cos(2 * fai0) - Crc * math.sin(2 * fai0));
    d_ii = 2 * d_fai0 * (Cis * math.cos(2 * fai0) - Cic * math.sin(2 * fai0));
    d_fai = d_fai0 + d_u;
    d_r = a * e * math.sin(E) * d_E + d_r;
    d_i = delt_i + d_ii;
    d_na = delt_omega - wie;
    v0 = d_r * np.matrix([[math.cos(fai), math.sin(fai), 0]]).transpose() + r * d_fai * np.matrix([
        -math.sin(fai), math.cos(fai), 0]).transpose();

    C2 = np.matrix([[-d_na * math.sin(na), -d_na * math.cos(na) * math.cos(i) + d_i * math.sin(na) * sin(i),
                     d_na * math.cos(na) * math.sin(i) + d_i * math.sin(na) * math.cos(i)], [
                        d_na * math.cos(na), -d_na * math.sin(na) * math.cos(i) - d_i * math.cos(na) * math.sin(i),
                        d_na * math.sin(na) * math.sin(i) - d_i * math.cos(na) * math.cos(i)
                    ], [0, d_i * math.cos(i), -d_i * math.sin(i)]])
    satv = C * v0 + C2 * X0;

    return satp, satv
