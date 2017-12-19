import math
import numpy as np

#已经验证没有问题
def sate_posivelo(t, orbit):
    GM = 3.986005e14  # 地球引力常数
    wie = 7.2921151467e-5  # 地球自转角速度
    Crs = orbit[0]
    delt_n = orbit[1]
    M0 = orbit[2]
    Cuc = orbit[3]
    e = orbit[4]
    Cus = orbit[5]
    a = pow(orbit[6], 2)
    toe = orbit[7]
    Cic = orbit[8]
    omega0 = orbit[9]
    Cis = orbit[10]
    i0 = orbit[11]
    Crc = orbit[12]
    w0 = orbit[13]
    delt_omega = orbit[14]
    d_i = orbit[15]
    af0 = orbit[16]
    af1 = orbit[17]
    af2 = orbit[18]
    toc = orbit[19]

    n = math.sqrt(GM / (a * a * a)) + delt_n
    delt_t = t - toe
    M = M0 + n * delt_t

    E0 = M
    E = E0
    delt_M = 1
    count = 0
    while abs(delt_M) > (1e-12):
        M1 = E0 - e * math.sin(E0)
        delt_M = M - M1
        delt_E = delt_M / (1 - e * math.cos(E0))
        E = E0 + delt_E
        E0 = E
        print(count, " ", (abs(delt_M)))
        count += 1

    r = a * (1 - e * math.cos(E))
    f = 2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(E / 2))
    fai0 = w0 + f
    delt_u = Cus * math.sin(2 * fai0) + Cuc * math.cos(2 * fai0)
    delt_r = Crs * math.sin(2 * fai0) + Crc * math.cos(2 * fai0)
    delt_i = Cis * math.sin(2 * fai0) + Cic * math.cos(2 * fai0)
    fai = fai0 + delt_u
    r = r + delt_r
    i = i0 + delt_i + d_i * delt_t
    na = omega0 + (delt_omega - wie) * delt_t - wie * toe
    X0 = r * np.mat([[math.cos(fai), math.sin(fai), 0]]).T

    C = np.matrix([[math.cos(na), -math.sin(na) * math.cos(i), math.sin(na) * math.sin(i)],
                   [math.sin(na), math.cos(na) * math.cos(i), -math.cos(na) * math.sin(i)],
                   [0, math.sin(i), math.cos(i)]])

    satp = np.dot(C, X0)
    d_E = n / (1 - e * math.cos(E))
    d_fai0 = d_E * math.sqrt((1 + e) / (1 - e)) * pow((math.cos(f / 2)), 2) / pow((math.cos(E / 2)), 2)
    d_u = 2 * d_fai0 * (Cus * math.cos(2 * fai0) - Cuc * math.sin(2 * fai0))
    d_r = 2 * d_fai0 * (Crs * math.cos(2 * fai0) - Crc * math.sin(2 * fai0))
    d_ii = 2 * d_fai0 * (Cis * math.cos(2 * fai0) - Cic * math.sin(2 * fai0))
    d_fai = d_fai0 + d_u
    d_r = a * e * math.sin(E) * d_E + d_r
    d_i = delt_i + d_ii
    d_na = delt_omega - wie
    v0 = d_r * np.matrix([[math.cos(fai), math.sin(fai), 0]]).transpose() + r * d_fai * np.matrix([
        -math.sin(fai), math.cos(fai), 0]).transpose()

    C2 = np.matrix([[-d_na * math.sin(na), -d_na * math.cos(na) * math.cos(i) + d_i * math.sin(na) * math.sin(i),
                     d_na * math.cos(na) * math.sin(i) + d_i * math.sin(na) * math.cos(i)], [
                        d_na * math.cos(na), -d_na * math.sin(na) * math.cos(i) - d_i * math.cos(na) * math.sin(i),
                        d_na * math.sin(na) * math.sin(i) - d_i * math.cos(na) * math.cos(i)
                    ], [0, d_i * math.cos(i), -d_i * math.sin(i)]])
    satv = C * v0 + C2 * X0

    return satp, satv


if __name__ == '__main__':

    #ans = 9.907832e+06 1.585084e+07 1.884309e+07

    #ans = -7.174652e+02 2.295145e+03  - 1.552965e+03
    time = 86400
    para = np.array([3.568750e+01,
                     4.572690e-09,
                     1.155528e-01,
                     1.844019e-06,
                     6.636133e-04,
                     1.122989e-05,
                     5.153600e+03,
                     86400,
                     -9.313226e-09,
                     -1.313863e+00,
                     7.450581e-09,
                     9.530541e-01,
                     1.576563e+02,
                     1.968896e+00,
                     -8.057121e-09,
                     -7.607460e-11,
                     4.069614e-04,
                     7.617018e-12,
                     0,
                     0])

    satp, satv = sate_posivelo(time, para)
    print(satp)
    print(satv)
